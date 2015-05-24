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
              rVd(1) = root_3*L*u_V_SDS
              rVd(2) = -(root_3*M*u_V_SDS)
              rVd(3) = 0
              rVd_t = (root_3*(L - M)*(L + M)*u_Vd_SDS)/2.
            end select
        end select
      case (ORB_F)
        select case (orb_dir_i)
          case (1)
            select case (orb_dir_j)
              case (1)
              rVd(1) = 3*Sqrt(2.5)*L*M*u_V_SFS
              rVd(2) = (-3*Sqrt(2.5)*(-1 + 2*Msq + Nsq)*u_V_SFS)/2.
              rVd(3) = 0
              rVd_t = -(Sqrt(2.5)*M*(-3*Lsq + Msq)*u_Vd_SFS)/2.
              case (2)
              rVd(1) = sqrt(15.)*M*N*u_V_SFS
              rVd(2) = sqrt(15.)*L*N*u_V_SFS
              rVd(3) = sqrt(15.)*L*M*u_V_SFS
              rVd_t = sqrt(15.)*L*M*N*u_Vd_SFS
              case (3)
              rVd(1) = 0
              rVd(2) = (Sqrt(1.5)*(-1 + 5*Nsq)*u_V_SFS)/2.
              rVd(3) = 5*Sqrt(1.5)*M*N*u_V_SFS
              rVd_t = (Sqrt(1.5)*M*(-1 + 5*Nsq)*u_Vd_SFS)/2.
              case (4)
              rVd(1) = 0
              rVd(2) = 0
              rVd(3) = (3*(-1 + 5*Nsq)*u_V_SFS)/2.
              rVd_t = (N*(-3 + 5*Nsq)*u_Vd_SFS)/2.
              case (5)
              rVd(1) = (Sqrt(1.5)*(-1 + 5*Nsq)*u_V_SFS)/2.
              rVd(2) = 0
              rVd(3) = 5*Sqrt(1.5)*L*N*u_V_SFS
              rVd_t = (Sqrt(1.5)*L*(-1 + 5*Nsq)*u_Vd_SFS)/2.
              case (6)
              rVd(1) = sqrt(15.)*L*N*u_V_SFS
              rVd(2) = -(sqrt(15.)*M*N*u_V_SFS)
              rVd(3) = (sqrt(15.)*(L - M)*(L + M)*u_V_SFS)/2.
              rVd_t = (sqrt(15.)*(L - M)*(L + M)*N*u_Vd_SFS)/2.
              case (7)
              rVd(1) = (3*Sqrt(2.5)*(Lsq - Msq)*u_V_SFS)/2.
              rVd(2) = -3*Sqrt(2.5)*L*M*u_V_SFS
              rVd(3) = 0
              rVd_t = (Sqrt(2.5)*(L**3 - 3*L*Msq)*u_Vd_SFS)/2.
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
              rVd(1) = 2*L*u_V_PPS
              rVd(2) = 2*M*u_V_PPP
              rVd(3) = 2*N*u_V_PPP
              rVd_t = Lsq*u_Vd_PPS + (Msq + Nsq)*u_Vd_PPP
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
              rVd_t = (M*((-1 + 3*Nsq)*u_Vd_PDS - 2*root_3*Nsq*u_Vd_PDP))/2.
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
              rVd(1) = L*N*(root_3*u_V_PDS - 2*u_V_PDP)
              rVd(2) = M*N*(-(root_3*u_V_PDS) + 2*u_V_PDP)
              rVd(3) = ((L - M)*(L + M)*(root_3*u_V_PDS - 2*u_V_PDP))/2.
              rVd_t = ((L - M)*(L + M)*N*(root_3*u_Vd_PDS - 2*u_Vd_PDP))/2.
            end select
          case (3)
            select case (orb_dir_j)
              case (1)
              rVd(1) = 2*root_3*L*M*u_V_PDS
              rVd(2) = root_3*Lsq*u_V_PDS + (-1 + 6*Msq + 2*Nsq)*u_V_PDP
              rVd(3) = 4*M*N*u_V_PDP
              rVd_t = root_3*Lsq*M*u_Vd_PDS + M*(-1 + 2*Msq + 2*Nsq)*u_Vd_PDP
              case (2)
              rVd(1) = M*N*(root_3*u_V_PDS - 2*u_V_PDP)
              rVd(2) = L*N*(root_3*u_V_PDS - 2*u_V_PDP)
              rVd(3) = L*M*(root_3*u_V_PDS - 2*u_V_PDP)
              rVd_t = L*M*N*(root_3*u_Vd_PDS - 2*u_Vd_PDP)
              case (3)
              rVd(1) = ((-1 + 3*Nsq)*u_V_PDS - 2*root_3*Nsq*u_V_PDP)/2.
              rVd(2) = 0
              rVd(3) = L*N*(3*u_V_PDS - 2*root_3*u_V_PDP)
              rVd_t = (L*((-1 + 3*Nsq)*u_Vd_PDS - 2*root_3*Nsq*u_Vd_PDP))/2.
              case (4)
              rVd(1) = 2*root_3*L*N*u_V_PDS
              rVd(2) = 4*M*N*u_V_PDP
              rVd(3) = root_3*Lsq*u_V_PDS + (-1 + 2*Msq + 6*Nsq)*u_V_PDP
              rVd_t = N*(root_3*Lsq*u_Vd_PDS + (-1 + 2*Msq + 2*Nsq)*u_Vd_PDP)
              case (5)
              rVd(1) = (-(root_3*(-1 + 2*Msq + Nsq)*u_V_PDS) +  &
      2*(2*Msq + Nsq)*u_V_PDP)/2.
              rVd(2) = L*M*(-2*root_3*u_V_PDS + 4*u_V_PDP)
              rVd(3) = L*N*(-(root_3*u_V_PDS) + 2*u_V_PDP)
              rVd_t = -(L*(root_3*(-1 + 2*Msq + Nsq)*u_Vd_PDS -  &
         2*(2*Msq + Nsq)*u_Vd_PDP))/2.
            end select
        end select
      case (ORB_F)
        select case (orb_dir_i)
          case (1)
            select case (orb_dir_j)
              case (1)
              rVd(1) = 0
              rVd(2) = (M*(sqrt(10.)*(3 - 8*Msq - 3*Nsq)*u_V_PFS +  &
        sqrt(15.)*(-5 + 8*Msq + 3*Nsq)*u_V_PFP))/2.
              rVd(3) = (-(sqrt(15.)*N*u_V_PFP) +  &
      Msq*N*(-3*sqrt(10.)*u_V_PFS + 3*sqrt(15.)*u_V_PFP))/2.
              rVd_t = (sqrt(5.)*(-(root_3*(-1 + Nsq)*u_Vd_PFP) +  &
        M**4*(-4*sqrt(2.)*u_Vd_PFS + 4*root_3*u_Vd_PFP) +  &
        Msq*(-3*sqrt(2.)*(-1 + Nsq)*u_Vd_PFS +  &
           root_3*(-5 + 3*Nsq)*u_Vd_PFP)))/4.
              case (2)
              rVd(1) = Sqrt(2.5)*N*u_V_PFP +  &
    Msq*N*(sqrt(15.)*u_V_PFS - 3*Sqrt(2.5)*u_V_PFP)
              rVd(2) = L*M*N*(2*sqrt(15.)*u_V_PFS - 3*sqrt(10.)*u_V_PFP)
              rVd(3) = L*(Sqrt(2.5)*u_V_PFP +  &
      Msq*(sqrt(15.)*u_V_PFS - 3*Sqrt(2.5)*u_V_PFP))
              rVd_t = (L*N*(sqrt(10.)*u_Vd_PFP +  &
        Msq*(2*sqrt(15.)*u_Vd_PFS - 3*sqrt(10.)*u_Vd_PFP)))/2.
              case (3)
              rVd(1) = 0
              rVd(2) = (M*(sqrt(6.)*(-1 + 5*Nsq)*u_V_PFS +  &
        (1 - 15*Nsq)*u_V_PFP))/2.
              rVd(3) = (5*N*(Msq*(sqrt(6.)*u_V_PFS - 3*u_V_PFP) + u_V_PFP))/2.
              rVd_t = ((-1 + 5*Nsq)*u_Vd_PFP +  &
      Msq*(sqrt(6.)*(-1 + 5*Nsq)*u_Vd_PFS + (1 - 15*Nsq)*u_Vd_PFP))/4.
              case (4)
              rVd(1) = 0
              rVd(2) = (N*(2*(-3 + 5*Nsq)*u_V_PFS +  &
        sqrt(6.)*(1 - 5*Nsq)*u_V_PFP))/4.
              rVd(3) = (M*(6*(-1 + 5*Nsq)*u_V_PFS +  &
        sqrt(6.)*(1 - 15*Nsq)*u_V_PFP))/4.
              rVd_t = (M*N*(2*(-3 + 5*Nsq)*u_Vd_PFS +  &
        sqrt(6.)*(1 - 5*Nsq)*u_Vd_PFP))/4.
              case (5)
              rVd(1) = (M*(sqrt(6.)*(-1 + 5*Nsq)*u_V_PFS +  &
        (1 - 15*Nsq)*u_V_PFP))/4.
              rVd(2) = (L*(sqrt(6.)*(-1 + 5*Nsq)*u_V_PFS +  &
        (1 - 15*Nsq)*u_V_PFP))/4.
              rVd(3) = (5*L*M*N*(sqrt(6.)*u_V_PFS - 3*u_V_PFP))/2.
              rVd_t = (L*M*(sqrt(6.)*(-1 + 5*Nsq)*u_Vd_PFS + (1 - 15*Nsq)*u_Vd_PFP))/ &
    4.
              case (6)
              rVd(1) = 0
              rVd(2) = (N*(-2*sqrt(15.)*(-1 + 6*Msq + Nsq)*u_V_PFS +  &
        sqrt(10.)*(-5 + 18*Msq + 3*Nsq)*u_V_PFP))/4.
              rVd(3) = (M*(-2*sqrt(15.)*(-1 + 2*Msq + 3*Nsq)*u_V_PFS +  &
        sqrt(10.)*(-5 + 6*Msq + 9*Nsq)*u_V_PFP))/4.
              rVd_t = (M*N*(-2*sqrt(15.)*(-1 + 2*Msq + Nsq)*u_Vd_PFS +  &
        sqrt(10.)*(-5 + 6*Msq + 3*Nsq)*u_Vd_PFP))/4.
              case (7)
              rVd(1) = (M*(-(sqrt(10.)*(-1 + 4*Msq + Nsq)*u_V_PFS) +  &
        sqrt(15.)*(-3 + 4*Msq + Nsq)*u_V_PFP))/4.
              rVd(2) = (L*(-(sqrt(10.)*(-1 + 12*Msq + Nsq)*u_V_PFS) +  &
        sqrt(15.)*(-3 + 12*Msq + Nsq)*u_V_PFP))/4.
              rVd(3) = (L*M*N*(-(sqrt(10.)*u_V_PFS) + sqrt(15.)*u_V_PFP))/2.
              rVd_t = (L*M*(-(sqrt(10.)*(-1 + 4*Msq + Nsq)*u_Vd_PFS) +  &
        sqrt(15.)*(-3 + 4*Msq + Nsq)*u_Vd_PFP))/4.
            end select
          case (2)
            select case (orb_dir_j)
              case (1)
              rVd(1) = 0
              rVd(2) = (-3*N*(-1 + 4*Msq + Nsq)* &
      (sqrt(10.)*u_V_PFS - sqrt(15.)*u_V_PFP))/4.
              rVd(3) = -(M*(-3 + 4*Msq + 9*Nsq)* &
       (sqrt(10.)*u_V_PFS - sqrt(15.)*u_V_PFP))/4.
              rVd_t = -(M*N*(-3 + 4*Msq + 3*Nsq)* &
       (sqrt(10.)*u_Vd_PFS - sqrt(15.)*u_Vd_PFP))/4.
              case (2)
              rVd(1) = (M*(2*sqrt(15.)*Nsq*u_V_PFS +  &
        sqrt(10.)*(1 - 3*Nsq)*u_V_PFP))/2.
              rVd(2) = (L*(2*sqrt(15.)*Nsq*u_V_PFS +  &
        sqrt(10.)*(1 - 3*Nsq)*u_V_PFP))/2.
              rVd(3) = L*M*N*(2*sqrt(15.)*u_V_PFS - 3*sqrt(10.)*u_V_PFP)
              rVd_t = (L*M*(2*sqrt(15.)*Nsq*u_Vd_PFS +  &
        sqrt(10.)*(1 - 3*Nsq)*u_Vd_PFP))/2.
              case (3)
              rVd(1) = 0
              rVd(2) = (N*(sqrt(6.)*(-1 + 5*Nsq)*u_V_PFS +  &
        (11 - 15*Nsq)*u_V_PFP))/4.
              rVd(3) = (M*(sqrt(6.)*(-1 + 15*Nsq)*u_V_PFS +  &
        (11 - 45*Nsq)*u_V_PFP))/4.
              rVd_t = (M*N*(sqrt(6.)*(-1 + 5*Nsq)*u_Vd_PFS +  &
        (11 - 15*Nsq)*u_Vd_PFP))/4.
              case (4)
              rVd(1) = 0
              rVd(2) = 0
              rVd(3) = N*((-3 + 10*Nsq)*u_V_PFS + sqrt(6.)*(3 - 5*Nsq)*u_V_PFP)
              rVd_t = (2*Nsq*(-3 + 5*Nsq)*u_Vd_PFS +  &
      sqrt(6.)*(-1 + 6*Nsq - 5*N**4)*u_Vd_PFP)/4.
              case (5)
              rVd(1) = (N*(sqrt(6.)*(-1 + 5*Nsq)*u_V_PFS +  &
        (11 - 15*Nsq)*u_V_PFP))/4.
              rVd(2) = 0
              rVd(3) = (L*(sqrt(6.)*(-1 + 15*Nsq)*u_V_PFS +  &
        (11 - 45*Nsq)*u_V_PFP))/4.
              rVd_t = (L*N*(sqrt(6.)*(-1 + 5*Nsq)*u_Vd_PFS +  &
        (11 - 15*Nsq)*u_Vd_PFP))/4.
              case (6)
              rVd(1) = L*(Sqrt(2.5)*u_V_PFP +  &
      Nsq*(sqrt(15.)*u_V_PFS - 3*Sqrt(2.5)*u_V_PFP))
              rVd(2) = -(M*(sqrt(10.)*u_V_PFP +  &
         Nsq*(2*sqrt(15.)*u_V_PFS - 3*sqrt(10.)*u_V_PFP)))/2.
              rVd(3) = -(N*(-1 + 2*Msq + Nsq)* &
       (2*sqrt(15.)*u_V_PFS - 3*sqrt(10.)*u_V_PFP))/2.
              rVd_t = (sqrt(5.)*(L - M)*(L + M)* &
      (2*root_3*Nsq*u_Vd_PFS + sqrt(2.)*(1 - 3*Nsq)*u_Vd_PFP))/4.
              case (7)
              rVd(1) = -(N*(-1 + 4*Msq + Nsq)* &
       (sqrt(10.)*u_V_PFS - sqrt(15.)*u_V_PFP))/4.
              rVd(2) = L*M*N*(-2*sqrt(10.)*u_V_PFS + 2*sqrt(15.)*u_V_PFP)
              rVd(3) = -(L*(-1 + 4*Msq + 3*Nsq)* &
       (sqrt(10.)*u_V_PFS - sqrt(15.)*u_V_PFP))/4.
              rVd_t = -(L*N*(-1 + 4*Msq + Nsq)* &
       (sqrt(10.)*u_Vd_PFS - sqrt(15.)*u_Vd_PFP))/4.
            end select
          case (3)
            select case (orb_dir_j)
              case (1)
              rVd(1) = (M*(sqrt(10.)*(3 - 4*Msq - 3*Nsq)*u_V_PFS +  &
        sqrt(15.)*(-1 + 4*Msq + 3*Nsq)*u_V_PFP))/4.
              rVd(2) = (L*(-3*sqrt(10.)*(-1 + 4*Msq + Nsq)*u_V_PFS +  &
        sqrt(15.)*(-1 + 12*Msq + 3*Nsq)*u_V_PFP))/4.
              rVd(3) = (-3*L*M*N*(sqrt(10.)*u_V_PFS - sqrt(15.)*u_V_PFP))/2.
              rVd_t = (L*M*(sqrt(10.)*(3 - 4*Msq - 3*Nsq)*u_Vd_PFS +  &
        sqrt(15.)*(-1 + 4*Msq + 3*Nsq)*u_Vd_PFP))/4.
              case (2)
              rVd(1) = (L*M*N*(2*sqrt(15.)*(-1 + Nsq)*u_V_PFS +  &
        sqrt(10.)*(2 - 3*Nsq)*u_V_PFP))/(-1 + Nsq)
              rVd(2) = (sqrt(5.)*N*(-3*sqrt(2.)*Msq*u_V_PFP +  &
        Lsq*(2*root_3*(-1 + Nsq)*u_V_PFS +  &
           sqrt(2.)*(2 - 3*Nsq)*u_V_PFP)))/(2.*(-1 + Nsq))
              rVd(3) = (sqrt(5.)*M*(sqrt(2.)*Msq*(1 + Nsq)*u_V_PFP +  &
        Lsq*(2*root_3*(-1 + Nsq)**2*u_V_PFS +  &
           sqrt(2.)*(-2 + 7*Nsq - 3*N**4)*u_V_PFP)))/(2.*(-1 + Nsq)**2)
              rVd_t = -(sqrt(5.)*M*N*(sqrt(2.)*Msq*u_Vd_PFP +  &
         Lsq*(-2*root_3*(-1 + Nsq)*u_Vd_PFS +  &
            sqrt(2.)*(-2 + 3*Nsq)*u_Vd_PFP)))/(2.*(-1 + Nsq))
              case (3)
              rVd(1) = (M*(sqrt(6.)*(-1 + 5*Nsq)*u_V_PFS +  &
        (1 - 15*Nsq)*u_V_PFP))/4.
              rVd(2) = (L*(sqrt(6.)*(-1 + 5*Nsq)*u_V_PFS +  &
        (1 - 15*Nsq)*u_V_PFP))/4.
              rVd(3) = (5*L*M*N*(sqrt(6.)*u_V_PFS - 3*u_V_PFP))/2.
              rVd_t = (L*M*(sqrt(6.)*(-1 + 5*Nsq)*u_Vd_PFS + (1 - 15*Nsq)*u_Vd_PFP))/ &
    4.
              case (4)
              rVd(1) = (N*(2*(-3 + 5*Nsq)*u_V_PFS +  &
        sqrt(6.)*(1 - 5*Nsq)*u_V_PFP))/4.
              rVd(2) = 0
              rVd(3) = (L*(6*(-1 + 5*Nsq)*u_V_PFS +  &
        sqrt(6.)*(1 - 15*Nsq)*u_V_PFP))/4.
              rVd_t = (L*N*(2*(-3 + 5*Nsq)*u_Vd_PFS +  &
        sqrt(6.)*(1 - 5*Nsq)*u_Vd_PFP))/4.
              case (5)
              rVd(1) = (L*(sqrt(6.)*(1 - 6*Nsq + 5*N**4)*u_V_PFS +  &
        Nsq*(11 - 15*Nsq)*u_V_PFP))/(2.*(-1 + Nsq))
              rVd(2) = (M*(1 - 5*Nsq)*u_V_PFP)/(2.*(-1 + Nsq))
              rVd(3) = (N*(4*Msq*u_V_PFP +  &
        Lsq*(5*sqrt(6.)*(-1 + Nsq)**2*u_V_PFS +  &
           (-11 + 30*Nsq - 15*N**4)*u_V_PFP)))/(2.*(-1 + Nsq)**2)
              rVd_t = (sqrt(6.)*Lsq*(1 - 6*Nsq + 5*N**4)*u_Vd_PFS +  &
      (Lsq*Nsq*(11 - 15*Nsq) + Msq*(1 - 5*Nsq))*u_Vd_PFP)/ &
    (4.*(-1 + Nsq))
              case (6)
              rVd(1) = (N*(-2*sqrt(15.)*(-1 + 2*Msq + Nsq)*u_V_PFS +  &
        sqrt(10.)*(-1 + 6*Msq + 3*Nsq)*u_V_PFP))/4.
              rVd(2) = L*M*N*(-2*sqrt(15.)*u_V_PFS + 3*sqrt(10.)*u_V_PFP)
              rVd(3) = (L*(-2*sqrt(15.)*(-1 + 2*Msq + 3*Nsq)*u_V_PFS +  &
        sqrt(10.)*(-1 + 6*Msq + 9*Nsq)*u_V_PFP))/4.
              rVd_t = (L*N*(-2*sqrt(15.)*(-1 + 2*Msq + Nsq)*u_Vd_PFS +  &
        sqrt(10.)*(-1 + 6*Msq + 3*Nsq)*u_Vd_PFP))/4.
              case (7)
              rVd(1) = (sqrt(5.)*(6*L*Msq*(-1 + Nsq)* &
         (-(sqrt(2.)*u_V_PFS) + root_3*u_V_PFP) +  &
        4*L**3*(sqrt(2.)*(-1 + Nsq)*u_V_PFS - root_3*Nsq*u_V_PFP)))/ &
    (4.*(-1 + Nsq))
              rVd(2) = (sqrt(5.)*(4*root_3*M**3*u_V_PFP +  &
        6*Lsq*M*(-1 + Nsq)*(-(sqrt(2.)*u_V_PFS) + root_3*u_V_PFP)))/ &
    (4.*(-1 + Nsq))
              rVd(3) = (sqrt(15.)*(L**4 - M**4)*N*u_V_PFP)/(2.*(-1 + Nsq)**2)
              rVd_t = (sqrt(5.)*(root_3*M**4*u_Vd_PFP +  &
        3*Lsq*Msq*(-1 + Nsq)* &
         (-(sqrt(2.)*u_Vd_PFS) + root_3*u_Vd_PFP) +  &
        L**4*(sqrt(2.)*(-1 + Nsq)*u_Vd_PFS - root_3*Nsq*u_Vd_PFP)))/ &
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
              rVd(1) = root_3*L*u_V_SDS
              rVd(2) = -(root_3*M*u_V_SDS)
              rVd(3) = 0
              rVd_t = (root_3*(L - M)*(L + M)*u_Vd_SDS)/2.
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
              rVd(1) = -2*root_3*L*M*u_V_PDS
              rVd(2) = -(root_3*Lsq*u_V_PDS) + (1 - 6*Msq - 2*Nsq)*u_V_PDP
              rVd(3) = -4*M*N*u_V_PDP
              rVd_t = -(root_3*Lsq*M*u_Vd_PDS) + M*(1 - 2*Msq - 2*Nsq)*u_Vd_PDP
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
              rVd(1) = (-2*L*N*(root_3*(-1 + Nsq)*u_V_PDS +  &
        (1 - 2*Nsq)*u_V_PDP))/(-1 + Nsq)
              rVd(2) = (2*M*N*u_V_PDP)/(-1 + Nsq)
              rVd(3) = (-(Msq*(1 + Nsq)*u_V_PDP) +  &
      Lsq*(-(root_3*(-1 + Nsq)**2*u_V_PDS) +  &
         (1 - 5*Nsq + 2*N**4)*u_V_PDP))/(-1 + Nsq)**2
              rVd_t = (N*(Msq*u_Vd_PDP -  &
        Lsq*(root_3*(-1 + Nsq)*u_Vd_PDS + (1 - 2*Nsq)*u_Vd_PDP)))/ &
    (-1 + Nsq)
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
              rVd(1) = L*N*(-(root_3*u_V_PDS) + 2*u_V_PDP)
              rVd(2) = M*N*(root_3*u_V_PDS - 2*u_V_PDP)
              rVd(3) = -((L - M)*(L + M)*(root_3*u_V_PDS - 2*u_V_PDP))/2.
              rVd_t = -((L - M)*(L + M)*N*(root_3*u_Vd_PDS - 2*u_Vd_PDP))/2.
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
              rVd(3) = 2*N*(-u_V_DDP + u_V_DDD -  &
      Msq*(3*u_V_DDS - 4*u_V_DDP + u_V_DDD))
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
              rVd_t = -(M*N*(3*(-1 + Msq + Nsq)*u_Vd_DDS +  &
        (3 - 4*Msq - 4*Nsq)*u_Vd_DDP + (Msq + Nsq)*u_Vd_DDD))
              case (5)
              rVd(1) = -(M*(-3 + 4*Msq + 3*Nsq)* &
       (3*u_V_DDS - 4*u_V_DDP + u_V_DDD))/2.
              rVd(2) = -(L*(-1 + 4*Msq + Nsq)* &
       (3*u_V_DDS - 4*u_V_DDP + u_V_DDD))/2.
              rVd(3) = 0
              rVd_t = (L*(L - M)*M*(L + M)*(3*u_Vd_DDS - 4*u_Vd_DDP + u_Vd_DDD))/2.
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
        2*(-3 + 4*Msq + 6*Nsq)*u_V_DDP + (3 - 2*Msq - 3*Nsq)*u_V_DDD &
 ))/2.
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
              rVd_t = ((1 - 3*Nsq)**2*u_Vd_DDS +  &
      3*(-1 + Nsq)*(-4*Nsq*u_Vd_DDP + (-1 + Nsq)*u_Vd_DDD))/4.
              case (4)
              rVd(1) = (root_3*N*((-1 + 3*Nsq)*u_V_DDS +  &
        (2 - 4*Nsq)*u_V_DDP + (-1 + Nsq)*u_V_DDD))/2.
              rVd(2) = 0
              rVd(3) = (root_3*L*((-1 + 9*Nsq)*u_V_DDS +  &
        (2 - 12*Nsq)*u_V_DDP + (-1 + 3*Nsq)*u_V_DDD))/2.
              rVd_t = (root_3*L*N*((-1 + 3*Nsq)*u_Vd_DDS + (2 - 4*Nsq)*u_Vd_DDP +  &
        (-1 + Nsq)*u_Vd_DDD))/2.
              case (5)
              rVd(1) = (root_3*L*((-1 + 3*Nsq)*u_V_DDS + u_V_DDD +  &
        Nsq*(-4*u_V_DDP + u_V_DDD)))/2.
              rVd(2) = -(root_3*M*((-1 + 3*Nsq)*u_V_DDS + u_V_DDD +  &
         Nsq*(-4*u_V_DDP + u_V_DDD)))/2.
              rVd(3) = -(root_3*N*(-1 + 2*Msq + Nsq)* &
       (3*u_V_DDS - 4*u_V_DDP + u_V_DDD))/2.
              rVd_t = (root_3*(L - M)*(L + M)* &
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
              rVd_t = -(M*N*(3*(-1 + Msq + Nsq)*u_Vd_DDS +  &
        (3 - 4*Msq - 4*Nsq)*u_Vd_DDP + (Msq + Nsq)*u_Vd_DDD))
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
              rVd(2) = 2*M*(-u_V_DDP + u_V_DDD -  &
      Nsq*(3*u_V_DDS - 4*u_V_DDP + u_V_DDD))
              rVd(3) = -2*N*(-1 + Msq + 2*Nsq)* &
    (3*u_V_DDS - 4*u_V_DDP + u_V_DDD)
              rVd_t = u_Vd_DDP - (-1 + Msq)*Nsq* &
     (3*u_Vd_DDS - 4*u_Vd_DDP + u_Vd_DDD) -  &
    N**4*(3*u_Vd_DDS - 4*u_Vd_DDP + u_Vd_DDD) +  &
    Msq*(-u_Vd_DDP + u_Vd_DDD)
              case (5)
              rVd(1) = -(N*(3*(-1 + 2*Msq + Nsq)*u_V_DDS +  &
         (2 - 8*Msq - 4*Nsq)*u_V_DDP + (1 + 2*Msq + Nsq)*u_V_DDD))/2.
              rVd(2) = -2*L*M*N*(3*u_V_DDS - 4*u_V_DDP + u_V_DDD)
              rVd(3) = (L*((3 - 6*Msq - 9*Nsq)*u_V_DDS +  &
        2*(-1 + 4*Msq + 6*Nsq)*u_V_DDP - (1 + 2*Msq + 3*Nsq)*u_V_DDD &
 ))/2.
              rVd_t = -(L*N*(3*(-1 + 2*Msq + Nsq)*u_Vd_DDS +  &
         (2 - 8*Msq - 4*Nsq)*u_Vd_DDP + (1 + 2*Msq + Nsq)*u_Vd_DDD))/2.
            end select
          case (5)
            select case (orb_dir_j)
              case (1)
              rVd(1) = -(M*(-3 + 4*Msq + 3*Nsq)* &
       (3*u_V_DDS - 4*u_V_DDP + u_V_DDD))/2.
              rVd(2) = -(L*(-1 + 4*Msq + Nsq)* &
       (3*u_V_DDS - 4*u_V_DDP + u_V_DDD))/2.
              rVd(3) = 0
              rVd_t = (L*(L - M)*M*(L + M)*(3*u_Vd_DDS - 4*u_Vd_DDP + u_Vd_DDD))/2.
              case (2)
              rVd(1) = 0
              rVd(2) = (N*(-3*(-1 + 6*Msq + Nsq)*u_V_DDS +  &
        (-6 + 24*Msq + 4*Nsq)*u_V_DDP - (-3 + 6*Msq + Nsq)*u_V_DDD))/ &
    2.
              rVd(3) = (M*((3 - 6*Msq - 9*Nsq)*u_V_DDS +  &
        2*(-3 + 4*Msq + 6*Nsq)*u_V_DDP + (3 - 2*Msq - 3*Nsq)*u_V_DDD &
 ))/2.
              rVd_t = -(M*N*(3*(-1 + 2*Msq + Nsq)*u_Vd_DDS +  &
         (6 - 8*Msq - 4*Nsq)*u_Vd_DDP + (-3 + 2*Msq + Nsq)*u_Vd_DDD))/ &
    2.
              case (3)
              rVd(1) = (root_3*L*((-1 + 3*Nsq)*u_V_DDS + u_V_DDD +  &
        Nsq*(-4*u_V_DDP + u_V_DDD)))/2.
              rVd(2) = -(root_3*M*((-1 + 3*Nsq)*u_V_DDS + u_V_DDD +  &
         Nsq*(-4*u_V_DDP + u_V_DDD)))/2.
              rVd(3) = -(root_3*N*(-1 + 2*Msq + Nsq)* &
       (3*u_V_DDS - 4*u_V_DDP + u_V_DDD))/2.
              rVd_t = (root_3*(L - M)*(L + M)* &
      ((-1 + 3*Nsq)*u_Vd_DDS + u_Vd_DDD +  &
        Nsq*(-4*u_Vd_DDP + u_Vd_DDD)))/4.
              case (4)
              rVd(1) = -(N*(3*(-1 + 2*Msq + Nsq)*u_V_DDS +  &
         (2 - 8*Msq - 4*Nsq)*u_V_DDP + (1 + 2*Msq + Nsq)*u_V_DDD))/2.
              rVd(2) = -2*L*M*N*(3*u_V_DDS - 4*u_V_DDP + u_V_DDD)
              rVd(3) = (L*((3 - 6*Msq - 9*Nsq)*u_V_DDS +  &
        2*(-1 + 4*Msq + 6*Nsq)*u_V_DDP - (1 + 2*Msq + 3*Nsq)*u_V_DDD &
 ))/2.
              rVd_t = -(L*N*(3*(-1 + 2*Msq + Nsq)*u_Vd_DDS +  &
         (2 - 8*Msq - 4*Nsq)*u_Vd_DDP + (1 + 2*Msq + Nsq)*u_Vd_DDD))/2.
              case (5)
              rVd(1) = 3*L*(Lsq - Msq)*u_V_DDS
              rVd(2) = M*(-1 + 2*Msq + Nsq)* &
    (3*u_V_DDS - 8*u_V_DDP + 2*u_V_DDD)
              rVd(3) = N*((2 - 8*Msq - 4*Nsq)*u_V_DDP +  &
      (1 + 2*Msq + Nsq)*u_V_DDD)
              rVd_t = (3*(Lsq - Msq)**2*u_Vd_DDS -  &
      4*(4*M**4 + 4*Msq*(-1 + Nsq) + Nsq*(-1 + Nsq))*u_Vd_DDP +  &
      (4*M**4 + 4*Msq*(-1 + Nsq) + (1 + Nsq)**2)*u_Vd_DDD)/4.
            end select
        end select
      case (ORB_F)
        select case (orb_dir_i)
          case (1)
            select case (orb_dir_j)
              case (1)
              rVd(1) = (root_3*(6*sqrt(10.)*Lsq*Msq*u_V_DFS -  &
        sqrt(10.)*Msq*(-3*Lsq + Msq)*u_V_DFS +  &
        sqrt(5.)*(1 + 8*M**4 - Nsq + 6*Msq*(-1 + Nsq))*u_V_DFP -  &
        sqrt(2.)*(4*M**4 - 2*Nsq + 3*Msq*(-1 + Nsq))*u_V_DFD))/4.
              rVd(2) = (L*M*(sqrt(30.)*(3 - 5*Msq - 3*Nsq)*u_V_DFS +  &
        (-3 + 8*Msq + 3*Nsq)*(2*sqrt(15.)*u_V_DFP - sqrt(6.)*u_V_DFD)))/ &
    2.
              rVd(3) = (L*N*(sqrt(15.)*(-1 + 6*Msq)*u_V_DFP +  &
        sqrt(6.)*(2 - 3*Msq)*u_V_DFD))/2.
              rVd_t = (root_3*L*(-(sqrt(10.)*Msq*(-3*Lsq + Msq)*u_Vd_DFS) +  &
        sqrt(5.)*(1 + 8*M**4 - Nsq + 6*Msq*(-1 + Nsq))*u_Vd_DFP -  &
        sqrt(2.)*(4*M**4 - 2*Nsq + 3*Msq*(-1 + Nsq))*u_Vd_DFD))/4.
              case (2)
              rVd(1) = 0
              rVd(2) = -6*M*N*(-1 + 2*Msq + Nsq)* &
    (sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD)
              rVd(3) = Sqrt(2.5)*(1 - 3*Nsq)*u_V_DFP + (-1 + 6*Nsq)*u_V_DFD -  &
    3*M**4*(sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD) -  &
    3*Msq*(-1 + 3*Nsq)*(sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP +  &
       u_V_DFD)
              rVd_t = (N*(-(sqrt(10.)*(-1 + Nsq)*u_Vd_DFP) +  &
        2*(-1 + 2*Nsq)*u_Vd_DFD -  &
        6*M**4*(sqrt(5.)*u_Vd_DFS - sqrt(10.)*u_Vd_DFP + u_Vd_DFD) -  &
        6*Msq*(-1 + Nsq)*(sqrt(5.)*u_Vd_DFS - sqrt(10.)*u_Vd_DFP +  &
           u_Vd_DFD)))/2.
              case (3)
              rVd(1) = ((-1 + 5*Nsq)*u_V_DFP - 2*sqrt(10.)*Nsq*u_V_DFD +  &
      Msq*(3*sqrt(2.)*(-1 + 5*Nsq)*u_V_DFS + (2 - 30*Nsq)*u_V_DFP +  &
         sqrt(10.)*(1 + 3*Nsq)*u_V_DFD))/4.
              rVd(2) = (L*M*(3*sqrt(2.)*(-1 + 5*Nsq)*u_V_DFS +  &
        (2 - 30*Nsq)*u_V_DFP + sqrt(10.)*(1 + 3*Nsq)*u_V_DFD))/2.
              rVd(3) = (L*N*(5*u_V_DFP - 2*sqrt(10.)*u_V_DFD +  &
        3*Msq*(5*sqrt(2.)*u_V_DFS - 10*u_V_DFP + sqrt(10.)*u_V_DFD)))/ &
    2.
              rVd_t = (L*((-1 + 5*Nsq)*u_Vd_DFP - 2*sqrt(10.)*Nsq*u_Vd_DFD +  &
        Msq*(3*sqrt(2.)*(-1 + 5*Nsq)*u_Vd_DFS + (2 - 30*Nsq)*u_Vd_DFP +  &
           sqrt(10.)*(1 + 3*Nsq)*u_Vd_DFD)))/4.
              case (4)
              rVd(1) = (root_3*M*N*((-3 + 5*Nsq)*u_V_DFS +  &
        sqrt(2.)*(1 - 5*Nsq)*u_V_DFP + sqrt(5.)*(1 + Nsq)*u_V_DFD))/2.
              rVd(2) = (root_3*L*N*((-3 + 5*Nsq)*u_V_DFS +  &
        sqrt(2.)*(1 - 5*Nsq)*u_V_DFP + sqrt(5.)*(1 + Nsq)*u_V_DFD))/2.
              rVd(3) = (root_3*L*M*(3*(-1 + 5*Nsq)*u_V_DFS +  &
        sqrt(2.)*(1 - 15*Nsq)*u_V_DFP + sqrt(5.)*(1 + 3*Nsq)*u_V_DFD))/2.
              rVd_t = (root_3*L*M*N*((-3 + 5*Nsq)*u_Vd_DFS +  &
        sqrt(2.)*(1 - 5*Nsq)*u_Vd_DFP + sqrt(5.)*(1 + Nsq)*u_Vd_DFD))/2.
              case (5)
              rVd(1) = (L*M*(3*sqrt(2.)*(1 - 6*Nsq + 5*N**4)*u_V_DFS +  &
        (-1 + 27*Nsq - 30*N**4)*u_V_DFP +  &
        sqrt(10.)*(-1 + 3*N**4)*u_V_DFD))/(2.*(-1 + Nsq))
              rVd(2) = (3*Msq*((1 - 5*Nsq)*u_V_DFP +  &
         2*sqrt(10.)*Nsq*u_V_DFD) +  &
      Lsq*(3*sqrt(2.)*(1 - 6*Nsq + 5*N**4)*u_V_DFS +  &
         (-1 + 27*Nsq - 30*N**4)*u_V_DFP +  &
         sqrt(10.)*(-1 + 3*N**4)*u_V_DFD))/(4.*(-1 + Nsq))
              rVd(3) = (M*N*(-15*sqrt(2.)*(-1 + Nsq)*(-1 + Msq + Nsq)*u_V_DFS +  &
        (26 - 60*Nsq + 30*N**4 + 30*Msq*(-1 + Nsq))*u_V_DFP +  &
        sqrt(10.)*(-1 + 6*Nsq - 3*N**4 - 3*Msq*(-1 + Nsq))*u_V_DFD))/ &
    (2.*(-1 + Nsq))
              rVd_t = (M*(Msq*((1 - 5*Nsq)*u_Vd_DFP + 2*sqrt(10.)*Nsq*u_Vd_DFD) +  &
        Lsq*(3*sqrt(2.)*(1 - 6*Nsq + 5*N**4)*u_Vd_DFS +  &
           (-1 + 27*Nsq - 30*N**4)*u_Vd_DFP +  &
           sqrt(10.)*(-1 + 3*N**4)*u_Vd_DFD)))/(4.*(-1 + Nsq))
              case (6)
              rVd(1) = (-3*M*(-3*Lsq + Msq)*N* &
      (sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD))/2.
              rVd(2) = (3*L*(Lsq - 3*Msq)*N* &
      (sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD))/2.
              rVd(3) = (3*L*(L - M)*M*(L + M)* &
      (sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD))/2.
              rVd_t = (3*L*(L - M)*M*(L + M)*N* &
      (sqrt(5.)*u_Vd_DFS - sqrt(10.)*u_Vd_DFP + u_Vd_DFD))/2.
              case (7)
              rVd(1) = Sqrt(7.5)*L*M*(2*Lsq - 3*Msq)*u_V_DFS
              rVd(2) = (root_3*(sqrt(10.)*L**4*u_V_DFS -  &
        9*sqrt(10.)*Lsq*Msq*u_V_DFS -  &
        sqrt(5.)*(3 + 40*M**4 - 5*Nsq + 2*N**4 + 30*Msq*(-1 + Nsq))* &
         u_V_DFP + sqrt(2.)*(1 + 20*M**4 - 4*Nsq + N**4 +  &
           15*Msq*(-1 + Nsq))*u_V_DFD))/4.
              rVd(3) = (M*N*(sqrt(15.)*(5 - 10*Msq - 4*Nsq)*u_V_DFP +  &
        sqrt(6.)*(-4 + 5*Msq + 2*Nsq)*u_V_DFD))/2.
              rVd_t = (root_3*M*(sqrt(10.)*Lsq*(Lsq - 3*Msq)*u_Vd_DFS -  &
        sqrt(5.)*(3 + 8*M**4 - 5*Nsq + 2*N**4 + 10*Msq*(-1 + Nsq))* &
         u_Vd_DFP + sqrt(2.)*(1 + 4*M**4 - 4*Nsq + N**4 +  &
           5*Msq*(-1 + Nsq))*u_Vd_DFD))/4.
            end select
          case (2)
            select case (orb_dir_j)
              case (1)
              rVd(1) = 0
              rVd(2) = (M*N*(sqrt(30.)*(3 - 8*Msq - 3*Nsq)*u_V_DFS +  &
        2*sqrt(15.)*(-4 + 8*Msq + 3*Nsq)*u_V_DFP +  &
        sqrt(6.)*(7 - 8*Msq - 3*Nsq)*u_V_DFD))/2.
              rVd(3) = (-((-1 + 3*Nsq)* &
         (sqrt(15.)*u_V_DFP - 2*sqrt(6.)*u_V_DFD)) -  &
      4*M**4*(sqrt(30.)*u_V_DFS - 2*sqrt(15.)*u_V_DFP +  &
         sqrt(6.)*u_V_DFD) + Msq* &
       (3*sqrt(30.)*(1 - 3*Nsq)*u_V_DFS +  &
         2*sqrt(15.)*(-4 + 9*Nsq)*u_V_DFP + sqrt(6.)*(7 - 9*Nsq)*u_V_DFD &
 ))/4.
              rVd_t = (N*(-((-1 + Nsq)*(sqrt(15.)*u_Vd_DFP - 2*sqrt(6.)*u_Vd_DFD)) -  &
        4*M**4*(sqrt(30.)*u_Vd_DFS - 2*sqrt(15.)*u_Vd_DFP +  &
           sqrt(6.)*u_Vd_DFD) +  &
        Msq*(-3*sqrt(30.)*(-1 + Nsq)*u_Vd_DFS +  &
           2*sqrt(15.)*(-4 + 3*Nsq)*u_Vd_DFP +  &
           sqrt(6.)*(7 - 3*Nsq)*u_Vd_DFD)))/4.
              case (2)
              rVd(1) = (Nsq*(sqrt(10.)*u_V_DFP - 4*u_V_DFD) + 2*u_V_DFD +  &
      Msq*(sqrt(10.)*u_V_DFP - 4*u_V_DFD +  &
         6*Nsq*(sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD)))/2.
              rVd(2) = L*M*(sqrt(10.)*u_V_DFP - 4*u_V_DFD +  &
      6*Nsq*(sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD))
              rVd(3) = L*N*(sqrt(10.)*u_V_DFP - 4*u_V_DFD +  &
      6*Msq*(sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD))
              rVd_t = (L*(Nsq*(sqrt(10.)*u_Vd_DFP - 4*u_Vd_DFD) + 2*u_Vd_DFD +  &
        Msq*(sqrt(10.)*u_Vd_DFP - 4*u_Vd_DFD +  &
           6*Nsq*(sqrt(5.)*u_Vd_DFS - sqrt(10.)*u_Vd_DFP + u_Vd_DFD))))/2.
              case (3)
              rVd(1) = 0
              rVd(2) = (3*M*N*(sqrt(2.)*(-1 + 5*Nsq)*u_V_DFS +  &
        (4 - 10*Nsq)*u_V_DFP + sqrt(10.)*(-1 + Nsq)*u_V_DFD))/2.
              rVd(3) = ((-1 + 15*Nsq)*u_V_DFP +  &
      2*sqrt(10.)*(1 - 3*Nsq)*u_V_DFD +  &
      Msq*(3*sqrt(2.)*(-1 + 15*Nsq)*u_V_DFS + (12 - 90*Nsq)*u_V_DFP +  &
         3*sqrt(10.)*(-1 + 3*Nsq)*u_V_DFD))/4.
              rVd_t = (N*((-1 + 5*Nsq)*u_Vd_DFP - 2*sqrt(10.)*(-1 + Nsq)*u_Vd_DFD +  &
        3*Msq*(sqrt(2.)*(-1 + 5*Nsq)*u_Vd_DFS + (4 - 10*Nsq)*u_Vd_DFP +  &
           sqrt(10.)*(-1 + Nsq)*u_Vd_DFD)))/4.
              case (4)
              rVd(1) = 0
              rVd(2) = (root_3*(sqrt(2.)*(-1 + 7*Nsq - 10*N**4)*u_V_DFP +  &
        2*Nsq*((-3 + 5*Nsq)*u_V_DFS + sqrt(5.)*(-1 + Nsq)*u_V_DFD)))/4.
              rVd(3) = (root_3*M*N*((-6 + 20*Nsq)*u_V_DFS +  &
        sqrt(2.)*(7 - 20*Nsq)*u_V_DFP + 2*sqrt(5.)*(-1 + 2*Nsq)*u_V_DFD) &
 )/2.
              rVd_t = (root_3*M*(sqrt(2.)*(-1 + 7*Nsq - 10*N**4)*u_Vd_DFP +  &
        2*Nsq*((-3 + 5*Nsq)*u_Vd_DFS + sqrt(5.)*(-1 + Nsq)*u_Vd_DFD)))/4.
              case (5)
              rVd(1) = (3*M*N*(sqrt(2.)*(-1 + 5*Nsq)*u_V_DFS +  &
        (4 - 10*Nsq)*u_V_DFP + sqrt(10.)*(-1 + Nsq)*u_V_DFD))/4.
              rVd(2) = (3*L*N*(sqrt(2.)*(-1 + 5*Nsq)*u_V_DFS +  &
        (4 - 10*Nsq)*u_V_DFP + sqrt(10.)*(-1 + Nsq)*u_V_DFD))/4.
              rVd(3) = (L*M*(3*sqrt(2.)*(-1 + 15*Nsq)*u_V_DFS +  &
        (12 - 90*Nsq)*u_V_DFP + 3*sqrt(10.)*(-1 + 3*Nsq)*u_V_DFD))/4.
              rVd_t = (3*L*M*N*(sqrt(2.)*(-1 + 5*Nsq)*u_Vd_DFS +  &
        (4 - 10*Nsq)*u_Vd_DFP + sqrt(10.)*(-1 + Nsq)*u_Vd_DFD))/4.
              case (6)
              rVd(1) = 0
              rVd(2) = (sqrt(10.)*(1 - 6*Msq)*u_V_DFP +  &
      8*(-1 + 3*Msq)*u_V_DFD -  &
      6*N**4*(sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD) +  &
      Nsq*(6*sqrt(5.)*(1 - 6*Msq)*u_V_DFS +  &
         9*sqrt(10.)*(-1 + 4*Msq)*u_V_DFP + 18*(1 - 2*Msq)*u_V_DFD))/4.
              rVd(3) = (3*M*N*(-2*sqrt(5.)*(-1 + 2*Msq + 2*Nsq)*u_V_DFS +  &
        sqrt(10.)*(-3 + 4*Msq + 4*Nsq)*u_V_DFP -  &
        2*(-3 + 2*Msq + 2*Nsq)*u_V_DFD))/2.
              rVd_t = (M*(sqrt(10.)*(1 - 2*Msq)*u_Vd_DFP + 8*(-1 + Msq)*u_Vd_DFD -  &
        6*N**4*(sqrt(5.)*u_Vd_DFS - sqrt(10.)*u_Vd_DFP + u_Vd_DFD) -  &
        3*Nsq*(2*sqrt(5.)*(-1 + 2*Msq)*u_Vd_DFS +  &
           sqrt(10.)*(3 - 4*Msq)*u_Vd_DFP + 2*(-3 + 2*Msq)*u_Vd_DFD)))/4.
              case (7)
              rVd(1) = (M*N*(-(sqrt(30.)*(-1 + 4*Msq + Nsq)*u_V_DFS) +  &
        2*sqrt(15.)*(-2 + 4*Msq + Nsq)*u_V_DFP -  &
        sqrt(6.)*(-5 + 4*Msq + Nsq)*u_V_DFD))/4.
              rVd(2) = (L*N*(-(sqrt(30.)*(-1 + 12*Msq + Nsq)*u_V_DFS) +  &
        2*sqrt(15.)*(-2 + 12*Msq + Nsq)*u_V_DFP -  &
        sqrt(6.)*(-5 + 12*Msq + Nsq)*u_V_DFD))/4.
              rVd(3) = (L*M*(sqrt(30.)*(1 - 4*Msq - 3*Nsq)*u_V_DFS +  &
        2*sqrt(15.)*(-2 + 4*Msq + 3*Nsq)*u_V_DFP +  &
        sqrt(6.)*(5 - 4*Msq - 3*Nsq)*u_V_DFD))/4.
              rVd_t = (L*M*N*(-(sqrt(30.)*(-1 + 4*Msq + Nsq)*u_Vd_DFS) +  &
        2*sqrt(15.)*(-2 + 4*Msq + Nsq)*u_Vd_DFP -  &
        sqrt(6.)*(-5 + 4*Msq + Nsq)*u_Vd_DFD))/4.
            end select
          case (3)
            select case (orb_dir_j)
              case (1)
              rVd(1) = (3*L*M*(sqrt(10.)*(-1 + 3*Nsq)*u_V_DFS -  &
        6*sqrt(5.)*Nsq*u_V_DFP + 3*sqrt(2.)*(1 + Nsq)*u_V_DFD))/4.
              rVd(2) = (3*(Lsq - Msq)*(sqrt(10.)*(-1 + 3*Nsq)*u_V_DFS -  &
        6*sqrt(5.)*Nsq*u_V_DFP + 3*sqrt(2.)*(1 + Nsq)*u_V_DFD))/8.
              rVd(3) = (-3*M*N*(-3 + 4*Msq + 3*Nsq)* &
      (sqrt(10.)*u_V_DFS - 2*sqrt(5.)*u_V_DFP + sqrt(2.)*u_V_DFD))/4.
              rVd_t = ((3*Lsq*M - M**3)*(sqrt(10.)*(-1 + 3*Nsq)*u_Vd_DFS -  &
        6*sqrt(5.)*Nsq*u_Vd_DFP + 3*sqrt(2.)*(1 + Nsq)*u_Vd_DFD))/8.
              case (2)
              rVd(1) = (root_3*M*N*(-1 + 3*Nsq)* &
      (sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD))/2.
              rVd(2) = (root_3*L*N*(-1 + 3*Nsq)* &
      (sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD))/2.
              rVd(3) = (root_3*L*M*(-1 + 9*Nsq)* &
      (sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD))/2.
              rVd_t = (root_3*L*M*N*(-1 + 3*Nsq)* &
      (sqrt(5.)*u_Vd_DFS - sqrt(10.)*u_Vd_DFP + u_Vd_DFD))/2.
              case (3)
              rVd(1) = 0
              rVd(2) = (root_3*(sqrt(2.)*(1 - 8*Nsq + 15*N**4)*u_V_DFS +  &
        2*Nsq*(11 - 15*Nsq)*u_V_DFP +  &
        sqrt(10.)*(1 - 4*Nsq + 3*N**4)*u_V_DFD))/8.
              rVd(3) = (root_3*M*N*(sqrt(2.)*(-4 + 15*Nsq)*u_V_DFS +  &
        (11 - 30*Nsq)*u_V_DFP + sqrt(10.)*(-2 + 3*Nsq)*u_V_DFD))/2.
              rVd_t = (root_3*M*(sqrt(2.)*(1 - 8*Nsq + 15*N**4)*u_Vd_DFS +  &
        2*Nsq*(11 - 15*Nsq)*u_Vd_DFP +  &
        sqrt(10.)*(1 - 4*Nsq + 3*N**4)*u_Vd_DFD))/8.
              case (4)
              rVd(1) = 0
              rVd(2) = 0
              rVd(3) = (3*((1 - 14*Nsq + 25*N**4)*u_V_DFS -  &
        sqrt(2.)*(1 - 18*Nsq + 25*N**4)*u_V_DFP +  &
        sqrt(5.)*(1 - 6*Nsq + 5*N**4)*u_V_DFD))/4.
              rVd_t = (N*((3 - 14*Nsq + 15*N**4)*u_Vd_DFS -  &
        3*sqrt(2.)*(1 - 6*Nsq + 5*N**4)*u_Vd_DFP +  &
        3*sqrt(5.)*(-1 + Nsq)**2*u_Vd_DFD))/4.
              case (5)
              rVd(1) = (root_3*(sqrt(2.)*(1 - 8*Nsq + 15*N**4)*u_V_DFS +  &
        2*Nsq*(11 - 15*Nsq)*u_V_DFP +  &
        sqrt(10.)*(1 - 4*Nsq + 3*N**4)*u_V_DFD))/8.
              rVd(2) = 0
              rVd(3) = (root_3*L*N*(sqrt(2.)*(-4 + 15*Nsq)*u_V_DFS +  &
        (11 - 30*Nsq)*u_V_DFP + sqrt(10.)*(-2 + 3*Nsq)*u_V_DFD))/2.
              rVd_t = (root_3*L*(sqrt(2.)*(1 - 8*Nsq + 15*N**4)*u_Vd_DFS +  &
        2*Nsq*(11 - 15*Nsq)*u_Vd_DFP +  &
        sqrt(10.)*(1 - 4*Nsq + 3*N**4)*u_Vd_DFD))/8.
              case (6)
              rVd(1) = (root_3*L*N*(-1 + 3*Nsq)* &
      (sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD))/2.
              rVd(2) = -(root_3*M*N*(-1 + 3*Nsq)* &
       (sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD))/2.
              rVd(3) = (root_3*(L - M)*(L + M)*(-1 + 9*Nsq)* &
      (sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD))/4.
              rVd_t = (root_3*(L - M)*(L + M)*N*(-1 + 3*Nsq)* &
      (sqrt(5.)*u_Vd_DFS - sqrt(10.)*u_Vd_DFP + u_Vd_DFD))/4.
              case (7)
              rVd(1) = (3*(Lsq - Msq)*(sqrt(10.)*(-1 + 3*Nsq)*u_V_DFS -  &
        6*sqrt(5.)*Nsq*u_V_DFP + 3*sqrt(2.)*(1 + Nsq)*u_V_DFD))/8.
              rVd(2) = (-3*L*M*(sqrt(10.)*(-1 + 3*Nsq)*u_V_DFS -  &
        6*sqrt(5.)*Nsq*u_V_DFP + 3*sqrt(2.)*(1 + Nsq)*u_V_DFD))/4.
              rVd(3) = (-3*L*N*(-1 + 4*Msq + Nsq)* &
      (sqrt(10.)*u_V_DFS - 2*sqrt(5.)*u_V_DFP + sqrt(2.)*u_V_DFD))/4.
              rVd_t = ((L**3 - 3*L*Msq)*(sqrt(10.)*(-1 + 3*Nsq)*u_Vd_DFS -  &
        6*sqrt(5.)*Nsq*u_Vd_DFP + 3*sqrt(2.)*(1 + Nsq)*u_Vd_DFD))/8.
            end select
          case (4)
            select case (orb_dir_j)
              case (1)
              rVd(1) = (M*N*(sqrt(30.)*(3 - 4*Msq - 3*Nsq)*u_V_DFS +  &
        2*sqrt(15.)*(-2 + 4*Msq + 3*Nsq)*u_V_DFP -  &
        sqrt(6.)*(1 + 4*Msq + 3*Nsq)*u_V_DFD))/4.
              rVd(2) = (L*N*(-3*sqrt(30.)*(-1 + 4*Msq + Nsq)*u_V_DFS +  &
        2*sqrt(15.)*(-2 + 12*Msq + 3*Nsq)*u_V_DFP -  &
        sqrt(6.)*(1 + 12*Msq + 3*Nsq)*u_V_DFD))/4.
              rVd(3) = (L*M*(sqrt(30.)*(3 - 4*Msq - 9*Nsq)*u_V_DFS +  &
        2*sqrt(15.)*(-2 + 4*Msq + 9*Nsq)*u_V_DFP -  &
        sqrt(6.)*(1 + 4*Msq + 9*Nsq)*u_V_DFD))/4.
              rVd_t = (L*M*N*(sqrt(30.)*(3 - 4*Msq - 3*Nsq)*u_Vd_DFS +  &
        2*sqrt(15.)*(-2 + 4*Msq + 3*Nsq)*u_Vd_DFP -  &
        sqrt(6.)*(1 + 4*Msq + 3*Nsq)*u_Vd_DFD))/4.
              case (2)
              rVd(1) = (L*M*(sqrt(10.)*(-1 + 6*Nsq - 6*N**4)*u_V_DFP +  &
        2*u_V_DFD + 6*Nsq*(-1 + Nsq)*(sqrt(5.)*u_V_DFS + u_V_DFD)))/ &
    (-1 + Nsq)
              rVd(2) = (Lsq*(sqrt(10.)*(-1 + 6*Nsq - 6*N**4)*u_V_DFP +  &
         2*u_V_DFD + 6*Nsq*(-1 + Nsq)*(sqrt(5.)*u_V_DFS + u_V_DFD))  &
 + Msq*(-6*u_V_DFD + Nsq*(-3*sqrt(10.)*u_V_DFP + 12*u_V_DFD)))/ &
    (2.*(-1 + Nsq))
              rVd(3) = (M*N*(-6*sqrt(5.)*(-1 + Nsq)*(-1 + Msq + Nsq)*u_V_DFS +  &
        sqrt(10.)*(5 - 12*Nsq + 6*N**4 + 6*Msq*(-1 + Nsq))*u_V_DFP -  &
        2*(2 - 6*Nsq + 3*N**4 + 3*Msq*(-1 + Nsq))*u_V_DFD))/(-1 + Nsq)
              rVd_t = (Lsq*M*(sqrt(10.)*(-1 + 6*Nsq - 6*N**4)*u_Vd_DFP +  &
         2*u_Vd_DFD + 6*Nsq*(-1 + Nsq)*(sqrt(5.)*u_Vd_DFS + u_Vd_DFD))  &
 + M**3*(-2*u_Vd_DFD + Nsq*(-(sqrt(10.)*u_Vd_DFP) + 4*u_Vd_DFD)))/ &
    (2.*(-1 + Nsq))
              case (3)
              rVd(1) = (3*M*N*(sqrt(2.)*(-1 + 5*Nsq)*u_V_DFS +  &
        (4 - 10*Nsq)*u_V_DFP + sqrt(10.)*(-1 + Nsq)*u_V_DFD))/4.
              rVd(2) = (3*L*N*(sqrt(2.)*(-1 + 5*Nsq)*u_V_DFS +  &
        (4 - 10*Nsq)*u_V_DFP + sqrt(10.)*(-1 + Nsq)*u_V_DFD))/4.
              rVd(3) = (L*M*(3*sqrt(2.)*(-1 + 15*Nsq)*u_V_DFS +  &
        (12 - 90*Nsq)*u_V_DFP + 3*sqrt(10.)*(-1 + 3*Nsq)*u_V_DFD))/4.
              rVd_t = (3*L*M*N*(sqrt(2.)*(-1 + 5*Nsq)*u_Vd_DFS +  &
        (4 - 10*Nsq)*u_Vd_DFP + sqrt(10.)*(-1 + Nsq)*u_Vd_DFD))/4.
              case (4)
              rVd(1) = (root_3*(sqrt(2.)*(-1 + 7*Nsq - 10*N**4)*u_V_DFP +  &
        2*Nsq*((-3 + 5*Nsq)*u_V_DFS + sqrt(5.)*(-1 + Nsq)*u_V_DFD)))/4.
              rVd(2) = 0
              rVd(3) = (root_3*L*N*((-6 + 20*Nsq)*u_V_DFS +  &
        sqrt(2.)*(7 - 20*Nsq)*u_V_DFP + 2*sqrt(5.)*(-1 + 2*Nsq)*u_V_DFD) &
 )/2.
              rVd_t = (root_3*L*(sqrt(2.)*(-1 + 7*Nsq - 10*N**4)*u_Vd_DFP +  &
        2*Nsq*((-3 + 5*Nsq)*u_Vd_DFS + sqrt(5.)*(-1 + Nsq)*u_Vd_DFD)))/4.
              case (5)
              rVd(1) = 0
              rVd(2) = (3*M*N*(sqrt(2.)*(1 - 5*Nsq)*u_V_DFS +  &
        2*(-2 + 5*Nsq)*u_V_DFP - sqrt(10.)*(-1 + Nsq)*u_V_DFD))/2.
              rVd(3) = (-3*sqrt(2.)*(1 - 18*Nsq + 25*N**4 + Msq*(-1 + 15*Nsq))* &
       u_V_DFS + (11 - 111*Nsq + 150*N**4 + 6*Msq*(-2 + 15*Nsq))* &
       u_V_DFP + sqrt(10.)*(-1 + 12*Nsq - 15*N**4 + Msq*(3 - 9*Nsq))* &
       u_V_DFD)/4.
              rVd_t = -(N*(3*sqrt(2.)*(-1 + Msq + Nsq)*(-1 + 5*Nsq)*u_Vd_DFS +  &
         (-11 + 37*Nsq - 30*N**4 - 6*Msq*(-2 + 5*Nsq))*u_Vd_DFP +  &
         sqrt(10.)*(-1 + Nsq)*(-1 + 3*Msq + 3*Nsq)*u_Vd_DFD))/4.
              case (6)
              rVd(1) = (sqrt(10.)*(1 - 2*Msq)*u_V_DFP + 8*Msq*u_V_DFD -  &
      6*N**4*(sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD) +  &
      Nsq*(6*sqrt(5.)*(1 - 2*Msq)*u_V_DFS +  &
         sqrt(10.)*(-5 + 12*Msq)*u_V_DFP + 2*(1 - 6*Msq)*u_V_DFD))/4.
              rVd(2) = L*M*(-(sqrt(10.)*u_V_DFP) + 4*u_V_DFD -  &
      6*Nsq*(sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD))
              rVd(3) = (L*N*(-6*sqrt(5.)*(-1 + 2*Msq + 2*Nsq)*u_V_DFS +  &
        sqrt(10.)*(-5 + 12*Msq + 12*Nsq)*u_V_DFP -  &
        2*(-1 + 6*Msq + 6*Nsq)*u_V_DFD))/2.
              rVd_t = (L*(sqrt(10.)*(1 - 2*Msq)*u_Vd_DFP + 8*Msq*u_Vd_DFD -  &
        6*N**4*(sqrt(5.)*u_Vd_DFS - sqrt(10.)*u_Vd_DFP + u_Vd_DFD) +  &
        Nsq*(6*sqrt(5.)*(1 - 2*Msq)*u_Vd_DFS +  &
           sqrt(10.)*(-5 + 12*Msq)*u_Vd_DFP + 2*(1 - 6*Msq)*u_Vd_DFD)))/4.
              case (7)
              rVd(1) = 0
              rVd(2) = (M*N*(sqrt(30.)*(-5 + 8*Msq + 5*Nsq)*u_V_DFS -  &
        2*sqrt(15.)*(-4 + 8*Msq + 5*Nsq)*u_V_DFP +  &
        sqrt(6.)*(-1 + 8*Msq + 5*Nsq)*u_V_DFD))/2.
              rVd(3) = (sqrt(30.)*(1 + 4*M**4 - 6*Nsq + 5*N**4 +  &
         5*Msq*(-1 + 3*Nsq))*u_V_DFS -  &
      sqrt(15.)*(1 + 8*M**4 - 9*Nsq + 10*N**4 + Msq*(-8 + 30*Nsq))* &
       u_V_DFP + sqrt(6.)*(-1 + 4*M**4 + 5*N**4 + Msq*(-1 + 15*Nsq))* &
       u_V_DFD)/4.
              rVd_t = (root_3*N*(sqrt(10.)* &
         (4*M**4 + 5*Msq*(-1 + Nsq) + (-1 + Nsq)**2)*u_Vd_DFS -  &
        sqrt(5.)*(1 + 8*M**4 - 3*Nsq + 2*N**4 + 2*Msq*(-4 + 5*Nsq))* &
         u_Vd_DFP + sqrt(2.)*(-1 + 4*M**4 + N**4 + Msq*(-1 + 5*Nsq))* &
         u_Vd_DFD))/4.
            end select
          case (5)
            select case (orb_dir_j)
              case (1)
              rVd(1) = Sqrt(7.5)*L*M*(3 - 5*Msq - 3*Nsq)*u_V_DFS
              rVd(2) = (root_3*(3*sqrt(10.)*L**4*u_V_DFS -  &
        12*sqrt(10.)*Lsq*Msq*u_V_DFS - 4*sqrt(5.)*u_V_DFP +  &
        3*sqrt(2.)*u_V_DFD + Nsq* &
         (10*sqrt(5.)*u_V_DFP - 2*sqrt(2.)*u_V_DFD) -  &
        30*Msq*(-1 + Nsq)*(2*sqrt(5.)*u_V_DFP - sqrt(2.)*u_V_DFD) +  &
        N**4*(-6*sqrt(5.)*u_V_DFP + 3*sqrt(2.)*u_V_DFD) +  &
        5*M**4*(sqrt(10.)*u_V_DFS - 16*sqrt(5.)*u_V_DFP +  &
           8*sqrt(2.)*u_V_DFD)))/8.
              rVd(3) = (M*N*(sqrt(15.)*(5 - 10*Msq - 6*Nsq)*u_V_DFP +  &
        sqrt(6.)*(-1 + 5*Msq + 3*Nsq)*u_V_DFD))/2.
              rVd_t = (root_3*M*(sqrt(10.)*(3*L**4 - 4*Lsq*Msq + M**4)*u_Vd_DFS -  &
        2*sqrt(5.)*(2 + 8*M**4 - 5*Nsq + 3*N**4 + 10*Msq*(-1 + Nsq))* &
         u_Vd_DFP + sqrt(2.)*(3 + 8*M**4 - 2*Nsq + 3*N**4 +  &
           10*Msq*(-1 + Nsq))*u_Vd_DFD))/8.
              case (2)
              rVd(1) = (-3*M*(-3*Lsq + Msq)*N* &
      (sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD))/2.
              rVd(2) = (3*L*(Lsq - 3*Msq)*N* &
      (sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD))/2.
              rVd(3) = (3*L*(L - M)*M*(L + M)* &
      (sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD))/2.
              rVd_t = (3*L*(L - M)*M*(L + M)*N* &
      (sqrt(5.)*u_Vd_DFS - sqrt(10.)*u_Vd_DFP + u_Vd_DFD))/2.
              case (3)
              rVd(1) = 0
              rVd(2) = (-3*sqrt(2.)*(-1 + 6*Msq + Nsq)*(-1 + 5*Nsq)*u_V_DFS +  &
      2*(2 - 21*Nsq + 15*N**4 + Msq*(-6 + 90*Nsq))*u_V_DFP +  &
      sqrt(10.)*(1 + 6*Nsq - 3*N**4 - 6*Msq*(1 + 3*Nsq))*u_V_DFD)/8.
              rVd(3) = (M*N*(-3*sqrt(2.)*(-3 + 5*Msq + 5*Nsq)*u_V_DFS +  &
        3*(-7 + 10*Msq + 10*Nsq)*u_V_DFP -  &
        3*sqrt(10.)*(-1 + Msq + Nsq)*u_V_DFD))/2.
              rVd_t = -(M*(3*sqrt(2.)*(-1 + 2*Msq + Nsq)*(-1 + 5*Nsq)*u_Vd_DFS +  &
         (-4 + 42*Nsq - 30*N**4 + Msq*(4 - 60*Nsq))*u_Vd_DFP +  &
         sqrt(10.)*(-1 - 6*Nsq + 3*N**4 + Msq*(2 + 6*Nsq))*u_Vd_DFD))/8.
              case (4)
              rVd(1) = (root_3*L*N*((-3 + 5*Nsq)*u_V_DFS +  &
        sqrt(2.)*(1 - 5*Nsq)*u_V_DFP + sqrt(5.)*(1 + Nsq)*u_V_DFD))/2.
              rVd(2) = -(root_3*M*N*((-3 + 5*Nsq)*u_V_DFS +  &
         sqrt(2.)*(1 - 5*Nsq)*u_V_DFP + sqrt(5.)*(1 + Nsq)*u_V_DFD))/2.
              rVd(3) = -(root_3*(-1 + 2*Msq + Nsq)* &
       (3*(-1 + 5*Nsq)*u_V_DFS + sqrt(2.)*(1 - 15*Nsq)*u_V_DFP +  &
         sqrt(5.)*(1 + 3*Nsq)*u_V_DFD))/4.
              rVd_t = (root_3*(L - M)*(L + M)*N* &
      ((-3 + 5*Nsq)*u_Vd_DFS + sqrt(2.)*(1 - 5*Nsq)*u_Vd_DFP +  &
        sqrt(5.)*(1 + Nsq)*u_Vd_DFD))/4.
              case (5)
              rVd(1) = (-3*sqrt(2.)*(-3 + 4*Msq + 3*Nsq)*(-1 + 5*Nsq)*u_V_DFS +  &
      sqrt(10.)*u_V_DFD + N**4*(30*u_V_DFP - 3*sqrt(10.)*u_V_DFD) -  &
      2*Nsq*(11*u_V_DFP + sqrt(10.)*u_V_DFD) +  &
      2*Msq*((-2 + 30*Nsq)*u_V_DFP - sqrt(10.)*(1 + 3*Nsq)*u_V_DFD))/8.
              rVd(2) = (L*M*(3*sqrt(2.)*(1 - 5*Nsq)*u_V_DFS +  &
        (-4 + 60*Nsq)*u_V_DFP - 2*sqrt(10.)*(1 + 3*Nsq)*u_V_DFD))/4.
              rVd(3) = (L*N*(-15*sqrt(2.)*(-1 + 2*Msq + Nsq)*u_V_DFS +  &
        (-22 + 60*Msq + 60*Nsq)*u_V_DFP -  &
        2*sqrt(10.)*(1 + 3*Msq + 3*Nsq)*u_V_DFD))/4.
              rVd_t = (L*(3*sqrt(2.)*(L - M)*(L + M)*(-1 + 5*Nsq)*u_Vd_DFS +  &
        2*(Nsq*(-11 + 15*Nsq) + Msq*(-2 + 30*Nsq))*u_Vd_DFP -  &
        sqrt(10.)*(-1 + 2*Nsq + 3*N**4 + Msq*(2 + 6*Nsq))*u_Vd_DFD))/8.
              case (6)
              rVd(1) = 3*sqrt(5.)*L*(Lsq - Msq)*N*u_V_DFS
              rVd(2) = 3*M*N*(-1 + 2*Msq + Nsq)* &
    (sqrt(5.)*u_V_DFS - 2*sqrt(10.)*u_V_DFP + 2*u_V_DFD)
              rVd(3) = (3*sqrt(5.)*(Lsq - Msq)**2*u_V_DFS -  &
      sqrt(10.)*(1 + 12*M**4 - 4*Nsq + 3*N**4 + 12*Msq*(-1 + Nsq))* &
       u_V_DFP + (-1 + 12*M**4 + 2*Nsq + 3*N**4 + 12*Msq*(-1 + Nsq))* &
       u_V_DFD + N*(-4*sqrt(10.)*N*(-2 + 6*Msq + 3*Nsq)*u_V_DFP +  &
         4*(N + 6*Msq*N + 3*N**3)*u_V_DFD))/4.
              rVd_t = (N*(3*sqrt(5.)*(Lsq - Msq)**2*u_Vd_DFS -  &
        sqrt(10.)*(1 + 12*M**4 - 4*Nsq + 3*N**4 + 12*Msq*(-1 + Nsq))* &
         u_Vd_DFP + (-1 + 12*M**4 + 2*Nsq + 3*N**4 + 12*Msq*(-1 + Nsq))* &
         u_Vd_DFD))/4.
              case (7)
              rVd(1) = (root_3*(4*sqrt(10.)*L*(L**3 - 2*L*Msq)*u_V_DFS +  &
        sqrt(10.)*(L**4 - 4*Lsq*Msq + 3*M**4)*u_V_DFS -  &
        2*sqrt(5.)*(8*M**4 + 6*Msq*(-1 + Nsq) + Nsq*(-1 + Nsq))* &
         u_V_DFP + sqrt(2.)*(8*M**4 + 6*Msq*(-1 + Nsq) + (1 + Nsq)**2)* &
         u_V_DFD))/8.
              rVd(2) = (L*M*(sqrt(30.)*(-2 + 5*Msq + 2*Nsq)*u_V_DFS -  &
        (-3 + 8*Msq + 3*Nsq)*(2*sqrt(15.)*u_V_DFP - sqrt(6.)*u_V_DFD)))/ &
    2.
              rVd(3) = (L*N*(sqrt(15.)*(1 - 6*Msq - 2*Nsq)*u_V_DFP +  &
        sqrt(6.)*(1 + 3*Msq + Nsq)*u_V_DFD))/2.
              rVd_t = (root_3*L*(sqrt(10.)*(L**4 - 4*Lsq*Msq + 3*M**4)*u_Vd_DFS -  &
        2*sqrt(5.)*(8*M**4 + 6*Msq*(-1 + Nsq) + Nsq*(-1 + Nsq))* &
         u_Vd_DFP + sqrt(2.)*(8*M**4 + 6*Msq*(-1 + Nsq) + (1 + Nsq)**2)* &
         u_Vd_DFD))/8.
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
              rVd(1) = -3*Sqrt(2.5)*L*M*u_V_SFS
              rVd(2) = (3*Sqrt(2.5)*(-Lsq + Msq)*u_V_SFS)/2.
              rVd(3) = 0
              rVd_t = (Sqrt(2.5)*M*(-3*Lsq + Msq)*u_Vd_SFS)/2.
            end select
          case (2)
            select case (orb_dir_j)
              case (1)
              rVd(1) = -(sqrt(15.)*M*N*u_V_SFS)
              rVd(2) = -(sqrt(15.)*L*N*u_V_SFS)
              rVd(3) = -(sqrt(15.)*L*M*u_V_SFS)
              rVd_t = -(sqrt(15.)*L*M*N*u_Vd_SFS)
            end select
          case (3)
            select case (orb_dir_j)
              case (1)
              rVd(1) = 0
              rVd(2) = (Sqrt(1.5)*(1 - 5*Nsq)*u_V_SFS)/2.
              rVd(3) = -5*Sqrt(1.5)*M*N*u_V_SFS
              rVd_t = (Sqrt(1.5)*M*(1 - 5*Nsq)*u_Vd_SFS)/2.
            end select
          case (4)
            select case (orb_dir_j)
              case (1)
              rVd(1) = 0
              rVd(2) = 0
              rVd(3) = (-3*(-1 + 5*Nsq)*u_V_SFS)/2.
              rVd_t = (N*(3 - 5*Nsq)*u_Vd_SFS)/2.
            end select
          case (5)
            select case (orb_dir_j)
              case (1)
              rVd(1) = (Sqrt(1.5)*(1 - 5*Nsq)*u_V_SFS)/2.
              rVd(2) = 0
              rVd(3) = -5*Sqrt(1.5)*L*N*u_V_SFS
              rVd_t = (Sqrt(1.5)*L*(1 - 5*Nsq)*u_Vd_SFS)/2.
            end select
          case (6)
            select case (orb_dir_j)
              case (1)
              rVd(1) = -(sqrt(15.)*L*N*u_V_SFS)
              rVd(2) = sqrt(15.)*M*N*u_V_SFS
              rVd(3) = -(sqrt(15.)*(L - M)*(L + M)*u_V_SFS)/2.
              rVd_t = -(sqrt(15.)*(L - M)*(L + M)*N*u_Vd_SFS)/2.
            end select
          case (7)
            select case (orb_dir_j)
              case (1)
              rVd(1) = (Sqrt(2.5)*(-1 + 4*Msq + Nsq)*u_V_SFS)/2.
              rVd(2) = 2*sqrt(10.)*L*M*u_V_SFS
              rVd(3) = Sqrt(2.5)*L*N*u_V_SFS
              rVd_t = (Sqrt(2.5)*L*(-1 + 4*Msq + Nsq)*u_Vd_SFS)/2.
            end select
        end select
      case (ORB_P)
        select case (orb_dir_i)
          case (1)
            select case (orb_dir_j)
              case (1)
              rVd(1) = 0
              rVd(2) = (M*(sqrt(10.)*(3 - 8*Msq - 3*Nsq)*u_V_PFS +  &
        sqrt(15.)*(-5 + 8*Msq + 3*Nsq)*u_V_PFP))/2.
              rVd(3) = (-(sqrt(15.)*N*u_V_PFP) +  &
      Msq*N*(-3*sqrt(10.)*u_V_PFS + 3*sqrt(15.)*u_V_PFP))/2.
              rVd_t = (sqrt(5.)*(-(root_3*(-1 + Nsq)*u_Vd_PFP) +  &
        M**4*(-4*sqrt(2.)*u_Vd_PFS + 4*root_3*u_Vd_PFP) +  &
        Msq*(-3*sqrt(2.)*(-1 + Nsq)*u_Vd_PFS +  &
           root_3*(-5 + 3*Nsq)*u_Vd_PFP)))/4.
              case (2)
              rVd(1) = 0
              rVd(2) = (-3*N*(-1 + 4*Msq + Nsq)* &
      (sqrt(10.)*u_V_PFS - sqrt(15.)*u_V_PFP))/4.
              rVd(3) = -(M*(-3 + 4*Msq + 9*Nsq)* &
       (sqrt(10.)*u_V_PFS - sqrt(15.)*u_V_PFP))/4.
              rVd_t = -(M*N*(-3 + 4*Msq + 3*Nsq)* &
       (sqrt(10.)*u_Vd_PFS - sqrt(15.)*u_Vd_PFP))/4.
              case (3)
              rVd(1) = (M*(sqrt(10.)*(3 - 4*Msq - 3*Nsq)*u_V_PFS +  &
        sqrt(15.)*(-1 + 4*Msq + 3*Nsq)*u_V_PFP))/4.
              rVd(2) = (L*(-3*sqrt(10.)*(-1 + 4*Msq + Nsq)*u_V_PFS +  &
        sqrt(15.)*(-1 + 12*Msq + 3*Nsq)*u_V_PFP))/4.
              rVd(3) = (-3*L*M*N*(sqrt(10.)*u_V_PFS - sqrt(15.)*u_V_PFP))/2.
              rVd_t = (L*M*(sqrt(10.)*(3 - 4*Msq - 3*Nsq)*u_Vd_PFS +  &
        sqrt(15.)*(-1 + 4*Msq + 3*Nsq)*u_Vd_PFP))/4.
            end select
          case (2)
            select case (orb_dir_j)
              case (1)
              rVd(1) = Sqrt(2.5)*N*u_V_PFP +  &
    Msq*N*(sqrt(15.)*u_V_PFS - 3*Sqrt(2.5)*u_V_PFP)
              rVd(2) = L*M*N*(2*sqrt(15.)*u_V_PFS - 3*sqrt(10.)*u_V_PFP)
              rVd(3) = L*(Sqrt(2.5)*u_V_PFP +  &
      Msq*(sqrt(15.)*u_V_PFS - 3*Sqrt(2.5)*u_V_PFP))
              rVd_t = (L*N*(sqrt(10.)*u_Vd_PFP +  &
        Msq*(2*sqrt(15.)*u_Vd_PFS - 3*sqrt(10.)*u_Vd_PFP)))/2.
              case (2)
              rVd(1) = (M*(2*sqrt(15.)*Nsq*u_V_PFS +  &
        sqrt(10.)*(1 - 3*Nsq)*u_V_PFP))/2.
              rVd(2) = (L*(2*sqrt(15.)*Nsq*u_V_PFS +  &
        sqrt(10.)*(1 - 3*Nsq)*u_V_PFP))/2.
              rVd(3) = L*M*N*(2*sqrt(15.)*u_V_PFS - 3*sqrt(10.)*u_V_PFP)
              rVd_t = (L*M*(2*sqrt(15.)*Nsq*u_Vd_PFS +  &
        sqrt(10.)*(1 - 3*Nsq)*u_Vd_PFP))/2.
              case (3)
              rVd(1) = (L*M*N*(2*sqrt(15.)*(-1 + Nsq)*u_V_PFS +  &
        sqrt(10.)*(2 - 3*Nsq)*u_V_PFP))/(-1 + Nsq)
              rVd(2) = (sqrt(5.)*N*(-3*sqrt(2.)*Msq*u_V_PFP +  &
        Lsq*(2*root_3*(-1 + Nsq)*u_V_PFS +  &
           sqrt(2.)*(2 - 3*Nsq)*u_V_PFP)))/(2.*(-1 + Nsq))
              rVd(3) = (sqrt(5.)*M*(sqrt(2.)*Msq*(1 + Nsq)*u_V_PFP +  &
        Lsq*(2*root_3*(-1 + Nsq)**2*u_V_PFS +  &
           sqrt(2.)*(-2 + 7*Nsq - 3*N**4)*u_V_PFP)))/(2.*(-1 + Nsq)**2)
              rVd_t = -(sqrt(5.)*M*N*(sqrt(2.)*Msq*u_Vd_PFP +  &
         Lsq*(-2*root_3*(-1 + Nsq)*u_Vd_PFS +  &
            sqrt(2.)*(-2 + 3*Nsq)*u_Vd_PFP)))/(2.*(-1 + Nsq))
            end select
          case (3)
            select case (orb_dir_j)
              case (1)
              rVd(1) = 0
              rVd(2) = (M*(sqrt(6.)*(-1 + 5*Nsq)*u_V_PFS +  &
        (1 - 15*Nsq)*u_V_PFP))/2.
              rVd(3) = (5*N*(Msq*(sqrt(6.)*u_V_PFS - 3*u_V_PFP) + u_V_PFP))/2.
              rVd_t = ((-1 + 5*Nsq)*u_Vd_PFP +  &
      Msq*(sqrt(6.)*(-1 + 5*Nsq)*u_Vd_PFS + (1 - 15*Nsq)*u_Vd_PFP))/4.
              case (2)
              rVd(1) = 0
              rVd(2) = (N*(sqrt(6.)*(-1 + 5*Nsq)*u_V_PFS +  &
        (11 - 15*Nsq)*u_V_PFP))/4.
              rVd(3) = (M*(sqrt(6.)*(-1 + 15*Nsq)*u_V_PFS +  &
        (11 - 45*Nsq)*u_V_PFP))/4.
              rVd_t = (M*N*(sqrt(6.)*(-1 + 5*Nsq)*u_Vd_PFS +  &
        (11 - 15*Nsq)*u_Vd_PFP))/4.
              case (3)
              rVd(1) = (M*(sqrt(6.)*(-1 + 5*Nsq)*u_V_PFS +  &
        (1 - 15*Nsq)*u_V_PFP))/4.
              rVd(2) = (L*(sqrt(6.)*(-1 + 5*Nsq)*u_V_PFS +  &
        (1 - 15*Nsq)*u_V_PFP))/4.
              rVd(3) = (5*L*M*N*(sqrt(6.)*u_V_PFS - 3*u_V_PFP))/2.
              rVd_t = (L*M*(sqrt(6.)*(-1 + 5*Nsq)*u_Vd_PFS + (1 - 15*Nsq)*u_Vd_PFP))/ &
    4.
            end select
          case (4)
            select case (orb_dir_j)
              case (1)
              rVd(1) = 0
              rVd(2) = (N*(2*(-3 + 5*Nsq)*u_V_PFS +  &
        sqrt(6.)*(1 - 5*Nsq)*u_V_PFP))/4.
              rVd(3) = (M*(6*(-1 + 5*Nsq)*u_V_PFS +  &
        sqrt(6.)*(1 - 15*Nsq)*u_V_PFP))/4.
              rVd_t = (M*N*(2*(-3 + 5*Nsq)*u_Vd_PFS +  &
        sqrt(6.)*(1 - 5*Nsq)*u_Vd_PFP))/4.
              case (2)
              rVd(1) = 0
              rVd(2) = 0
              rVd(3) = N*((-3 + 10*Nsq)*u_V_PFS + sqrt(6.)*(3 - 5*Nsq)*u_V_PFP)
              rVd_t = (2*Nsq*(-3 + 5*Nsq)*u_Vd_PFS +  &
      sqrt(6.)*(-1 + 6*Nsq - 5*N**4)*u_Vd_PFP)/4.
              case (3)
              rVd(1) = (N*(2*(-3 + 5*Nsq)*u_V_PFS +  &
        sqrt(6.)*(1 - 5*Nsq)*u_V_PFP))/4.
              rVd(2) = 0
              rVd(3) = (L*(6*(-1 + 5*Nsq)*u_V_PFS +  &
        sqrt(6.)*(1 - 15*Nsq)*u_V_PFP))/4.
              rVd_t = (L*N*(2*(-3 + 5*Nsq)*u_Vd_PFS +  &
        sqrt(6.)*(1 - 5*Nsq)*u_Vd_PFP))/4.
            end select
          case (5)
            select case (orb_dir_j)
              case (1)
              rVd(1) = (M*(sqrt(6.)*(-1 + 5*Nsq)*u_V_PFS +  &
        (1 - 15*Nsq)*u_V_PFP))/4.
              rVd(2) = (L*(sqrt(6.)*(-1 + 5*Nsq)*u_V_PFS +  &
        (1 - 15*Nsq)*u_V_PFP))/4.
              rVd(3) = (5*L*M*N*(sqrt(6.)*u_V_PFS - 3*u_V_PFP))/2.
              rVd_t = (L*M*(sqrt(6.)*(-1 + 5*Nsq)*u_Vd_PFS + (1 - 15*Nsq)*u_Vd_PFP))/ &
    4.
              case (2)
              rVd(1) = (N*(sqrt(6.)*(-1 + 5*Nsq)*u_V_PFS +  &
        (11 - 15*Nsq)*u_V_PFP))/4.
              rVd(2) = 0
              rVd(3) = (L*(sqrt(6.)*(-1 + 15*Nsq)*u_V_PFS +  &
        (11 - 45*Nsq)*u_V_PFP))/4.
              rVd_t = (L*N*(sqrt(6.)*(-1 + 5*Nsq)*u_Vd_PFS +  &
        (11 - 15*Nsq)*u_Vd_PFP))/4.
              case (3)
              rVd(1) = (L*(sqrt(6.)*(1 - 6*Nsq + 5*N**4)*u_V_PFS +  &
        Nsq*(11 - 15*Nsq)*u_V_PFP))/(2.*(-1 + Nsq))
              rVd(2) = (M*(1 - 5*Nsq)*u_V_PFP)/(2.*(-1 + Nsq))
              rVd(3) = (N*(4*Msq*u_V_PFP +  &
        Lsq*(5*sqrt(6.)*(-1 + Nsq)**2*u_V_PFS +  &
           (-11 + 30*Nsq - 15*N**4)*u_V_PFP)))/(2.*(-1 + Nsq)**2)
              rVd_t = (sqrt(6.)*Lsq*(1 - 6*Nsq + 5*N**4)*u_Vd_PFS +  &
      (Lsq*Nsq*(11 - 15*Nsq) + Msq*(1 - 5*Nsq))*u_Vd_PFP)/ &
    (4.*(-1 + Nsq))
            end select
          case (6)
            select case (orb_dir_j)
              case (1)
              rVd(1) = 0
              rVd(2) = (N*(-2*sqrt(15.)*(-1 + 6*Msq + Nsq)*u_V_PFS +  &
        sqrt(10.)*(-5 + 18*Msq + 3*Nsq)*u_V_PFP))/4.
              rVd(3) = (M*(-2*sqrt(15.)*(-1 + 2*Msq + 3*Nsq)*u_V_PFS +  &
        sqrt(10.)*(-5 + 6*Msq + 9*Nsq)*u_V_PFP))/4.
              rVd_t = (M*N*(-2*sqrt(15.)*(-1 + 2*Msq + Nsq)*u_Vd_PFS +  &
        sqrt(10.)*(-5 + 6*Msq + 3*Nsq)*u_Vd_PFP))/4.
              case (2)
              rVd(1) = L*(Sqrt(2.5)*u_V_PFP +  &
      Nsq*(sqrt(15.)*u_V_PFS - 3*Sqrt(2.5)*u_V_PFP))
              rVd(2) = -(M*(sqrt(10.)*u_V_PFP +  &
         Nsq*(2*sqrt(15.)*u_V_PFS - 3*sqrt(10.)*u_V_PFP)))/2.
              rVd(3) = -(N*(-1 + 2*Msq + Nsq)* &
       (2*sqrt(15.)*u_V_PFS - 3*sqrt(10.)*u_V_PFP))/2.
              rVd_t = (sqrt(5.)*(L - M)*(L + M)* &
      (2*root_3*Nsq*u_Vd_PFS + sqrt(2.)*(1 - 3*Nsq)*u_Vd_PFP))/4.
              case (3)
              rVd(1) = (N*(-2*sqrt(15.)*(-1 + 2*Msq + Nsq)*u_V_PFS +  &
        sqrt(10.)*(-1 + 6*Msq + 3*Nsq)*u_V_PFP))/4.
              rVd(2) = L*M*N*(-2*sqrt(15.)*u_V_PFS + 3*sqrt(10.)*u_V_PFP)
              rVd(3) = (L*(-2*sqrt(15.)*(-1 + 2*Msq + 3*Nsq)*u_V_PFS +  &
        sqrt(10.)*(-1 + 6*Msq + 9*Nsq)*u_V_PFP))/4.
              rVd_t = (L*N*(-2*sqrt(15.)*(-1 + 2*Msq + Nsq)*u_Vd_PFS +  &
        sqrt(10.)*(-1 + 6*Msq + 3*Nsq)*u_Vd_PFP))/4.
            end select
          case (7)
            select case (orb_dir_j)
              case (1)
              rVd(1) = (M*(-(sqrt(10.)*(-1 + 4*Msq + Nsq)*u_V_PFS) +  &
        sqrt(15.)*(-3 + 4*Msq + Nsq)*u_V_PFP))/4.
              rVd(2) = (L*(-(sqrt(10.)*(-1 + 12*Msq + Nsq)*u_V_PFS) +  &
        sqrt(15.)*(-3 + 12*Msq + Nsq)*u_V_PFP))/4.
              rVd(3) = (L*M*N*(-(sqrt(10.)*u_V_PFS) + sqrt(15.)*u_V_PFP))/2.
              rVd_t = (L*M*(-(sqrt(10.)*(-1 + 4*Msq + Nsq)*u_Vd_PFS) +  &
        sqrt(15.)*(-3 + 4*Msq + Nsq)*u_Vd_PFP))/4.
              case (2)
              rVd(1) = -(N*(-1 + 4*Msq + Nsq)* &
       (sqrt(10.)*u_V_PFS - sqrt(15.)*u_V_PFP))/4.
              rVd(2) = L*M*N*(-2*sqrt(10.)*u_V_PFS + 2*sqrt(15.)*u_V_PFP)
              rVd(3) = -(L*(-1 + 4*Msq + 3*Nsq)* &
       (sqrt(10.)*u_V_PFS - sqrt(15.)*u_V_PFP))/4.
              rVd_t = -(L*N*(-1 + 4*Msq + Nsq)* &
       (sqrt(10.)*u_Vd_PFS - sqrt(15.)*u_Vd_PFP))/4.
              case (3)
              rVd(1) = (sqrt(5.)*(6*L*Msq*(-1 + Nsq)* &
         (-(sqrt(2.)*u_V_PFS) + root_3*u_V_PFP) +  &
        4*L**3*(sqrt(2.)*(-1 + Nsq)*u_V_PFS - root_3*Nsq*u_V_PFP)))/ &
    (4.*(-1 + Nsq))
              rVd(2) = (sqrt(5.)*(4*root_3*M**3*u_V_PFP +  &
        6*Lsq*M*(-1 + Nsq)*(-(sqrt(2.)*u_V_PFS) + root_3*u_V_PFP)))/ &
    (4.*(-1 + Nsq))
              rVd(3) = (sqrt(15.)*(L**4 - M**4)*N*u_V_PFP)/(2.*(-1 + Nsq)**2)
              rVd_t = (sqrt(5.)*(root_3*M**4*u_Vd_PFP +  &
        3*Lsq*Msq*(-1 + Nsq)* &
         (-(sqrt(2.)*u_Vd_PFS) + root_3*u_Vd_PFP) +  &
        L**4*(sqrt(2.)*(-1 + Nsq)*u_Vd_PFS - root_3*Nsq*u_Vd_PFP)))/ &
    (4.*(-1 + Nsq))
            end select
        end select
      case (ORB_D)
        select case (orb_dir_i)
          case (1)
            select case (orb_dir_j)
              case (1)
              rVd(1) = (root_3*(-6*sqrt(10.)*Lsq*Msq*u_V_DFS +  &
        sqrt(10.)*Msq*(-3*Lsq + Msq)*u_V_DFS -  &
        sqrt(5.)*(1 + 8*M**4 - Nsq + 6*Msq*(-1 + Nsq))*u_V_DFP +  &
        sqrt(2.)*(4*M**4 - 2*Nsq + 3*Msq*(-1 + Nsq))*u_V_DFD))/4.
              rVd(2) = (L*M*(sqrt(30.)*(-3 + 5*Msq + 3*Nsq)*u_V_DFS -  &
        (-3 + 8*Msq + 3*Nsq)*(2*sqrt(15.)*u_V_DFP - sqrt(6.)*u_V_DFD)))/ &
    2.
              rVd(3) = (L*N*(sqrt(15.)*(1 - 6*Msq)*u_V_DFP +  &
        sqrt(6.)*(-2 + 3*Msq)*u_V_DFD))/2.
              rVd_t = -(root_3*L*(-(sqrt(10.)*Msq*(-3*Lsq + Msq)*u_Vd_DFS) +  &
         sqrt(5.)*(1 + 8*M**4 - Nsq + 6*Msq*(-1 + Nsq))*u_Vd_DFP -  &
         sqrt(2.)*(4*M**4 - 2*Nsq + 3*Msq*(-1 + Nsq))*u_Vd_DFD))/4.
              case (2)
              rVd(1) = 0
              rVd(2) = (M*N*(sqrt(30.)*(-3 + 8*Msq + 3*Nsq)*u_V_DFS -  &
        2*sqrt(15.)*(-4 + 8*Msq + 3*Nsq)*u_V_DFP +  &
        sqrt(6.)*(-7 + 8*Msq + 3*Nsq)*u_V_DFD))/2.
              rVd(3) = ((-1 + 3*Nsq)*(sqrt(15.)*u_V_DFP - 2*sqrt(6.)*u_V_DFD) +  &
      4*M**4*(sqrt(30.)*u_V_DFS - 2*sqrt(15.)*u_V_DFP +  &
         sqrt(6.)*u_V_DFD) + Msq* &
       (3*sqrt(30.)*(-1 + 3*Nsq)*u_V_DFS +  &
         2*sqrt(15.)*(4 - 9*Nsq)*u_V_DFP + sqrt(6.)*(-7 + 9*Nsq)*u_V_DFD &
 ))/4.
              rVd_t = (N*((-1 + Nsq)*(sqrt(15.)*u_Vd_DFP - 2*sqrt(6.)*u_Vd_DFD) +  &
        4*M**4*(sqrt(30.)*u_Vd_DFS - 2*sqrt(15.)*u_Vd_DFP +  &
           sqrt(6.)*u_Vd_DFD) +  &
        Msq*(3*sqrt(30.)*(-1 + Nsq)*u_Vd_DFS +  &
           2*sqrt(15.)*(4 - 3*Nsq)*u_Vd_DFP +  &
           sqrt(6.)*(-7 + 3*Nsq)*u_Vd_DFD)))/4.
              case (3)
              rVd(1) = (-3*L*M*(sqrt(10.)*(-1 + 3*Nsq)*u_V_DFS -  &
        6*sqrt(5.)*Nsq*u_V_DFP + 3*sqrt(2.)*(1 + Nsq)*u_V_DFD))/4.
              rVd(2) = (-3*(Lsq - Msq)* &
      (sqrt(10.)*(-1 + 3*Nsq)*u_V_DFS - 6*sqrt(5.)*Nsq*u_V_DFP +  &
        3*sqrt(2.)*(1 + Nsq)*u_V_DFD))/8.
              rVd(3) = (3*M*N*(-3 + 4*Msq + 3*Nsq)* &
      (sqrt(10.)*u_V_DFS - 2*sqrt(5.)*u_V_DFP + sqrt(2.)*u_V_DFD))/4.
              rVd_t = -((3*Lsq*M - M**3)*(sqrt(10.)*(-1 + 3*Nsq)*u_Vd_DFS -  &
         6*sqrt(5.)*Nsq*u_Vd_DFP + 3*sqrt(2.)*(1 + Nsq)*u_Vd_DFD))/8.
              case (4)
              rVd(1) = (M*N*(sqrt(30.)*(-3 + 4*Msq + 3*Nsq)*u_V_DFS -  &
        2*sqrt(15.)*(-2 + 4*Msq + 3*Nsq)*u_V_DFP +  &
        sqrt(6.)*(1 + 4*Msq + 3*Nsq)*u_V_DFD))/4.
              rVd(2) = (L*N*(3*sqrt(30.)*(-1 + 4*Msq + Nsq)*u_V_DFS -  &
        2*sqrt(15.)*(-2 + 12*Msq + 3*Nsq)*u_V_DFP +  &
        sqrt(6.)*(1 + 12*Msq + 3*Nsq)*u_V_DFD))/4.
              rVd(3) = (L*M*(sqrt(30.)*(-3 + 4*Msq + 9*Nsq)*u_V_DFS -  &
        2*sqrt(15.)*(-2 + 4*Msq + 9*Nsq)*u_V_DFP +  &
        sqrt(6.)*(1 + 4*Msq + 9*Nsq)*u_V_DFD))/4.
              rVd_t = (L*M*N*(sqrt(30.)*(-3 + 4*Msq + 3*Nsq)*u_Vd_DFS -  &
        2*sqrt(15.)*(-2 + 4*Msq + 3*Nsq)*u_Vd_DFP +  &
        sqrt(6.)*(1 + 4*Msq + 3*Nsq)*u_Vd_DFD))/4.
              case (5)
              rVd(1) = Sqrt(7.5)*M*(-3*L**3 + 2*L*Msq)*u_V_DFS
              rVd(2) = (-3*sqrt(30.)*L**4*u_V_DFS +  &
      12*sqrt(30.)*Lsq*Msq*u_V_DFS + 4*sqrt(15.)*u_V_DFP -  &
      3*sqrt(6.)*u_V_DFD + N**4* &
       (6*sqrt(15.)*u_V_DFP - 3*sqrt(6.)*u_V_DFD) +  &
      30*Msq*(-1 + Nsq)*(2*sqrt(15.)*u_V_DFP - sqrt(6.)*u_V_DFD) +  &
      Nsq*(-10*sqrt(15.)*u_V_DFP + 2*sqrt(6.)*u_V_DFD) -  &
      5*M**4*(sqrt(30.)*u_V_DFS - 16*sqrt(15.)*u_V_DFP +  &
         8*sqrt(6.)*u_V_DFD))/8.
              rVd(3) = (M*N*(sqrt(15.)*(-5 + 10*Msq + 6*Nsq)*u_V_DFP +  &
        sqrt(6.)*(1 - 5*Msq - 3*Nsq)*u_V_DFD))/2.
              rVd_t = -(root_3*M*(sqrt(10.)*(3*L**4 - 4*Lsq*Msq + M**4)*u_Vd_DFS -  &
         2*sqrt(5.)*(2 + 8*M**4 - 5*Nsq + 3*N**4 + 10*Msq*(-1 + Nsq))* &
          u_Vd_DFP + sqrt(2.)*(3 + 8*M**4 - 2*Nsq + 3*N**4 +  &
            10*Msq*(-1 + Nsq))*u_Vd_DFD))/8.
            end select
          case (2)
            select case (orb_dir_j)
              case (1)
              rVd(1) = 0
              rVd(2) = 6*M*N*(-1 + 2*Msq + Nsq)* &
    (sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD)
              rVd(3) = (sqrt(10.)*(-1 + 3*Nsq)*u_V_DFP +  &
      2*(1 - 6*Nsq)*u_V_DFD +  &
      6*M**4*(sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD) +  &
      6*Msq*(-1 + 3*Nsq)*(sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP +  &
         u_V_DFD))/2.
              rVd_t = (N*(sqrt(10.)*(-1 + Nsq)*u_Vd_DFP + 2*(1 - 2*Nsq)*u_Vd_DFD +  &
        6*M**4*(sqrt(5.)*u_Vd_DFS - sqrt(10.)*u_Vd_DFP + u_Vd_DFD) +  &
        6*Msq*(-1 + Nsq)*(sqrt(5.)*u_Vd_DFS - sqrt(10.)*u_Vd_DFP +  &
           u_Vd_DFD)))/2.
              case (2)
              rVd(1) = (-2*u_V_DFD + Nsq*(-(sqrt(10.)*u_V_DFP) + 4*u_V_DFD) +  &
      Msq*(-(sqrt(10.)*u_V_DFP) + 4*u_V_DFD -  &
         6*Nsq*(sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD)))/2.
              rVd(2) = L*M*(-(sqrt(10.)*u_V_DFP) + 4*u_V_DFD -  &
      6*Nsq*(sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD))
              rVd(3) = L*N*(-(sqrt(10.)*u_V_DFP) + 4*u_V_DFD -  &
      6*Msq*(sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD))
              rVd_t = -(L*(Nsq*(sqrt(10.)*u_Vd_DFP - 4*u_Vd_DFD) + 2*u_Vd_DFD +  &
         Msq*(sqrt(10.)*u_Vd_DFP - 4*u_Vd_DFD +  &
            6*Nsq*(sqrt(5.)*u_Vd_DFS - sqrt(10.)*u_Vd_DFP + u_Vd_DFD))))/ &
    2.
              case (3)
              rVd(1) = -(root_3*M*N*(-1 + 3*Nsq)* &
       (sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD))/2.
              rVd(2) = -(root_3*L*N*(-1 + 3*Nsq)* &
       (sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD))/2.
              rVd(3) = -(root_3*L*M*(-1 + 9*Nsq)* &
       (sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD))/2.
              rVd_t = -(root_3*L*M*N*(-1 + 3*Nsq)* &
       (sqrt(5.)*u_Vd_DFS - sqrt(10.)*u_Vd_DFP + u_Vd_DFD))/2.
              case (4)
              rVd(1) = 0
              rVd(2) = (sqrt(10.)*(-1 + 3*Msq)*u_V_DFP +  &
      2*(1 - 6*Msq)*u_V_DFD +  &
      6*(-1 + 3*Msq)*Nsq*(sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP +  &
         u_V_DFD) + 6*N**4*(sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP +  &
         u_V_DFD))/2.
              rVd(3) = 6*M*N*(-1 + Msq + 2*Nsq)* &
    (sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD)
              rVd_t = (M*(sqrt(10.)*(-1 + Msq)*u_Vd_DFP + 2*(1 - 2*Msq)*u_Vd_DFD +  &
        6*(-1 + Msq)*Nsq*(sqrt(5.)*u_Vd_DFS - sqrt(10.)*u_Vd_DFP +  &
           u_Vd_DFD) + 6*N**4* &
         (sqrt(5.)*u_Vd_DFS - sqrt(10.)*u_Vd_DFP + u_Vd_DFD)))/2.
              case (5)
              rVd(1) = (3*M*(-3*Lsq + Msq)*N* &
      (sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD))/2.
              rVd(2) = (3*L*N*(-1 + 4*Msq + Nsq)* &
      (sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD))/2.
              rVd(3) = (-3*L*(L - M)*M*(L + M)* &
      (sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD))/2.
              rVd_t = (-3*L*(L - M)*M*(L + M)*N* &
      (sqrt(5.)*u_Vd_DFS - sqrt(10.)*u_Vd_DFP + u_Vd_DFD))/2.
            end select
          case (3)
            select case (orb_dir_j)
              case (1)
              rVd(1) = ((1 - 5*Nsq)*u_V_DFP + 2*sqrt(10.)*Nsq*u_V_DFD -  &
      Msq*(3*sqrt(2.)*(-1 + 5*Nsq)*u_V_DFS + (2 - 30*Nsq)*u_V_DFP +  &
         sqrt(10.)*(1 + 3*Nsq)*u_V_DFD))/4.
              rVd(2) = -(L*M*(3*sqrt(2.)*(-1 + 5*Nsq)*u_V_DFS +  &
         (2 - 30*Nsq)*u_V_DFP + sqrt(10.)*(1 + 3*Nsq)*u_V_DFD))/2.
              rVd(3) = -(L*N*(5*u_V_DFP - 2*sqrt(10.)*u_V_DFD +  &
         3*Msq*(5*sqrt(2.)*u_V_DFS - 10*u_V_DFP + sqrt(10.)*u_V_DFD)))/ &
    2.
              rVd_t = (L*((1 - 5*Nsq)*u_Vd_DFP + 2*sqrt(10.)*Nsq*u_Vd_DFD -  &
        Msq*(3*sqrt(2.)*(-1 + 5*Nsq)*u_Vd_DFS + (2 - 30*Nsq)*u_Vd_DFP +  &
           sqrt(10.)*(1 + 3*Nsq)*u_Vd_DFD)))/4.
              case (2)
              rVd(1) = 0
              rVd(2) = (-3*M*N*(sqrt(2.)*(-1 + 5*Nsq)*u_V_DFS +  &
        (4 - 10*Nsq)*u_V_DFP + sqrt(10.)*(-1 + Nsq)*u_V_DFD))/2.
              rVd(3) = ((1 - 15*Nsq)*u_V_DFP +  &
      2*sqrt(10.)*(-1 + 3*Nsq)*u_V_DFD +  &
      3*Msq*(sqrt(2.)*(1 - 15*Nsq)*u_V_DFS + (-4 + 30*Nsq)*u_V_DFP +  &
         sqrt(10.)*(1 - 3*Nsq)*u_V_DFD))/4.
              rVd_t = (N*((1 - 5*Nsq)*u_Vd_DFP + 2*sqrt(10.)*(-1 + Nsq)*u_Vd_DFD -  &
        3*Msq*(sqrt(2.)*(-1 + 5*Nsq)*u_Vd_DFS + (4 - 10*Nsq)*u_Vd_DFP +  &
           sqrt(10.)*(-1 + Nsq)*u_Vd_DFD)))/4.
              case (3)
              rVd(1) = 0
              rVd(2) = -(root_3*(sqrt(2.)*(1 - 8*Nsq + 15*N**4)*u_V_DFS +  &
         2*Nsq*(11 - 15*Nsq)*u_V_DFP +  &
         sqrt(10.)*(1 - 4*Nsq + 3*N**4)*u_V_DFD))/8.
              rVd(3) = (root_3*M*N*(sqrt(2.)*(4 - 15*Nsq)*u_V_DFS +  &
        (-11 + 30*Nsq)*u_V_DFP + sqrt(10.)*(2 - 3*Nsq)*u_V_DFD))/2.
              rVd_t = -(root_3*M*(sqrt(2.)*(1 - 8*Nsq + 15*N**4)*u_Vd_DFS +  &
         2*Nsq*(11 - 15*Nsq)*u_Vd_DFP +  &
         sqrt(10.)*(1 - 4*Nsq + 3*N**4)*u_Vd_DFD))/8.
              case (4)
              rVd(1) = (-3*M*N*(sqrt(2.)*(-1 + 5*Nsq)*u_V_DFS +  &
        (4 - 10*Nsq)*u_V_DFP + sqrt(10.)*(-1 + Nsq)*u_V_DFD))/4.
              rVd(2) = (-3*L*N*(sqrt(2.)*(-1 + 5*Nsq)*u_V_DFS +  &
        (4 - 10*Nsq)*u_V_DFP + sqrt(10.)*(-1 + Nsq)*u_V_DFD))/4.
              rVd(3) = (3*L*M*(sqrt(2.)*(1 - 15*Nsq)*u_V_DFS +  &
        (-4 + 30*Nsq)*u_V_DFP + sqrt(10.)*(1 - 3*Nsq)*u_V_DFD))/4.
              rVd_t = (-3*L*M*N*(sqrt(2.)*(-1 + 5*Nsq)*u_Vd_DFS +  &
        (4 - 10*Nsq)*u_Vd_DFP + sqrt(10.)*(-1 + Nsq)*u_Vd_DFD))/4.
              case (5)
              rVd(1) = 0
              rVd(2) = (3*sqrt(2.)*(-1 + 6*Msq + Nsq)*(-1 + 5*Nsq)*u_V_DFS -  &
      2*(2 - 21*Nsq + 15*N**4 + Msq*(-6 + 90*Nsq))*u_V_DFP +  &
      sqrt(10.)*(-1 - 6*Nsq + 3*N**4 + 6*Msq*(1 + 3*Nsq))*u_V_DFD)/8.
              rVd(3) = (3*M*N*(sqrt(2.)*(-3 + 5*Msq + 5*Nsq)*u_V_DFS +  &
        (7 - 10*Msq - 10*Nsq)*u_V_DFP +  &
        sqrt(10.)*(-1 + Msq + Nsq)*u_V_DFD))/2.
              rVd_t = (M*(3*sqrt(2.)*(-1 + 2*Msq + Nsq)*(-1 + 5*Nsq)*u_Vd_DFS +  &
        (-4 + 42*Nsq - 30*N**4 + Msq*(4 - 60*Nsq))*u_Vd_DFP +  &
        sqrt(10.)*(-1 - 6*Nsq + 3*N**4 + Msq*(2 + 6*Nsq))*u_Vd_DFD))/8.
            end select
          case (4)
            select case (orb_dir_j)
              case (1)
              rVd(1) = -(root_3*M*N*((-3 + 5*Nsq)*u_V_DFS +  &
         sqrt(2.)*(1 - 5*Nsq)*u_V_DFP + sqrt(5.)*(1 + Nsq)*u_V_DFD))/2.
              rVd(2) = -(root_3*L*N*((-3 + 5*Nsq)*u_V_DFS +  &
         sqrt(2.)*(1 - 5*Nsq)*u_V_DFP + sqrt(5.)*(1 + Nsq)*u_V_DFD))/2.
              rVd(3) = -(root_3*L*M*(3*(-1 + 5*Nsq)*u_V_DFS +  &
         sqrt(2.)*(1 - 15*Nsq)*u_V_DFP + sqrt(5.)*(1 + 3*Nsq)*u_V_DFD))/ &
    2.
              rVd_t = -(root_3*L*M*N*((-3 + 5*Nsq)*u_Vd_DFS +  &
         sqrt(2.)*(1 - 5*Nsq)*u_Vd_DFP + sqrt(5.)*(1 + Nsq)*u_Vd_DFD))/2.
              case (2)
              rVd(1) = 0
              rVd(2) = -(root_3*(sqrt(2.)*(-1 + 7*Nsq - 10*N**4)*u_V_DFP +  &
         2*Nsq*((-3 + 5*Nsq)*u_V_DFS + sqrt(5.)*(-1 + Nsq)*u_V_DFD)))/ &
    4.
              rVd(3) = (root_3*M*N*((6 - 20*Nsq)*u_V_DFS +  &
        sqrt(2.)*(-7 + 20*Nsq)*u_V_DFP + 2*sqrt(5.)*(1 - 2*Nsq)*u_V_DFD) &
 )/2.
              rVd_t = -(root_3*M*(sqrt(2.)*(-1 + 7*Nsq - 10*N**4)*u_Vd_DFP +  &
         2*Nsq*((-3 + 5*Nsq)*u_Vd_DFS + sqrt(5.)*(-1 + Nsq)*u_Vd_DFD)))/ &
    4.
              case (3)
              rVd(1) = 0
              rVd(2) = 0
              rVd(3) = ((-3 + 42*Nsq - 75*N**4)*u_V_DFS +  &
      3*sqrt(2.)*(1 - 18*Nsq + 25*N**4)*u_V_DFP -  &
      3*sqrt(5.)*(1 - 6*Nsq + 5*N**4)*u_V_DFD)/4.
              rVd_t = -(N*((3 - 14*Nsq + 15*N**4)*u_Vd_DFS -  &
         3*sqrt(2.)*(1 - 6*Nsq + 5*N**4)*u_Vd_DFP +  &
         3*sqrt(5.)*(-1 + Nsq)**2*u_Vd_DFD))/4.
              case (4)
              rVd(1) = -(root_3*(sqrt(2.)*(-1 + 7*Nsq - 10*N**4)*u_V_DFP +  &
         2*Nsq*((-3 + 5*Nsq)*u_V_DFS + sqrt(5.)*(-1 + Nsq)*u_V_DFD)))/ &
    4.
              rVd(2) = 0
              rVd(3) = (root_3*L*N*((6 - 20*Nsq)*u_V_DFS +  &
        sqrt(2.)*(-7 + 20*Nsq)*u_V_DFP + 2*sqrt(5.)*(1 - 2*Nsq)*u_V_DFD) &
 )/2.
              rVd_t = -(root_3*L*(sqrt(2.)*(-1 + 7*Nsq - 10*N**4)*u_Vd_DFP +  &
         2*Nsq*((-3 + 5*Nsq)*u_Vd_DFS + sqrt(5.)*(-1 + Nsq)*u_Vd_DFD)))/ &
    4.
              case (5)
              rVd(1) = -(root_3*L*N*((-3 + 5*Nsq)*u_V_DFS +  &
         sqrt(2.)*(1 - 5*Nsq)*u_V_DFP + sqrt(5.)*(1 + Nsq)*u_V_DFD))/2.
              rVd(2) = (root_3*M*N*((-3 + 5*Nsq)*u_V_DFS +  &
        sqrt(2.)*(1 - 5*Nsq)*u_V_DFP + sqrt(5.)*(1 + Nsq)*u_V_DFD))/2.
              rVd(3) = (root_3*(-1 + 2*Msq + Nsq)* &
      (3*(-1 + 5*Nsq)*u_V_DFS + sqrt(2.)*(1 - 15*Nsq)*u_V_DFP +  &
        sqrt(5.)*(1 + 3*Nsq)*u_V_DFD))/4.
              rVd_t = -(root_3*(L - M)*(L + M)*N* &
       ((-3 + 5*Nsq)*u_Vd_DFS + sqrt(2.)*(1 - 5*Nsq)*u_Vd_DFP +  &
         sqrt(5.)*(1 + Nsq)*u_Vd_DFD))/4.
            end select
          case (5)
            select case (orb_dir_j)
              case (1)
              rVd(1) = (L*M*(-3*sqrt(2.)*(1 - 6*Nsq + 5*N**4)*u_V_DFS +  &
        (1 - 27*Nsq + 30*N**4)*u_V_DFP + sqrt(10.)*(1 - 3*N**4)*u_V_DFD) &
 )/(2.*(-1 + Nsq))
              rVd(2) = (Msq*(3*(-1 + 5*Nsq)*u_V_DFP -  &
         6*sqrt(10.)*Nsq*u_V_DFD) +  &
      Lsq*(-3*sqrt(2.)*(1 - 6*Nsq + 5*N**4)*u_V_DFS +  &
         (1 - 27*Nsq + 30*N**4)*u_V_DFP + sqrt(10.)*(1 - 3*N**4)*u_V_DFD &
 ))/(4.*(-1 + Nsq))
              rVd(3) = (M*N*(15*sqrt(2.)*(-1 + Nsq)*(-1 + Msq + Nsq)*u_V_DFS +  &
        (-26 + 60*Nsq - 30*N**4 - 30*Msq*(-1 + Nsq))*u_V_DFP +  &
        sqrt(10.)*(1 - 6*Nsq + 3*N**4 + 3*Msq*(-1 + Nsq))*u_V_DFD))/ &
    (2.*(-1 + Nsq))
              rVd_t = -(M*(Msq*((1 - 5*Nsq)*u_Vd_DFP + 2*sqrt(10.)*Nsq*u_Vd_DFD) +  &
         Lsq*(3*sqrt(2.)*(1 - 6*Nsq + 5*N**4)*u_Vd_DFS +  &
            (-1 + 27*Nsq - 30*N**4)*u_Vd_DFP +  &
            sqrt(10.)*(-1 + 3*N**4)*u_Vd_DFD)))/(4.*(-1 + Nsq))
              case (2)
              rVd(1) = (-3*M*N*(sqrt(2.)*(-1 + 5*Nsq)*u_V_DFS +  &
        (4 - 10*Nsq)*u_V_DFP + sqrt(10.)*(-1 + Nsq)*u_V_DFD))/4.
              rVd(2) = (-3*L*N*(sqrt(2.)*(-1 + 5*Nsq)*u_V_DFS +  &
        (4 - 10*Nsq)*u_V_DFP + sqrt(10.)*(-1 + Nsq)*u_V_DFD))/4.
              rVd(3) = (3*L*M*(sqrt(2.)*(1 - 15*Nsq)*u_V_DFS +  &
        (-4 + 30*Nsq)*u_V_DFP + sqrt(10.)*(1 - 3*Nsq)*u_V_DFD))/4.
              rVd_t = (-3*L*M*N*(sqrt(2.)*(-1 + 5*Nsq)*u_Vd_DFS +  &
        (4 - 10*Nsq)*u_Vd_DFP + sqrt(10.)*(-1 + Nsq)*u_Vd_DFD))/4.
              case (3)
              rVd(1) = -(root_3*(sqrt(2.)*(1 - 8*Nsq + 15*N**4)*u_V_DFS +  &
         2*Nsq*(11 - 15*Nsq)*u_V_DFP +  &
         sqrt(10.)*(1 - 4*Nsq + 3*N**4)*u_V_DFD))/8.
              rVd(2) = 0
              rVd(3) = (root_3*L*N*(sqrt(2.)*(4 - 15*Nsq)*u_V_DFS +  &
        (-11 + 30*Nsq)*u_V_DFP + sqrt(10.)*(2 - 3*Nsq)*u_V_DFD))/2.
              rVd_t = -(root_3*L*(sqrt(2.)*(1 - 8*Nsq + 15*N**4)*u_Vd_DFS +  &
         2*Nsq*(11 - 15*Nsq)*u_Vd_DFP +  &
         sqrt(10.)*(1 - 4*Nsq + 3*N**4)*u_Vd_DFD))/8.
              case (4)
              rVd(1) = 0
              rVd(2) = (M*N*(3*sqrt(2.)*(-1 + 5*Nsq)*u_V_DFS -  &
        6*(-2 + 5*Nsq)*u_V_DFP + 3*sqrt(10.)*(-1 + Nsq)*u_V_DFD))/2.
              rVd(3) = (3*sqrt(2.)*(1 - 18*Nsq + 25*N**4 + Msq*(-1 + 15*Nsq))* &
       u_V_DFS + (-11 + 111*Nsq - 150*N**4 + Msq*(12 - 90*Nsq))* &
       u_V_DFP + sqrt(10.)*(1 - 12*Nsq + 15*N**4 + Msq*(-3 + 9*Nsq))* &
       u_V_DFD)/4.
              rVd_t = (N*(3*sqrt(2.)*(-1 + Msq + Nsq)*(-1 + 5*Nsq)*u_Vd_DFS +  &
        (-11 + 37*Nsq - 30*N**4 - 6*Msq*(-2 + 5*Nsq))*u_Vd_DFP +  &
        sqrt(10.)*(-1 + Nsq)*(-1 + 3*Msq + 3*Nsq)*u_Vd_DFD))/4.
              case (5)
              rVd(1) = (3*sqrt(2.)*(-3 + 4*Msq + 3*Nsq)*(-1 + 5*Nsq)*u_V_DFS -  &
      sqrt(10.)*u_V_DFD + Nsq*(22*u_V_DFP + 2*sqrt(10.)*u_V_DFD) +  &
      N**4*(-30*u_V_DFP + 3*sqrt(10.)*u_V_DFD) +  &
      Msq*((4 - 60*Nsq)*u_V_DFP + 2*sqrt(10.)*(1 + 3*Nsq)*u_V_DFD))/8.
              rVd(2) = (L*M*(3*sqrt(2.)*(-1 + 5*Nsq)*u_V_DFS +  &
        (4 - 60*Nsq)*u_V_DFP + 2*sqrt(10.)*(1 + 3*Nsq)*u_V_DFD))/4.
              rVd(3) = (L*N*(15*sqrt(2.)*(-1 + 2*Msq + Nsq)*u_V_DFS +  &
        (22 - 60*Msq - 60*Nsq)*u_V_DFP +  &
        2*sqrt(10.)*(1 + 3*Msq + 3*Nsq)*u_V_DFD))/4.
              rVd_t = -(L*(3*sqrt(2.)*(L - M)*(L + M)*(-1 + 5*Nsq)*u_Vd_DFS +  &
         2*(Nsq*(-11 + 15*Nsq) + Msq*(-2 + 30*Nsq))*u_Vd_DFP -  &
         sqrt(10.)*(-1 + 2*Nsq + 3*N**4 + Msq*(2 + 6*Nsq))*u_Vd_DFD))/8.
            end select
          case (6)
            select case (orb_dir_j)
              case (1)
              rVd(1) = (3*M*(-3*Lsq + Msq)*N* &
      (sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD))/2.
              rVd(2) = (3*L*N*(-1 + 4*Msq + Nsq)* &
      (sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD))/2.
              rVd(3) = (-3*L*(L - M)*M*(L + M)* &
      (sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD))/2.
              rVd_t = (-3*L*(L - M)*M*(L + M)*N* &
      (sqrt(5.)*u_Vd_DFS - sqrt(10.)*u_Vd_DFP + u_Vd_DFD))/2.
              case (2)
              rVd(1) = 0
              rVd(2) = (sqrt(10.)*(-1 + 6*Msq)*u_V_DFP +  &
      8*(1 - 3*Msq)*u_V_DFD +  &
      6*N**4*(sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD) +  &
      Nsq*(6*sqrt(5.)*(-1 + 6*Msq)*u_V_DFS +  &
         9*(sqrt(10.)*(1 - 4*Msq)*u_V_DFP + 2*(-1 + 2*Msq)*u_V_DFD)))/4.
              rVd(3) = (3*M*N*(2*sqrt(5.)*(-1 + 2*Msq + 2*Nsq)*u_V_DFS +  &
        sqrt(10.)*(3 - 4*Msq - 4*Nsq)*u_V_DFP +  &
        2*(-3 + 2*Msq + 2*Nsq)*u_V_DFD))/2.
              rVd_t = (M*(sqrt(10.)*(-1 + 2*Msq)*u_Vd_DFP - 8*(-1 + Msq)*u_Vd_DFD +  &
        6*N**4*(sqrt(5.)*u_Vd_DFS - sqrt(10.)*u_Vd_DFP + u_Vd_DFD) +  &
        3*Nsq*(2*sqrt(5.)*(-1 + 2*Msq)*u_Vd_DFS +  &
           sqrt(10.)*(3 - 4*Msq)*u_Vd_DFP + 2*(-3 + 2*Msq)*u_Vd_DFD)))/4.
              case (3)
              rVd(1) = -(root_3*L*N*(-1 + 3*Nsq)* &
       (sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD))/2.
              rVd(2) = (root_3*M*N*(-1 + 3*Nsq)* &
      (sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD))/2.
              rVd(3) = -(root_3*(L - M)*(L + M)*(-1 + 9*Nsq)* &
       (sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD))/4.
              rVd_t = -(root_3*(L - M)*(L + M)*N*(-1 + 3*Nsq)* &
       (sqrt(5.)*u_Vd_DFS - sqrt(10.)*u_Vd_DFP + u_Vd_DFD))/4.
              case (4)
              rVd(1) = (sqrt(10.)*(-1 + 2*Msq)*u_V_DFP - 8*Msq*u_V_DFD +  &
      6*N**4*(sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD) +  &
      Nsq*(6*sqrt(5.)*(-1 + 2*Msq)*u_V_DFS +  &
         sqrt(10.)*(5 - 12*Msq)*u_V_DFP + 2*(-1 + 6*Msq)*u_V_DFD))/4.
              rVd(2) = L*M*(sqrt(10.)*u_V_DFP - 4*u_V_DFD +  &
      6*Nsq*(sqrt(5.)*u_V_DFS - sqrt(10.)*u_V_DFP + u_V_DFD))
              rVd(3) = (L*N*(6*sqrt(5.)*(-1 + 2*Msq + 2*Nsq)*u_V_DFS +  &
        sqrt(10.)*(5 - 12*Msq - 12*Nsq)*u_V_DFP +  &
        2*(-1 + 6*Msq + 6*Nsq)*u_V_DFD))/2.
              rVd_t = (L*(sqrt(10.)*(-1 + 2*Msq)*u_Vd_DFP - 8*Msq*u_Vd_DFD +  &
        6*N**4*(sqrt(5.)*u_Vd_DFS - sqrt(10.)*u_Vd_DFP + u_Vd_DFD) +  &
        Nsq*(6*sqrt(5.)*(-1 + 2*Msq)*u_Vd_DFS +  &
           sqrt(10.)*(5 - 12*Msq)*u_Vd_DFP + 2*(-1 + 6*Msq)*u_Vd_DFD)))/4.
              case (5)
              rVd(1) = 3*sqrt(5.)*L*N*(-1 + 2*Msq + Nsq)*u_V_DFS
              rVd(2) = -3*M*N*(-1 + 2*Msq + Nsq)* &
    (sqrt(5.)*u_V_DFS - 2*sqrt(10.)*u_V_DFP + 2*u_V_DFD)
              rVd(3) = (-3*sqrt(5.)*(Lsq - Msq)**2*u_V_DFS +  &
      sqrt(10.)*(1 + 12*M**4 - 4*Nsq + 3*N**4 + 12*Msq*(-1 + Nsq))* &
       u_V_DFP - (-1 + 12*M**4 + 2*Nsq + 3*N**4 + 12*Msq*(-1 + Nsq))* &
       u_V_DFD - N*(-4*sqrt(10.)*N*(-2 + 6*Msq + 3*Nsq)*u_V_DFP +  &
         4*(N + 6*Msq*N + 3*N**3)*u_V_DFD))/4.
              rVd_t = -(N*(3*sqrt(5.)*(Lsq - Msq)**2*u_Vd_DFS -  &
         sqrt(10.)*(1 + 12*M**4 - 4*Nsq + 3*N**4 + 12*Msq*(-1 + Nsq))* &
          u_Vd_DFP + (-1 + 12*M**4 + 2*Nsq + 3*N**4 + 12*Msq*(-1 + Nsq))* &
          u_Vd_DFD))/4.
            end select
          case (7)
            select case (orb_dir_j)
              case (1)
              rVd(1) = Sqrt(7.5)*L*M*(-2*Lsq + 3*Msq)*u_V_DFS
              rVd(2) = (root_3*(-(sqrt(10.)*L**4*u_V_DFS) +  &
        9*sqrt(10.)*Lsq*Msq*u_V_DFS +  &
        sqrt(5.)*(3 + 40*M**4 - 5*Nsq + 2*N**4 + 30*Msq*(-1 + Nsq))* &
         u_V_DFP - sqrt(2.)*(1 + 20*M**4 - 4*Nsq + N**4 +  &
           15*Msq*(-1 + Nsq))*u_V_DFD))/4.
              rVd(3) = (M*N*(sqrt(15.)*(-5 + 10*Msq + 4*Nsq)*u_V_DFP +  &
        sqrt(6.)*(4 - 5*Msq - 2*Nsq)*u_V_DFD))/2.
              rVd_t = -(root_3*M*(sqrt(10.)*Lsq*(Lsq - 3*Msq)*u_Vd_DFS -  &
         sqrt(5.)*(3 + 8*M**4 - 5*Nsq + 2*N**4 + 10*Msq*(-1 + Nsq))* &
          u_Vd_DFP + sqrt(2.)*(1 + 4*M**4 - 4*Nsq + N**4 +  &
            5*Msq*(-1 + Nsq))*u_Vd_DFD))/4.
              case (2)
              rVd(1) = (M*N*(sqrt(30.)*(-1 + 4*Msq + Nsq)*u_V_DFS -  &
        2*sqrt(15.)*(-2 + 4*Msq + Nsq)*u_V_DFP +  &
        sqrt(6.)*(-5 + 4*Msq + Nsq)*u_V_DFD))/4.
              rVd(2) = (L*N*(sqrt(30.)*(-1 + 12*Msq + Nsq)*u_V_DFS -  &
        2*sqrt(15.)*(-2 + 12*Msq + Nsq)*u_V_DFP +  &
        sqrt(6.)*(-5 + 12*Msq + Nsq)*u_V_DFD))/4.
              rVd(3) = (L*M*(sqrt(30.)*(-1 + 4*Msq + 3*Nsq)*u_V_DFS -  &
        2*sqrt(15.)*(-2 + 4*Msq + 3*Nsq)*u_V_DFP +  &
        sqrt(6.)*(-5 + 4*Msq + 3*Nsq)*u_V_DFD))/4.
              rVd_t = (L*M*N*(sqrt(30.)*(-1 + 4*Msq + Nsq)*u_Vd_DFS -  &
        2*sqrt(15.)*(-2 + 4*Msq + Nsq)*u_Vd_DFP +  &
        sqrt(6.)*(-5 + 4*Msq + Nsq)*u_Vd_DFD))/4.
              case (3)
              rVd(1) = (-3*(Lsq - Msq)* &
      (sqrt(10.)*(-1 + 3*Nsq)*u_V_DFS - 6*sqrt(5.)*Nsq*u_V_DFP +  &
        3*sqrt(2.)*(1 + Nsq)*u_V_DFD))/8.
              rVd(2) = (3*L*M*(sqrt(10.)*(-1 + 3*Nsq)*u_V_DFS -  &
        6*sqrt(5.)*Nsq*u_V_DFP + 3*sqrt(2.)*(1 + Nsq)*u_V_DFD))/4.
              rVd(3) = (3*L*N*(-1 + 4*Msq + Nsq)* &
      (sqrt(10.)*u_V_DFS - 2*sqrt(5.)*u_V_DFP + sqrt(2.)*u_V_DFD))/4.
              rVd_t = -((L**3 - 3*L*Msq)*(sqrt(10.)*(-1 + 3*Nsq)*u_Vd_DFS -  &
         6*sqrt(5.)*Nsq*u_Vd_DFP + 3*sqrt(2.)*(1 + Nsq)*u_Vd_DFD))/8.
              case (4)
              rVd(1) = 0
              rVd(2) = (M*N*(sqrt(30.)*(5 - 8*Msq - 5*Nsq)*u_V_DFS +  &
        2*sqrt(15.)*(-4 + 8*Msq + 5*Nsq)*u_V_DFP +  &
        sqrt(6.)*(1 - 8*Msq - 5*Nsq)*u_V_DFD))/2.
              rVd(3) = (-(sqrt(30.)*(1 + 4*M**4 - 6*Nsq + 5*N**4 +  &
           5*Msq*(-1 + 3*Nsq))*u_V_DFS) +  &
      sqrt(15.)*(1 + 8*M**4 - 9*Nsq + 10*N**4 + Msq*(-8 + 30*Nsq))* &
       u_V_DFP + sqrt(6.)*(1 - 4*M**4 - 5*N**4 + Msq*(1 - 15*Nsq))* &
       u_V_DFD)/4.
              rVd_t = -(root_3*N*(sqrt(10.)* &
          (4*M**4 + 5*Msq*(-1 + Nsq) + (-1 + Nsq)**2)*u_Vd_DFS -  &
         sqrt(5.)*(1 + 8*M**4 - 3*Nsq + 2*N**4 + 2*Msq*(-4 + 5*Nsq))* &
          u_Vd_DFP + sqrt(2.)*(-1 + 4*M**4 + N**4 + Msq*(-1 + 5*Nsq))* &
          u_Vd_DFD))/4.
              case (5)
              rVd(1) = (root_3*(-4*sqrt(10.)*L*(L**3 - 2*L*Msq)*u_V_DFS -  &
        sqrt(10.)*(L**4 - 4*Lsq*Msq + 3*M**4)*u_V_DFS +  &
        2*sqrt(5.)*(8*M**4 + 6*Msq*(-1 + Nsq) + Nsq*(-1 + Nsq))* &
         u_V_DFP - sqrt(2.)*(8*M**4 + 6*Msq*(-1 + Nsq) + (1 + Nsq)**2)* &
         u_V_DFD))/8.
              rVd(2) = (L*M*(sqrt(30.)*(2 - 5*Msq - 2*Nsq)*u_V_DFS +  &
        (-3 + 8*Msq + 3*Nsq)*(2*sqrt(15.)*u_V_DFP - sqrt(6.)*u_V_DFD)))/ &
    2.
              rVd(3) = (L*N*(sqrt(15.)*(-1 + 6*Msq + 2*Nsq)*u_V_DFP -  &
        sqrt(6.)*(1 + 3*Msq + Nsq)*u_V_DFD))/2.
              rVd_t = -(root_3*L*(sqrt(10.)*(L**4 - 4*Lsq*Msq + 3*M**4)*u_Vd_DFS -  &
         2*sqrt(5.)*(8*M**4 + 6*Msq*(-1 + Nsq) + Nsq*(-1 + Nsq))* &
          u_Vd_DFP + sqrt(2.)*(8*M**4 + 6*Msq*(-1 + Nsq) + (1 + Nsq)**2)* &
          u_Vd_DFD))/8.
            end select
        end select
      case (ORB_F)
        select case (orb_dir_i)
          case (1)
            select case (orb_dir_j)
              case (1)
              rVd(1) = (-15*L*(8*M**4*(2*u_V_FFS - 3*u_V_FFP) +  &
        6*Msq*(-1 + Nsq)*(2*u_V_FFS - 3*u_V_FFP) +  &
        3*(-1 + Nsq)*u_V_FFP))/8.
              rVd(2) = (3*M*(10*(8*M**4 + 10*Msq*(-1 + Nsq) + 3*(-1 + Nsq)**2)* &
         u_V_FFS - 15*(2 + 8*M**4 - 5*Nsq + 3*N**4 +  &
           10*Msq*(-1 + Nsq))*u_V_FFP +  &
        (16*M**4 + 16*Msq*(-1 + Nsq) + 3*(-1 + Nsq)**2)* &
         (6*u_V_FFD - u_V_FFF)))/8.
              rVd(3) = (-3*N*(5*(-1 + Nsq)*u_V_FFP -  &
        4*(2 + 12*M**4 - 4*Nsq + 9*Msq*(-1 + Nsq))*u_V_FFD +  &
        2*(-1 + Msq)*(1 + 4*Msq + 3*Nsq)*u_V_FFF))/8.
              rVd_t = (10*(-3*Lsq*M + M**3)**2*u_Vd_FFS -  &
      (15*((L**3 - 3*L*Msq)**2 + Msq*(-3*Lsq + Msq)**2*Nsq)*u_Vd_FFP)/ &
       (-1 + Nsq) + 6*(16*M**6 + 24*M**4*(-1 + Nsq) -  &
         4*Nsq*(-1 + Nsq) + 9*Msq*(-1 + Nsq)**2)*u_Vd_FFD -  &
      (-1 + Msq)*(16*M**4 + 8*Msq*(-1 + 3*Nsq) + (1 + 3*Nsq)**2)* &
       u_Vd_FFF)/16.
              case (2)
              rVd(1) = -(Sqrt(1.5)*N*(5*(-1 + Nsq)*u_V_FFP + 4*u_V_FFD +  &
         3*Msq*(-1 + Nsq)*(10*u_V_FFS - 15*u_V_FFP + 6*u_V_FFD -  &
            u_V_FFF) + u_V_FFF + Nsq*(-8*u_V_FFD + 3*u_V_FFF) +  &
         M**4*(40*u_V_FFS - 4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF))) &
 )/4.
              rVd(2) = -(Sqrt(1.5)*L*M*N*(-3 + 8*Msq + 3*Nsq)* &
       (10*u_V_FFS - 15*u_V_FFP + 6*u_V_FFD - u_V_FFF))/2.
              rVd(3) = -(Sqrt(1.5)*L*(-5*u_V_FFP + 4*u_V_FFD +  &
         3*Msq*(-1 + 3*Nsq)*(10*u_V_FFS - 15*u_V_FFP +  &
            6*u_V_FFD - u_V_FFF) + u_V_FFF +  &
         3*Nsq*(5*u_V_FFP - 8*u_V_FFD + 3*u_V_FFF) +  &
         M**4*(40*u_V_FFS - 4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF))) &
 )/4.
              rVd_t = -(Sqrt(1.5)*L*N*(5*(-1 + Nsq)*u_Vd_FFP + 4*u_Vd_FFD +  &
         3*Msq*(-1 + Nsq)*(10*u_Vd_FFS - 15*u_Vd_FFP + 6*u_Vd_FFD -  &
            u_Vd_FFF) + u_Vd_FFF + Nsq*(-8*u_Vd_FFD + 3*u_Vd_FFF) +  &
         M**4*(40*u_Vd_FFS - 4*(15*u_Vd_FFP - 6*u_Vd_FFD + u_Vd_FFF))) &
 )/4.
              case (3)
              rVd(1) = (sqrt(15.)*(4*L**3* &
         (u_V_FFP + Nsq*(-5*u_V_FFP + 8*u_V_FFD - 3*u_V_FFF) -  &
           u_V_FFF) + 6*L*Msq*(-1 + Nsq)* &
         (2*(-1 + 5*Nsq)*u_V_FFS + (1 - 15*Nsq)*u_V_FFP +  &
           2*(1 + 3*Nsq)*u_V_FFD - (1 + Nsq)*u_V_FFF)))/ &
    (16.*(-1 + Nsq))
              rVd(2) = (sqrt(15.)*(6*Lsq*M*(-1 + Nsq)* &
         (2*(-1 + 5*Nsq)*u_V_FFS + (1 - 15*Nsq)*u_V_FFP +  &
           2*(1 + 3*Nsq)*u_V_FFD - (1 + Nsq)*u_V_FFF) +  &
        4*M**3*(-2*(1 - 6*Nsq + 5*N**4)*u_V_FFS + 2*u_V_FFD +  &
           Nsq*((-11 + 15*Nsq)*u_V_FFP - 2*(2 + 3*Nsq)*u_V_FFD +  &
              (3 + Nsq)*u_V_FFF))))/(16.*(-1 + Nsq))
              rVd(3) = (sqrt(15.)*N*(-4*M**4*(-1 + Nsq)* &
         (10*u_V_FFS - 15*u_V_FFP + 6*u_V_FFD - u_V_FFF) +  &
        4*(-1 + Nsq)*(u_V_FFP - 2*u_V_FFD + u_V_FFF) +  &
        Msq*(-30*(-1 + Nsq)**2*u_V_FFS + 53*u_V_FFP - 34*u_V_FFD +  &
           11*u_V_FFF - 6*Nsq* &
            (15*u_V_FFP - 6*u_V_FFD + u_V_FFF) +  &
           3*N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF))))/ &
    (8.*(-1 + Nsq))
              rVd_t = (sqrt(15.)*(L**4*(u_Vd_FFP +  &
           Nsq*(-5*u_Vd_FFP + 8*u_Vd_FFD - 3*u_Vd_FFF) - u_Vd_FFF) +  &
        3*Lsq*Msq*(-1 + Nsq)* &
         (2*(-1 + 5*Nsq)*u_Vd_FFS + (1 - 15*Nsq)*u_Vd_FFP +  &
           2*(1 + 3*Nsq)*u_Vd_FFD - (1 + Nsq)*u_Vd_FFF) +  &
        M**4*(-2*(1 - 6*Nsq + 5*N**4)*u_Vd_FFS + 2*u_Vd_FFD +  &
           Nsq*((-11 + 15*Nsq)*u_Vd_FFP - 2*(2 + 3*Nsq)*u_Vd_FFD +  &
              (3 + Nsq)*u_Vd_FFF))))/(16.*(-1 + Nsq))
              case (4)
              rVd(1) = (3*Sqrt(2.5)*L*M*N* &
      (2*(-3 + 5*Nsq)*u_V_FFS + (3 - 15*Nsq)*u_V_FFP +  &
        6*(1 + Nsq)*u_V_FFD - (3 + Nsq)*u_V_FFF))/4.
              rVd(2) = (-3*Sqrt(2.5)*N*(-1 + 2*Msq + Nsq)* &
      (2*(-3 + 5*Nsq)*u_V_FFS + 3*u_V_FFP + 6*u_V_FFD -  &
        3*u_V_FFF - Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF)))/8.
              rVd(3) = (-3*Sqrt(2.5)*M*(-3 + 4*Msq + 3*Nsq)* &
      (2*(-1 + 5*Nsq)*u_V_FFS + u_V_FFP + 2*u_V_FFD - u_V_FFF -  &
        Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF)))/8.
              rVd_t = -(Sqrt(2.5)*M*(-3*Lsq + Msq)*N* &
       (2*(-3 + 5*Nsq)*u_Vd_FFS + (3 - 15*Nsq)*u_Vd_FFP +  &
         6*(1 + Nsq)*u_Vd_FFD - (3 + Nsq)*u_Vd_FFF))/8.
              case (5)
              rVd(1) = -(sqrt(15.)*M*(2*(-3 + 4*Msq + 3*Nsq)*(-1 + 5*Nsq)* &
          u_V_FFS - u_V_FFP - 6*u_V_FFD +  &
         Nsq*(38*u_V_FFP + 4*u_V_FFD - 6*u_V_FFF) +  &
         4*Msq*(u_V_FFP + 2*u_V_FFD - u_V_FFF) + u_V_FFF -  &
         4*Msq*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
         3*N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF)))/16.
              rVd(2) = -(sqrt(15.)*L*(6*(-1 + 4*Msq + Nsq)*(-1 + 5*Nsq)*u_V_FFS -  &
         u_V_FFP - 6*u_V_FFD +  &
         Nsq*(38*u_V_FFP + 4*u_V_FFD - 6*u_V_FFF) +  &
         12*Msq*(u_V_FFP + 2*u_V_FFD - u_V_FFF) + u_V_FFF -  &
         12*Msq*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
         3*N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF)))/16.
              rVd(3) = -(sqrt(15.)*L*M*N*(2*(-9 + 10*Msq + 15*Nsq)*u_V_FFS +  &
         19*u_V_FFP + 2*u_V_FFD - 3*u_V_FFF -  &
         2*Msq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
         3*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF)))/4.
              rVd_t = -(sqrt(15.)*L*M*(2*(-3 + 4*Msq + 3*Nsq)*(-1 + 5*Nsq)* &
          u_Vd_FFS - u_Vd_FFP - 6*u_Vd_FFD +  &
         Nsq*(38*u_Vd_FFP + 4*u_Vd_FFD - 6*u_Vd_FFF) +  &
         4*Msq*(u_Vd_FFP + 2*u_Vd_FFD - u_Vd_FFF) + u_Vd_FFF -  &
         4*Msq*Nsq*(15*u_Vd_FFP - 6*u_Vd_FFD + u_Vd_FFF) -  &
         3*N**4*(15*u_Vd_FFP - 6*u_Vd_FFD + u_Vd_FFF)))/16.
              case (6)
              rVd(1) = 0
              rVd(2) = (Sqrt(1.5)*N*(10*(40*M**4 + 30*Msq*(-1 + Nsq) +  &
           3*(-1 + Nsq)**2)*u_V_FFS - 35*u_V_FFP +  &
        Nsq*(80*u_V_FFP - 20*u_V_FFD) + 10*u_V_FFD - 5*u_V_FFF +  &
        30*Msq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        40*M**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        30*Msq*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        3*N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF)))/8.
              rVd(3) = (Sqrt(1.5)*M*(10*(8*M**4 + 10*Msq*(-1 + 3*Nsq) +  &
           3*(1 - 6*Nsq + 5*N**4))*u_V_FFS - 35*u_V_FFP +  &
        60*Nsq*(4*u_V_FFP - u_V_FFD) + 10*u_V_FFD - 5*u_V_FFF +  &
        10*Msq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        8*M**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        30*Msq*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        15*N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF)))/8.
              rVd_t = (Sqrt(1.5)*M*N*(10*(8*M**4 + 10*Msq*(-1 + Nsq) +  &
           3*(-1 + Nsq)**2)*u_Vd_FFS - 35*u_Vd_FFP +  &
        Nsq*(80*u_Vd_FFP - 20*u_Vd_FFD) + 10*u_Vd_FFD - 5*u_Vd_FFF +  &
        10*Msq*(15*u_Vd_FFP - 6*u_Vd_FFD + u_Vd_FFF) -  &
        8*M**4*(15*u_Vd_FFP - 6*u_Vd_FFD + u_Vd_FFF) -  &
        10*Msq*Nsq*(15*u_Vd_FFP - 6*u_Vd_FFD + u_Vd_FFF) -  &
        3*N**4*(15*u_Vd_FFP - 6*u_Vd_FFD + u_Vd_FFF)))/8.
              case (7)
              rVd(1) = (3*M*(5*L**4 - 10*Lsq*Msq + M**4)* &
      (10*u_V_FFS - 15*u_V_FFP + 6*u_V_FFD - u_V_FFF))/16.
              rVd(2) = (3*L*(L**4 - 10*Lsq*Msq + 5*M**4)* &
      (10*u_V_FFS - 15*u_V_FFP + 6*u_V_FFD - u_V_FFF))/16.
              rVd(3) = 0
              rVd_t = (L*M*(3*L**4 - 10*Lsq*Msq + 3*M**4)* &
      (10*u_Vd_FFS - 15*u_Vd_FFP + 6*u_Vd_FFD - u_Vd_FFF))/16.
            end select
          case (2)
            select case (orb_dir_j)
              case (1)
              rVd(1) = -(Sqrt(1.5)*N*(5*(-1 + Nsq)*u_V_FFP + 4*u_V_FFD +  &
         3*Msq*(-1 + Nsq)*(10*u_V_FFS - 15*u_V_FFP + 6*u_V_FFD -  &
            u_V_FFF) + u_V_FFF + Nsq*(-8*u_V_FFD + 3*u_V_FFF) +  &
         M**4*(40*u_V_FFS - 4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF))) &
 )/4.
              rVd(2) = -(Sqrt(1.5)*L*M*N*(-3 + 8*Msq + 3*Nsq)* &
       (10*u_V_FFS - 15*u_V_FFP + 6*u_V_FFD - u_V_FFF))/2.
              rVd(3) = -(Sqrt(1.5)*L*(-5*u_V_FFP + 4*u_V_FFD +  &
         3*Msq*(-1 + 3*Nsq)*(10*u_V_FFS - 15*u_V_FFP +  &
            6*u_V_FFD - u_V_FFF) + u_V_FFF +  &
         3*Nsq*(5*u_V_FFP - 8*u_V_FFD + 3*u_V_FFF) +  &
         M**4*(40*u_V_FFS - 4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF))) &
 )/4.
              rVd_t = -(Sqrt(1.5)*L*N*(5*(-1 + Nsq)*u_Vd_FFP + 4*u_Vd_FFD +  &
         3*Msq*(-1 + Nsq)*(10*u_Vd_FFS - 15*u_Vd_FFP + 6*u_Vd_FFD -  &
            u_Vd_FFF) + u_Vd_FFF + Nsq*(-8*u_Vd_FFD + 3*u_Vd_FFF) +  &
         M**4*(40*u_Vd_FFS - 4*(15*u_Vd_FFP - 6*u_Vd_FFD + u_Vd_FFF))) &
 )/4.
              case (2)
              rVd(1) = (L*(-4*u_V_FFD -  &
        2*Nsq*(5*u_V_FFP - 8*u_V_FFD + 3*u_V_FFF) +  &
        2*N**4*(5*u_V_FFP - 8*u_V_FFD + 3*u_V_FFF) +  &
        Msq*(-1 + Nsq)*(5*u_V_FFP - 8*u_V_FFD +  &
           Nsq*(30*u_V_FFS - 45*u_V_FFP + 18*u_V_FFD -  &
              3*u_V_FFF) + 3*u_V_FFF)))/(-1 + Nsq)
              rVd(2) = (4*M**3*(2*u_V_FFD +  &
         Nsq*(-1 + Nsq)*(-5*u_V_FFP + 8*u_V_FFD - 3*u_V_FFF)) +  &
      2*Lsq*M*(5*u_V_FFP - 4*u_V_FFD + 3*u_V_FFF +  &
         Nsq*(30*(-1 + Nsq)**2*u_V_FFS -  &
            5*(9 - 17*Nsq + 9*N**4)*u_V_FFP +  &
            2*(9 - 14*Nsq + 9*N**4)*u_V_FFD -  &
            3*(1 - Nsq + N**4)*u_V_FFF)))/(2.*(-1 + Nsq)**2)
              rVd(3) = (N*(5*(-1 + Nsq)*u_V_FFP + 4*u_V_FFD -  &
        3*M**4*(-1 + Nsq)*(10*u_V_FFS - 15*u_V_FFP + 6*u_V_FFD -  &
           u_V_FFF) - 3*Msq*(-1 + Nsq)**2* &
         (10*u_V_FFS - 15*u_V_FFP + 6*u_V_FFD - u_V_FFF) -  &
        3*u_V_FFF + Nsq*(-8*u_V_FFD + 3*u_V_FFF)))/(-1 + Nsq)
              rVd_t = (L**4*(2*u_Vd_FFD +  &
         Nsq*(-1 + Nsq)*(-5*u_Vd_FFP + 8*u_Vd_FFD - 3*u_Vd_FFF)) +  &
      M**4*(2*u_Vd_FFD + Nsq*(-1 + Nsq)* &
          (-5*u_Vd_FFP + 8*u_Vd_FFD - 3*u_Vd_FFF)) +  &
      Lsq*Msq*(5*u_Vd_FFP - 4*u_Vd_FFD + 3*u_Vd_FFF +  &
         Nsq*(30*(-1 + Nsq)**2*u_Vd_FFS -  &
            5*(9 - 17*Nsq + 9*N**4)*u_Vd_FFP +  &
            2*(9 - 14*Nsq + 9*N**4)*u_Vd_FFD -  &
            3*(1 - Nsq + N**4)*u_Vd_FFF)))/(2.*(-1 + Nsq)**2)
              case (3)
              rVd(1) = (Sqrt(2.5)*N*((-1 + 5*Nsq)*u_V_FFP + 4*u_V_FFD -  &
        3*u_V_FFF + Nsq*(-8*u_V_FFD + 3*u_V_FFF) +  &
        Msq*(6*(-1 + 5*Nsq)*u_V_FFS + 13*u_V_FFP - 10*u_V_FFD +  &
           3*u_V_FFF - 3*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF))) &
 )/4.
              rVd(2) = (Sqrt(2.5)*L*M*N*(6*(-1 + 5*Nsq)*u_V_FFS + 13*u_V_FFP -  &
        10*u_V_FFD + 3*u_V_FFF -  &
        3*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF)))/2.
              rVd(3) = (Sqrt(2.5)*L*(-u_V_FFP + 4*u_V_FFD - 3*u_V_FFF +  &
        3*Nsq*(5*u_V_FFP - 8*u_V_FFD + 3*u_V_FFF) +  &
        Msq*((-6 + 90*Nsq)*u_V_FFS + 13*u_V_FFP - 10*u_V_FFD +  &
           3*u_V_FFF - 9*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF))) &
 )/4.
              rVd_t = (Sqrt(2.5)*L*N*((-1 + 5*Nsq)*u_Vd_FFP + 4*u_Vd_FFD -  &
        3*u_Vd_FFF + Nsq*(-8*u_Vd_FFD + 3*u_Vd_FFF) +  &
        Msq*(6*(-1 + 5*Nsq)*u_Vd_FFS + 13*u_Vd_FFP - 10*u_Vd_FFD +  &
           3*u_Vd_FFF - 3*Nsq*(15*u_Vd_FFP - 6*u_Vd_FFD + u_Vd_FFF))) &
 )/4.
              case (4)
              rVd(1) = (sqrt(15.)*M*(-u_V_FFP -  &
        2*Nsq*(3*u_V_FFS - 4*u_V_FFP + u_V_FFD) +  &
        N**4*(10*u_V_FFS - 15*u_V_FFP + 6*u_V_FFD - u_V_FFF) +  &
        u_V_FFF))/4.
              rVd(2) = (sqrt(15.)*L*(-u_V_FFP -  &
        2*Nsq*(3*u_V_FFS - 4*u_V_FFP + u_V_FFD) +  &
        N**4*(10*u_V_FFS - 15*u_V_FFP + 6*u_V_FFD - u_V_FFF) +  &
        u_V_FFF))/4.
              rVd(3) = sqrt(15.)*L*M*N*((-3 + 10*Nsq)*u_V_FFS + 4*u_V_FFP -  &
      u_V_FFD - Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF))
              rVd_t = (sqrt(15.)*L*M*(-u_Vd_FFP -  &
        2*Nsq*(3*u_Vd_FFS - 4*u_Vd_FFP + u_Vd_FFD) +  &
        N**4*(10*u_Vd_FFS - 15*u_Vd_FFP + 6*u_Vd_FFD - u_Vd_FFF) +  &
        u_Vd_FFF))/4.
              case (5)
              rVd(1) = (Sqrt(2.5)*L*M*N*(6*(1 - 6*Nsq + 5*N**4)*u_V_FFS +  &
        (-12 + 53*Nsq - 45*N**4)*u_V_FFP +  &
        2*(3 - 10*Nsq + 9*N**4)*u_V_FFD - 3*Nsq*(-1 + Nsq)*u_V_FFF))/ &
    (2.*(-1 + Nsq))
              rVd(2) = (Sqrt(2.5)*N*(3*Msq* &
         ((1 - 5*Nsq)*u_V_FFP + (-4 + 8*Nsq)*u_V_FFD -  &
           3*(-1 + Nsq)*u_V_FFF) -  &
        (-1 + Msq + Nsq)*(6*(1 - 6*Nsq + 5*N**4)*u_V_FFS -  &
           12*u_V_FFP + 6*u_V_FFD -  &
           3*N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) +  &
           Nsq*(53*u_V_FFP - 20*u_V_FFD + 3*u_V_FFF))))/ &
    (4.*(-1 + Nsq))
              rVd(3) = (Sqrt(2.5)*M*(-6*(-1 + Msq + Nsq)*(1 - 16*Nsq + 15*N**4)* &
         u_V_FFS - 12*u_V_FFP + 6*u_V_FFD +  &
        9*Msq*N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) +  &
        9*N**6*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        4*Msq*Nsq*(37*u_V_FFP - 16*u_V_FFD + 3*u_V_FFF) +  &
        Msq*(13*u_V_FFP - 10*u_V_FFD + 3*u_V_FFF) -  &
        2*N**4*(139*u_V_FFP - 55*u_V_FFD + 9*u_V_FFF) +  &
        Nsq*(147*u_V_FFP + 9*(-6*u_V_FFD + u_V_FFF))))/ &
    (4.*(-1 + Nsq))
              rVd_t = (Sqrt(2.5)*M*N*(Msq* &
         ((1 - 5*Nsq)*u_Vd_FFP + (-4 + 8*Nsq)*u_Vd_FFD -  &
           3*(-1 + Nsq)*u_Vd_FFF) +  &
        Lsq*(6*(1 - 6*Nsq + 5*N**4)*u_Vd_FFS +  &
           (-12 + 53*Nsq - 45*N**4)*u_Vd_FFP +  &
           2*(3 - 10*Nsq + 9*N**4)*u_Vd_FFD - 3*Nsq*(-1 + Nsq)*u_Vd_FFF &
 )))/(4.*(-1 + Nsq))
              case (6)
              rVd(1) = (M*(-3 + 4*Msq + 3*Nsq)* &
      (-5*u_V_FFP + 8*u_V_FFD - 3*u_V_FFF +  &
        3*Nsq*(-10*u_V_FFS + 15*u_V_FFP - 6*u_V_FFD + u_V_FFF)))/ &
    4.
              rVd(2) = (L*(-1 + 4*Msq + Nsq)* &
      (-5*u_V_FFP + 8*u_V_FFD - 3*u_V_FFF +  &
        3*Nsq*(-10*u_V_FFS + 15*u_V_FFP - 6*u_V_FFD + u_V_FFF)))/ &
    4.
              rVd(3) = (3*L*(L - M)*M*(L + M)*N* &
      (10*u_V_FFS - 15*u_V_FFP + 6*u_V_FFD - u_V_FFF))/2.
              rVd_t = (L*(L - M)*M*(L + M)* &
      (5*u_Vd_FFP - 8*u_Vd_FFD + 3*u_Vd_FFF +  &
        Nsq*(30*u_Vd_FFS - 3*(15*u_Vd_FFP - 6*u_Vd_FFD + u_Vd_FFF))))/ &
    4.
              case (7)
              rVd(1) = 0
              rVd(2) = (Sqrt(1.5)*N*(10*(20*M**4 + 15*Msq*(-1 + Nsq) +  &
           (-1 + Nsq)**2)*u_V_FFS - 20*u_V_FFP + 10*u_V_FFD +  &
        15*Msq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        20*M**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        15*Msq*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) +  &
        5*Nsq*(7*u_V_FFP - 4*u_V_FFD + u_V_FFF)))/4.
              rVd(3) = (Sqrt(1.5)*M*(10*(1 + 4*M**4 - 6*Nsq + 5*N**4 +  &
           5*Msq*(-1 + 3*Nsq))*u_V_FFS - 20*u_V_FFP + 10*u_V_FFD +  &
        5*Msq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        4*M**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        15*Msq*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        5*N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) +  &
        15*Nsq*(7*u_V_FFP - 4*u_V_FFD + u_V_FFF)))/4.
              rVd_t = (Sqrt(1.5)*M*N*(10*(4*M**4 + 5*Msq*(-1 + Nsq) + (-1 + Nsq)**2)* &
         u_Vd_FFS - 20*u_Vd_FFP + 10*u_Vd_FFD +  &
        5*Msq*(15*u_Vd_FFP - 6*u_Vd_FFD + u_Vd_FFF) -  &
        4*M**4*(15*u_Vd_FFP - 6*u_Vd_FFD + u_Vd_FFF) -  &
        5*Msq*Nsq*(15*u_Vd_FFP - 6*u_Vd_FFD + u_Vd_FFF) -  &
        N**4*(15*u_Vd_FFP - 6*u_Vd_FFD + u_Vd_FFF) +  &
        5*Nsq*(7*u_Vd_FFP - 4*u_Vd_FFD + u_Vd_FFF)))/4.
            end select
          case (3)
            select case (orb_dir_j)
              case (1)
              rVd(1) = (sqrt(15.)*(4*L**3* &
         (u_V_FFP + Nsq*(-5*u_V_FFP + 8*u_V_FFD - 3*u_V_FFF) -  &
           u_V_FFF) + 6*L*Msq*(-1 + Nsq)* &
         (2*(-1 + 5*Nsq)*u_V_FFS + (1 - 15*Nsq)*u_V_FFP +  &
           2*(1 + 3*Nsq)*u_V_FFD - (1 + Nsq)*u_V_FFF)))/ &
    (16.*(-1 + Nsq))
              rVd(2) = (sqrt(15.)*(6*Lsq*M*(-1 + Nsq)* &
         (2*(-1 + 5*Nsq)*u_V_FFS + (1 - 15*Nsq)*u_V_FFP +  &
           2*(1 + 3*Nsq)*u_V_FFD - (1 + Nsq)*u_V_FFF) +  &
        4*M**3*(-2*(1 - 6*Nsq + 5*N**4)*u_V_FFS + 2*u_V_FFD +  &
           Nsq*((-11 + 15*Nsq)*u_V_FFP - 2*(2 + 3*Nsq)*u_V_FFD +  &
              (3 + Nsq)*u_V_FFF))))/(16.*(-1 + Nsq))
              rVd(3) = (sqrt(15.)*N*(-4*M**4*(-1 + Nsq)* &
         (10*u_V_FFS - 15*u_V_FFP + 6*u_V_FFD - u_V_FFF) +  &
        4*(-1 + Nsq)*(u_V_FFP - 2*u_V_FFD + u_V_FFF) +  &
        Msq*(-30*(-1 + Nsq)**2*u_V_FFS + 53*u_V_FFP - 34*u_V_FFD +  &
           11*u_V_FFF - 6*Nsq* &
            (15*u_V_FFP - 6*u_V_FFD + u_V_FFF) +  &
           3*N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF))))/ &
    (8.*(-1 + Nsq))
              rVd_t = (sqrt(15.)*(L**4*(u_Vd_FFP +  &
           Nsq*(-5*u_Vd_FFP + 8*u_Vd_FFD - 3*u_Vd_FFF) - u_Vd_FFF) +  &
        3*Lsq*Msq*(-1 + Nsq)* &
         (2*(-1 + 5*Nsq)*u_Vd_FFS + (1 - 15*Nsq)*u_Vd_FFP +  &
           2*(1 + 3*Nsq)*u_Vd_FFD - (1 + Nsq)*u_Vd_FFF) +  &
        M**4*(-2*(1 - 6*Nsq + 5*N**4)*u_Vd_FFS + 2*u_Vd_FFD +  &
           Nsq*((-11 + 15*Nsq)*u_Vd_FFP - 2*(2 + 3*Nsq)*u_Vd_FFD +  &
              (3 + Nsq)*u_Vd_FFF))))/(16.*(-1 + Nsq))
              case (2)
              rVd(1) = (Sqrt(2.5)*N*((-1 + 5*Nsq)*u_V_FFP + 4*u_V_FFD -  &
        3*u_V_FFF + Nsq*(-8*u_V_FFD + 3*u_V_FFF) +  &
        Msq*(6*(-1 + 5*Nsq)*u_V_FFS + 13*u_V_FFP - 10*u_V_FFD +  &
           3*u_V_FFF - 3*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF))) &
 )/4.
              rVd(2) = (Sqrt(2.5)*L*M*N*(6*(-1 + 5*Nsq)*u_V_FFS + 13*u_V_FFP -  &
        10*u_V_FFD + 3*u_V_FFF -  &
        3*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF)))/2.
              rVd(3) = (Sqrt(2.5)*L*(-u_V_FFP + 4*u_V_FFD - 3*u_V_FFF +  &
        3*Nsq*(5*u_V_FFP - 8*u_V_FFD + 3*u_V_FFF) +  &
        Msq*((-6 + 90*Nsq)*u_V_FFS + 13*u_V_FFP - 10*u_V_FFD +  &
           3*u_V_FFF - 9*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF))) &
 )/4.
              rVd_t = (Sqrt(2.5)*L*N*((-1 + 5*Nsq)*u_Vd_FFP + 4*u_Vd_FFD -  &
        3*u_Vd_FFF + Nsq*(-8*u_Vd_FFD + 3*u_Vd_FFF) +  &
        Msq*(6*(-1 + 5*Nsq)*u_Vd_FFS + 13*u_Vd_FFP - 10*u_Vd_FFD +  &
           3*u_Vd_FFF - 3*Nsq*(15*u_Vd_FFP - 6*u_Vd_FFD + u_Vd_FFF))) &
 )/4.
              case (3)
              rVd(1) = 0
              rVd(2) = (M*(6*(1 - 5*Nsq)**2*u_V_FFS +  &
        (-1 + 130*Nsq - 225*N**4)*u_V_FFP +  &
        5*(-1 + Nsq)*(2*(-1 + 9*Nsq)*u_V_FFD - 3*(-1 + Nsq)*u_V_FFF)) &
 )/8.
              rVd(3) = (5*N*((-1 + 5*Nsq)*u_V_FFP + 4*u_V_FFD - 3*u_V_FFF +  &
        Nsq*(-8*u_V_FFD + 3*u_V_FFF) +  &
        Msq*(6*(-1 + 5*Nsq)*u_V_FFS + 13*u_V_FFP - 10*u_V_FFD +  &
           3*u_V_FFF - 3*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF))) &
 )/4.
              rVd_t = ((1 - 5*Nsq)**2*u_Vd_FFP -  &
      5*(-1 + Nsq)*(Nsq*(8*u_Vd_FFD - 3*u_Vd_FFF) + 3*u_Vd_FFF) +  &
      Msq*(6*(1 - 5*Nsq)**2*u_Vd_FFS +  &
         (-1 + 130*Nsq - 225*N**4)*u_Vd_FFP +  &
         5*(-1 + Nsq)*(2*(-1 + 9*Nsq)*u_Vd_FFD - 3*(-1 + Nsq)*u_Vd_FFF) &
 ))/16.
              case (4)
              rVd(1) = 0
              rVd(2) = (Sqrt(1.5)*N*((6 - 40*Nsq + 50*N**4)*u_V_FFS +  &
        (-11 + 70*Nsq - 75*N**4)*u_V_FFP +  &
        5*(-1 + Nsq)*((-2 + 6*Nsq)*u_V_FFD - (-1 + Nsq)*u_V_FFF)))/8.
              rVd(3) = (Sqrt(1.5)*M*(2*(3 - 60*Nsq + 125*N**4)*u_V_FFS -  &
        11*u_V_FFP + 10*u_V_FFD - 5*u_V_FFF -  &
        25*N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) +  &
        30*Nsq*(7*u_V_FFP - 4*u_V_FFD + u_V_FFF)))/8.
              rVd_t = (Sqrt(1.5)*M*N*((6 - 40*Nsq + 50*N**4)*u_Vd_FFS +  &
        (-11 + 70*Nsq - 75*N**4)*u_Vd_FFP +  &
        5*(-1 + Nsq)*((-2 + 6*Nsq)*u_Vd_FFD - (-1 + Nsq)*u_Vd_FFF)))/8.
              case (5)
              rVd(1) = (M*(6*(1 - 5*Nsq)**2*u_V_FFS +  &
        (-1 + 130*Nsq - 225*N**4)*u_V_FFP +  &
        5*(-1 + Nsq)*(2*(-1 + 9*Nsq)*u_V_FFD - 3*(-1 + Nsq)*u_V_FFF)) &
 )/16.
              rVd(2) = (L*(6*(1 - 5*Nsq)**2*u_V_FFS +  &
        (-1 + 130*Nsq - 225*N**4)*u_V_FFP +  &
        5*(-1 + Nsq)*(2*(-1 + 9*Nsq)*u_V_FFD - 3*(-1 + Nsq)*u_V_FFF)) &
 )/16.
              rVd(3) = (5*L*M*N*(6*(-1 + 5*Nsq)*u_V_FFS + 13*u_V_FFP -  &
        10*u_V_FFD + 3*u_V_FFF -  &
        3*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF)))/4.
              rVd_t = (L*M*(6*(1 - 5*Nsq)**2*u_Vd_FFS +  &
        (-1 + 130*Nsq - 225*N**4)*u_Vd_FFP +  &
        5*(-1 + Nsq)*(2*(-1 + 9*Nsq)*u_Vd_FFD - 3*(-1 + Nsq)*u_Vd_FFF)) &
 )/16.
              case (6)
              rVd(1) = 0
              rVd(2) = (Sqrt(2.5)*N*(-6*(-1 + 6*Msq + Nsq)*(-1 + 5*Nsq)* &
         u_V_FFS + 15*u_V_FFP - 18*u_V_FFD +  &
        Msq*(-78*u_V_FFP + 60*u_V_FFD - 18*u_V_FFF) +  &
        Nsq*(-68*u_V_FFP + 44*u_V_FFD - 12*u_V_FFF) +  &
        9*u_V_FFF + 18*Msq*Nsq* &
         (15*u_V_FFP - 6*u_V_FFD + u_V_FFF) +  &
        3*N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF)))/8.
              rVd(3) = (Sqrt(2.5)*M*(-6*(1 - 18*Nsq + 25*N**4 + Msq*(-2 + 30*Nsq))* &
         u_V_FFS + 15*u_V_FFP - 18*u_V_FFD +  &
        Msq*(-26*u_V_FFP + 20*u_V_FFD - 6*u_V_FFF) + 9*u_V_FFF +  &
        18*Msq*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) +  &
        15*N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        12*Nsq*(17*u_V_FFP - 11*u_V_FFD + 3*u_V_FFF)))/8.
              rVd_t = (Sqrt(2.5)*M*N*(-6*(-1 + 2*Msq + Nsq)*(-1 + 5*Nsq)*u_Vd_FFS +  &
        (15 - 68*Nsq + 45*N**4 + Msq*(-26 + 90*Nsq))*u_Vd_FFP -  &
        18*u_Vd_FFD + Nsq*(44*u_Vd_FFD - 12*u_Vd_FFF) +  &
        Msq*(20*u_Vd_FFD - 6*u_Vd_FFF) + 9*u_Vd_FFF +  &
        6*Msq*Nsq*(-6*u_Vd_FFD + u_Vd_FFF) +  &
        3*N**4*(-6*u_Vd_FFD + u_Vd_FFF)))/8.
              case (7)
              rVd(1) = -(sqrt(15.)*M*(2*(-1 + 4*Msq + Nsq)*(-1 + 5*Nsq)*u_V_FFS -  &
         3*u_V_FFP - 2*u_V_FFD +  &
         4*Msq*(u_V_FFP + 2*u_V_FFD - u_V_FFF) + 3*u_V_FFF -  &
         4*Msq*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
         N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) +  &
         Nsq*(26*u_V_FFP - 20*u_V_FFD + 6*u_V_FFF)))/16.
              rVd(2) = -(sqrt(15.)*L*(2*(-1 + 12*Msq + Nsq)*(-1 + 5*Nsq)* &
          u_V_FFS - 3*u_V_FFP - 2*u_V_FFD +  &
         12*Msq*(u_V_FFP + 2*u_V_FFD - u_V_FFF) + 3*u_V_FFF -  &
         12*Msq*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
         N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) +  &
         Nsq*(26*u_V_FFP - 20*u_V_FFD + 6*u_V_FFF)))/16.
              rVd(3) = -(sqrt(15.)*L*M*N*(2*(-3 + 10*Msq + 5*Nsq)*u_V_FFS +  &
         13*u_V_FFP - 10*u_V_FFD + 3*u_V_FFF -  &
         2*Msq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
         Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF)))/4.
              rVd_t = -(sqrt(15.)*L*M*(2*(-1 + 4*Msq + Nsq)*(-1 + 5*Nsq)*u_Vd_FFS -  &
         3*u_Vd_FFP - 2*u_Vd_FFD +  &
         4*Msq*(u_Vd_FFP + 2*u_Vd_FFD - u_Vd_FFF) + 3*u_Vd_FFF -  &
         4*Msq*Nsq*(15*u_Vd_FFP - 6*u_Vd_FFD + u_Vd_FFF) -  &
         N**4*(15*u_Vd_FFP - 6*u_Vd_FFD + u_Vd_FFF) +  &
         Nsq*(26*u_Vd_FFP - 20*u_Vd_FFD + 6*u_Vd_FFF)))/16.
            end select
          case (4)
            select case (orb_dir_j)
              case (1)
              rVd(1) = (3*Sqrt(2.5)*L*M*N* &
      (2*(-3 + 5*Nsq)*u_V_FFS + (3 - 15*Nsq)*u_V_FFP +  &
        6*(1 + Nsq)*u_V_FFD - (3 + Nsq)*u_V_FFF))/4.
              rVd(2) = (-3*Sqrt(2.5)*N*(-1 + 2*Msq + Nsq)* &
      (2*(-3 + 5*Nsq)*u_V_FFS + 3*u_V_FFP + 6*u_V_FFD -  &
        3*u_V_FFF - Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF)))/8.
              rVd(3) = (-3*Sqrt(2.5)*M*(-3 + 4*Msq + 3*Nsq)* &
      (2*(-1 + 5*Nsq)*u_V_FFS + u_V_FFP + 2*u_V_FFD - u_V_FFF -  &
        Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF)))/8.
              rVd_t = -(Sqrt(2.5)*M*(-3*Lsq + Msq)*N* &
       (2*(-3 + 5*Nsq)*u_Vd_FFS + (3 - 15*Nsq)*u_Vd_FFP +  &
         6*(1 + Nsq)*u_Vd_FFD - (3 + Nsq)*u_Vd_FFF))/8.
              case (2)
              rVd(1) = (sqrt(15.)*M*(-u_V_FFP -  &
        2*Nsq*(3*u_V_FFS - 4*u_V_FFP + u_V_FFD) +  &
        N**4*(10*u_V_FFS - 15*u_V_FFP + 6*u_V_FFD - u_V_FFF) +  &
        u_V_FFF))/4.
              rVd(2) = (sqrt(15.)*L*(-u_V_FFP -  &
        2*Nsq*(3*u_V_FFS - 4*u_V_FFP + u_V_FFD) +  &
        N**4*(10*u_V_FFS - 15*u_V_FFP + 6*u_V_FFD - u_V_FFF) +  &
        u_V_FFF))/4.
              rVd(3) = sqrt(15.)*L*M*N*((-3 + 10*Nsq)*u_V_FFS + 4*u_V_FFP -  &
      u_V_FFD - Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF))
              rVd_t = (sqrt(15.)*L*M*(-u_Vd_FFP -  &
        2*Nsq*(3*u_Vd_FFS - 4*u_Vd_FFP + u_Vd_FFD) +  &
        N**4*(10*u_Vd_FFS - 15*u_Vd_FFP + 6*u_Vd_FFD - u_Vd_FFF) +  &
        u_Vd_FFF))/4.
              case (3)
              rVd(1) = 0
              rVd(2) = (Sqrt(1.5)*N*((6 - 40*Nsq + 50*N**4)*u_V_FFS +  &
        (-11 + 70*Nsq - 75*N**4)*u_V_FFP +  &
        5*(-1 + Nsq)*((-2 + 6*Nsq)*u_V_FFD - (-1 + Nsq)*u_V_FFF)))/8.
              rVd(3) = (Sqrt(1.5)*M*(2*(3 - 60*Nsq + 125*N**4)*u_V_FFS -  &
        11*u_V_FFP + 10*u_V_FFD - 5*u_V_FFF -  &
        25*N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) +  &
        30*Nsq*(7*u_V_FFP - 4*u_V_FFD + u_V_FFF)))/8.
              rVd_t = (Sqrt(1.5)*M*N*((6 - 40*Nsq + 50*N**4)*u_Vd_FFS +  &
        (-11 + 70*Nsq - 75*N**4)*u_Vd_FFP +  &
        5*(-1 + Nsq)*((-2 + 6*Nsq)*u_Vd_FFD - (-1 + Nsq)*u_Vd_FFF)))/8.
              case (4)
              rVd(1) = 0
              rVd(2) = 0
              rVd(3) = (3*N*((6 - 40*Nsq + 50*N**4)*u_V_FFS +  &
        (-11 + 70*Nsq - 75*N**4)*u_V_FFP +  &
        5*(-1 + Nsq)*((-2 + 6*Nsq)*u_V_FFD - (-1 + Nsq)*u_V_FFF)))/4.
              rVd_t = (2*Nsq*(3 - 5*Nsq)**2*u_Vd_FFS -  &
      3*(1 - 5*Nsq)**2*(-1 + Nsq)*u_Vd_FFP +  &
      30*Nsq*(-1 + Nsq)**2*u_Vd_FFD - 5*(-1 + Nsq)**3*u_Vd_FFF)/8.
              case (5)
              rVd(1) = (Sqrt(1.5)*N*((6 - 40*Nsq + 50*N**4)*u_V_FFS +  &
        (-11 + 70*Nsq - 75*N**4)*u_V_FFP +  &
        5*(-1 + Nsq)*((-2 + 6*Nsq)*u_V_FFD - (-1 + Nsq)*u_V_FFF)))/8.
              rVd(2) = 0
              rVd(3) = (Sqrt(1.5)*L*(2*(3 - 60*Nsq + 125*N**4)*u_V_FFS -  &
        11*u_V_FFP + 10*u_V_FFD - 5*u_V_FFF -  &
        25*N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) +  &
        30*Nsq*(7*u_V_FFP - 4*u_V_FFD + u_V_FFF)))/8.
              rVd_t = (Sqrt(1.5)*L*N*((6 - 40*Nsq + 50*N**4)*u_Vd_FFS +  &
        (-11 + 70*Nsq - 75*N**4)*u_Vd_FFP +  &
        5*(-1 + Nsq)*((-2 + 6*Nsq)*u_Vd_FFD - (-1 + Nsq)*u_Vd_FFF)))/8.
              case (6)
              rVd(1) = 0
              rVd(2) = (sqrt(15.)*M*(u_V_FFP +  &
        2*Nsq*(3*u_V_FFS - 4*u_V_FFP + u_V_FFD) - u_V_FFF +  &
        N**4*(-10*u_V_FFS + 15*u_V_FFP - 6*u_V_FFD + u_V_FFF)))/2.
              rVd(3) = (sqrt(15.)*N*(u_V_FFP +  &
        2*Nsq*(3*u_V_FFS - 4*u_V_FFP + u_V_FFD) - u_V_FFF +  &
        N**4*(-10*u_V_FFS + 15*u_V_FFP - 6*u_V_FFD + u_V_FFF) +  &
        2*(-1 + 2*Msq + Nsq)*((3 - 10*Nsq)*u_V_FFS +  &
           (-4 + 15*Nsq)*u_V_FFP + u_V_FFD +  &
           Nsq*(-6*u_V_FFD + u_V_FFF))))/4.
              rVd_t = (sqrt(15.)*(-1 + 2*Msq + Nsq)* &
      (u_Vd_FFP + 2*Nsq*(3*u_Vd_FFS - 4*u_Vd_FFP + u_Vd_FFD) -  &
        u_Vd_FFF + N**4*(-10*u_Vd_FFS + 15*u_Vd_FFP - 6*u_Vd_FFD +  &
           u_Vd_FFF)))/8.
              case (7)
              rVd(1) = (-3*Sqrt(2.5)*N*(-1 + 2*Msq + Nsq)* &
      (2*(-3 + 5*Nsq)*u_V_FFS + 3*u_V_FFP + 6*u_V_FFD -  &
        3*u_V_FFF - Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF)))/8.
              rVd(2) = (-3*Sqrt(2.5)*L*M*N* &
      (2*(-3 + 5*Nsq)*u_V_FFS + (3 - 15*Nsq)*u_V_FFP +  &
        6*(1 + Nsq)*u_V_FFD - (3 + Nsq)*u_V_FFF))/4.
              rVd(3) = (-3*Sqrt(2.5)*L*(-1 + 4*Msq + Nsq)* &
      (2*(-1 + 5*Nsq)*u_V_FFS + u_V_FFP + 2*u_V_FFD - u_V_FFF -  &
        Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF)))/8.
              rVd_t = (Sqrt(2.5)*L*(Lsq - 3*Msq)*N* &
      (2*(-3 + 5*Nsq)*u_Vd_FFS + (3 - 15*Nsq)*u_Vd_FFP +  &
        6*(1 + Nsq)*u_Vd_FFD - (3 + Nsq)*u_Vd_FFF))/8.
            end select
          case (5)
            select case (orb_dir_j)
              case (1)
              rVd(1) = -(sqrt(15.)*M*(2*(-3 + 4*Msq + 3*Nsq)*(-1 + 5*Nsq)* &
          u_V_FFS - u_V_FFP - 6*u_V_FFD +  &
         Nsq*(38*u_V_FFP + 4*u_V_FFD - 6*u_V_FFF) +  &
         4*Msq*(u_V_FFP + 2*u_V_FFD - u_V_FFF) + u_V_FFF -  &
         4*Msq*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
         3*N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF)))/16.
              rVd(2) = -(sqrt(15.)*L*(6*(-1 + 4*Msq + Nsq)*(-1 + 5*Nsq)*u_V_FFS -  &
         u_V_FFP - 6*u_V_FFD +  &
         Nsq*(38*u_V_FFP + 4*u_V_FFD - 6*u_V_FFF) +  &
         12*Msq*(u_V_FFP + 2*u_V_FFD - u_V_FFF) + u_V_FFF -  &
         12*Msq*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
         3*N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF)))/16.
              rVd(3) = -(sqrt(15.)*L*M*N*(2*(-9 + 10*Msq + 15*Nsq)*u_V_FFS +  &
         19*u_V_FFP + 2*u_V_FFD - 3*u_V_FFF -  &
         2*Msq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
         3*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF)))/4.
              rVd_t = -(sqrt(15.)*L*M*(2*(-3 + 4*Msq + 3*Nsq)*(-1 + 5*Nsq)* &
          u_Vd_FFS - u_Vd_FFP - 6*u_Vd_FFD +  &
         Nsq*(38*u_Vd_FFP + 4*u_Vd_FFD - 6*u_Vd_FFF) +  &
         4*Msq*(u_Vd_FFP + 2*u_Vd_FFD - u_Vd_FFF) + u_Vd_FFF -  &
         4*Msq*Nsq*(15*u_Vd_FFP - 6*u_Vd_FFD + u_Vd_FFF) -  &
         3*N**4*(15*u_Vd_FFP - 6*u_Vd_FFD + u_Vd_FFF)))/16.
              case (2)
              rVd(1) = (Sqrt(2.5)*L*M*N*(6*(1 - 6*Nsq + 5*N**4)*u_V_FFS +  &
        (-12 + 53*Nsq - 45*N**4)*u_V_FFP +  &
        2*(3 - 10*Nsq + 9*N**4)*u_V_FFD - 3*Nsq*(-1 + Nsq)*u_V_FFF))/ &
    (2.*(-1 + Nsq))
              rVd(2) = (Sqrt(2.5)*N*(3*Msq* &
         ((1 - 5*Nsq)*u_V_FFP + (-4 + 8*Nsq)*u_V_FFD -  &
           3*(-1 + Nsq)*u_V_FFF) -  &
        (-1 + Msq + Nsq)*(6*(1 - 6*Nsq + 5*N**4)*u_V_FFS -  &
           12*u_V_FFP + 6*u_V_FFD -  &
           3*N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) +  &
           Nsq*(53*u_V_FFP - 20*u_V_FFD + 3*u_V_FFF))))/ &
    (4.*(-1 + Nsq))
              rVd(3) = (Sqrt(2.5)*M*(-6*(-1 + Msq + Nsq)*(1 - 16*Nsq + 15*N**4)* &
         u_V_FFS - 12*u_V_FFP + 6*u_V_FFD +  &
        9*Msq*N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) +  &
        9*N**6*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        4*Msq*Nsq*(37*u_V_FFP - 16*u_V_FFD + 3*u_V_FFF) +  &
        Msq*(13*u_V_FFP - 10*u_V_FFD + 3*u_V_FFF) -  &
        2*N**4*(139*u_V_FFP - 55*u_V_FFD + 9*u_V_FFF) +  &
        Nsq*(147*u_V_FFP + 9*(-6*u_V_FFD + u_V_FFF))))/ &
    (4.*(-1 + Nsq))
              rVd_t = (Sqrt(2.5)*M*N*(Msq* &
         ((1 - 5*Nsq)*u_Vd_FFP + (-4 + 8*Nsq)*u_Vd_FFD -  &
           3*(-1 + Nsq)*u_Vd_FFF) +  &
        Lsq*(6*(1 - 6*Nsq + 5*N**4)*u_Vd_FFS +  &
           (-12 + 53*Nsq - 45*N**4)*u_Vd_FFP +  &
           2*(3 - 10*Nsq + 9*N**4)*u_Vd_FFD - 3*Nsq*(-1 + Nsq)*u_Vd_FFF &
 )))/(4.*(-1 + Nsq))
              case (3)
              rVd(1) = (M*(6*(1 - 5*Nsq)**2*u_V_FFS +  &
        (-1 + 130*Nsq - 225*N**4)*u_V_FFP +  &
        5*(-1 + Nsq)*(2*(-1 + 9*Nsq)*u_V_FFD - 3*(-1 + Nsq)*u_V_FFF)) &
 )/16.
              rVd(2) = (L*(6*(1 - 5*Nsq)**2*u_V_FFS +  &
        (-1 + 130*Nsq - 225*N**4)*u_V_FFP +  &
        5*(-1 + Nsq)*(2*(-1 + 9*Nsq)*u_V_FFD - 3*(-1 + Nsq)*u_V_FFF)) &
 )/16.
              rVd(3) = (5*L*M*N*(6*(-1 + 5*Nsq)*u_V_FFS + 13*u_V_FFP -  &
        10*u_V_FFD + 3*u_V_FFF -  &
        3*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF)))/4.
              rVd_t = (L*M*(6*(1 - 5*Nsq)**2*u_Vd_FFS +  &
        (-1 + 130*Nsq - 225*N**4)*u_Vd_FFP +  &
        5*(-1 + Nsq)*(2*(-1 + 9*Nsq)*u_Vd_FFD - 3*(-1 + Nsq)*u_Vd_FFF)) &
 )/16.
              case (4)
              rVd(1) = (Sqrt(1.5)*N*((6 - 40*Nsq + 50*N**4)*u_V_FFS +  &
        (-11 + 70*Nsq - 75*N**4)*u_V_FFP +  &
        5*(-1 + Nsq)*((-2 + 6*Nsq)*u_V_FFD - (-1 + Nsq)*u_V_FFF)))/8.
              rVd(2) = 0
              rVd(3) = (Sqrt(1.5)*L*(2*(3 - 60*Nsq + 125*N**4)*u_V_FFS -  &
        11*u_V_FFP + 10*u_V_FFD - 5*u_V_FFF -  &
        25*N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) +  &
        30*Nsq*(7*u_V_FFP - 4*u_V_FFD + u_V_FFF)))/8.
              rVd_t = (Sqrt(1.5)*L*N*((6 - 40*Nsq + 50*N**4)*u_Vd_FFS +  &
        (-11 + 70*Nsq - 75*N**4)*u_Vd_FFP +  &
        5*(-1 + Nsq)*((-2 + 6*Nsq)*u_Vd_FFD - (-1 + Nsq)*u_Vd_FFF)))/8.
              case (5)
              rVd(1) = (L*(6*(1 - 5*Nsq)**2*(-1 + Nsq)*u_V_FFS - 10*u_V_FFD +  &
        Nsq*(-((11 - 15*Nsq)**2*u_V_FFP) +  &
           10*(7 - 15*Nsq + 9*N**4)*u_V_FFD - 15*(-1 + Nsq)**2*u_V_FFF &
 )))/(8.*(-1 + Nsq))
              rVd(2) = (M*(-((1 - 5*Nsq)**2*u_V_FFP) +  &
        5*(-1 + Nsq)*(8*Nsq*u_V_FFD - 3*(-1 + Nsq)*u_V_FFF)))/ &
    (8.*(-1 + Nsq))
              rVd(3) = -(N*(60*(-1 + Msq + Nsq)*(1 - 6*Nsq + 5*N**4)*u_V_FFS +  &
         (121 - 660*Nsq + 1005*N**4 - 450*N**6 -  &
            10*Msq*(13 - 58*Nsq + 45*N**4))*u_V_FFP +  &
         5*(-1 + Nsq)*(4*(3 - 12*Nsq + 9*N**4 + Msq*(-5 + 9*Nsq))* &
             u_V_FFD - 3*(-1 + Nsq)*(-1 + 2*Msq + 2*Nsq)*u_V_FFF)))/ &
    (8.*(-1 + Nsq))
              rVd_t = (Msq*(-((1 - 5*Nsq)**2*u_Vd_FFP) +  &
         5*(-1 + Nsq)*(8*Nsq*u_Vd_FFD - 3*(-1 + Nsq)*u_Vd_FFF)) +  &
      Lsq*(6*(1 - 5*Nsq)**2*(-1 + Nsq)*u_Vd_FFS - 10*u_Vd_FFD +  &
         Nsq*(-((11 - 15*Nsq)**2*u_Vd_FFP) +  &
            10*(7 - 15*Nsq + 9*N**4)*u_Vd_FFD -  &
            15*(-1 + Nsq)**2*u_Vd_FFF)))/(16.*(-1 + Nsq))
              case (6)
              rVd(1) = (Sqrt(2.5)*N*(-6*(-1 + 2*Msq + Nsq)*(-1 + 5*Nsq)* &
         u_V_FFS + (11 - 48*Nsq + 45*N**4 + Msq*(-26 + 90*Nsq))* &
         u_V_FFP - 2*u_V_FFD + 12*Nsq*u_V_FFD +  &
        Msq*(20*u_V_FFD - 6*u_V_FFF) - 3*u_V_FFF +  &
        6*Msq*Nsq*(-6*u_V_FFD + u_V_FFF) +  &
        3*N**4*(-6*u_V_FFD + u_V_FFF)))/8.
              rVd(2) = (Sqrt(2.5)*L*M*N*((6 - 30*Nsq)*u_V_FFS +  &
        (-13 + 45*Nsq)*u_V_FFP + 10*u_V_FFD - 3*u_V_FFF +  &
        3*Nsq*(-6*u_V_FFD + u_V_FFF)))/2.
              rVd(3) = (Sqrt(2.5)*L*(-6*(1 - 18*Nsq + 25*N**4 + Msq*(-2 + 30*Nsq))* &
         u_V_FFS + 11*u_V_FFP - 2*u_V_FFD +  &
        36*Nsq*(-4*u_V_FFP + u_V_FFD) +  &
        Msq*(-26*u_V_FFP + 20*u_V_FFD - 6*u_V_FFF) - 3*u_V_FFF +  &
        18*Msq*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) +  &
        15*N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF)))/8.
              rVd_t = (Sqrt(2.5)*L*N*(-6*(-1 + 2*Msq + Nsq)*(-1 + 5*Nsq)*u_Vd_FFS +  &
        (11 - 48*Nsq + 45*N**4 + Msq*(-26 + 90*Nsq))*u_Vd_FFP -  &
        2*u_Vd_FFD + 12*Nsq*u_Vd_FFD +  &
        Msq*(20*u_Vd_FFD - 6*u_Vd_FFF) - 3*u_Vd_FFF +  &
        6*Msq*Nsq*(-6*u_Vd_FFD + u_Vd_FFF) +  &
        3*N**4*(-6*u_Vd_FFD + u_Vd_FFF)))/8.
              case (7)
              rVd(1) = (sqrt(15.)*(-6*L*Msq*(-1 + Nsq)* &
         (2*(-1 + 5*Nsq)*u_V_FFS + (1 - 15*Nsq)*u_V_FFP +  &
           2*(1 + 3*Nsq)*u_V_FFD - (1 + Nsq)*u_V_FFF) +  &
        4*L**3*(2*(1 - 6*Nsq + 5*N**4)*u_V_FFS - 2*u_V_FFD +  &
           Nsq*((11 - 15*Nsq)*u_V_FFP + (4 + 6*Nsq)*u_V_FFD -  &
              (3 + Nsq)*u_V_FFF))))/(16.*(-1 + Nsq))
              rVd(2) = (sqrt(15.)*(-6*Lsq*M*(-1 + Nsq)* &
         (2*(-1 + 5*Nsq)*u_V_FFS + (1 - 15*Nsq)*u_V_FFP +  &
           2*(1 + 3*Nsq)*u_V_FFD - (1 + Nsq)*u_V_FFF) +  &
        4*M**3*((-1 + 5*Nsq)*u_V_FFP + u_V_FFF +  &
           Nsq*(-8*u_V_FFD + 3*u_V_FFF))))/(16.*(-1 + Nsq))
              rVd(3) = (sqrt(15.)*N*(10*(-1 + Nsq)* &
         (4*M**4 + 5*Msq*(-1 + Nsq) + (-1 + Nsq)**2)*u_V_FFS +  &
        11*u_V_FFP + 2*u_V_FFD - 3*u_V_FFF +  &
        4*M**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) +  &
        10*Msq*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        4*M**4*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) +  &
        3*N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        5*Msq*N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        N**6*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) +  &
        Nsq*(-41*u_V_FFP + 10*u_V_FFD + u_V_FFF) +  &
        Msq*(-67*u_V_FFP + 14*u_V_FFD + 3*u_V_FFF)))/(8.*(-1 + Nsq))
              rVd_t = (sqrt(15.)*(-3*Lsq*Msq*(-1 + Nsq)* &
         (2*(-1 + 5*Nsq)*u_Vd_FFS + (1 - 15*Nsq)*u_Vd_FFP +  &
           2*(1 + 3*Nsq)*u_Vd_FFD - (1 + Nsq)*u_Vd_FFF) +  &
        M**4*((-1 + 5*Nsq)*u_Vd_FFP + u_Vd_FFF +  &
           Nsq*(-8*u_Vd_FFD + 3*u_Vd_FFF)) +  &
        L**4*(2*(1 - 6*Nsq + 5*N**4)*u_Vd_FFS - 2*u_Vd_FFD +  &
           Nsq*((11 - 15*Nsq)*u_Vd_FFP + (4 + 6*Nsq)*u_Vd_FFD -  &
              (3 + Nsq)*u_Vd_FFF))))/(16.*(-1 + Nsq))
            end select
          case (6)
            select case (orb_dir_j)
              case (1)
              rVd(1) = 0
              rVd(2) = (Sqrt(1.5)*N*(10*(40*M**4 + 30*Msq*(-1 + Nsq) +  &
           3*(-1 + Nsq)**2)*u_V_FFS - 35*u_V_FFP +  &
        Nsq*(80*u_V_FFP - 20*u_V_FFD) + 10*u_V_FFD - 5*u_V_FFF +  &
        30*Msq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        40*M**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        30*Msq*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        3*N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF)))/8.
              rVd(3) = (Sqrt(1.5)*M*(10*(8*M**4 + 10*Msq*(-1 + 3*Nsq) +  &
           3*(1 - 6*Nsq + 5*N**4))*u_V_FFS - 35*u_V_FFP +  &
        60*Nsq*(4*u_V_FFP - u_V_FFD) + 10*u_V_FFD - 5*u_V_FFF +  &
        10*Msq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        8*M**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        30*Msq*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        15*N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF)))/8.
              rVd_t = (Sqrt(1.5)*M*N*(10*(8*M**4 + 10*Msq*(-1 + Nsq) +  &
           3*(-1 + Nsq)**2)*u_Vd_FFS - 35*u_Vd_FFP +  &
        Nsq*(80*u_Vd_FFP - 20*u_Vd_FFD) + 10*u_Vd_FFD - 5*u_Vd_FFF +  &
        10*Msq*(15*u_Vd_FFP - 6*u_Vd_FFD + u_Vd_FFF) -  &
        8*M**4*(15*u_Vd_FFP - 6*u_Vd_FFD + u_Vd_FFF) -  &
        10*Msq*Nsq*(15*u_Vd_FFP - 6*u_Vd_FFD + u_Vd_FFF) -  &
        3*N**4*(15*u_Vd_FFP - 6*u_Vd_FFD + u_Vd_FFF)))/8.
              case (2)
              rVd(1) = (M*(-3 + 4*Msq + 3*Nsq)* &
      (-5*u_V_FFP + 8*u_V_FFD - 3*u_V_FFF +  &
        3*Nsq*(-10*u_V_FFS + 15*u_V_FFP - 6*u_V_FFD + u_V_FFF)))/ &
    4.
              rVd(2) = (L*(-1 + 4*Msq + Nsq)* &
      (-5*u_V_FFP + 8*u_V_FFD - 3*u_V_FFF +  &
        3*Nsq*(-10*u_V_FFS + 15*u_V_FFP - 6*u_V_FFD + u_V_FFF)))/ &
    4.
              rVd(3) = (3*L*(L - M)*M*(L + M)*N* &
      (10*u_V_FFS - 15*u_V_FFP + 6*u_V_FFD - u_V_FFF))/2.
              rVd_t = (L*(L - M)*M*(L + M)* &
      (5*u_Vd_FFP - 8*u_Vd_FFD + 3*u_Vd_FFF +  &
        Nsq*(30*u_Vd_FFS - 3*(15*u_Vd_FFP - 6*u_Vd_FFD + u_Vd_FFF))))/ &
    4.
              case (3)
              rVd(1) = 0
              rVd(2) = (Sqrt(2.5)*N*(-6*(-1 + 6*Msq + Nsq)*(-1 + 5*Nsq)* &
         u_V_FFS + 15*u_V_FFP - 18*u_V_FFD +  &
        Msq*(-78*u_V_FFP + 60*u_V_FFD - 18*u_V_FFF) +  &
        Nsq*(-68*u_V_FFP + 44*u_V_FFD - 12*u_V_FFF) +  &
        9*u_V_FFF + 18*Msq*Nsq* &
         (15*u_V_FFP - 6*u_V_FFD + u_V_FFF) +  &
        3*N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF)))/8.
              rVd(3) = (Sqrt(2.5)*M*(-6*(1 - 18*Nsq + 25*N**4 + Msq*(-2 + 30*Nsq))* &
         u_V_FFS + 15*u_V_FFP - 18*u_V_FFD +  &
        Msq*(-26*u_V_FFP + 20*u_V_FFD - 6*u_V_FFF) + 9*u_V_FFF +  &
        18*Msq*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) +  &
        15*N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        12*Nsq*(17*u_V_FFP - 11*u_V_FFD + 3*u_V_FFF)))/8.
              rVd_t = (Sqrt(2.5)*M*N*(-6*(-1 + 2*Msq + Nsq)*(-1 + 5*Nsq)*u_Vd_FFS +  &
        (15 - 68*Nsq + 45*N**4 + Msq*(-26 + 90*Nsq))*u_Vd_FFP -  &
        18*u_Vd_FFD + Nsq*(44*u_Vd_FFD - 12*u_Vd_FFF) +  &
        Msq*(20*u_Vd_FFD - 6*u_Vd_FFF) + 9*u_Vd_FFF +  &
        6*Msq*Nsq*(-6*u_Vd_FFD + u_Vd_FFF) +  &
        3*N**4*(-6*u_Vd_FFD + u_Vd_FFF)))/8.
              case (4)
              rVd(1) = 0
              rVd(2) = (sqrt(15.)*M*(u_V_FFP +  &
        2*Nsq*(3*u_V_FFS - 4*u_V_FFP + u_V_FFD) - u_V_FFF +  &
        N**4*(-10*u_V_FFS + 15*u_V_FFP - 6*u_V_FFD + u_V_FFF)))/2.
              rVd(3) = (sqrt(15.)*N*(u_V_FFP +  &
        2*Nsq*(3*u_V_FFS - 4*u_V_FFP + u_V_FFD) - u_V_FFF +  &
        N**4*(-10*u_V_FFS + 15*u_V_FFP - 6*u_V_FFD + u_V_FFF) +  &
        2*(-1 + 2*Msq + Nsq)*((3 - 10*Nsq)*u_V_FFS +  &
           (-4 + 15*Nsq)*u_V_FFP + u_V_FFD +  &
           Nsq*(-6*u_V_FFD + u_V_FFF))))/4.
              rVd_t = (sqrt(15.)*(-1 + 2*Msq + Nsq)* &
      (u_Vd_FFP + 2*Nsq*(3*u_Vd_FFS - 4*u_Vd_FFP + u_Vd_FFD) -  &
        u_Vd_FFF + N**4*(-10*u_Vd_FFS + 15*u_Vd_FFP - 6*u_Vd_FFD +  &
           u_Vd_FFF)))/8.
              case (5)
              rVd(1) = (Sqrt(2.5)*N*(-6*(-1 + 2*Msq + Nsq)*(-1 + 5*Nsq)* &
         u_V_FFS + (11 - 48*Nsq + 45*N**4 + Msq*(-26 + 90*Nsq))* &
         u_V_FFP - 2*u_V_FFD + 12*Nsq*u_V_FFD +  &
        Msq*(20*u_V_FFD - 6*u_V_FFF) - 3*u_V_FFF +  &
        6*Msq*Nsq*(-6*u_V_FFD + u_V_FFF) +  &
        3*N**4*(-6*u_V_FFD + u_V_FFF)))/8.
              rVd(2) = (Sqrt(2.5)*L*M*N*((6 - 30*Nsq)*u_V_FFS +  &
        (-13 + 45*Nsq)*u_V_FFP + 10*u_V_FFD - 3*u_V_FFF +  &
        3*Nsq*(-6*u_V_FFD + u_V_FFF)))/2.
              rVd(3) = (Sqrt(2.5)*L*(-6*(1 - 18*Nsq + 25*N**4 + Msq*(-2 + 30*Nsq))* &
         u_V_FFS + 11*u_V_FFP - 2*u_V_FFD +  &
        36*Nsq*(-4*u_V_FFP + u_V_FFD) +  &
        Msq*(-26*u_V_FFP + 20*u_V_FFD - 6*u_V_FFF) - 3*u_V_FFF +  &
        18*Msq*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) +  &
        15*N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF)))/8.
              rVd_t = (Sqrt(2.5)*L*N*(-6*(-1 + 2*Msq + Nsq)*(-1 + 5*Nsq)*u_Vd_FFS +  &
        (11 - 48*Nsq + 45*N**4 + Msq*(-26 + 90*Nsq))*u_Vd_FFP -  &
        2*u_Vd_FFD + 12*Nsq*u_Vd_FFD +  &
        Msq*(20*u_Vd_FFD - 6*u_Vd_FFF) - 3*u_Vd_FFF +  &
        6*Msq*Nsq*(-6*u_Vd_FFD + u_Vd_FFF) +  &
        3*N**4*(-6*u_Vd_FFD + u_Vd_FFF)))/8.
              case (6)
              rVd(1) = 15*L*(Lsq - Msq)*Nsq*u_V_FFS
              rVd(2) = M*(-1 + 2*Msq + Nsq)* &
    (5*u_V_FFP - 8*u_V_FFD +  &
      3*Nsq*(5*u_V_FFS - 15*u_V_FFP + 6*u_V_FFD - u_V_FFF) +  &
      3*u_V_FFF)
              rVd(3) = (N*(30*(-1 + 2*Msq + Nsq)**2*u_V_FFS - 35*u_V_FFP +  &
        2*u_V_FFD + 6*Nsq*(25*u_V_FFP - 4*u_V_FFD - u_V_FFF) +  &
        3*u_V_FFF - 12*M**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        24*Msq*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        9*N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) +  &
        8*Msq*(25*u_V_FFP - 13*u_V_FFD + 3*u_V_FFF)))/4.
              rVd_t = (30*(Lsq - Msq)**2*Nsq*u_Vd_FFS -  &
      5*((1 - 3*Nsq)**2*(-1 + Nsq) + 4*M**4*(-1 + 9*Nsq) +  &
         4*Msq*(1 - 10*Nsq + 9*N**4))*u_Vd_FFP +  &
      2*(Nsq*(1 - 3*Nsq)**2 + 4*M**4*(-4 + 9*Nsq) +  &
         4*Msq*(4 - 13*Nsq + 9*N**4))*u_Vd_FFD -  &
      3*(-1 + Nsq)*(4*M**4 + 4*Msq*(-1 + Nsq) + (1 + Nsq)**2)*u_Vd_FFF)/ &
    8.
              case (7)
              rVd(1) = (Sqrt(1.5)*N*(10*(8*M**4 + 6*Msq*(-1 + Nsq) +  &
           (-1 + Nsq)**2)*u_V_FFS - 5*u_V_FFP - 2*u_V_FFD +  &
        4*Nsq*(5*u_V_FFP + u_V_FFD - u_V_FFF) - 3*u_V_FFF +  &
        6*Msq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        8*M**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        6*Msq*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF)))/8.
              rVd(2) = (Sqrt(1.5)*L*M*N*(-3 + 8*Msq + 3*Nsq)* &
      (10*u_V_FFS - 15*u_V_FFP + 6*u_V_FFD - u_V_FFF))/2.
              rVd(3) = (Sqrt(1.5)*L*(10*(1 + 8*M**4 - 6*Nsq + 5*N**4 +  &
           6*Msq*(-1 + 3*Nsq))*u_V_FFS - 5*u_V_FFP - 2*u_V_FFD +  &
        12*Nsq*(5*u_V_FFP + u_V_FFD - u_V_FFF) - 3*u_V_FFF +  &
        6*Msq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        8*M**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        18*Msq*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        5*N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF)))/8.
              rVd_t = (Sqrt(1.5)*L*N*(10*(8*M**4 + 6*Msq*(-1 + Nsq) + (-1 + Nsq)**2)* &
         u_Vd_FFS - 5*u_Vd_FFP - 2*u_Vd_FFD +  &
        4*Nsq*(5*u_Vd_FFP + u_Vd_FFD - u_Vd_FFF) - 3*u_Vd_FFF +  &
        6*Msq*(15*u_Vd_FFP - 6*u_Vd_FFD + u_Vd_FFF) -  &
        8*M**4*(15*u_Vd_FFP - 6*u_Vd_FFD + u_Vd_FFF) -  &
        6*Msq*Nsq*(15*u_Vd_FFP - 6*u_Vd_FFD + u_Vd_FFF) -  &
        N**4*(15*u_Vd_FFP - 6*u_Vd_FFD + u_Vd_FFF)))/8.
            end select
          case (7)
            select case (orb_dir_j)
              case (1)
              rVd(1) = (3*M*(5*L**4 - 10*Lsq*Msq + M**4)* &
      (10*u_V_FFS - 15*u_V_FFP + 6*u_V_FFD - u_V_FFF))/16.
              rVd(2) = (3*L*(L**4 - 10*Lsq*Msq + 5*M**4)* &
      (10*u_V_FFS - 15*u_V_FFP + 6*u_V_FFD - u_V_FFF))/16.
              rVd(3) = 0
              rVd_t = (L*M*(3*L**4 - 10*Lsq*Msq + 3*M**4)* &
      (10*u_Vd_FFS - 15*u_Vd_FFP + 6*u_Vd_FFD - u_Vd_FFF))/16.
              case (2)
              rVd(1) = 0
              rVd(2) = (Sqrt(1.5)*N*(10*(20*M**4 + 15*Msq*(-1 + Nsq) +  &
           (-1 + Nsq)**2)*u_V_FFS - 20*u_V_FFP + 10*u_V_FFD +  &
        15*Msq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        20*M**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        15*Msq*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) +  &
        5*Nsq*(7*u_V_FFP - 4*u_V_FFD + u_V_FFF)))/4.
              rVd(3) = (Sqrt(1.5)*M*(10*(1 + 4*M**4 - 6*Nsq + 5*N**4 +  &
           5*Msq*(-1 + 3*Nsq))*u_V_FFS - 20*u_V_FFP + 10*u_V_FFD +  &
        5*Msq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        4*M**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        15*Msq*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        5*N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) +  &
        15*Nsq*(7*u_V_FFP - 4*u_V_FFD + u_V_FFF)))/4.
              rVd_t = (Sqrt(1.5)*M*N*(10*(4*M**4 + 5*Msq*(-1 + Nsq) + (-1 + Nsq)**2)* &
         u_Vd_FFS - 20*u_Vd_FFP + 10*u_Vd_FFD +  &
        5*Msq*(15*u_Vd_FFP - 6*u_Vd_FFD + u_Vd_FFF) -  &
        4*M**4*(15*u_Vd_FFP - 6*u_Vd_FFD + u_Vd_FFF) -  &
        5*Msq*Nsq*(15*u_Vd_FFP - 6*u_Vd_FFD + u_Vd_FFF) -  &
        N**4*(15*u_Vd_FFP - 6*u_Vd_FFD + u_Vd_FFF) +  &
        5*Nsq*(7*u_Vd_FFP - 4*u_Vd_FFD + u_Vd_FFF)))/4.
              case (3)
              rVd(1) = -(sqrt(15.)*M*(2*(-1 + 4*Msq + Nsq)*(-1 + 5*Nsq)*u_V_FFS -  &
         3*u_V_FFP - 2*u_V_FFD +  &
         4*Msq*(u_V_FFP + 2*u_V_FFD - u_V_FFF) + 3*u_V_FFF -  &
         4*Msq*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
         N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) +  &
         Nsq*(26*u_V_FFP - 20*u_V_FFD + 6*u_V_FFF)))/16.
              rVd(2) = -(sqrt(15.)*L*(2*(-1 + 12*Msq + Nsq)*(-1 + 5*Nsq)* &
          u_V_FFS - 3*u_V_FFP - 2*u_V_FFD +  &
         12*Msq*(u_V_FFP + 2*u_V_FFD - u_V_FFF) + 3*u_V_FFF -  &
         12*Msq*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
         N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) +  &
         Nsq*(26*u_V_FFP - 20*u_V_FFD + 6*u_V_FFF)))/16.
              rVd(3) = -(sqrt(15.)*L*M*N*(2*(-3 + 10*Msq + 5*Nsq)*u_V_FFS +  &
         13*u_V_FFP - 10*u_V_FFD + 3*u_V_FFF -  &
         2*Msq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
         Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF)))/4.
              rVd_t = -(sqrt(15.)*L*M*(2*(-1 + 4*Msq + Nsq)*(-1 + 5*Nsq)*u_Vd_FFS -  &
         3*u_Vd_FFP - 2*u_Vd_FFD +  &
         4*Msq*(u_Vd_FFP + 2*u_Vd_FFD - u_Vd_FFF) + 3*u_Vd_FFF -  &
         4*Msq*Nsq*(15*u_Vd_FFP - 6*u_Vd_FFD + u_Vd_FFF) -  &
         N**4*(15*u_Vd_FFP - 6*u_Vd_FFD + u_Vd_FFF) +  &
         Nsq*(26*u_Vd_FFP - 20*u_Vd_FFD + 6*u_Vd_FFF)))/16.
              case (4)
              rVd(1) = (-3*Sqrt(2.5)*N*(-1 + 2*Msq + Nsq)* &
      (2*(-3 + 5*Nsq)*u_V_FFS + 3*u_V_FFP + 6*u_V_FFD -  &
        3*u_V_FFF - Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF)))/8.
              rVd(2) = (-3*Sqrt(2.5)*L*M*N* &
      (2*(-3 + 5*Nsq)*u_V_FFS + (3 - 15*Nsq)*u_V_FFP +  &
        6*(1 + Nsq)*u_V_FFD - (3 + Nsq)*u_V_FFF))/4.
              rVd(3) = (-3*Sqrt(2.5)*L*(-1 + 4*Msq + Nsq)* &
      (2*(-1 + 5*Nsq)*u_V_FFS + u_V_FFP + 2*u_V_FFD - u_V_FFF -  &
        Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF)))/8.
              rVd_t = (Sqrt(2.5)*L*(Lsq - 3*Msq)*N* &
      (2*(-3 + 5*Nsq)*u_Vd_FFS + (3 - 15*Nsq)*u_Vd_FFP +  &
        6*(1 + Nsq)*u_Vd_FFD - (3 + Nsq)*u_Vd_FFF))/8.
              case (5)
              rVd(1) = (sqrt(15.)*(-6*L*Msq*(-1 + Nsq)* &
         (2*(-1 + 5*Nsq)*u_V_FFS + (1 - 15*Nsq)*u_V_FFP +  &
           2*(1 + 3*Nsq)*u_V_FFD - (1 + Nsq)*u_V_FFF) +  &
        4*L**3*(2*(1 - 6*Nsq + 5*N**4)*u_V_FFS - 2*u_V_FFD +  &
           Nsq*((11 - 15*Nsq)*u_V_FFP + (4 + 6*Nsq)*u_V_FFD -  &
              (3 + Nsq)*u_V_FFF))))/(16.*(-1 + Nsq))
              rVd(2) = (sqrt(15.)*(-6*Lsq*M*(-1 + Nsq)* &
         (2*(-1 + 5*Nsq)*u_V_FFS + (1 - 15*Nsq)*u_V_FFP +  &
           2*(1 + 3*Nsq)*u_V_FFD - (1 + Nsq)*u_V_FFF) +  &
        4*M**3*((-1 + 5*Nsq)*u_V_FFP + u_V_FFF +  &
           Nsq*(-8*u_V_FFD + 3*u_V_FFF))))/(16.*(-1 + Nsq))
              rVd(3) = (sqrt(15.)*N*(10*(-1 + Nsq)* &
         (4*M**4 + 5*Msq*(-1 + Nsq) + (-1 + Nsq)**2)*u_V_FFS +  &
        11*u_V_FFP + 2*u_V_FFD - 3*u_V_FFF +  &
        4*M**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) +  &
        10*Msq*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        4*M**4*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) +  &
        3*N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        5*Msq*N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        N**6*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) +  &
        Nsq*(-41*u_V_FFP + 10*u_V_FFD + u_V_FFF) +  &
        Msq*(-67*u_V_FFP + 14*u_V_FFD + 3*u_V_FFF)))/(8.*(-1 + Nsq))
              rVd_t = (sqrt(15.)*(-3*Lsq*Msq*(-1 + Nsq)* &
         (2*(-1 + 5*Nsq)*u_Vd_FFS + (1 - 15*Nsq)*u_Vd_FFP +  &
           2*(1 + 3*Nsq)*u_Vd_FFD - (1 + Nsq)*u_Vd_FFF) +  &
        M**4*((-1 + 5*Nsq)*u_Vd_FFP + u_Vd_FFF +  &
           Nsq*(-8*u_Vd_FFD + 3*u_Vd_FFF)) +  &
        L**4*(2*(1 - 6*Nsq + 5*N**4)*u_Vd_FFS - 2*u_Vd_FFD +  &
           Nsq*((11 - 15*Nsq)*u_Vd_FFP + (4 + 6*Nsq)*u_Vd_FFD -  &
              (3 + Nsq)*u_Vd_FFF))))/(16.*(-1 + Nsq))
              case (6)
              rVd(1) = (Sqrt(1.5)*N*(10*(8*M**4 + 6*Msq*(-1 + Nsq) +  &
           (-1 + Nsq)**2)*u_V_FFS - 5*u_V_FFP - 2*u_V_FFD +  &
        4*Nsq*(5*u_V_FFP + u_V_FFD - u_V_FFF) - 3*u_V_FFF +  &
        6*Msq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        8*M**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        6*Msq*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF)))/8.
              rVd(2) = (Sqrt(1.5)*L*M*N*(-3 + 8*Msq + 3*Nsq)* &
      (10*u_V_FFS - 15*u_V_FFP + 6*u_V_FFD - u_V_FFF))/2.
              rVd(3) = (Sqrt(1.5)*L*(10*(1 + 8*M**4 - 6*Nsq + 5*N**4 +  &
           6*Msq*(-1 + 3*Nsq))*u_V_FFS - 5*u_V_FFP - 2*u_V_FFD +  &
        12*Nsq*(5*u_V_FFP + u_V_FFD - u_V_FFF) - 3*u_V_FFF +  &
        6*Msq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        8*M**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        18*Msq*Nsq*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF) -  &
        5*N**4*(15*u_V_FFP - 6*u_V_FFD + u_V_FFF)))/8.
              rVd_t = (Sqrt(1.5)*L*N*(10*(8*M**4 + 6*Msq*(-1 + Nsq) + (-1 + Nsq)**2)* &
         u_Vd_FFS - 5*u_Vd_FFP - 2*u_Vd_FFD +  &
        4*Nsq*(5*u_Vd_FFP + u_Vd_FFD - u_Vd_FFF) - 3*u_Vd_FFF +  &
        6*Msq*(15*u_Vd_FFP - 6*u_Vd_FFD + u_Vd_FFF) -  &
        8*M**4*(15*u_Vd_FFP - 6*u_Vd_FFD + u_Vd_FFF) -  &
        6*Msq*Nsq*(15*u_Vd_FFP - 6*u_Vd_FFD + u_Vd_FFF) -  &
        N**4*(15*u_Vd_FFP - 6*u_Vd_FFD + u_Vd_FFF)))/8.
              case (7)
              rVd(1) = (15*L*(2*(8*M**4 + 6*Msq*(-1 + Nsq) + (-1 + Nsq)**2)* &
         u_V_FFS - 3*(8*M**4 + 6*Msq*(-1 + Nsq) + Nsq*(-1 + Nsq))* &
         u_V_FFP))/8.
              rVd(2) = (3*M*(-20*L**4*u_V_FFS + 60*Lsq*Msq*u_V_FFS +  &
        15*(3 + 8*M**4 - 5*Nsq + 2*N**4 + 10*Msq*(-1 + Nsq))*u_V_FFP -  &
        (16*M**4 + 16*Msq*(-1 + Nsq) + 3*(-1 + Nsq)**2)* &
         (6*u_V_FFD - u_V_FFF)))/8.
              rVd(3) = (-3*N*(5*(-1 + Nsq)*u_V_FFP +  &
        (-2 + 48*M**4 + 4*Nsq + 6*N**4 + 36*Msq*(-1 + Nsq))*u_V_FFD -  &
        (3 + 8*M**4 + 4*Nsq + N**4 + 6*Msq*(-1 + Nsq))*u_V_FFF))/8.
              rVd_t = (10*(L**3 - 3*L*Msq)**2*u_Vd_FFS -  &
      (15*((-3*Lsq*M + M**3)**2 + Lsq*(Lsq - 3*Msq)**2*Nsq)*u_Vd_FFP)/ &
       (-1 + Nsq) - 6*(16*M**6 + 24*M**4*(-1 + Nsq) +  &
         9*Msq*(-1 + Nsq)**2 + (-1 + Nsq)*(1 + Nsq)**2)*u_Vd_FFD +  &
      (16*M**6 + 24*M**4*(-1 + Nsq) + 9*Msq*(-1 + Nsq)**2 +  &
         Nsq*(3 + Nsq)**2)*u_Vd_FFF)/16.
            end select
        end select
    end select
end select
