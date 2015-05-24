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
  case (ORB_F)
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
          V(1,2) = Sqrt(1.5)/2.
          V(2,1) = -Sqrt(1.5)/2.
          V(2,2) = 0
          case (3)
          V(1,1) = 0
          V(1,2) = 0
          V(2,1) = 0
          V(2,2) = 0
          case (4)
          V(1,1) = 0
          V(1,2) = 0
          V(2,1) = 0
          V(2,2) = 0
          case (5)
          V(1,1) = 0
          V(1,2) = 0
          V(2,1) = 0
          V(2,2) = 0
          case (6)
          V(1,1) = 0
          V(1,2) = cmplx(0,-0.5,dp)*Sqrt(1.5)
          V(2,1) = cmplx(0,-0.5,dp)*Sqrt(1.5)
          V(2,2) = 0
          case (7)
          V(1,1) = cmplx(0,1.5,dp)
          V(1,2) = 0
          V(2,1) = 0
          V(2,2) = cmplx(0,-1.5,dp)
        end select
      case (2)
        select case (orb_dir_j)
          case (1)
          V(1,1) = 0
          V(1,2) = -Sqrt(1.5)/2.
          V(2,1) = Sqrt(1.5)/2.
          V(2,2) = 0
          case (2)
          V(1,1) = 0
          V(1,2) = 0
          V(2,1) = 0
          V(2,2) = 0
          case (3)
          V(1,1) = 0
          V(1,2) = Sqrt(2.5)/2.
          V(2,1) = -Sqrt(2.5)/2.
          V(2,2) = 0
          case (4)
          V(1,1) = 0
          V(1,2) = 0
          V(2,1) = 0
          V(2,2) = 0
          case (5)
          V(1,1) = 0
          V(1,2) = cmplx(0,-0.5,dp)*Sqrt(2.5)
          V(2,1) = cmplx(0,-0.5,dp)*Sqrt(2.5)
          V(2,2) = 0
          case (6)
          V(1,1) = cmplx(0,1,dp)
          V(1,2) = 0
          V(2,1) = 0
          V(2,2) = cmplx(0,-1,dp)
          case (7)
          V(1,1) = 0
          V(1,2) = cmplx(0,-0.5,dp)*Sqrt(1.5)
          V(2,1) = cmplx(0,-0.5,dp)*Sqrt(1.5)
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
          V(1,2) = -Sqrt(2.5)/2.
          V(2,1) = Sqrt(2.5)/2.
          V(2,2) = 0
          case (3)
          V(1,1) = 0
          V(1,2) = 0
          V(2,1) = 0
          V(2,2) = 0
          case (4)
          V(1,1) = 0
          V(1,2) = cmplx(0,-1,dp)*Sqrt(1.5)
          V(2,1) = cmplx(0,-1,dp)*Sqrt(1.5)
          V(2,2) = 0
          case (5)
          V(1,1) = cmplx(0,0.5,dp)
          V(1,2) = 0
          V(2,1) = 0
          V(2,2) = cmplx(0,-0.5,dp)
          case (6)
          V(1,1) = 0
          V(1,2) = cmplx(0,-0.5,dp)*Sqrt(2.5)
          V(2,1) = cmplx(0,-0.5,dp)*Sqrt(2.5)
          V(2,2) = 0
          case (7)
          V(1,1) = 0
          V(1,2) = 0
          V(2,1) = 0
          V(2,2) = 0
        end select
      case (4)
        select case (orb_dir_j)
          case (1)
          V(1,1) = 0
          V(1,2) = 0
          V(2,1) = 0
          V(2,2) = 0
          case (2)
          V(1,1) = 0
          V(1,2) = 0
          V(2,1) = 0
          V(2,2) = 0
          case (3)
          V(1,1) = 0
          V(1,2) = cmplx(0,1,dp)*Sqrt(1.5)
          V(2,1) = cmplx(0,1,dp)*Sqrt(1.5)
          V(2,2) = 0
          case (4)
          V(1,1) = 0
          V(1,2) = 0
          V(2,1) = 0
          V(2,2) = 0
          case (5)
          V(1,1) = 0
          V(1,2) = -Sqrt(1.5)
          V(2,1) = Sqrt(1.5)
          V(2,2) = 0
          case (6)
          V(1,1) = 0
          V(1,2) = 0
          V(2,1) = 0
          V(2,2) = 0
          case (7)
          V(1,1) = 0
          V(1,2) = 0
          V(2,1) = 0
          V(2,2) = 0
        end select
      case (5)
        select case (orb_dir_j)
          case (1)
          V(1,1) = 0
          V(1,2) = 0
          V(2,1) = 0
          V(2,2) = 0
          case (2)
          V(1,1) = 0
          V(1,2) = cmplx(0,0.5,dp)*Sqrt(2.5)
          V(2,1) = cmplx(0,0.5,dp)*Sqrt(2.5)
          V(2,2) = 0
          case (3)
          V(1,1) = cmplx(0,-0.5,dp)
          V(1,2) = 0
          V(2,1) = 0
          V(2,2) = cmplx(0,0.5,dp)
          case (4)
          V(1,1) = 0
          V(1,2) = Sqrt(1.5)
          V(2,1) = -Sqrt(1.5)
          V(2,2) = 0
          case (5)
          V(1,1) = 0
          V(1,2) = 0
          V(2,1) = 0
          V(2,2) = 0
          case (6)
          V(1,1) = 0
          V(1,2) = -Sqrt(2.5)/2.
          V(2,1) = Sqrt(2.5)/2.
          V(2,2) = 0
          case (7)
          V(1,1) = 0
          V(1,2) = 0
          V(2,1) = 0
          V(2,2) = 0
        end select
      case (6)
        select case (orb_dir_j)
          case (1)
          V(1,1) = 0
          V(1,2) = cmplx(0,0.5,dp)*Sqrt(1.5)
          V(2,1) = cmplx(0,0.5,dp)*Sqrt(1.5)
          V(2,2) = 0
          case (2)
          V(1,1) = cmplx(0,-1,dp)
          V(1,2) = 0
          V(2,1) = 0
          V(2,2) = cmplx(0,1,dp)
          case (3)
          V(1,1) = 0
          V(1,2) = cmplx(0,0.5,dp)*Sqrt(2.5)
          V(2,1) = cmplx(0,0.5,dp)*Sqrt(2.5)
          V(2,2) = 0
          case (4)
          V(1,1) = 0
          V(1,2) = 0
          V(2,1) = 0
          V(2,2) = 0
          case (5)
          V(1,1) = 0
          V(1,2) = Sqrt(2.5)/2.
          V(2,1) = -Sqrt(2.5)/2.
          V(2,2) = 0
          case (6)
          V(1,1) = 0
          V(1,2) = 0
          V(2,1) = 0
          V(2,2) = 0
          case (7)
          V(1,1) = 0
          V(1,2) = -Sqrt(1.5)/2.
          V(2,1) = Sqrt(1.5)/2.
          V(2,2) = 0
        end select
      case (7)
        select case (orb_dir_j)
          case (1)
          V(1,1) = cmplx(0,-1.5,dp)
          V(1,2) = 0
          V(2,1) = 0
          V(2,2) = cmplx(0,1.5,dp)
          case (2)
          V(1,1) = 0
          V(1,2) = cmplx(0,0.5,dp)*Sqrt(1.5)
          V(2,1) = cmplx(0,0.5,dp)*Sqrt(1.5)
          V(2,2) = 0
          case (3)
          V(1,1) = 0
          V(1,2) = 0
          V(2,1) = 0
          V(2,2) = 0
          case (4)
          V(1,1) = 0
          V(1,2) = 0
          V(2,1) = 0
          V(2,2) = 0
          case (5)
          V(1,1) = 0
          V(1,2) = 0
          V(2,1) = 0
          V(2,2) = 0
          case (6)
          V(1,1) = 0
          V(1,2) = Sqrt(1.5)/2.
          V(2,1) = -Sqrt(1.5)/2.
          V(2,2) = 0
          case (7)
          V(1,1) = 0
          V(1,2) = 0
          V(2,1) = 0
          V(2,2) = 0
        end select
    end select
end select
