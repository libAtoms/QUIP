module angular_functions_module

use system_module
use units_module
use linearalgebra_module

implicit none
private


   real(dp), dimension(:,:,:,:,:,:), allocatable, save :: cg_array
   integer, save :: cg_j1_max=0, cg_m1_max=0, cg_j2_max=0, cg_m2_max=0, cg_j_max=0, cg_m_max=0 
   logical, save :: cg_initialised = .false.

   public :: SphericalYCartesian, GradSphericalYCartesian
   public :: SphericalYCartesian_all, GradSphericalYCartesian_all
   public :: SolidRCartesian

   public :: wigner3j
   public :: cg_initialise, cg_finalise, cg_array

contains

   !#################################################################################
   !#
   !% Solid Harmonic function using Cartesian coordinates
   !%
   !% $ R_{l m} = \sqrt{\frac{4 \pi}{2 l + 1}} r^l Y_{l m} $
   !#
   !#################################################################################

   function SolidRCartesian(l, m, x)

     complex(dp) :: SolidRCartesian
     integer, intent(in) :: l, m
     real(dp), intent(in) :: x(3)
     integer :: p, q, s

     SolidRCartesian = CPLX_ZERO

     do p = 0, l
        q = p - m
        s = l - p - q

        if ((q >= 0) .and. (s >= 0)) then
           SolidRCartesian = SolidRCartesian + ((cmplx(-0.5_dp * x(1), -0.5_dp * x(2), dp)**p) &
                                             * (cmplx(0.5_dp * x(1), -0.5_dp * x(2), dp)**q) &
                                             * (x(3)**s) &
                                             / (factorial(p) * factorial(q) * factorial(s)))
        end if
     end do

     SolidRCartesian = SolidRCartesian * sqrt(factorial(l + m) * factorial(l - m))

   end function SolidRCartesian

   function SolidRCartesian_all(l_max, x)

     integer, intent(in) :: l_max
     real(dp), intent(in) :: x(3)
     complex(dp) :: SolidRCartesian_all(0:l_max, -l_max:l_max)

     integer :: l, m
     integer :: p, q, s
     complex(kind=dp) :: cm, cp, cm_term(0:l_max), cp_term(0:2*l_max), x3_term(0:2*l_max)
     real(dp) :: factorials(0:2*l_max)

     SolidRCartesian_all = CPLX_ZERO

     cm = cmplx(-0.5_dp * x(1), -0.5_dp * x(2), dp)
     cp = cmplx(0.5_dp * x(1), -0.5_dp * x(2), dp)

     ! p = 0 .. l
     do p=0, l_max
         cm_term(p) = cm**p
     end do
     ! q, s = 0 .. 2l
     do q=0, 2*l_max
        cp_term(q) = cp**q
        x3_term(q) = x(3)**q
        factorials(q) = factorial(q)
     end do

     do l=0, l_max
     do m=-l, l
         do p = 0, l
            q = p - m
            s = l - p - q

            if ((q >= 0) .and. (s >= 0)) then
               SolidRCartesian_all(l,m) = SolidRCartesian_all(l,m) + &
                ( cm_term(p) * cp_term(q) * x3_term(s) / &
                  ( factorials(p) * factorials(q) * factorials(s) ) )
            end if
         end do

         SolidRCartesian_all(l,m) = SolidRCartesian_all(l,m) * sqrt(factorials(l + m) * factorials(l - m))
     end do
     end do

   end function SolidRCartesian_all

   !#################################################################################
   !#
   !% Spherical Harmonic function using Cartesian coordinates
   !#
   !#################################################################################

   function SphericalYCartesian(l, m, x)

     complex(dp) :: SphericalYCartesian
     integer, intent(in) :: l, m
     real(dp), intent(in) :: x(3)

     SphericalYCartesian = SolidRCartesian(l, m, x) * sqrt(((2.0_dp * l) + 1) / (4.0_dp * PI)) &
                                                    * (normsq(x)**(-0.5_dp * l))

   end function SphericalYCartesian

   function SphericalYCartesian_all(l_max, x)

     integer, intent(in) :: l_max
     real(dp), intent(in) :: x(3)
     complex(dp) :: SphericalYCartesian_all(0:l_max, -l_max:l_max)

     real(dp) :: normsq_x
     integer :: l, m

     normsq_x = normsq(x)
     SphericalYCartesian_all = SolidRCartesian_all(l_max, x) 
     do l=0, l_max
         SphericalYCartesian_all(l,-l:l) = SphericalYCartesian_all(l,-l:l) * &
            sqrt(((2.0_dp * l) + 1) / (4.0_dp * PI)) * (normsq_x**(-0.5_dp * l))
     end do

   end function SphericalYCartesian_all

    !#################################################################################
    !#
    !% Derivative of Spherical Harmonic function using Cartesian coordinates
    !#
    !#################################################################################

    function GradSphericalYCartesian(l, m, x)

      complex(dp) :: GradSphericalYCartesian(3)
      integer, intent(in) :: l, m
      real(dp), intent(in) :: x(3)
      integer :: p, q, s

      GradSphericalYCartesian = CPLX_ZERO

      do p = 0, l
         q = p - m
         s = l - p - q

         if ((p >= 1) .and. (q >= 0) .and. (s >= 0)) then
            GradSphericalYCartesian(1) = GradSphericalYCartesian(1) &
                                       - ((cmplx(-0.5_dp * x(1), -0.5_dp * x(2), dp)**(p - 1)) &
                                       * (cmplx(0.5_dp * x(1), -0.5_dp * x(2), dp)**q) &
                                       * (x(3)**s) &
                                       * 0.5_dp &
                                       / (factorial(p - 1) * factorial(q) * factorial(s)))
            GradSphericalYCartesian(2) = GradSphericalYCartesian(2) &
                                       - ((cmplx(-0.5_dp * x(1), -0.5_dp * x(2), dp)**(p - 1)) &
                                       * (cmplx(0.5_dp * x(1), -0.5_dp * x(2), dp)**q) &
                                       * (x(3)**s) &
                                       * 0.5_dp * cmplx(0.0_dp, 1.0_dp, dp) &
                                       / (factorial(p - 1) * factorial(q) * factorial(s)))
         end if

         if ((p >= 0) .and. (q >= 1) .and. (s >= 0)) then
            GradSphericalYCartesian(1) = GradSphericalYCartesian(1) &
                                       + ((cmplx(-0.5_dp * x(1), -0.5_dp * x(2), dp)**p) &
                                       * (cmplx(0.5_dp * x(1), -0.5_dp * x(2), dp)**(q - 1)) &
                                       * (x(3)**s) &
                                       * 0.5_dp &
                                       / (factorial(p) * factorial(q - 1) * factorial(s)))
            GradSphericalYCartesian(2) = GradSphericalYCartesian(2) &
                                       - ((cmplx(-0.5_dp * x(1), -0.5_dp * x(2), dp)**p) &
                                       * (cmplx(0.5_dp * x(1), -0.5_dp * x(2), dp)**(q - 1)) &
                                       * (x(3)**s) &
                                       * 0.5_dp * cmplx(0.0_dp, 1.0_dp, dp) &
                                       / (factorial(p) * factorial(q - 1) * factorial(s)))
         end if

         if ((p >= 0) .and. (q >= 0) .and. (s >= 1)) then
            GradSphericalYCartesian(3) = GradSphericalYCartesian(3) &
                                       + ((cmplx(-0.5_dp * x(1), -0.5_dp * x(2), dp)**p) &
                                       * (cmplx(0.5_dp * x(1), -0.5_dp * x(2), dp)**q) &
                                       * (x(3)**(s - 1)) &
                                       / (factorial(p) * factorial(q) * factorial(s - 1)))
         end if
      end do

      GradSphericalYCartesian = GradSphericalYCartesian &
                              * sqrt(factorial(l + m) * factorial(l - m) * ((2.0_dp * l) + 1) / (4.0_dp * PI)) &
                              * (normsq(x)**(-0.5_dp * l))

      GradSphericalYCartesian = GradSphericalYCartesian &
                              - (l * x * SphericalYCartesian(l, m, x) / normsq(x))

    end function GradSphericalYCartesian

    function GradSphericalYCartesian_all(l_max, x)

      integer, intent(in) :: l_max
      real(dp), intent(in) :: x(3)
      complex(dp) :: GradSphericalYCartesian_all(0:l_max, -l_max:l_max, 3)

      integer :: l, m
      integer :: p, q, s
      complex(kind=dp) :: cm, cp, cm_term(0:l_max), cp_term(0:2*l_max), x3_term(0:2*l_max)
      real(dp) :: factorials(0:2*l_max)
      real(dp) :: normsq_x
      complex(dp) :: tt

      GradSphericalYCartesian_all = CPLX_ZERO
 
      cm = cmplx(-0.5_dp * x(1), -0.5_dp * x(2), dp)
      cp = cmplx(0.5_dp * x(1), -0.5_dp * x(2), dp)
 
      ! p = 0 .. l
      do p=0, l_max
          cm_term(p) = cm**p
      end do
      ! q, s = 0 .. 2l
      do q=0, 2*l_max
         cp_term(q) = cp**q
         x3_term(q) = x(3)**q
         factorials(q) = factorial(q)
      end do

      normsq_x = normsq(x)
 
      do l=0, l_max
      do m=-l, l
          do p = 0, l
             q = p - m
             s = l - p - q

             if ((p >= 1) .and. (q >= 0) .and. (s >= 0)) then
                tt = (cm_term(p-1) * cp_term(q) * x3_term(s) &
                      * 0.5_dp &
                      / (factorials(p - 1) * factorials(q) * factorials(s)))
                GradSphericalYCartesian_all(l,m,1) = GradSphericalYCartesian_all(l,m,1) &
                                           - tt
                GradSphericalYCartesian_all(l,m,2) = GradSphericalYCartesian_all(l,m,2) &
                                           - tt * cmplx(0.0_dp, 1.0_dp, dp)
             end if

             if ((p >= 0) .and. (q >= 1) .and. (s >= 0)) then
                tt = (cm_term(p) * cp_term(q-1) * x3_term(s) &
                      * 0.5_dp &
                      / (factorials(p) * factorials(q - 1) * factorials(s)))
                GradSphericalYCartesian_all(l,m,1) = GradSphericalYCartesian_all(l,m,1) &
                                           + tt
                GradSphericalYCartesian_all(l,m,2) = GradSphericalYCartesian_all(l,m,2) &
                                           - tt * cmplx(0.0_dp, 1.0_dp, dp)
             end if

             if ((p >= 0) .and. (q >= 0) .and. (s >= 1)) then
                GradSphericalYCartesian_all(l,m,3) = GradSphericalYCartesian_all(l,m,3) &
                                           + (cm_term(p) * cp_term(q) * x3_term(s-1) &
                                           / (factorials(p) * factorials(q) * factorials(s - 1)))
             end if
          end do

          GradSphericalYCartesian_all(l,m,:) = GradSphericalYCartesian_all(l,m,:) &
                                  * sqrt(factorials(l + m) * factorials(l - m) * ((2.0_dp * l) + 1) / (4.0_dp * PI)) &
                                  * (normsq_x**(-0.5_dp * l))

          GradSphericalYCartesian_all(l,m,:) = GradSphericalYCartesian_all(l,m,:) &
                                  - (l * x * SphericalYCartesian(l, m, x) / normsq_x)
        end do
        end do

    end function GradSphericalYCartesian_all

   subroutine cg_initialise(j,denom)

      integer, intent(in) :: j
      integer :: i_j1,i_m1,i_j2,i_m2,i_j,i_m
      integer, intent(in), optional :: denom

      integer :: my_denom

      if (cg_initialised .and. j > cg_j_max) then ! need to reinitialise
	 call cg_finalise()
	 cg_initialised = .false.
      endif
      if (cg_initialised) return

      my_denom = optional_default(1,denom)

      cg_j1_max = j
      cg_m1_max = j
      cg_j2_max = j
      cg_m2_max = j
      cg_j_max = j !(j1_max+j2_max)
      cg_m_max = j !(j1_max+j2_max)

      allocate( cg_array(0:cg_j1_max,-cg_m1_max:cg_m1_max,0:cg_j2_max,-cg_m2_max:cg_m2_max,&
      & 0:cg_j_max,-cg_j_max:cg_j_max) )
 
      cg_array = 0.0_dp

      do i_j1 = 0, cg_j1_max
      do i_m1 = -i_j1, i_j1, my_denom
      do i_j2 = 0, cg_j2_max
      do i_m2 = -i_j2, i_j2, my_denom
      do i_j = abs(i_j1-i_j2), min(cg_j_max,i_j1+i_j2)
      do i_m = -i_j, i_j, my_denom


         cg_array(i_j1,i_m1,i_j2,i_m2,i_j,i_m) = cg_calculate(i_j1,i_m1,i_j2,i_m2,i_j,i_m,denom)

      enddo
      enddo
      enddo
      enddo
      enddo
      enddo

      cg_initialised = .true.
    
   endsubroutine cg_initialise

   !#################################################################################
   !#
   !% Finalise global CG arrays
   !#
   !#################################################################################

   subroutine cg_finalise

      cg_j1_max = 0
      cg_m1_max = 0
      cg_j2_max = 0
      cg_m2_max = 0
      cg_j_max = 0
      cg_m_max = 0

      if(allocated(cg_array)) deallocate( cg_array )
      cg_initialised = .false.

   endsubroutine cg_finalise

   !#################################################################################
   !#
   !% Look up CG coefficient from CG array, previously calculated.
   !#
   !#################################################################################

   function cg_lookup(j1,m1,j2,m2,j,m,denom) result(cg)

     real(dp)            :: cg
     integer, intent(in) :: j1,m1,j2,m2,j,m
     integer, intent(in), optional :: denom

     cg=0.0_dp

     if ( .not. cg_check(j1,m1,j2,m2,j,m,denom) ) then
        return
     endif

     if( j1<=cg_j1_max .and. j2<=cg_j2_max .and. j<=cg_j_max .and. &
       abs(m1)<=cg_m1_max .and. abs(m2)<=cg_m2_max .and. abs(m) <= cg_m_max .and. cg_initialised ) then
         cg = cg_array(j1,m1,j2,m2,j,m)
     else
         cg = cg_calculate(j1,m1,j2,m2,j,m,denom)
     endif

   endfunction cg_lookup

   !#################################################################################
   !#
   !% Check if input variables for CG make sense.
   !% Source: http://mathworld.wolfram.com/Clebsch-GordanCoefficient.html \\
   !%                                                                    
   !% $ j_1 + j_2 \ge j $ \\
   !% $ j_1 - j_2 \ge -j $ \\
   !% $ j_1 - j_2 \le j $ \\
   !% $ j_1 \ge m_1 \ge j_1 $ \\
   !% $ j_2 \ge m_2 \ge j_2 $ \\
   !% $ j \ge m \ge j $ \\
   !% $ m_1 + m_2 = m $ \\
   !#
   !#################################################################################

   function cg_check(j1,m1,j2,m2,j,m,denom)

      logical :: cg_check
      integer, intent(in) :: j1,m1,j2,m2,j,m
      integer, intent(in), optional :: denom

      integer :: my_denom

      my_denom = optional_default(1,denom)
      cg_check = (j1>=0) .and. (j2>=0) .and. (j>=0) .and. &
               (abs(m1)<=j1) .and. (abs(m2)<=j2) .and. (abs(m)<=j) &
               .and. (m1+m2==m) .and. (j1+j2 >= j) .and. (abs(j1-j2) <= j) &
               .and. (mod(j1+j2+j,my_denom)==0)
      
   endfunction cg_check

   !#################################################################################
   !#
   !% Calculate a Clebsch-Gordan coefficient $\left< j_1 m_1 j_2 m_2 | j m \right>$
   !% Source: http://mathworld.wolfram.com/Clebsch-GordanCoefficient.html \\
   !% $ \left< j_1 m_1 j_2 m_2 | j m \right> = (-1)^{m+j_1-j_2) 
   !% \sqrt{2j+1} \left( \begin{array}{ccc}
   !% j_1 & j_2 & j \\
   !% m_1 & m_2 & -m \\
   !% \end{array} \right) $
   !% where the thing on the right-hand side is the Wigner 3J symbol.
   !#
   !#################################################################################

   function cg_calculate(j1,m1,j2,m2,j,m,denom) result(cg)

     real(dp)            :: cg
     integer, intent(in) :: j1,m1,j2,m2,j,m
     integer, intent(in), optional :: denom

     integer :: my_denom

     my_denom = optional_default(1,denom)

     cg=0.0_dp
     if ( .not. cg_check(j1,m1,j2,m2,j,m,denom) ) return

     cg = oscillate((m+j1-j2)/my_denom) * sqrt(2.0_dp*real(j,dp)/real(my_denom,dp)+1.0_dp) * &
     wigner3j(j1,m1,j2,m2,j,-m,denom)

   end function cg_calculate

   !#################################################################################
   !#
   !% Triangle coefficient
   !% Source: http://mathworld.wolfram.com/TriangleCoefficient.html
   !% 
   !% $ \Delta(a,b,c) = \frac{ (a+b-c)! (a-b+c)! (-a+b+c)! }{ (a+b+c+1)! } $
   !#
   !#################################################################################

   function tc(a,b,c,denom)

      real(dp) :: tc
      integer, intent(in) :: a, b, c
      integer, intent(in), optional :: denom
      integer :: my_denom

      my_denom = optional_default(1,denom)

      tc = factorial((a+b-c)/my_denom) * factorial((a-b+c)/my_denom) * &
      factorial((-a+b+c)/my_denom) / factorial((a+b+c)/my_denom+1)

   endfunction tc

   !#################################################################################
   !#
   !% Wigner 3J symbol
   !% Source: http://mathworld.wolfram.com/Wigner3j-Symbol.html
   !% 
   !% \[
   !% \left( \begin{array}{ccc}
   !% j_1 & j_2 & j \\
   !% m_1 & m_2 & m \\
   !% \end{array} \right) = (-1)^{j_1-j_2-m) \sqrt{ \Delta (j_1,j_2,j) }
   !% \sqrt{ (j_1+m_1)! (j_1-m_1)! (j_2+m_2)! (j_2-m_2)! (j+m)! (j-m)! }
   !% \sum_k \frac{ (-1)^k }{k! (j-j_2+k+m_1)! (j-j_1+k-m_2)! (j_1+j_2-j-k)!
   !% (j_1-k-m_1)! (j_2-k+m_2)! }
   !% \]
   !% the summation index k runs on all integers where none of the argument of
   !% factorials are negative
   !% $\Delta(a,b,c)$ is the triangle coefficient.
   !#
   !#################################################################################

   function wigner3j(j1,m1,j2,m2,j,m,denom)

       real(dp)            :: wigner3j
       integer, intent(in) :: j1,m1,j2,m2,j,m
       integer, intent(in), optional :: denom

       real(dp) :: pre_fac, triang_coeff, main_coeff, sum_coeff, sum_term
       integer  :: k, kmin, kmax, my_denom

       my_denom = optional_default(1,denom)

       pre_fac = oscillate((j1-j2-m)/my_denom)

       triang_coeff = sqrt( tc(j1,j2,j,denom) )

       main_coeff = sqrt( &
                  factorial((j1+m1)/my_denom) * factorial((j1-m1)/my_denom) * &
                  factorial((j2+m2)/my_denom) * factorial((j2-m2)/my_denom) * &
                  factorial((j+m)/my_denom) * factorial((j-m)/my_denom) )
                  
       sum_coeff = 0.0_dp

       kmin = max( j2-j-m1, j1+m2-j, 0 ) / my_denom
       kmax = min( j1+j2-j, j1-m1, j2+m2) / my_denom

       do k = kmin, kmax

          sum_term = 1.0_dp / ( factorial(k) * factorial((j-j2+m1)/my_denom+k) * &
                   factorial((j-j1-m2)/my_denom+k) * factorial((j1+j2-j)/my_denom-k) * &
                   factorial((j1-m1)/my_denom-k) * factorial((j2+m2)/my_denom-k) )

          sum_term = oscillate(k) * sum_term
          
          sum_coeff = sum_coeff + sum_term

       enddo

       wigner3j = pre_fac * triang_coeff * main_coeff * sum_coeff

   endfunction wigner3j

end module angular_functions_module
