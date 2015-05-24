program test_gamma

   use gamma_module

   implicit none
   integer, parameter :: n_a = 5
   integer, parameter :: n_x = 9
   integer :: i, j, k
   real(dp) :: a(n_a), x(n_x), gamma_mathematica(n_a*n_x)

   a = (/ 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp /)
   x = (/ 0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp, 7.0_dp, 8.0_dp /)
   gamma_mathematica = (/ &
   1.00000000000000e0_dp, &
   3.67879441171442e-1_dp, &
   1.35335283236613e-1_dp, &
   4.9787068367864e-2_dp, &
   1.83156388887342e-2_dp, &
   6.73794699908547e-3_dp, &
   2.47875217666636e-3_dp, &
   9.11881965554516e-4_dp, &
   3.35462627902512e-4_dp, &
   1.00000000000000e0_dp, &
   7.35758882342885e-1_dp, &
   4.06005849709838e-1_dp, &
   1.99148273471456e-1_dp, &
   9.15781944436709e-2_dp, &
   4.04276819945128e-2_dp, &
   1.73512652366645e-2_dp, &
   7.29505572443613e-3_dp, &
   3.01916365112261e-3_dp, &
   2.00000000000000e0_dp, &
   1.83939720585721e0_dp, &
   1.35335283236613e0_dp, &
   8.46380162253687e-1_dp, &
   4.76206611107089e-1_dp, &
   2.49304038966162e-1_dp, &
   1.23937608833318e-1_dp, &
   5.92723277610436e-2_dp, &
   2.7507935488006e-2_dp, &
   6.0000000000000e0_dp, &
   5.88607105874308e0_dp, &
   5.14274076299128e0_dp, &
   3.88339133269339e0_dp, &
   2.60082072220025e0_dp, &
   1.59015549178417e0_dp, &
   9.07223296659887e-1_dp, &
   4.9059249746833e-1_dp, &
   2.54280671950104e-1_dp, &
   2.40000000000000e1_dp, &
   2.39121636761438e1_dp, &
   2.27363275837509e1_dp, &
   1.95663178685705e1_dp, &
   1.5092086444317e1_dp, &
   1.05718388415651e1_dp, &
   6.84135600759915e0_dp, &
   4.15179858916971e0_dp, &
   2.39117761168911e0_dp/)

   !print*,'TESTING GAMMA FUNCTION - GFORTRAN ONLY'
   !do i = 1, n_a
   !   print*, a(i), exp(ln_gamma(a(i))), gamma(a(i)), abs(exp(ln_gamma(a(i)))-gamma(a(i)))
   !enddo

   print*,'TESTING INCOMPLETE GAMMA FUNCTION - Mathematica results'
   print'(2a10,2a25,a14)','a','x','Gamma(a,x) Mathematica','Gamma(a,x) Current','Delta'
   k = 0
   do i = 1, n_a
      do j = 1, n_x
         k = k + 1
         print'(2f10.6,2e25.15,e14.5)',a(i),x(j),gamma_mathematica(k),gamma_incomplete_upper(a(i),x(j)), &
            abs(gamma_mathematica(k)-gamma_incomplete_upper(a(i),x(j)))
      enddo
   enddo

endprogram test_gamma
