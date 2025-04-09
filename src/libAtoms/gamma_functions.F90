module gamma_module

   use system_module,only : dp

   ! Incomplete GAMMA function. Adapted from Numerical Recipes in Fortran 77.

   implicit none
!   integer, parameter :: dp = kind(1.0d0)
   integer, parameter :: max_iteration = 100
   real(dp), parameter :: eps = 1.0e-12_dp
   real(dp), parameter :: min_float = tiny(1.0_dp) * 10.0_dp
   real(dp), parameter ::  log_sqrt_2_pi = log(sqrt(2.0_dp*3.14159265358979323846264338327950288_dp) )

   interface gamma_incomplete_upper
      module procedure gamma_incomplete_upper, gamma_incomplete_upper_xarray1d, &
                       gamma_incomplete_upper_ainteger
   endinterface gamma_incomplete_upper

   contains

   function gamma_incomplete_upper_xarray1d(a,x)
      real(dp), intent(in) :: a
      real(dp), intent(in), dimension(:) :: x
      real(dp), dimension(size(x)) :: gamma_incomplete_upper_xarray1d

      integer :: i

      do i = 1, size(x)
         gamma_incomplete_upper_xarray1d(i) = gamma_incomplete_upper(a,x(i))
      enddo

   endfunction gamma_incomplete_upper_xarray1d

   function gamma_incomplete_upper_ainteger(a,x)
      integer, intent(in) :: a
      real(dp), intent(in) :: x
      real(dp) :: gamma_incomplete_upper_ainteger

      gamma_incomplete_upper_ainteger = gamma_incomplete_upper(real(a,dp),x)

   endfunction gamma_incomplete_upper_ainteger

   function gamma_incomplete_upper(a,x)
      
      real(dp), intent(in) :: a, x
      real(dp) :: gamma_incomplete_upper
      ! USES gcf,gser
      ! Returns the incomplete gamma function Q(a, x) ≡ 1 − P(a, x).
      real(dp) :: gammcf, gamser, gln
      
      if(x < 0.0_dp .or. a <= 0.0_dp ) stop 'bad arguments in gamma_incomplete_upper'

      if(x < a+1.0_dp)then                         !Use the series representation
         call gser(gamser,a,x,gln)
         gamma_incomplete_upper = 1.0_dp - gamser  !and take its complement.
      else                                         !Use the continued fraction representation.
         call gcf(gammcf,a,x,gln)
         gamma_incomplete_upper = gammcf
      endif

      gamma_incomplete_upper = gamma_incomplete_upper * exp(gln)

   endfunction gamma_incomplete_upper

   function gammp(a,x)

      real(dp), intent(in) :: a, x
      real(dp) :: gammp
      ! USES gcf,gser
      ! Returns the incomplete gamma function P(a, x).
      real(dp) :: gammcf,gamser,gln

      if(x < 0.0_dp .or. a <= 0.0_dp) stop 'bad arguments in gammp'
      
      if(x < a + 1.0_dp )then              !Use the series representation.
         call gser(gamser,a,x,gln)
         gammp=gamser
      else                                 !Use the continued fraction representation
         call gcf(gammcf,a,x,gln)
         gammp =1.0_dp - gammcf            !and take its complement.
      endif

   endfunction gammp

   function gammq(a,x)
      
      real(dp), intent(in) :: a, x
      real(dp) :: gammq
      ! USES gcf,gser
      ! Returns the incomplete gamma function Q(a, x) ≡ 1 − P(a, x).
      real(dp) :: gammcf, gamser, gln
      
      if(x < 0.0_dp .or. a <= 0.0_dp ) stop 'bad arguments in gammq'

      if(x < a+1.0_dp)then                  !Use the series representation
         call gser(gamser,a,x,gln)
         gammq=1.0_dp - gamser              !and take its complement.
      else                                  !Use the continued fraction representation.
         call gcf(gammcf,a,x,gln)
         gammq=gammcf
      endif

   endfunction gammq

   subroutine gser(gamser,a,x,gln)
      real(dp), intent(in) :: a, x
      real(dp), intent(out) :: gamser, gln
      
      ! USES ln_gamma
      !Returns the incomplete gamma function P(a, x) evaluated by its series
      !representation as
      !gamser. Also returns ln \Gamma(a) as gln.

      integer :: i
      real(dp) :: ap, del, sum

      gln = ln_gamma(a)

      if(x <= 0.0_dp)then
         if(x < 0.0_dp) stop 'x < 0 in gser'
         gamser=0.0_dp
      else
         ap = a
         sum = 1.0_dp / a
         del = sum
         i = 0
         do
            i = i + 1
            ap = ap + 1.0_dp
            del = del * x / ap
            sum = sum + del
            if(abs(del) < abs(sum)*eps) exit
            if(i == max_iteration) stop 'Number of iterations reached maximum &
               number of iterations, stopping. Consider increasing max_iteration.'
         enddo
         gamser = sum * exp(-x + a * log(x) - gln)
      endif

   endsubroutine gser

   subroutine gcf(gammcf,a,x,gln)
      real(dp), intent(in) :: a, x
      real(dp), intent(out) :: gammcf, gln

      ! USES ln_gamma
      ! Returns the incomplete gamma function Q(a, x) evaluated by its continued
      ! fraction representation as gammcf. Also returns ln \Gamma(a) as gln.

      integer :: i
      real(dp) :: an,b,c,d,del,h

      gln = ln_gamma(a)

      b = x + 1.0_dp - a               ! Set up for evaluating continued fraction by modiﬁed
      c = 1.0 / min_float              ! Lentz’s method (5.2) with b0 = 0.
      d = 1.0_dp / b
      h = d
      i = 0
      do                               ! Iterate to convergence.
         i = i + 1
         an = -i * (i-a)
         b = b + 2.0_dp
         d = an * d + b
         
         if(abs(d) < min_float) d = min_float

         c = b + an / c

         if(abs(c) < min_float) c = min_float

         d = 1.0_dp / d
         del = d * c
         h = h * del
         if(abs(del-1.0_dp) < eps) exit
         if(i == max_iteration) stop 'Number of iterations reached maximum &
            number of iterations, stopping. Consider increasing max_iteration.'
      enddo
      gammcf = exp(-x + a*log(x) - gln)*h ! Put factors in front.

   endsubroutine gcf

!   function ln_gamma(xx)
!      real(dp), intent(in) :: xx
!      real(dp) :: ln_gamma
!
!      ! Returns the value ln[ \Gamma(xx)] for xx > 0.
!      integer :: j
!      real(dp) :: ser, tmp, x, y
!
!      real(dp), dimension(6), parameter :: cof = (/76.18009172947146_dp,-86.50532032941677_dp, &
!      24.01409824083091_dp,-1.231739572450155_dp,0.1208650973866179e-2_dp,-0.5395239384953e-5_dp/)
!      real(dp), parameter :: stp = 2.5066282746310005_dp
!
!      x = xx
!      y = x
!      tmp = x + 5.5_dp
!      tmp = (x + 0.5_dp) * log(tmp) - tmp
!      ser = 1.000000000190015_dp
!      do j = 1, 6
!         y = y + 1.0_dp
!         ser = ser + cof(j) / y
!      enddo
!      ln_gamma = tmp + log(stp*ser/x)
!   endfunction ln_gamma

   function ln_gamma(xx)
      real(dp), intent(in) :: xx
      real(dp) :: ln_gamma

      integer :: i
      real(dp) :: z, x, t

      real(dp), parameter :: g = 7.0_dp
      integer, parameter :: n = 9
      real(dp), dimension(n), parameter :: p = &
        (/0.99999999999980993227684700473478_dp, 676.520368121885098567009190444019_dp, &
          -1259.13921672240287047156078755283_dp, 771.3234287776530788486528258894_dp, &
          -176.61502916214059906584551354_dp, 12.507343278686904814458936853_dp, &
          -0.13857109526572011689554707_dp, 9.984369578019570859563e-6_dp, 1.50563273514931155834e-7_dp /)

      !real(dp), parameter :: g = 4.74218750_dp
      !integer, parameter :: n = 15
      !real(dp), dimension(n), parameter :: p = &
      !(/0.99999999999999709182_dp,  57.156235665862923517_dp, -59.597960355475491248_dp, &
      !  14.136097974741747174_dp,  -0.49191381609762019978_dp,  .33994649984811888699e-4_dp, &
      !  0.46523628927048575665e-4_dp,  -0.98374475304879564677e-4_dp,  0.15808870322491248884e-3_dp, &
      !  -0.21026444172410488319e-3_dp,  0.21743961811521264320e-3_dp,  -0.16431810653676389022e-3_dp, &
      !  0.84418223983852743293e-4_dp,   -0.26190838401581408670e-4_dp,     0.36899182659531622704e-5_dp/)
 
      z = xx - 1.0_dp

      x = 0.0_dp
      do i = n, 2, -1
         x = x + p(i)/(z+real(i-1,dp))
      enddo
      x = x + p(1)

      t = z + g + 0.5_dp

      ln_gamma =  log_sqrt_2_pi - t + log(x) + (z+0.5_dp)*log(t)

   endfunction ln_gamma

endmodule gamma_module
