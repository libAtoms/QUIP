program test_UI_k
use libatoms_module
implicit none
   integer :: n, i, j
   real(dp), allocatable :: Vx(:), Vy(:), dF_dl_x(:), dF_dl_y(:), integrand_x(:), integrand_y(:)
   real(dp) :: min_x, max_x, x, l, dx
   real(dp) :: k, beta, T
   type(spline) :: V0_spline, dF_dl_spline
   integer :: n_samples = 500
   real(dp) :: min_V0, min_V0_x, max_V0_x, F, min_F_1, min_F_2, max_F
   real(dp) :: normalization, prob_left, prob_right

   call system_initialise()

   read *, n
   allocate(Vx(n), Vy(n))
   do i=1, n
      read *, Vx(i), Vy(i)
   end do
   read *, min_x, max_x
   read *, k, T
   beta = 1.0_dp/T

   ! calculate potential spline V0
   call initialise(V0_spline, Vx, Vy, (Vy(2)-Vy(1))/(Vx(2)-Vx(1)), (Vy(n)-Vy(n-1))/(Vx(n)-Vx(n-1)))
   min_V0 = 1e38
   dx = (max_x-min_x)/real(n_samples-1,dp)
   do i=1, n_samples
      x = min_x + real(i-1,dp)*dx
      if (spline_value(V0_spline,x) < min_V0) then
	 min_V0_x = x
	 min_V0 = spline_value(V0_spline, x)
      endif
      if (i > 1 .and. i < n_samples) then
	 if (spline_value(V0_spline,x) > spline_value(V0_spline,x-dx) .and. spline_value(V0_spline,x) > spline_value(V0_spline,x+dx)) then
	    max_V0_x = x
	 endif
      endif
      call print ("V0 "// x//" " // spline_value(V0_spline, x) //" " // spline_deriv(V0_spline, x))
   end do

   call print("V0 min max x " // min_V0_x // " " // max_V0_x)

   ! dF/dl(l) = \integral_{-\infty}^{\infty} dx dVr/dl(x,l,k) \exp(-beta (V0(x) + Vr(x,l,k)))
   allocate(dF_dl_x(n_samples), dF_dl_y(n_samples))
   allocate(integrand_x(3*n_samples), integrand_y(3*n_samples))
   do i=1, n_samples
      l = min_x + (max_x-min_x)*real(i-1,dp)/real(n_samples-1,dp)
      dF_dl_x(i) = l

      do j=1, 3*n_samples
	 x = min_x-(max_x-min_x) + 3*(max_x-min_x)*real(j-1,dp)/real(3*n_samples-1,dp)
	 integrand_x(j) = x
	 integrand_y(j) = dVr_dl(x,l,k) * exp(-beta*(spline_value(V0_spline,x)+Vr(x,l,k)))
! call print("dF_dl integrand_y " // (-beta*(spline_value(V0_spline,x)+Vr(x,l,k))) // " " // integrand_y(j))
      end do
      dF_dl_y(i) = TrapezoidIntegral(integrand_x, integrand_y)

      do j=1, 3*n_samples
	 x = min_x-(max_x-min_x) + 3*(max_x-min_x)*real(j-1,dp)/real(3*n_samples-1,dp)
	 integrand_x(j) = x
	 integrand_y(j) = exp(-beta*(spline_value(V0_spline,x)+Vr(x,l,k)))
      end do
      dF_dl_y(i) = dF_dl_y(i) / TrapezoidIntegral(integrand_x, integrand_y)

   end do
   call initialise(dF_dl_spline, dF_dl_x, dF_dl_y, (dF_dl_y(2)-dF_dl_y(1))/(dF_dl_x(2)-dF_dl_x(1)), (dF_dl_y(n)-dF_dl_y(n-1))/(dF_dl_x(n)-dF_dl_x(n-1)))
   do i=1, n_samples
      x = min_x + (max_x-min_x)*real(i-1,dp)/real(n_samples-1,dp)
      call print ("dF_dl "// x//" " // spline_value(dF_dl_spline, x) // " " // spline_nintegrate(dF_dl_spline, min_V0_x, x))
   end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   dx = 3*(max_x-min_x)/real(3*n_samples-1,dp)
   do j=1, 3*n_samples
      x = min_x-(max_x-min_x) + real(j-1,dp)*dx
      integrand_x(j) = x
      F = spline_value(V0_spline, x)
      integrand_y(j) = exp(-beta*F)
   end do
   normalization = TrapezoidIntegral(integrand_x, integrand_y)

   do j=1, 3*n_samples
      x = min_x-(max_x-min_x) + 3*(max_x-min_x)*real(j-1,dp)/real(3*n_samples-1,dp)
      integrand_x(j) = x
      if (x <= max_V0_x) then
	 integrand_y(j) = exp(-beta*spline_value(V0_spline,x))
      else
	 integrand_y(j) = 0.0_dp
      endif
   end do
   prob_left = TrapezoidIntegral(integrand_x, integrand_y)/normalization
   
   do j=1, 3*n_samples
      x = min_x-(max_x-min_x) + 3*(max_x-min_x)*real(j-1,dp)/real(3*n_samples-1,dp)
      integrand_x(j) = x
      if (x > max_V0_x) then
	 integrand_y(j) = exp(-beta*spline_value(V0_spline,x))
      else
	 integrand_y(j) = 0.0_dp
      endif
   end do
   prob_right = TrapezoidIntegral(integrand_x, integrand_y)/normalization
   call print("V0 integrals: left basin probability " // prob_left // " right basin probability " // prob_right // &
      " ratio " // (prob_right/prob_left))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   dx = 3*(max_x-min_x)/real(3*n_samples-1,dp)
   min_F_1 = 1e38_dp
   min_F_2 = 1e38_dp
   do j=1, 3*n_samples
      x = min_x-(max_x-min_x) + real(j-1,dp)*dx
      integrand_x(j) = x
      F = spline_nintegrate(dF_dl_spline, min_V0_x, x)
      if (j > 1 .and. j < 3*n_samples) then
	 if (F > spline_nintegrate(dF_dl_spline, min_V0_x, x-dx) .and. F > spline_nintegrate(dF_dl_spline, min_V0_x, x+dx)) then
	    max_F = F
	 endif
	 if (F < spline_nintegrate(dF_dl_spline, min_V0_x, x-dx) .and. F < spline_nintegrate(dF_dl_spline, min_V0_x, x+dx)) then
	    if (min_F_1 < 1e37_dp) then
	       min_F_2 = F
	    else
	       min_F_1 = F
	    endif
	 endif
      endif
      integrand_y(j) = exp(-beta*F)
   end do
   normalization = TrapezoidIntegral(integrand_x, integrand_y)

   do j=1, 3*n_samples
      x = min_x-(max_x-min_x) + 3*(max_x-min_x)*real(j-1,dp)/real(3*n_samples-1,dp)
      integrand_x(j) = x
      if (x <= max_V0_x) then
	 integrand_y(j) = exp(-beta*spline_nintegrate(dF_dl_spline, min_V0_x, x))
      else
	 integrand_y(j) = 0.0_dp
      endif
   end do
   prob_left = TrapezoidIntegral(integrand_x, integrand_y)/normalization
   
   do j=1, 3*n_samples
      x = min_x-(max_x-min_x) + 3*(max_x-min_x)*real(j-1,dp)/real(3*n_samples-1,dp)
      integrand_x(j) = x
      if (x > max_V0_x) then
	 integrand_y(j) = exp(-beta*spline_nintegrate(dF_dl_spline, min_V0_x, x))
      else
	 integrand_y(j) = 0.0_dp
      endif
   end do
   prob_right = TrapezoidIntegral(integrand_x, integrand_y)/normalization
   call print("F integrals: left basin probability " // prob_left // " right basin probability " // prob_right // &
      " ratio " // (prob_right/prob_left))
   call print("F minima: " // min_F_1 // " " // min_F_2 // " prob ratio " // exp(-beta*(min_F_2-min_F_1)) // " max " // max_F)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   call finalise(V0_spline)
   call finalise(dF_dl_spline)
   call system_finalise()

contains
   function Vr(x,l,k)
      real(dp), intent(in) :: x, l, k
      real(dp) :: Vr

      Vr = k*(x-l)**2
   end function Vr

   function dVr_dl(x,l,k)
      real(dp), intent(in) :: x, l, k
      real(dp) :: dVr_dl

      dVr_dl = -2.0_dp*k*(x-l)
   end function dVr_dl

end program
