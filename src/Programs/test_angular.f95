program test_angular

	use system_module
	use angular_functions_module

	implicit none

	integer :: l_max, n, i, j
	real(dp) :: x(3, 1), b( -3:3, 0:3, 1)

	call system_initialise()
	call enable_timing()

	!x = transpose(reshape((/ 0.1_dp, 0.1_dp, 0.1_dp, 0.2_dp, 0.2_dp, 0.2_dp, 0.3_dp, 0.3_dp, 0.3_dp /), shape(x)))
	x(1,1) = 1.0
	x(2,1) = 1.0
	x(3,1) = 1.0

	l_max = 3
	n = 1

	b = IterativeHarmonics(l_max, x)

	do i=0, l_max
		do j=-i, i
			print *, "l=", i, "m=", j, "Q:", b(j,i,n)
		end do
	end do

	! This is what I'll need to use for benchmarking
	!CALL system_timer("timer")
	!b = IterativeHarmonics(l_max, n, x)
	!CALL system_timer("timer")

	!print *, b

	call system_finalise()

end program test_angular
