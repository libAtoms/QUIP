program test_angular

	use system_module
	use angular_functions_module

	implicit none

	integer :: l_max, n, i, j
	real(dp) :: x(3, 1), b(-4:4, 0:4, 1, 3)

	CALL system_initialise()
	!CALL enable_timing()

	l_max = 4
	x(1,1) = 8.0
	x(2,1) = 5.0
	x(3,1) = 3.0	

	!CALL RANDOM_NUMBER(x)
	
	!CALL system_timer("timer")
	b = GradSphericalIterative(l_max, x)
	!CALL system_timer("timer")

	call system_finalise()

end program test_angular
