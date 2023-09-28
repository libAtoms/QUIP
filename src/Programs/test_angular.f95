program test_angular

	use system_module
	use angular_functions_module

	implicit none

	integer :: l_max, n, i
	real(dp) :: x(3, 100000000), b(-5:5, 0:5, 100000000, 3)

	CALL system_initialise()
	CALL enable_timing()

	l_max = 8
	
	CALL RANDOM_NUMBER(x)
	
	CALL system_timer("timer")
	b = GradSphericalIterative(l_max, x)
	CALL system_timer("timer")

	call system_finalise()

end program test_angular
