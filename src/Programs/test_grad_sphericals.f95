program test_grad_sphericals

	use system_module
	use angular_functions_module

	implicit none

	integer :: l_max, n, i, j
	real(dp) :: x(3, 1), b(-10:10, 0:10, 1:1, 1:3), y(3), c(0:10, -10:10)

	CALL system_initialise()
	CALL enable_timing()

	l_max = 10

	CALL RANDOM_NUMBER(x)
	
	CALL system_timer("timer")
	do i=1, SIZE(x,2)
		y(1) = x(1,i)
		y(2) = x(2,i)
		y(3) = x(3,i)
		c = SphericalYCartesian_all(l_max, y)
	end do
	CALL system_timer("timer")

	call system_finalise()

end program test_grad_sphericals