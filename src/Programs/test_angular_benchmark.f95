program test_angular_benchmark

    use system_module
    use angular_functions_module

    implicit none

    integer :: l_max, n, i, j
	real(dp) :: x(3, 1), y(3), c(0:4, -4:4)

    CALL system_initialise()
    CALL enable_timing()
    CALL random_number(x)

    l_max = 3

    CALL system_timer("timer")
    do i=1, 1
		y(1) = x(1,i)
		y(2) = x(2,i)
		y(3) = x(3,i)
		c = SphericalYCartesian_all(l_max, y)
		print *, c
	end do
	CALL system_timer("timer")

	call system_finalise()

end program test_angular_benchmark



    