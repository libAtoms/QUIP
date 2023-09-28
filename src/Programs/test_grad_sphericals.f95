
program test_grad_sphericals

	use system_module
	use angular_functions_module

	implicit none

    ! Variables
	integer :: l_max, m
	real(dp) :: x(3)
	complex(dp) :: a, b(3)

	CALL system_initialise()

	l_max = 3
    m = 0
    x(1) = 1.0_dp
    x(2) = 1.0_dp
    x(3) = 1.0_dp
	
	a = SolidRCartesian(l_max, m , x)
	b = GradSphericalYCartesian(l_max, m, x)
	print *, a
    print *, b

	call system_finalise()

end program test_grad_sphericals