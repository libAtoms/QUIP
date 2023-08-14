program test_angular

	use system_module
	use angular_functions_module

	implicit none

	integer :: l_max, n, i, j
	real(dp) :: x(4, 1), b(-4:4, 0:4, 1)

	CALL system_initialise()

	l_max = 4
	n = 1
	x(1, 1) = 1.0
	x(2, 1) = 5.0
	x(3, 1) = 1.0

	b = IterativeHarmonics(l_max, x)

	do i=0, l_max
		do j=-i, i
			print *, "l=", i, "m=", j, "Q:", b(j,i,n)
		end do
	end do

	call system_finalise()

end program test_angular
