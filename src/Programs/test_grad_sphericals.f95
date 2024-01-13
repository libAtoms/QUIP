program test_grad_sphericals

	use system_module
	use angular_functions_module

	implicit none

	integer :: l_max, i, j, k
	real(dp) :: x(3, 1)
	real(dp), allocatable :: b(:,:,:,:)

	l_max = 4
	allocate(b(-l_max:l_max, 0:l_max, SIZE(x,2), 3))
	x(1,1) = 0.0
	x(2,1) = 0.0
	x(3,1) = 0.0

	CALL system_initialise()

	b = GradSphericalIterative(l_max, x)

	do k=1, SIZE(x,2)
		do i=0, l_max
			do j=-i, i
				print *, b(j,i,k,1), b(j,i,k,2), b(j,i,k,3)
			end do
		end do
	end do

	deallocate(b)
	CALL system_finalise()

end program test_grad_sphericals