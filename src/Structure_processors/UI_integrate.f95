program UI_integrate
use libatoms_module
implicit none
   type(InOutput) :: io
   integer :: N, n_lines
   real(dp), allocatable :: v(:), v_err(:), f(:), f_err(:), f_int(:), f_int_err_sq(:)
   real(dp) :: t_err_sq
   character(len=256), allocatable :: line_array(:)
   integer :: i, i_zero

   call system_initialise(verbosity=PRINT_SILENT)
   call verbosity_push(PRINT_NORMAL)

   call initialise(io, "stdin")
   call read_file(io, line_array, n_lines)
   call finalise(io)
   allocate(v(n_lines+1))
   allocate(v_err(n_lines+1))
   allocate(f(n_lines+1))
   allocate(f_err(n_lines+1))
   allocate(f_int(n_lines+1))
   allocate(f_int_err_sq(n_lines+1))
   N=n_lines
   do i=1, N
      read (unit=line_array(i), fmt=*) v(i), v_err(i), f(i), f_err(i)
   end do
   i_zero = 1
   do i=2, N
      if (f(i-1)*f(i) < 0) then ! force changed sign, must have crossed 0 here
	 v(i+1:n+1) = v(i:n)
	 v_err(i+1:n+1) = v_err(i:n)
	 f(i+1:n+1) = f(i:n)
	 f_err(i+1:N+1) = f_err(i:N)

	 v(i) = v(i-1)*f(i+1)/(f(i+1)-f(i-1)) - v(i+1)*f(i-1)/(f(i+1)-f(i-1))
	 v_err(i) =  Sqrt((f(i+1)**2*v_err(i-1)**2)/(-f(i-1) + f(i+1))**2 + \
	    f_err(i+1)**2*(-((f(i+1)*v(i-1))/(-f(i-1) + f(i+1))**2) + v(i-1)/(-f(i-1) + f(i+1)) + \
	    (f(i-1)*v(i+1))/(-f(i-1) + f(i+1))**2)**2 + f_err(i-1)**2*((f(i+1)*v(i-1))/(-f(i-1) + f(i+1))**2 \
	    - (f(i-1)*v(i+1))/(-f(i-1) + f(i+1))**2 - v(i+1)/(-f(i-1) + f(i+1)))**2 + \
	    (f(i-1)**2*v_err(i+1)**2)/(-f(i-1) + f(i+1))**2)

	 f(i) = 0.0_dp
	 f_err(i) =  Sqrt((f(i+1)**2*f_err(i-1)**2)/(-f(i-1) + f(i+1))**2 + \
	    f_err(i+1)**2*(-((f(i+1)*f(i-1))/(-f(i-1) + f(i+1))**2) + f(i-1)/(-f(i-1) + f(i+1)) + \
	    (f(i-1)*f(i+1))/(-f(i-1) + f(i+1))**2)**2 + f_err(i-1)**2*((f(i+1)*f(i-1))/(-f(i-1) + f(i+1))**2 \
	    - (f(i-1)*f(i+1))/(-f(i-1) + f(i+1))**2 - f(i+1)/(-f(i-1) + f(i+1)))**2 + \
	    (f(i-1)**2*f_err(i+1)**2)/(-f(i-1) + f(i+1))**2)

	 i_zero = i
	 N = N+1
	 exit
      endif
   end do

   f_int(i_zero) = 0.0_dp
   f_int_err_sq(i_zero) = 0.0_dp
   do i=i_zero-1, 1, -1
      f_int(i) = f_int(i+1) - (v(i+1)-v(i))*(f(i+1)+f(i))/2.0_dp
      t_err_sq = ((f(i+1)+f(i))/2.0_dp)**2*v_err(i+1)**2 + ((f(i+1)+f(i))/2.0_dp)**2*v_err(i)**2 + &
	 ((v(i+1)-v(i))/2.0_dp)**2*f_err(i+1)**2 + ((v(i+1)-v(i))/2.0_dp)**2*f_err(i)**2
      f_int_err_sq(i) = f_int_err_sq(i+1) + t_err_sq
   end do
   do i=i_zero, N-1
      f_int(i+1) = f_int(i) + (v(i+1)-v(i))*(f(i+1)+f(i))/2.0_dp
      t_err_sq = ((f(i+1)+f(i))/2.0_dp)**2*v_err(i+1)**2 + ((f(i+1)+f(i))/2.0_dp)**2*v_err(i)**2 + &
	 ((v(i+1)-v(i))/2.0_dp)**2*f_err(i+1)**2 + ((v(i+1)-v(i))/2.0_dp)**2*f_err(i)**2
      f_int_err_sq(i+1) = f_int_err_sq(i) + t_err_sq
   end do

   call print("# col_val err  grad_val err int_val err")
   do i=1, N
      call print(v(i)//" "//v_err(i)//" "//f(i)//" "//f_err(i)//" "//f_int(i)//" "//sqrt(f_int_err_sq(i)))
   end do

   call verbosity_pop()
   call system_finalise()
end program UI_integrate
