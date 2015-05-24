module statistics_module
use system_module, only : dp, PRINT_VERBOSE, PRINT_NORMAL, print, current_verbosity, operator(//)
implicit none
private

public :: mean, variance

public :: mean_var_decorrelated_err
interface mean_var_decorrelated_err
   module procedure mean_var_decorrelated_err_r1
end interface mean_var_decorrelated_err

contains

subroutine mean_var_decorrelated_err_r1(A, m, v, decorrelated_err, decorrelation_t)
   real(dp) :: A(:)
   real(dp), optional :: m, v, decorrelated_err, decorrelation_t

   integer :: bin_i, bin_size, n_bins, N
   real(dp) :: A_mean, A_var, Nr, A_decorrelation_t
   real(dp), allocatable :: binned_means(:)
   real(dp) :: binned_means_var
   real(dp) :: p_d_t, pp_d_t
   logical :: found_max

   N = size(A)
   Nr = N

   A_mean = mean(A)
   if (present(m)) m = A_mean
   if (present(v) .or. present(decorrelation_t) .or. present(decorrelated_err)) A_var = variance(A, A_mean)
   if (present(v)) v = A_var

   found_max = .false.
   p_d_t = -HUGE(1.0_dp)
   pp_d_t = HUGE(1.0_dp)
   if (present(decorrelation_t) .or. present(decorrelated_err)) then
      bin_size = 1
      do while (bin_size <= N)
	 n_bins = N/bin_size
	 allocate(binned_means(n_bins))
	 do bin_i=1, n_bins
	    binned_means(bin_i) = mean(A((bin_i-1)*bin_size+1:(bin_i-1)*bin_size+bin_size))
	 end do
	 binned_means_var = variance(binned_means, mean(binned_means))
	 deallocate(binned_means)
	 A_decorrelation_t = binned_means_var*real(bin_size,dp)/A_var
	 call print("bin_size " // bin_size // " decorrelation_t " // A_decorrelation_t, PRINT_VERBOSE)
	 if (A_decorrelation_t < p_d_t .and. p_d_t > pp_d_t) then
	    if (.not. found_max) then
	       if (present(decorrelation_t)) decorrelation_t = A_decorrelation_t
	       if (present(decorrelated_err)) decorrelated_err = sqrt(binned_means_var/real(n_bins,dp))
	       found_max = .true.
	    endif
	    if (current_verbosity() <= PRINT_NORMAL) exit
	 endif
	 pp_d_t = p_d_t
	 p_d_t = A_decorrelation_t
	 if (bin_size == 1) then
	    bin_size = bin_size * 2
	 else
	    bin_size = bin_size * 1.5
	 endif
      end do
   end if

end subroutine mean_var_decorrelated_err_r1

function mean(A)
   real(dp), intent(in) :: A(:)
   real(dp) :: mean

   mean = sum(A)/real(size(A),dp)
end function mean

function variance(A, A_mean)
   real(dp), intent(in) :: A(:), A_mean
   real(dp) :: variance

   variance = sum((A-A_mean)**2)/real(size(A),dp)
end function variance

end module statistics_module
