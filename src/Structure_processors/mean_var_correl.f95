! H0 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! H0 X
! H0 X   libAtoms+QUIP: atomistic simulation library
! H0 X
! H0 X   Portions of this code were written by
! H0 X     Albert Bartok-Partay, Silvia Cereda, Gabor Csanyi, James Kermode,
! H0 X     Ivan Solt, Wojciech Szlachta, Csilla Varnai, Steven Winfield.
! H0 X
! H0 X   Copyright 2006-2010.
! H0 X
! H0 X   These portions of the source code are released under the GNU General
! H0 X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html
! H0 X
! H0 X   If you would like to license the source code under different terms,
! H0 X   please contact Gabor Csanyi, gabor@csanyi.net
! H0 X
! H0 X   Portions of this code were written by Noam Bernstein as part of
! H0 X   his employment for the U.S. Government, and are not subject
! H0 X   to copyright in the USA.
! H0 X
! H0 X
! H0 X   When using this software, please cite the following reference:
! H0 X
! H0 X   http://www.libatoms.org
! H0 X
! H0 X  Additional contributions by
! H0 X    Alessio Comisso, Chiara Gattinoni, and Gianpietro Moras
! H0 X
! H0 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

! compute mean and/or variance and/or correlation of a grid of data over one of the two indices
! row index is called "bins", column index is called "time"

module mean_var_correl_util_mod
use system_module
implicit none

contains

subroutine calc_correl(data, n_bins, n_data, over_bins, correlation_max_lag, correlation_subtract_mean, data_correl)
   real(dp) :: data(:,:)
   integer :: n_bins, n_data
   logical :: over_bins
   integer :: correlation_max_lag
   logical :: correlation_subtract_mean
   real(dp), allocatable :: data_correl(:,:)

   real(dp), allocatable :: data_mean(:)

   integer :: i, i_lag

    if (allocated(data_correl)) deallocate(data_correl)
    if (over_bins) then
      allocate(data_mean(size(data,2)))
      data_mean = sum(data,1)/real(size(data,1),dp)
      allocate(data_correl(min(n_bins, correlation_max_lag+1), n_data))
      data_correl = 0.0_dp
      if (correlation_subtract_mean) then
	do i_lag=0, min(n_bins-1, correlation_max_lag)
	  do i=1, n_bins-i_lag
	    data_correl(i_lag+1, :) = data_correl(i_lag+1, :) + (data(i,:)-data_mean(:))*(data(i+i_lag,:)-data_mean(:))
	  end do
	  data_correl(i_lag+1, :) = data_correl(i_lag+1, :) / real(n_bins-i_lag, dp)
	end do
      else
	do i_lag=0, min(n_bins-1, correlation_max_lag)
	  do i=1, n_bins-i_lag
	    data_correl(i_lag+1, :) = data_correl(i_lag+1, :) + (data(i,:))*(data(i+i_lag,:))
	  end do
	  data_correl(i_lag+1, :) = data_correl(i_lag+1, :) / real(n_bins-i_lag, dp)
	end do
      endif
    else ! over bins false
      allocate(data_mean(size(data,1)))
      data_mean = sum(data,2)/real(size(data,2),dp)
      allocate(data_correl(n_bins, min(correlation_max_lag+1,n_data)))
      data_correl = 0.0_dp
      if (correlation_subtract_mean) then
	do i_lag=0, min(n_data-1, correlation_max_lag)
	  do i=1, n_data-i_lag
	    data_correl(:, i_lag+1) = data_correl(:, i_lag+1) + (data(:,i)-data_mean(:))*(data(:,i+i_lag)-data_mean(:))
	  end do
	  data_correl(:, i_lag+1) = data_correl(:, i_lag+1) / real(n_data-i_lag, dp)
	end do
      else
	do i_lag=0, min(n_data-1, correlation_max_lag)
	  do i=1, n_data-i_lag
	    data_correl(:, i_lag+1) = data_correl(:, i_lag+1) + (data(:,i))*(data(:,i+i_lag))
	  end do
	  data_correl(:, i_lag+1) = data_correl(:, i_lag+1) / real(n_data-i_lag, dp)
	end do
      endif
    endif
    deallocate(data_mean)
end subroutine

! effective N_ind samples from autocorrelation variance
subroutine calc_correlation_var_effective_N(data_correl, over_bins, correlation_var_effective_N_long_lag, effective_N)
   real(dp) :: data_correl(:,:)
   logical :: over_bins
   integer :: correlation_var_effective_N_long_lag
   real(dp) :: effective_N(:)

   real(dp) :: correl_0, correl_mean, correl_std_dev
   integer :: i

   if (over_bins) then
     if (correlation_var_effective_N_long_lag+1 < size(data_correl,1)) then
       do i=1, size(data_correl,2)
	 correl_0 = data_correl(1, i)
	 correl_mean = sum(data_correl(correlation_var_effective_N_long_lag+1:,i)) / real(size(data_correl,1)-correlation_var_effective_N_long_lag, dp)
	 correl_std_dev = sqrt( sum((data_correl(correlation_var_effective_N_long_lag+1:,i)-correl_mean)**2) / real(size(data_correl,1)-correlation_var_effective_N_long_lag, dp) )
	 if (correl_std_dev > 1.0e-8_dp) then
	   effective_N(i) = (correl_0/correl_std_dev)**2
	 else
	   effective_N(i) = -1
	 endif
       end do
     endif ! long_lag short enough
   else ! over_bins false
     if (correlation_var_effective_N_long_lag+1 < size(data_correl,2)) then
       do i=1, size(data_correl,1)
	 correl_0 = data_correl(i, 1)
	 correl_mean = sum(data_correl(i,correlation_var_effective_N_long_lag+1:)) / real(size(data_correl,2)-correlation_var_effective_N_long_lag, dp)
	 correl_std_dev = sqrt( sum((data_correl(i,correlation_var_effective_N_long_lag+1:)-correl_mean)**2) / real(size(data_correl,2)-correlation_var_effective_N_long_lag, dp) )
	 if (correl_std_dev > 1.0e-8_dp) then
	   effective_N(i) = (correl_0/correl_std_dev)**2
	 else
	   effective_N(i) = -1
	 endif
	 call print(" i " // i // " correl_0 " // correl_0 // " mean " // correl_mean // " std_dev " // correl_std_dev // " eff_N " // effective_N(i), verbosity=PRINT_VERBOSE)
       end do
     endif ! long_lag short enough
   endif  ! over_bins
end subroutine

! effective N_ind samples from summed ac, a la Sokal
subroutine calc_summed_ac_effective_N(N, data_correl, over_bins, max_lag, effective_N)
   integer :: N
   real(dp) :: data_correl(:,:)
   logical :: over_bins
   integer :: max_lag
   real(dp) :: effective_N(:)

   real(dp) :: correl_0, two_tau_int
   integer :: i

   if (over_bins) then
    do i=1, size(data_correl,2)
      correl_0 = data_correl(1, i)
      two_tau_int = 1.0_dp + 2.0_dp*sum(data_correl(2:max_lag+1,i))/correl_0
      effective_N(i) = real(N,dp)/two_tau_int
      call print("summed_ac i " // i // " correl_0 " // correl_0 // " two_tau_int " // two_tau_int // " effective_N "// effective_N(i), verbosity=PRINT_VERBOSE)
    end do
   else ! over_bins false
    do i=1, size(data_correl,1)
      correl_0 = data_correl(i, 1)
      two_tau_int = 1.0_dp + 2.0_dp*sum(data_correl(i,2:max_lag+1))/correl_0
      effective_N(i) = real(N,dp)/two_tau_int
      call print("summed_ac i " // i // " correl_0 " // correl_0 // " two_tau_int " // two_tau_int // " effective_N "// effective_N(i), verbosity=PRINT_VERBOSE)
    end do
   endif  ! over_bins
end subroutine

function calc_c_hat(data, k) result(c_hat)
   real(dp) :: data(:)
   integer :: k
   real(dp) :: c_hat

   real(dp)  :: sum_0, sum_k, sum_sq_0, sum_sq_k, cross_covar
   real(dp) :: N_minus_k
   integer :: N

   N = size(data)
   sum_0 = sum(data(1:N-k))
   sum_k = sum(data(1+k:N))
   sum_sq_0 = sum(data(1:N-k)**2)
   sum_sq_k = sum(data(1+k:N)**2)
   cross_covar = sum(data(1:N-k)*data(1+k:N))
   N_minus_k = real(N-k,dp)

   c_hat = ( cross_covar/N_minus_k - sum_0*sum_k/N_minus_k**2 ) / &
           ( sqrt(sum_sq_0 - (sum_0**2)/N_minus_k)* &
	     sqrt(sum_sq_k - (sum_k**2)/N_minus_k) / (N_minus_k-1.0_dp) )
end function

! effective N_ind samples from summed ac, a la Kerl
subroutine calc_sliding_window_effective_N(data, over_bins, max_k, effective_N)
   real(dp) :: data(:,:)
   logical :: over_bins
   integer :: max_k
   real(dp) :: effective_N(:)

   real(dp), allocatable :: c_hat(:)
   real(dp) :: two_tau_int
   integer :: i, k

   ! John Kerl Ph.D. Thesis
   allocate(c_hat(1:max_k))

   if (over_bins) then
    do i=1, size(data,2)
      do k=1, max_k
	 c_hat(k) = calc_c_hat(data(:,i),k)
	 call print("c_hat " // k // " " // c_hat(k)//" "//(1.0_dp+2.0_dp*sum(c_hat(1:k))), verbosity=PRINT_VERBOSE)
      end do
      two_tau_int = 1.0_dp + 2.0_dp*sum(c_hat)
      effective_N(i) = real(size(data,1),dp)/two_tau_int
      call print(" i " // i // " two_tau_int " // two_tau_int // " effective_N "// effective_N(i), verbosity=PRINT_VERBOSE)
    end do
   else ! over_bins false
    do i=1, size(data,1)
      do k=1, max_k
	 c_hat(k) = calc_c_hat(data(i,:),k)
	 call print("c_hat " // k // " " // c_hat(k)//" "//(1.0_dp+2.0_dp*sum(c_hat(1:k))), verbosity=PRINT_VERBOSE)
      end do
      two_tau_int = 1.0_dp + 2.0_dp*sum(c_hat)
      effective_N(i) = real(size(data,2),dp)/two_tau_int
      call print(" i " // i // " two_tau_int " // two_tau_int // " effective_N "// effective_N(i), verbosity=PRINT_VERBOSE)
    end do
   endif  ! over_bins

   deallocate(c_hat)
end subroutine

! effective N_ind samples from binning (blocking), a la Hartman
subroutine calc_binning_effective_N(data, over_bins, effective_N)
   real(dp) :: data(:,:)
   logical :: over_bins
   real(dp) :: effective_N(:)

   real(dp) :: bin_size_d
   integer :: i, bin_size, N
   real(dp) :: err_est, prev_err_est, best_err_est, mean, means_variance_1

   if (over_bins) then
      do i=1, size(data,2)
	 N = size(data,1)
	 err_est=0.0_dp
	 best_err_est = -1.0_dp
	 bin_size_d = 1
	 do while (bin_size_d > 0 .and. bin_size_d <= N/5)
	    bin_size = bin_size_d

	    prev_err_est = err_est
	    err_est = binned_err_estimator(data(:,i),bin_size)

	    if (err_est < prev_err_est .and. abs(err_est-prev_err_est)/(0.5_dp*(err_est+prev_err_est)) < 0.1_dp) best_err_est = err_est
	    call print("binning_effective_N " // bin_size//" "//err_est, verbosity=PRINT_VERBOSE)
	    bin_size_d = bin_size_d * 4
	 end do

	 if (best_err_est < 0.0_dp) best_err_est = prev_err_est

	 mean = sum(data(:,i))/real(N,dp)
	 means_variance_1 = sum((data(:,i)-mean)**2)/real(N,dp)
	 effective_N(i) = means_variance_1/best_err_est**2
	 if (effective_N(i) > N) effective_N(i) = N
      end do
   else
      do i=1, size(data,1)
	 N = size(data,2)
	 err_est=0.0_dp
	 best_err_est = -1.0_dp
	 bin_size_d = 1
	 do while (bin_size_d > 0 .and. bin_size_d <= N/5)
	    bin_size = bin_size_d

	    prev_err_est = err_est
	    err_est = binned_err_estimator(data(i,:),bin_size)

	    if (err_est < prev_err_est .and. abs(err_est-prev_err_est)/(0.5_dp*(err_est+prev_err_est)) < 0.1_dp) best_err_est = err_est
	    call print("binning_effective_N " // bin_size//" "//err_est, verbosity=PRINT_VERBOSE)
	    bin_size_d = bin_size_d * 4
	 end do

	 if (best_err_est < 0.0_dp) best_err_est = prev_err_est

	 mean = sum(data(i,:))/real(N,dp)
	 means_variance_1 = sum((data(i,:)-mean)**2)/real(N,dp)
	 effective_N(i) = means_variance_1/best_err_est**2
	 if (effective_N(i) > N) effective_N(i) = N
      end do
   end if
end subroutine

function binned_err_estimator(data, bin_size) result(err_est)
   real(dp) :: data(:)
   integer :: bin_size
   real(dp) :: err_est

   integer :: n_bins, bin_i, N
   real(dp), allocatable :: means(:)
   real(dp) :: means_mean, means_variance

   N = size(data)
   n_bins = N/bin_size
   allocate(means(n_bins))

   do bin_i=1, n_bins
      means(bin_i) = sum(data((bin_i-1)*bin_size+1:bin_i*bin_size))/real(bin_size,dp)
   end do
   means_mean = sum(means)/real(n_bins,dp)
   means_variance = sum((means-means_mean)**2)/real(n_bins,dp)
   err_est = sqrt(means_variance/real(n_bins,dp))

   deallocate(means)
end function

end module mean_var_correl_util_mod

program mean_var_correl
use libatoms_module
use mean_var_correl_util_mod
implicit none
  integer :: n_bins, n_data, n_weights, i, j, bin_i, skip, max_frame
  logical :: do_weights
  real(dp), allocatable :: data(:,:), weights(:), data_line(:)
  character(len=128), allocatable :: bin_labels(:)
  type(Dictionary) :: cli_params, data_params
  logical :: do_mean, do_var, do_histogram, do_correl, correlation_subtract_mean, do_correlation_var_effective_N, do_summed_ac_effective_N, do_sliding_window_effective_N, do_binning_effective_N, do_exp_smoothing
  real(dp) :: exp_smoothing_time
  integer :: exp_smoothing_bin_i
  integer :: histogram_n_bins
  real(dp) :: histogram_min_v, histogram_max_v, histogram_cur_min_v, histogram_cur_max_v, histogram_extra_width, histogram_bin_width
  logical :: do_histogram_effective_N, do_histogram_correl
  character(len=STRING_LENGTH) :: infile_name, outfile_name
  character(len=102400) :: myline
  type(inoutput) :: infile, outfile
  real(dp), allocatable :: data_mean(:), data_var(:), data_correl(:,:), data_histogram(:,:,:), data_histogram_data(:,:,:), data_histogram_correl(:,:)
  real(dp), allocatable :: smooth_data_v(:)
  integer :: reduction_index, other_index, sz, r_sz, correlation_max_lag, n_correl_print, correlation_var_effective_N_long_lag, sliding_window_effective_N_max_k, summed_ac_effective_N_max_lag
  logical :: over_bins, over_time
  real(dp), allocatable :: correlation_var_effective_N(:), summed_ac_effective_N(:), sliding_window_effective_N(:), binning_effective_N(:), histogram_effective_N(:,:)
  character(len=STRING_LENGTH) :: verbosity_str

  call system_initialise()

  call initialise(cli_params)
  call param_register(cli_params, "infile", "stdin", infile_name, help_string="input filename")
  call param_register(cli_params, "outfile", "stdout", outfile_name, help_string="output filename")
  call param_register(cli_params, "skip", "0", skip, help_string="skip this many initial frames")
  call param_register(cli_params, "max_frame", "0", max_frame, help_string="stop at this frame. if 0, continue until end")
  call param_register(cli_params, "mean", "F", do_mean, help_string="calculate mean")
  call param_register(cli_params, "exp_smoothing", "F", do_exp_smoothing, help_string="calculate exponentially smoothed data")
  call param_register(cli_params, "exp_smoothing_time", "100", exp_smoothing_time, help_string="time constant for exponential smoothing (in units of frames). 0 => no smoothing")
  call param_register(cli_params, "exp_smoothing_bin_i", "0", exp_smoothing_bin_i, help_string="bin to evaluate for exponential smoothing")
  call param_register(cli_params, "variance", "F", do_var, help_string="calculate variance")
  call param_register(cli_params, "correlation_var_effective_N", "F", do_correlation_var_effective_N, help_string="calculate effective N from ratio of initial variance to long time autocorrelation variance")
  call param_register(cli_params, "correlation_var_effective_N_long_lag", "1001", &
    correlation_var_effective_N_long_lag, help_string="lag after which autocorrelation is assumed to be in long time regime (for effective_N calculation)")
  call param_register(cli_params, "summed_ac_effective_N", "F", do_summed_ac_effective_N, help_string="calculate Sokals's summed_ac \tau_{int} based effective N from integrated autocorrelation")
  call param_register(cli_params, "summed_ac_effective_N_max_lag", "1001", &
    summed_ac_effective_N_max_lag, help_string="max lag to calculate autocorrelation estimator for summed_ac effective N")
  call param_register(cli_params, "sliding_window_effective_N", "F", do_sliding_window_effective_N, help_string="calculate Kerl's sliding_window \tau_{int} based effective N from fancy integrated autocorrelation")
  call param_register(cli_params, "sliding_window_effective_N_max_k", "1001", &
    sliding_window_effective_N_max_k, help_string="max lag to calculate autocorrelation estimator for sliding window effective N")
  call param_register(cli_params, "binning_effective_N", "F", do_binning_effective_N, help_string="calculate binning effective N")
  call param_register(cli_params, "correlation", "F", do_correl, help_string="calculate autocorrelation")
  call param_register(cli_params, "correlation_subtract_mean", "T", correlation_subtract_mean, help_string="subtract mean before calculating autocorrelation")
  call param_register(cli_params, "correlation_max_lag", "2000", correlation_max_lag, help_string="maxmimum lag to compute autocorrelation for")
  call param_register(cli_params, "histogram", "F", do_histogram, help_string="compute histogram")
  call param_register(cli_params, "histogram_n_bins", "10", histogram_n_bins, help_string="number of bins for histogram calculation")
  call param_register(cli_params, "histogram_extra_width", "0.1", histogram_extra_width, help_string="extra range to do histogram over, in fractions of max_val-min_val (ignored if histogram_min,max_v are used)")
  call param_register(cli_params, "histogram_effective_N", "F", do_histogram_effective_N, help_string="if true calculate effective N for bins based on autocorrelation")
  call param_register(cli_params, "histogram_correl", "F", do_histogram_correl, help_string="if true calculate autocorrelation for histogram bins")
  call param_register(cli_params, "histogram_min_v", "0.0", histogram_min_v, help_string="minimum value to use for histogram range (if not specified, automatic)")
  call param_register(cli_params, "histogram_max_v", "-1.0", histogram_max_v, help_string="maximum value to use for histogram range (if not specified, automatic)")
  call param_register(cli_params, "over_bins", "F", over_bins, help_string="do mean/variance/correlation over bins")
  call param_register(cli_params, "over_time", "F", over_time, help_string="do mean/variance/correlation over time")
  call param_register(cli_params, "verbosity", "NORMAL", verbosity_str, help_string="verbosity level")
  if (.not.param_read_args(cli_params)) then
    call print("Usage: "//trim(EXEC_NAME)//" infile=stdin outfile=stdout mean=F variance=F correlation=F effective_N=F", PRINT_ALWAYS)
    call print("        correlation_max_lag=1000 correlation_var_effective_N_long_lag=1001 over_bins=F over_time=T", PRINT_ALWAYS)
    call system_abort("Unable to parse command line")
  endif
  call finalise(cli_params)

  if (over_bins .and. over_time) &
    call system_abort("specified both over_bins="//over_bins//" over_time="//over_time)

  call verbosity_push(verbosity_of_str(trim(verbosity_str)))

  if (.not. over_bins .and. .not. over_time) then
    over_bins = .false.
    over_time = .true.
  endif

  if (over_bins) then
    reduction_index = 1
    other_index = 2
  else
    reduction_index = 2
    other_index = 1
  endif

  call initialise(data_params)
  call param_register(data_params, "n_bins", param_mandatory, n_bins, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(data_params, "n_data", param_mandatory, n_data, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(data_params, "do_weights", 'F', do_weights, help_string="If true, do weighted mean")

  call initialise(infile, infile_name, INPUT)
  myline=read_line(infile) ! comment
  myline=read_line(infile) ! data description
  if (.not.param_read_line(data_params, trim(myline), ignore_unknown=.false.)) &
    call system_abort("Couldn't parse n_bins or n_data in 2nd line of file")

  if (do_weights) then
    if (over_bins) call system_abort("Can't do weighted mean over bins")
    n_weights = 1
  else
    n_weights = 0
  end if

  allocate(bin_labels(n_bins))
  do i=1, n_bins
    bin_labels(i) = read_line(infile)
  end do
  if(max_frame > 0 .and. max_frame <= n_data) then
     n_data = max_frame
  end if
  allocate(data_line(n_bins+n_weights))
  if(skip > 0) then
     do i=1,skip
        call read_ascii(infile, data_line(:))
     end do
     n_data = n_data-skip
  end if
  allocate(data(n_bins, n_data))
  if (n_weights > 0) then
    allocate(weights(n_data))
  endif
  do i=1, n_data
    call read_ascii(infile, data_line(:))
    if (n_weights > 0) then
        weights(i) = data_line(1)
        data(:,i) = data_line(2:)
    else
        data(:,i) = data_line(:)
    endif
  end do
  deallocate(data_line)


  sz = size(data, other_index)
  r_sz = size(data, reduction_index)

  allocate(data_mean(sz))
  if (n_weights > 0) then
      do i=1, n_data
         data(:,i) = data(:,i) * weights(i)
      end do
      data_mean = sum(data,reduction_index) / sum(weights)
  else
      data_mean = sum(data,reduction_index)/real(size(data,reduction_index),dp)
  endif

  call initialise(outfile, outfile_name, OUTPUT)

  if (do_var) then
    allocate(data_var(sz))
    if (over_bins) then
      do i=1, sz
        data_var(i) = sum((data(:,i)-data_mean(i))**2)/size(data,reduction_index)
      end do
    else 
      do i=1, sz
        data_var(i) = sum((data(i,:)-data_mean(i))**2)/size(data,reduction_index)
      end do
    endif
  endif


  if (do_correl .or. do_correlation_var_effective_N .or. do_summed_ac_effective_N) then
    call calc_correl(data, n_bins, n_data, over_bins, correlation_max_lag, correlation_subtract_mean, data_correl)

    if (do_correl) then
      call print("# correl", file=outfile)
      if (over_bins) then
	n_correl_print = n_data
      else
	n_correl_print = min(n_data, correlation_max_lag+1)
      endif
      do i=1, n_correl_print
	call print(""//data_correl(:,i), file=outfile)
      end do
    endif

    if (do_correlation_var_effective_N) then
      allocate(correlation_var_effective_N(sz))
      call calc_correlation_var_effective_N(data_correl, over_bins, correlation_var_effective_N_long_lag, correlation_var_effective_N)
    endif 
    if (do_summed_ac_effective_N) then
      allocate(summed_ac_effective_N(sz))
      call calc_summed_ac_effective_N(size(data,reduction_index), data_correl, over_bins, summed_ac_effective_N_max_lag, summed_ac_effective_N)
    endif 
  endif ! do_correl or do_effective_N
  if (do_sliding_window_effective_N) then
    allocate(sliding_window_effective_N(sz))
    call calc_sliding_window_effective_N(data, over_bins, sliding_window_effective_N_max_k, sliding_window_effective_N)
  endif 
  if (do_binning_effective_N) then
    allocate(binning_effective_N(sz))
    call calc_binning_effective_N(data, over_bins, binning_effective_N)
  endif 

  ! data(n_bins, n_data)
  if (do_exp_smoothing) then
    if (over_bins) then
      call system_abort("Can't actually do exponentially smoothed time trace over bins")
    else
      myline="# bin_label value"
      if (exp_smoothing_bin_i <= 0 .or. exp_smoothing_bin_i > n_bins) then
	 call system_abort("exp_smoothing_bin_i out of range "//exp_smoothing_bin_i)
      endif
    endif
    call print(trim(myline), file=outfile)
    if (over_bins) then
      call system_abort("Can't actually do exponentially smoothed time trace over bins")
    else
       allocate(smooth_data_v(n_bins))
       smooth_data_v(1:n_bins) = data(1:n_bins, 1)
       call print(trim(bin_labels(exp_smoothing_bin_i))//" "//smooth_data_v(exp_smoothing_bin_i), file=outfile)
       do i=2, n_data
	 if (exp_smoothing_time > 0.0_dp) then
	    smooth_data_v(1:n_bins) = (1.0_dp - 1.0_dp/exp_smoothing_time)*smooth_data_v(1:n_bins) + 1.0_dp/exp_smoothing_time*data(1:n_bins, i)
	 else
	    smooth_data_v(1:n_bins) = data(1:n_bins, i)
	 endif
	 call print(trim(bin_labels(exp_smoothing_bin_i))//" "//smooth_data_v(exp_smoothing_bin_i), file=outfile)
       end do
    endif
  end if ! do_exp_smoothing

  if (do_mean .or. do_var .or. do_correlation_var_effective_N .or. do_summed_ac_effective_N .or. do_sliding_window_effective_N .or. do_binning_effective_N) then
    if (over_bins) then
      myline = "#"
    else
      myline = "# bin_label"
    endif
    if (do_mean) myline = trim(myline) //" mean"
    if (do_var) myline = trim(myline) //" var N"
    if (do_correlation_var_effective_N) myline = trim(myline) // " correlation_var_effective_N correlation_var_effective_decorrel_time"
    if (do_summed_ac_effective_N) myline = trim(myline) // " summed_ac_effective_N summed_ac_effective_decorrel_time"
    if (do_sliding_window_effective_N) myline = trim(myline) // " sliding_window_effective_N sliding_window_effective_decorrel_time"
    if (do_binning_effective_N) myline = trim(myline) // " binning_effective_N binning_effective_decorrel_time"

    call print(trim(myline), file=outfile)

    do i=1, sz
      if (over_bins) then
        myline = ""
      else
        myline = trim(bin_labels(i))
      endif
      if (do_mean) myline = trim(myline) //" " // data_mean(i)
      if (do_var) myline = trim(myline) //" " // data_var(i)//" "//size(data,reduction_index)
      if (do_correlation_var_effective_N) then
	if (correlation_var_effective_N(i) > 0) then
	  myline = trim(myline) // " " // correlation_var_effective_N(i) // " " // (size(data,reduction_index)/correlation_var_effective_N(i))
	else
	  myline = trim(myline) // " " // correlation_var_effective_N(i) // " " // 0
	endif
      endif
      if (do_summed_ac_effective_N) then
	if (summed_ac_effective_N(i) > 0) then
	  myline = trim(myline) // " " // summed_ac_effective_N(i) // " " // (size(data,reduction_index)/summed_ac_effective_N(i))
	else
	  myline = trim(myline) // " " // summed_ac_effective_N(i) // " " // 0
	endif
      endif
      if (do_sliding_window_effective_N) then
	if (sliding_window_effective_N(i) > 0) then
	  myline = trim(myline) // " " // sliding_window_effective_N(i) // " " // (size(data,reduction_index)/sliding_window_effective_N(i))
	else
	  myline = trim(myline) // " " // sliding_window_effective_N(i) // " " // 0
	endif
      endif
      if (do_binning_effective_N) then
	if (binning_effective_N(i) > 0) then
	  myline = trim(myline) // " " // binning_effective_N(i) // " " // (size(data,reduction_index)/binning_effective_N(i))
	else
	  myline = trim(myline) // " " // binning_effective_N(i) // " " // 0
	endif
      endif
      call print(trim(myline), file=outfile)
    end do
  end if

  if (do_histogram) then
    ! data(n_bins, n_data)
    allocate(data_histogram(3, histogram_n_bins, sz))
    if (do_histogram_effective_N) then
      allocate(data_histogram_data(histogram_n_bins, r_sz, sz))
      data_histogram_data = 0.0_dp
    endif
    data_histogram = 0.0_dp
    do i=1, sz
      if (histogram_max_v > histogram_min_v) then
	 histogram_cur_min_v = histogram_min_v
	 histogram_cur_max_v = histogram_max_v
      else
	 if (over_bins) then
	   histogram_cur_min_v = minval(data(:,i))
	   histogram_cur_max_v = maxval(data(:,i))
	 else
	   histogram_cur_min_v = minval(data(i,:))
	   histogram_cur_max_v = maxval(data(i,:))
	 endif
	 histogram_extra_width = (histogram_cur_max_v-histogram_cur_min_v)*histogram_extra_width
	 histogram_cur_min_v = histogram_cur_min_v - histogram_extra_width
	 histogram_cur_max_v = histogram_cur_max_v + histogram_extra_width
      endif
      histogram_bin_width = (histogram_cur_max_v-histogram_cur_min_v)/histogram_n_bins
      if (histogram_bin_width == 0.0_dp) histogram_bin_width = 1.0e-8_dp
      do bin_i=1, histogram_n_bins
	data_histogram(1, bin_i, i) = histogram_cur_min_v+(real(bin_i,dp)-0.5_dp)*histogram_bin_width
      end do
      do j=1, r_sz
	if (over_bins) then
	  bin_i = (data(j,i) - histogram_cur_min_v)/histogram_bin_width+1
	else
	  bin_i = (data(i,j) - histogram_cur_min_v)/histogram_bin_width+1
	endif
	if (bin_i >= 1 .and. bin_i <= histogram_n_bins) then
	  data_histogram(2,bin_i,i) = data_histogram(2,bin_i,i) + 1.0_dp
	  if (do_histogram_effective_N) then
	     data_histogram_data(bin_i, j, i) = 1
	  endif
	endif
      end do ! r_sz
    end do ! i
    data_histogram(2,:,:) = data_histogram(2,:,:)/real(r_sz, dp)
    ! root variance
    data_histogram(3,:,:) = sqrt(data_histogram(2,:,:) - data_histogram(2,:,:)**2)
    ! normalize
    data_histogram(2:3,:,:) = data_histogram(2:3,:,:)/histogram_bin_width
    if (do_histogram_effective_N) then
      allocate(histogram_effective_N(histogram_n_bins,sz))
      do i=1, sz
	 call calc_correl(data_histogram_data(:,:,i), histogram_n_bins, n_data, .false., correlation_max_lag, .true., data_histogram_correl)
	 if (do_histogram_correl) then
	    do bin_i=1, histogram_n_bins
	       call print("# data_histogram_correl i="//i // " bin=" // bin_i)
	       do j=1, size(data_histogram_correl,2)
		  call print(data_histogram_correl(bin_i,j))
	       end do
	       call print("")
	       call print("")
	    end do
	 endif
	 mainlog%prefix="EFF_N value_i="//i
	 call calc_correlation_var_effective_N(data_histogram_correl, .false., correlation_var_effective_N_long_lag, histogram_effective_N(:,i))
	 mainlog%prefix=""
      end do
    endif

    if (do_histogram_effective_N) then
       call print("# v p(v) rms(p(v)) effective_N", file=outfile)
    else
       call print("# v p(v) rms(p(v))", file=outfile)
    endif

    do i=1, sz
      if (over_bins) then
	myline = "#"
      else
	myline = "# "//trim(bin_labels(i))
      endif
      call print(trim(myline), file=outfile)
      do bin_i=1, histogram_n_bins
	if (do_histogram_effective_N) then
	   call print(data_histogram(:,bin_i,i)//" "//histogram_effective_N(bin_i,i), file=outfile)
	else
	   call print(data_histogram(:,bin_i,i), file=outfile)
	endif
      end do
      call print("", file=outfile)
      call print("", file=outfile)
    end do
  end if ! do_histogram

  call finalise(outfile)
  call system_finalise()

end program

