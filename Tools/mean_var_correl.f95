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

program mean_var_correl
use libatoms_module
implicit none
  integer :: n_bins, n_data, i, j, bin_i
  real(dp), allocatable :: data(:,:)
  character(len=128), allocatable :: bin_labels(:)
  type(Dictionary) :: cli_params, data_params
  logical :: do_mean, do_var, do_histogram, do_correl, correlation_subtract_mean, do_effective_N
  integer :: histogram_n_bins
  real(dp) :: histogram_min_v, histogram_max_v, histogram_bin_width
  character(len=FIELD_LENGTH) :: infile_name, outfile_name
  character(len=102400) :: myline
  type(inoutput) :: infile, outfile
  real(dp), allocatable :: data_mean(:), data_var(:), data_correl(:,:), data_histogram(:,:,:)
  integer :: reduction_index, other_index, sz, correlation_max_lag, i_lag, n_correl_print, correlation_effective_N_long_lag
  logical :: over_bins, over_time
  real(dp) :: correl_0, correl_mean, correl_std_dev
  real(dp), allocatable :: effective_N(:)

  call system_initialise()

  call initialise(cli_params)
  call param_register(cli_params, "infile", "stdin", infile_name, help_string="input filename")
  call param_register(cli_params, "outfile", "stdout", outfile_name, help_string="output filename")
  call param_register(cli_params, "mean", "F", do_mean, help_string="calculate mean")
  call param_register(cli_params, "variance", "F", do_var, help_string="calculate variance")
  call param_register(cli_params, "effective_N", "F", do_effective_N, help_string="calculate effective N from ratio of initial to long time autocorrelation")
  call param_register(cli_params, "correlation", "F", do_correl, help_string="calculate autocorrelation")
  call param_register(cli_params, "correlation_subtract_mean", "T", correlation_subtract_mean, help_string="subtract mean before calculating autocorrelation")
  call param_register(cli_params, "correlation_max_lag", "1000", correlation_max_lag, help_string="maxmimum lag to compute autocorrelation for")
  call param_register(cli_params, "correlation_effective_N_long_lag", "1001", correlation_effective_N_long_lag, help_string="lag after which autocorrelation is assumed to be in long time regime (for effective_N calculation)")
  call param_register(cli_params, "histogram", "F", do_histogram, help_string="compute histogram")
  call param_register(cli_params, "histogram_n_bins", "10", histogram_n_bins, help_string="number of bins for histogram calculation")
  call param_register(cli_params, "over_bins", "F", over_bins, help_string="do mean/variance/correlation over bins")
  call param_register(cli_params, "over_time", "F", over_time, help_string="do mean/variance/correlation over time")
  if (.not.param_read_args(cli_params)) then
    call print("Usage: "//trim(EXEC_NAME)//" infile=stdin outfile=stdout mean=F variance=F correlation=F effective_N=F", PRINT_ALWAYS)
    call print("        correlation_max_lag=1000 correlation_effective_N_long_lag=1001 over_bins=F over_time=T", PRINT_ALWAYS)
    call system_abort("Unable to parse command line")
  endif
  call finalise(cli_params)

  if (over_bins .and. over_time) &
    call system_abort("specified both over_bins="//over_bins//" over_time="//over_time)

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

  call initialise(infile, infile_name, INPUT)
  myline=read_line(infile) ! comment
  myline=read_line(infile) ! data description
  if (.not.param_read_line(data_params, trim(myline), ignore_unknown=.false.)) &
    call system_abort("Couldn't parse n_bins or n_data in 2nd line of file")

  allocate(bin_labels(n_bins))
  do i=1, n_bins
    bin_labels(i) = read_line(infile)
  end do
  allocate(data(n_bins, n_data))
  do i=1, n_data
    call read_ascii(infile, data(:,i))
  end do

  sz = size(data, other_index)

  allocate(data_mean(sz))
  data_mean = sum(data,reduction_index)/real(size(data,reduction_index),dp)

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


  if (do_correl .or. do_effective_N) then
    if (over_bins) then
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
    else
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

    if (do_effective_N) then
      if (over_bins) then
	if (correlation_effective_N_long_lag+1 < size(data_correl,1)) then
	  allocate(effective_N(size(data_correl,2)))
	  call print("", file=outfile)
	  call print("", file=outfile)
	  call print("# effective_N", file=outfile)
	  do i=1, size(data_correl,2)
	    correl_0 = data_correl(1, i)
	    correl_mean = sum(data_correl(correlation_effective_N_long_lag+1:,i)) / real(size(data_correl,1)-correlation_effective_N_long_lag, dp)
	    correl_std_dev = sqrt( sum((data_correl(correlation_effective_N_long_lag+1:,i)-correl_mean)**2) / real(size(data_correl,1)-correlation_effective_N_long_lag, dp) )
	    if (correl_std_dev > 0.0_dp) then
	      effective_N(i) = (correl_0/correl_std_dev)**2
	    else
	      effective_N(i) = 0
	    endif
	  end do
	endif
      else
	if (correlation_effective_N_long_lag+1 < size(data_correl,2)) then
	  allocate(effective_N(size(data_correl,1)))
	  call print("", file=outfile)
	  call print("", file=outfile)
	  call print("# effective_N", file=outfile)
	  do i=1, size(data_correl,1)
	    correl_0 = data_correl(i, 1)
	    correl_mean = sum(data_correl(i,correlation_effective_N_long_lag+1:)) / real(size(data_correl,2)-correlation_effective_N_long_lag, dp)
	    correl_std_dev = sqrt( sum((data_correl(i,correlation_effective_N_long_lag+1:)-correl_mean)**2) / real(size(data_correl,2)-correlation_effective_N_long_lag, dp) )
	    if (correl_std_dev > 0.0_dp) then
	      effective_N(i) = (correl_0/correl_std_dev)**2
	    else
	      effective_N(i) = 0
	    endif
	  end do
	endif
      endif
    endif
  endif ! do_correl or do_effective_N

  if (do_mean .or. do_var .or. do_effective_N) then
    if (over_bins) then
      myline = "#"
    else
      myline = "# bin_label"
    endif
    if (do_mean) myline = trim(myline) //" mean"
    if (do_var) myline = trim(myline) //" var N"
    if (do_effective_N) myline = trim(myline) // " effective_N effective_decorrel_time"

    call print(trim(myline), file=outfile)

    do i=1, sz
      if (over_bins) then
        myline = ""
      else
        myline = trim(bin_labels(i))
      endif
      if (do_mean) myline = trim(myline) //" " // data_mean(i)
      if (do_var) myline = trim(myline) //" " // data_var(i)//" "//size(data,reduction_index)
      if (do_effective_N) then
	if (effective_N(i) > 0) then
	  myline = trim(myline) // " " // effective_N(i) // " " // (size(data,reduction_index)/effective_N(i))
	else
	  myline = trim(myline) // " " // effective_N(i) // " " // 0
	endif
      endif
      call print(trim(myline), file=outfile)
    end do
  end if

  if (do_histogram) then
    ! data(n_bins, n_data)
    allocate(data_histogram(2, histogram_n_bins, sz))
    data_histogram = 0.0_dp
    do i=1, sz
      if (over_bins) then
	histogram_min_v = minval(data(:,i))
	histogram_max_v = maxval(data(:,i))
      else
	histogram_min_v = minval(data(i,:))
	histogram_max_v = maxval(data(i,:))
      endif
      histogram_bin_width = (histogram_max_v-histogram_min_v)/histogram_n_bins
      if (histogram_bin_width == 0.0_dp) histogram_bin_width = 1.0e-8_dp
      do bin_i=1, histogram_n_bins
	data_histogram(1, bin_i, i) = histogram_min_v+(real(bin_i,dp)-0.5_dp)*histogram_bin_width
      end do
      do j=1, size(data, reduction_index)
	if (over_bins) then
	  bin_i = (data(j,i) - histogram_min_v)/histogram_bin_width
	else
	  bin_i = (data(i,j) - histogram_min_v)/histogram_bin_width
	endif
	bin_i = max(1, bin_i)
	bin_i = min(histogram_n_bins, bin_i)
	data_histogram(2,bin_i,i) = data_histogram(2,bin_i,i) + 1.0_dp/histogram_bin_width
      end do
    end do

    call print("# v p(v)", file=outfile)

    do i=1, sz
      if (over_bins) then
	myline = "#"
      else
	myline = "# "//trim(bin_labels(i))
      endif
      call print(trim(myline), file=outfile)
      do bin_i=1, histogram_n_bins
	call print(data_histogram(:,bin_i,i), file=outfile)
      end do
      call print("", file=outfile)
      call print("", file=outfile)
    end do
  end if

  call finalise(outfile)
  call system_finalise()

end program
