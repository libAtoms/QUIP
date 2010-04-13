! compute mean and/or variance and/or correlation of a grid of data over one of the two indices
! row index is called "bins", column index is called "time"

program mean_var_correl
use libatoms_module
implicit none
  integer :: n_bins, n_data, i
  real(dp), allocatable :: data(:,:)
  character(len=128), allocatable :: bin_labels(:)
  type(Dictionary) :: cli_params, data_params
  logical :: do_mean, do_var, do_correl, correlation_subtract_mean, do_effective_N
  character(len=FIELD_LENGTH) :: infile_name, outfile_name
  character(len=102400) :: myline
  type(inoutput) :: infile, outfile
  real(dp), allocatable :: data_mean(:), data_var(:), data_correl(:,:)
  integer :: reduction_index, other_index, sz, correlation_max_lag, i_lag, n_correl_print, correlation_effective_N_long_lag
  logical :: over_bins, over_time
  real(dp) :: correl_0, correl_mean, correl_std_dev
  real(dp), allocatable :: effective_N(:)

  call system_initialise()

  call initialise(cli_params)
  call param_register(cli_params, "infile", "stdin", infile_name)
  call param_register(cli_params, "outfile", "stdout", outfile_name)
  call param_register(cli_params, "mean", "F", do_mean)
  call param_register(cli_params, "variance", "F", do_var)
  call param_register(cli_params, "effective_N", "F", do_effective_N)
  call param_register(cli_params, "correlation", "F", do_correl)
  call param_register(cli_params, "correlation_subtract_mean", "T", correlation_subtract_mean)
  call param_register(cli_params, "correlation_max_lag", "1000", correlation_max_lag)
  call param_register(cli_params, "correlation_effective_N_long_lag", "1001", correlation_effective_N_long_lag)
  call param_register(cli_params, "over_bins", "F", over_bins)
  call param_register(cli_params, "over_time", "F", over_time)
  if (.not.param_read_args(cli_params, do_check=.true.)) then
    call print("Usage: "//trim(EXEC_NAME)//" infile=stdin outfile=stdout mean=F variance=F correlation=F effective_N=F", ERROR)
    call print("        correlation_max_lag=1000 correlation_effective_N_long_lag=1001 over_bins=F over_time=T", ERROR)
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
  call param_register(data_params, "n_bins", param_mandatory, n_bins)
  call param_register(data_params, "n_data", param_mandatory, n_data)

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
    myline=read_line(infile) ! data
    read (unit=myline,fmt=*) data(:,i)
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
      if (effective_N(i) > 0) then
	if (do_effective_N) myline = trim(myline) // " " // effective_N(i) // " " // (size(data,reduction_index)/effective_N(i))
      else
	if (do_effective_N) myline = trim(myline) // " " // effective_N(i) // " " // 0
      endif
      call print(trim(myline), file=outfile)
    end do
  end if

  call finalise(outfile)
  call system_finalise()

end program
