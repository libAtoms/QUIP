! Calculates density in a grid of boxes (possibly smoothed) with error bars
! still very incomplete

program density_new
use libatoms_module
implicit none

  type(Dictionary) :: cli_params
  character(len=FIELD_LENGTH) :: infilename
  logical :: infile_is_list
  character(len=FIELD_LENGTH) :: outfilename
  type(Inoutput) :: outfile
  character(len=FIELD_LENGTH) :: mask
  integer :: decimation
  real(dp) :: min_time, max_time
  logical :: Gaussian_smoothing
  real(dp) :: Gaussian_sigma
  logical :: sort_Time, no_Time_dups
  real(dp) :: min_p(3), bin_width(3)
  integer :: n_bins(3), min_bin(3), max_bin(3)
  logical :: quiet
  logical :: autocorrelation
  integer :: autocorrelation_max_lag

  integer :: n_histos
  type(Atoms_ll) :: structure_ll
  integer, allocatable :: histo_count(:,:,:,:)
  integer :: i, i1, i2, i3
  real(dp), allocatable :: v_autocorr(:,:)

  call system_initialise(NORMAL)

  call initialise(cli_params)
  call param_register(cli_params, 'infile', param_mandatory, infilename)
  call param_register(cli_params, 'infile_is_list', 'F', infile_is_list)
  call param_register(cli_params, 'outfile', 'data.den1', outfilename)

  mask = ""
  call param_register(cli_params, 'AtomMask', "", mask)
  call param_register(cli_params, 'min_p', PARAM_MANDATORY, min_p)
  call param_register(cli_params, 'bin_width', PARAM_MANDATORY, bin_width)
  call param_register(cli_params, 'n_bins', PARAM_MANDATORY, n_bins)
  call param_register(cli_params, 'min_bin', '0 0 0', min_bin)
  call param_register(cli_params, 'max_bin', '0 0 0', max_bin)
  call param_register(cli_params, 'decimation', '1', decimation)
  call param_register(cli_params, 'min_time', '-1.0', min_time)
  call param_register(cli_params, 'max_time', '-1.0', max_time)
  call param_register(cli_params, 'Gaussian', 'F', Gaussian_smoothing)
  call param_register(cli_params, 'sigma', '0.0', Gaussian_sigma)
  call param_register(cli_params, 'sort_Time', 'F', sort_Time)
  call param_register(cli_params, 'no_Time_dups', 'F', no_Time_dups)
  call param_register(cli_params, 'autocorrelation', 'F', autocorrelation)
  call param_register(cli_params, 'autocorrelation_max_lag', '10000', autocorrelation_max_lag)
  call param_register(cli_params, 'quiet', 'F', quiet)
  if (.not. param_read_args(cli_params, do_check = .true.)) then
      call print_usage()
    call system_abort('could not parse argument line')
  end if
  call finalise(cli_params)

  call print("infile " // trim(infilename) // " infile_is_list " // infile_is_list)
  call print("outfilename " // trim(outfilename))
  call print("AtomMask " // trim(mask))
  call print("min_p " // min_p // " bin_width " // bin_width // " n_bins " // n_bins)
  call print("min_bin " // min_bin // " max_bin " // max_bin)
  call print("decimation " // decimation // " min_time " // min_time // " max_time " // max_time)
  call print("Gaussian " // gaussian_smoothing // " sigma " // Gaussian_sigma)
  call print("sort_Time " // sort_Time // " no_Time_dups " // no_Time_dups)
  call print("autocorrelation " // autocorrelation // " autocorrelation_max_lag " // autocorrelation_max_lag)

  call print("Reading configurations")
  call read_xyz(structure_ll, infilename, infile_is_list, decimation, min_time, max_time, sort_Time, no_Time_dups, quiet)

  call print("Calculating densities")
  call calc_histos(histo_count, n_histos, min_p, bin_width, n_bins, structure_ll, 0.0_dp)

  call initialise(outfile, outfilename, OUTPUT)

  if (autocorrelation) then
    call print("Calculating autocorrelations")
    allocate(v_autocorr(n_bins(1),autocorrelation_max_lag))
    v_autocorr = autocorrelation_vector(real(histo_count(:,1,1,:),dp), max_lag=autocorrelation_max_lag)

    call print("Printing autocorrelations")
    call print("# Autocorrelation", file=outfile)
    do i=1, autocorrelation_max_lag
      call print(i // " " // v_autocorr(1:n_bins(1),i), file=outfile)
    end do
  endif

  if (any(min_bin /= 0) .or. any (max_bin /= 0)) then
    if (any(min_bin <= 0) .or. any(min_bin > n_bins) .or. any(max_bin <= 0) .or. any(max_bin > n_bins)) &
      call system_abort ("bin range  out of range min " // min_bin // " max " // max_bin)
    if (autocorrelation) then
      call print("", file=outfile)
      call print("", file=outfile)
    endif
    call print("# Density", file=outfile)
    do i1=min_bin(1), max_bin(1)
    do i2=min_bin(2), max_bin(2)
    do i3=min_bin(3), max_bin(3)
      call print( (min_p+bin_width*(/ i1-0.5_dp, i2-0.5_dp, i3-0.5_dp /)) // " " // sum(histo_count(i1, i2, i3, :))/real(n_histos,dp), file=outfile)
    end do
    end do
    end do
  endif

  call finalise(outfile)

  call system_finalise()

contains

  subroutine print_usage()

    character(len=1024) my_exec_name

    if (EXEC_NAME == "<UNKNOWN>") then
      my_exec_name=""
    else
      my_exec_name=EXEC_NAME
    endif

    call print("Usage: " // trim(my_exec_name)//" infile=filename [infile_is_list=logical(F)]", ERROR)
    call print("       outfile=filename [AtomMask=species()] min_p={x y z}", ERROR)
    call print("       bin_width={x y z} n_bins={nx ny nz} [decimation=n(1)]", ERROR)
    call print("       [min_time=t(-1.0)] [max_time=t(-1.0)] [Gaussian=logical(F)] [sigma=s(1.0)]", ERROR)
    call print("       [sort_Time(F)] [no_Time_dups(F)] autocorrelation_max_lag=N(10000)", ERROR)
  end subroutine print_usage

  function autocorrelation_scalar(v, max_lag) result(autocorrelation)
    real(dp), intent(in) :: v(:)
    integer, intent(in) :: max_lag
    real(dp) :: autocorrelation(max_lag)

    integer :: lag, i
    real(dp) :: t_sum, mean
    integer :: nt

    nt = size(v,1)

    mean = sum(v)/real(nt,dp)
    do lag=1, max_lag
      t_sum = 0.0_dp
      do i=1, nt-lag
	t_sum = t_sum + (v(i)-mean)*(v(i+lag)-mean)
      end do
      autocorrelation(lag) = t_sum/real(nt-lag,dp)
    end do

  end function autocorrelation_scalar

  function autocorrelation_vector(v, max_lag) result(autocorrelation)
    real(dp), intent(in) :: v(:,:)
    integer, intent(in) :: max_lag
    real(dp) :: autocorrelation(size(v,1),max_lag)

    integer :: lag, i
    real(dp), allocatable :: t_sum(:), mean(:)
    integer :: n1, nt

    n1 = size(v,1)
    nt = size(v,2)
    allocate(t_sum(n1), mean(n1))

    mean = sum(v,2)/real(nt,dp)
    do lag=1, max_lag
      t_sum = 0.0_dp
      do i=1, nt-lag
	t_sum(1:n1) = t_sum(1:n1) + (v(1:n1,i)-mean)*(v(1:n1,i+lag)-mean)
      end do
      autocorrelation(1:n1,lag) = t_sum(1:n1)/real(nt-lag,dp)
    end do

    deallocate(t_sum, mean)

  end function autocorrelation_vector

  subroutine calc_histos(histo_count, n_histos, min_p, bin_width, n_bins, structure_ll, lag)
    integer, intent(inout), allocatable :: histo_count(:,:,:,:)
    integer, intent(out) :: n_histos
    real(dp), intent(in) :: min_p(3), bin_width(3)
    integer, intent(in) :: n_bins(3)
    type(atoms_ll), intent(in) :: structure_ll
    real(dp), intent(in) :: lag
    
    real(dp) :: last_time, cur_time
    type(atoms_ll_entry), pointer :: entry
    logical :: do_this_histo

    n_histos = 0
    allocate(histo_count(n_bins(1),n_bins(2),n_bins(3),10))
    histo_count = 0.0_dp

    entry => structure_ll%first
    last_time = -1.0_dp
    do while (associated(entry))
      do_this_histo = .false.
      if (lag > 0.0_dp) then
	if (.not. get_value(entry%at%params, "Time", cur_time)) call system_abort("calc_histos called with lag="//lag//" > 0, but no Time value in entry")
	if ((cur_time .feq. last_time) .or. (cur_time >= last_time+lag)) then
	  do_this_histo = .true.
	  if (cur_time >= last_time+lag) last_time = cur_time
	endif
      else
	do_this_histo = .true.
      endif

      if (do_this_histo) then
	! if (get_value(entry%at%params, "Time", cur_time)) call print("doing histo for config with time " // cur_time)
	n_histos = n_histos + 1
	call reallocate_histos(histo_count, n_histos, n_bins)
	call add_histo_count(histo_count(:,:,:,n_histos), entry%at, min_p, bin_width, n_bins)
      endif
      entry => entry%next
    end do
  end subroutine calc_histos

  subroutine add_histo_count(histo_count, at, min_p, bin_width, n_bins)
    integer, intent(inout) :: histo_count(:,:,:)
    type(Atoms), intent(in) :: at
    real(dp) :: min_p(3), bin_width(3)
    integer :: n_bins(3)

    real(dp) :: min_p_lat(3), bin_width_lat(3), p_lat(3)
    integer :: i, bin(3)

    min_p_lat = at%g .mult. min_p
    bin_width_lat = at%g .mult. bin_width
    do i=1, at%N
      p_lat = at%g .mult. at%pos(:,i)
      bin = floor((p_lat-min_p_lat)/bin_width_lat)+1
      if (all(bin >= 1) .and. all (bin <= n_bins)) histo_count(bin(1),bin(2),bin(3)) = histo_count(bin(1),bin(2),bin(3)) + 1
    end do
  end subroutine add_histo_count

  subroutine reallocate_histos(histos, n, n_bins)
    integer, allocatable, intent(inout) :: histos(:,:,:,:)
    integer, intent(in) :: n, n_bins(3)

    integer, allocatable :: t_histos(:,:,:,:)

    if (allocated(histos)) then
      if (n <= size(histos,4)) return
      allocate(t_histos(size(histos,1),size(histos,2),size(histos,3),size(histos,4)))
      t_histos = histos
      deallocate(histos)
      allocate(histos(size(t_histos,1), size(t_histos,2), size(t_histos,3), 2*size(t_histos,4)))
      histos(1:size(t_histos,1),1:size(t_histos,2),1:size(t_histos,3),1:size(t_histos,4)) = t_histos(1:size(t_histos,1),1:size(t_histos,2),1:size(t_histos,3),1:size(t_histos,4))
      histos(1:size(t_histos,1),1:size(t_histos,2),1:size(t_histos,3),size(t_histos,4)+1:size(histos,4)) = 0.0_dp
      deallocate(t_histos)
    else
      allocate(histos(n_bins(1), n_bins(2), n_bins(3), n))
    endif
  end subroutine reallocate_histos

end program density_new
