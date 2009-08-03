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
  character(len=FIELD_LENGTH) :: commandfilename
  type(Inoutput) :: commandfile
  character(len=1024) :: args_str

  character(len=FIELD_LENGTH) :: mask_str
  integer :: decimation
  real(dp) :: min_time, max_time
  logical :: gaussian_smoothing
  real(dp) :: gaussian_sigma
  logical :: radial_histo
  logical :: sort_Time, no_Time_dups
  real(dp) :: min_p(3), bin_width(3)
  integer :: n_bins(3)
  logical :: quiet
  logical :: autocorrelation, mean
  integer :: autocorrelation_max_lag
  real(dp) :: mean_decorrelation_time
  logical :: have_params
  integer :: status, arg_line_no

  logical :: no_compute_index
  integer :: n_histos
  type(Atoms_ll) :: structure_ll
  real(dp), allocatable :: histo_raw(:,:,:,:), histo_mean(:,:,:), histo_var(:,:,:)
  integer :: i_lag, i1, i2, i3
  real(dp), allocatable :: autocorr(:,:,:,:)
 real(dp) :: density
 logical, allocatable :: mask_a(:)
 integer :: num_atoms

  call system_initialise(NORMAL)

  call initialise(cli_params)
  call register_cli_params(cli_params,.true., infilename, infile_is_list, outfilename, commandfilename, &
    mask_str, min_p, bin_width, n_bins, decimation, min_time, max_time, gaussian_smoothing, gaussian_sigma, radial_histo, &
    sort_Time, no_Time_dups, mean, mean_decorrelation_time, autocorrelation, autocorrelation_max_lag, quiet)
  if (.not. param_read_args(cli_params, do_check = .true.)) then
      call print_usage()
    call system_abort('could not parse argument line')
  end if

  call print("Reading configurations")
  call read_xyz(structure_ll, infilename, infile_is_list, decimation, min_time, max_time, sort_Time, no_Time_dups, quiet, no_compute_index)

  if (len_trim(commandfilename) == 0) then
    args_str = write_string(cli_params)
    call finalise(cli_params)
  else
    call finalise(cli_params)
    call print("reading commands from file '"//trim(commandfilename)//"'")
    call initialise(commandfile, trim(commandfilename), INPUT)
    arg_line_no = 1
    args_str = read_line(commandfile)
    call initialise(cli_params)
    call print("got arguments line '"//trim(args_str)//"'")
    call register_cli_params(cli_params,.false., infilename, infile_is_list, outfilename, commandfilename, &
      mask_str, min_p, bin_width, n_bins, decimation, min_time, max_time, gaussian_smoothing, gaussian_sigma, radial_histo, &
      sort_Time, no_Time_dups, mean, mean_decorrelation_time, autocorrelation, autocorrelation_max_lag, quiet)
    if (.not. param_read_line(cli_params, trim(args_str), ignore_unknown = .true.)) then
      call print_usage()
      call system_abort("could not parse argument line "//arg_line_no//" from file '"//trim(commandfilename//"'"))
    end if
    call finalise(cli_params)
  endif

  have_params = .true.
  do while (have_params)
    if (mean .and. autocorrelation) &
      call system_abort("You probably really don't want to do mean and decorrelation with the same parameters, so I won't let you")

    no_compute_index=.false.
    if ((.not. sort_Time) .and. (.not. no_Time_dups)) no_compute_index=.true.

    call print("infile " // trim(infilename) // " infile_is_list " // infile_is_list)
    call print("outfilename " // trim(outfilename))
    call print("AtomMask " // trim(mask_str))
    call print("min_p " // min_p)
    call print("bin_width " // bin_width // " n_bins " // n_bins)
    call print("decimation " // decimation // " min_time " // min_time // " max_time " // max_time)
    call print("gaussian " // gaussian_smoothing // " sigma " // gaussian_sigma)
    call print("sort_Time " // sort_Time // " no_Time_dups " // no_Time_dups)
    call print("mean " // mean // " mean_decorrelation_time " // mean_decorrelation_time)
    call print("autocorrelation " // autocorrelation // " autocorrelation_max_lag " // autocorrelation_max_lag)

    call print("Calculating densities")
    call calc_histos(histo_raw, n_histos, min_p, bin_width, n_bins, structure_ll, mean_decorrelation_time, gaussian_smoothing, gaussian_sigma, radial_histo, mask_str)

allocate(mask_a(structure_ll%first%at%N))
call is_in_mask(mask_a, structure_ll%first%at, mask_str)
num_atoms = count(mask_a(1:size(mask_a)))
density = real(num_atoms,dp)/cell_volume(structure_ll%first%at)
call print('density = '//density)
deallocate(mask_a)
histo_raw(:,:,:,:) = histo_raw(:,:,:,:) / density

    call initialise(outfile, outfilename, OUTPUT)

    if (autocorrelation) then
      call print("Calculating autocorrelations")
      allocate(autocorr(n_bins(1),n_bins(2),n_bins(3),autocorrelation_max_lag+1))
      autocorr = autocorrelation_3array(histo_raw(:,:,:,1:n_histos), max_lag=autocorrelation_max_lag)

      call print("Printing autocorrelations")
      call print("# Autocorrelation", file=outfile)
      do i_lag=0, autocorrelation_max_lag
	call print(i_lag // " " // reshape(autocorr(1:n_bins(1),1:n_bins(2),1:n_bins(3),i_lag+1), &
					      (/ n_bins(1)*n_bins(2)*n_bins(3) /)),  file=outfile)
      end do
    endif

    if (mean) then
      if (autocorrelation) then
	call print("", file=outfile)
	call print("", file=outfile)
      endif
      allocate(histo_mean(n_bins(1), n_bins(2), n_bins(3)))
      allocate(histo_var(n_bins(1), n_bins(2), n_bins(3)))
      call calc_mean_var_3array(histo_raw(:,:,:,1:n_histos), histo_mean, histo_var)

      if (radial_histo) then
	call print("# Density ", file=outfile)
	call print("# r mean var n_samples", file=outfile)
	do i1=1, n_bins(1)
	  call print((bin_width(1)*(real(i1,dp)-0.5)) // " " // histo_mean(i1, 1, 1) // " " // &
	  histo_var(i1, 1, 1) // " " // n_histos, file=outfile)
	end do
      else
	call print("# Density", file=outfile)
	call print("# x y z   mean   var  n_samples", file=outfile)
	do i1=1,n_bins(1)
	do i2=1,n_bins(2)
	do i3=1,n_bins(3)
	  call print( (min_p+bin_width*(/ i1-0.5_dp, i2-0.5_dp, i3-0.5_dp /)) // &
		      " " // histo_mean(i1, i2, i3) // &
		      " " // histo_var(i1, i2, i3) // " " // n_histos, file=outfile)
	end do
        call print('', file=outfile)
	end do
        call print('', file=outfile)
	end do
      endif
    end if

    call finalise(outfile)

    if (allocated(histo_raw)) deallocate(histo_raw)
    if (allocated(autocorr)) deallocate(autocorr)
    if (allocated(histo_mean)) deallocate(histo_mean)
    if (allocated(histo_var)) deallocate(histo_var)

    if (len_trim(commandfilename) == 0) then
      have_params = .false.
    else
      arg_line_no = arg_line_no + 1
      args_str = read_line(commandfile, status=status)
      if (status == 0) then
	call print("got arguments line '"//trim(args_str)//"'")
	call initialise(cli_params)
	call register_cli_params(cli_params,.false., infilename, infile_is_list, outfilename, commandfilename, &
	  mask_str, min_p, bin_width, n_bins, decimation, min_time, max_time, gaussian_smoothing, gaussian_sigma, radial_histo, &
	  sort_Time, no_Time_dups, mean, mean_decorrelation_time, autocorrelation, autocorrelation_max_lag, quiet)
	if (.not. param_read_line(cli_params, trim(args_str), ignore_unknown = .true.)) then
	  call print_usage()
	  call system_abort("could not parse argument line "//arg_line_no//" from file '"//trim(commandfilename//"'"))
	end if
	call finalise(cli_params)
      else
	have_params = .false.
      endif
    endif

  end do

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
    call print("       outfile=filename [commandfile=filename()] [AtomMask=species()] min_p={x y z}", ERROR)
    call print("       bin_width={x y z} n_bins={nx ny nz} [decimation=n(1)]", ERROR)
    call print("       [min_time=t(-1.0)] [max_time=t(-1.0)] [gaussian=logical(F)] [sigma=s(1.0)] [radial_histo=logical(F)]", ERROR)
    call print("       [sort_Time(F)] [no_Time_dups(F)]", ERROR)
    call print("       [mean=logical(T)] [mean_decorrelation_time=t(0.0)]", ERROR)
    call print("       [autocorrelation=logical(F)] [autocorrelation_max_lag=N(10000)]", ERROR)
    call print("       [quiet=logical(F)]", ERROR)
  end subroutine print_usage

  subroutine register_cli_params(cli_params, initial, infilename, infile_is_list, outfilename, commandfilename, &
    mask_str, min_p, bin_width, n_bins, decimation, min_time, max_time, gaussian_smoothing, gaussian_sigma, radial_histo, &
    sort_Time, no_Time_dups, mean, mean_decorrelation_time, autocorrelation, autocorrelation_max_lag, quiet)
    type(Dictionary), intent(inout) :: cli_params
    logical, intent(in) :: initial
    character(len=*), intent(inout) :: infilename
    logical, intent(inout) :: infile_is_list
    character(len=*), intent(inout) :: outfilename
    character(len=*), intent(inout) :: commandfilename
    character(len=*), intent(inout) :: mask_str
    real(dp), intent(inout) :: min_p(3), bin_width(3)
    integer, intent(inout) :: n_bins(3), decimation
    real(dp), intent(inout) :: min_time, max_time
    logical, intent(inout) :: gaussian_smoothing
    real(dp), intent(inout) :: gaussian_sigma
    logical, intent(inout) :: radial_histo
    logical, intent(inout) :: sort_Time, no_Time_dups
    logical, intent(inout) :: mean
    real(dp), intent(inout) :: mean_decorrelation_time
    logical, intent(inout) :: autocorrelation
    integer, intent(inout) :: autocorrelation_max_lag
    logical, intent(inout) :: quiet
    character(len=FIELD_LENGTH) :: t_field


    if (initial) then
      call param_register(cli_params, 'infile', param_mandatory, infilename)
      call param_register(cli_params, 'infile_is_list', 'F', infile_is_list)
      call param_register(cli_params, 'outfile', 'stdout', outfilename)
      mask_str = ""
      commandfilename = ""
      call param_register(cli_params, 'commandfile', "", commandfilename)
      call param_register(cli_params, 'AtomMask', "", mask_str)
      call param_register(cli_params, 'min_p', PARAM_MANDATORY, min_p)
      call param_register(cli_params, 'bin_width', PARAM_MANDATORY, bin_width)
      call param_register(cli_params, 'n_bins', PARAM_MANDATORY, n_bins)
      call param_register(cli_params, 'decimation', '1', decimation)
      call param_register(cli_params, 'min_time', '-1.0', min_time)
      call param_register(cli_params, 'max_time', '-1.0', max_time)
      call param_register(cli_params, 'gaussian', 'F', gaussian_smoothing)
      call param_register(cli_params, 'sigma', '0.0', gaussian_sigma)
      call param_register(cli_params, 'radial_histo', 'F', radial_histo)
      call param_register(cli_params, 'sort_Time', 'F', sort_Time)
      call param_register(cli_params, 'no_Time_dups', 'F', no_Time_dups)
      call param_register(cli_params, 'mean', 'T', mean)
      call param_register(cli_params, 'mean_decorrelation_time', '0.0', mean_decorrelation_time)
      call param_register(cli_params, 'autocorrelation', 'F', autocorrelation)
      call param_register(cli_params, 'autocorrelation_max_lag', '10000', autocorrelation_max_lag)
      call param_register(cli_params, 'quiet', 'F', quiet)
    else
      t_field=infilename
      call param_register(cli_params, 'infile', trim(t_field), infilename)
      call param_register(cli_params, 'infile_is_list', ""//infile_is_list, infile_is_list)
      t_field=outfilename
      call param_register(cli_params, 'outfile', trim(t_field), outfilename)
      t_field=commandfilename
      call param_register(cli_params, 'commandfile', trim(t_field), commandfilename)
      t_field=mask_str
      call param_register(cli_params, 'AtomMask', trim(t_field), mask_str)
      call param_register(cli_params, 'min_p', ""//min_p, min_p)
      call param_register(cli_params, 'bin_width', ""//bin_width, bin_width)
      call param_register(cli_params, 'n_bins', ""//n_bins, n_bins)
      call param_register(cli_params, 'min_time', ""//min_time, min_time)
      call param_register(cli_params, 'max_time', ""//max_time, max_time)
      call param_register(cli_params, 'gaussian', ""//gaussian_smoothing, gaussian_smoothing)
      call param_register(cli_params, 'sigma', ""//gaussian_sigma, gaussian_sigma)
      call param_register(cli_params, 'radial_histo', ""//radial_histo, radial_histo)
      call param_register(cli_params, 'mean', ""//mean, mean)
      call param_register(cli_params, 'mean_decorrelation_time', ""//mean_decorrelation_time, mean_decorrelation_time)
      call param_register(cli_params, 'autocorrelation', ""//autocorrelation, autocorrelation)
      call param_register(cli_params, 'autocorrelation_max_lag', ""//autocorrelation_max_lag, autocorrelation_max_lag)
    endif
  end subroutine register_cli_params

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

  function autocorrelation_3array(v, max_lag) result(autocorrelation)
    real(dp), intent(in) :: v(:,:,:,:)
    integer, intent(in) :: max_lag
    real(dp) :: autocorrelation(size(v,1),size(v,2),size(v,3),max_lag+1)

    integer :: lag, i
    real(dp), allocatable :: t_sum(:,:,:), mean(:,:,:)
    integer :: n1, n2, n3, nt

    n1 = size(v,1)
    n2 = size(v,2)
    n3 = size(v,3)
    nt = size(v,4)
    allocate(t_sum(n1,n2,n3), mean(n1,n2,n3))

    mean = sum(v,4)/real(nt,dp)
    do lag=0, max_lag
      t_sum = 0.0_dp
      do i=1, nt-lag
	t_sum(1:n1,1:n2,1:n3) = t_sum(1:n1,1:n2,1:n3) + (v(1:n1,1:n2,1:n3,i)-mean)*(v(1:n1,1:n2,1:n3,i+lag)-mean)
      end do
      autocorrelation(1:n1,1:n2,1:n3,lag+1) = t_sum(1:n1,1:n2,1:n3)/real(nt-lag,dp)
    end do

    deallocate(t_sum, mean)

  end function autocorrelation_3array

  subroutine calc_histos(histo_count, n_histos, min_p, bin_width, n_bins, structure_ll, interval, gaussian, gaussian_sigma, radial_histo, mask_str)
    real(dp), intent(inout), allocatable :: histo_count(:,:,:,:)
    integer, intent(out) :: n_histos
    real(dp), intent(in) :: min_p(3), bin_width(3)
    integer, intent(in) :: n_bins(3)
    type(atoms_ll), intent(in) :: structure_ll
    real(dp), intent(in) :: interval
    logical, intent(in) :: gaussian
    real(dp), intent(in) :: gaussian_sigma
    logical, intent(in) :: radial_histo
    character(len=*), optional, intent(in) :: mask_str
    
    real(dp) :: last_time, cur_time
    type(atoms_ll_entry), pointer :: entry
    logical :: do_this_histo

    n_histos = 0
    allocate(histo_count(n_bins(1),n_bins(2),n_bins(3),10))
    histo_count = 0.0_dp

    entry => structure_ll%first
    cur_time = -1.0_dp
    last_time = -1.0e38_dp
    do while (associated(entry))
      do_this_histo = .false.
      if (interval > 0.0_dp) then
	if (cur_time < last_time) &
	  call system_abort("calc_histos called with interval="//interval//" > 0, but unsorted sequence: cur_time="// &
			    cur_time//" < last_time="//last_time)
	if (.not. get_value(entry%at%params, "Time", cur_time)) &
	  call system_abort("calc_histos called with interval="//interval//" > 0, but no Time value in entry")
	if ((cur_time .feq. last_time) .or. (cur_time >= last_time+interval)) then
	  do_this_histo = .true.
	  if (cur_time >= last_time+interval) last_time = cur_time
	endif
      else
	do_this_histo = .true.
      endif
      if (do_this_histo) then
	! if (get_value(entry%at%params, "Time", cur_time)) call print("doing histo for config with time " // cur_time)
	n_histos = n_histos + 1
	if (mod(n_histos,10) == 0) write (mainlog%unit,'(I1,$)') mod(n_histos/10,10)
	if (mod(n_histos,1000) == 0) write (mainlog%unit,'(a)') " "
	call reallocate_histos(histo_count, n_histos, n_bins)
	if (radial_histo) then
	  call accumulate_radial_histo_count(histo_count(:,1,1,n_histos), entry%at, bin_width(1), n_bins(1), gaussian_sigma, mask_str)
	else
	  call accumulate_histo_count(histo_count(:,:,:,n_histos), entry%at, min_p, bin_width, n_bins, gaussian, gaussian_sigma, mask_str)
	endif
      endif
      entry => entry%next
    end do
  end subroutine calc_histos

  subroutine accumulate_radial_histo_count(histo_count, at, bin_width, n_bins, gaussian_sigma, mask_str)
    real(dp), intent(inout) :: histo_count(:)
    type(Atoms), intent(in) :: at
    real(dp), intent(in) :: bin_width
    integer, intent(in) :: n_bins
    real(dp), intent(in) :: gaussian_sigma
    character(len=*), optional, intent(in) :: mask_str

    real(dp) :: bin_r, p(3), dist, r
    logical, allocatable :: mask_a(:)
    integer at_i, t_i, p_i, bin_i
    integer :: n_t = 4, n_p = 2
    real(dp), allocatable :: theta(:), phi(:), w(:)

    allocate(mask_a(at%N))
    call is_in_mask(mask_a, at, mask_str)

    allocate(theta(n_t), phi(n_p), w(n_p+1))
    do t_i=1, n_t
      theta(t_i) = (t_i-1)*2.0_dp*PI/real(n_t,dp)
    end do
    do p_i=1, n_p
      phi(p_i) = (p_i-1-floor((n_p-1)/2.0_dp))*PI/real(n_p+1,dp)
!      w(p_i) = cos(phi(p_i))/(real(n_t*n_p,dp)*(gaussian_sigma*sqrt(PI))**3)
      w(p_i) = 1._dp/(real(n_t*n_p,dp)*(gaussian_sigma*sqrt(2.0_dp*PI))**3)
    end do

    do at_i=1, at%N
      if (.not. mask_a(at_i)) cycle
      r = norm(at%pos(:,at_i))
      do bin_i=1, n_bins
	bin_r = (real(bin_i,dp)-0.5_dp)*bin_width
!Include all the atoms, density fn won't curve down at the end of the plot
!	if (abs(r-bin_r) > 4.0_dp*gaussian_sigma) cycle
	do t_i=1, n_t
	do p_i=1, n_p
	  p(1) = bin_r*cos(theta(t_i))*cos(phi(p_i))
	  p(2) = bin_r*sin(theta(t_i))*cos(phi(p_i))
	  p(3) = bin_r*sin(phi(p_i))
	  dist = distance_min_image(at,p,at%pos(:,at_i))
!Include all the atoms, slow but minimises error
!	  if (dist > 4.0_dp*gaussian_sigma) cycle
!	  histo_count(bin_i) = histo_count(bin_i) + exp(-(dist/gaussian_sigma)**2)*w(p_i)
          histo_count(bin_i) = histo_count(bin_i) + exp(-0.5_dp*(dist/(gaussian_sigma))**2)*w(p_i)
	end do
	end do
      end do
    end do

    deallocate(mask_a, theta, phi, w)
  end subroutine accumulate_radial_histo_count

  subroutine accumulate_histo_count(histo_count, at, min_p, bin_width, n_bins, gaussian, gaussian_sigma, mask_str)
    real(dp), intent(inout) :: histo_count(:,:,:)
    type(Atoms), intent(in) :: at
    real(dp), intent(in) :: min_p(3), bin_width(3)
    integer, intent(in) :: n_bins(3)
    logical, intent(in) :: gaussian
    real(dp), intent(in) :: gaussian_sigma
    character(len=*), optional, intent(in) :: mask_str

    ! real(dp) :: min_p_lat(3), bin_width_lat(3), p_lat(3)
    integer :: i, bin(3), bin_1, bin_2, bin_3
    logical, allocatable :: mask_a(:)
    integer :: i1, i2, i3
    real(dp) :: r1, r2, r3, r1_sq, r2_sq, r3_sq
    integer :: n_samples = 20
    real(dp) :: p_center(3), p(3)
    real(dp) :: range
    real(dp) :: weight, normalization, n_samples_d, sigma_sq

    allocate(mask_a(at%N))
    call is_in_mask(mask_a, at, mask_str)

    range = 3.0_dp*gaussian_sigma
    normalization = ((2.0*range/real(2*n_samples,dp))**3)/(gaussian_sigma*sqrt(PI))**3

    n_samples_d = real(n_samples,dp)
    sigma_sq = gaussian_sigma**2

    ! min_p_lat = at%g .mult. min_p
    ! bin_width_lat = at%g .mult. bin_width
    do i=1, at%N
if (mod(i,10) == 0) call print("accumulate_histo_count i " // i,ERROR)
      if (.not. mask_a(i)) cycle
      if (gaussian) then
	p_center = at%pos(:,i)
	do i1=-n_samples, n_samples
	  r1 = real(i1,dp)/real(n_samples,dp)*range
	  r1_sq = r1*r1
	  p(1) = p_center(1) + r1
	  bin_1 = floor(p(1)-min_p(1)/bin_width(1))+1
	  if (bin_1 <= 0 .or. bin_1 > n_bins(1)) cycle
	  do i2=-n_samples, n_samples
	    r2 = real(i2,dp)/real(n_samples,dp)*range
	    r2_sq = r2*r2
	    p(2) = p_center(2) + r2
	    bin_2 = floor(p(2)-min_p(2)/bin_width(2))+1
	    if (bin_2 <= 0 .or. bin_2 > n_bins(2)) cycle
	    do i3=-n_samples, n_samples
	      r3 = real(i3,dp)/real(n_samples,dp)*range
	      r3_sq = r3*r3
	      p(3) = p_center(3) + r3
	      bin_3 = floor(p(3)-min_p(3)/bin_width(3))+1
	      if (bin_3 <= 0 .or. bin_3 > n_bins(3)) cycle
	      ! p_lat = at%g .mult. p
	      ! bin = floor((p_lat-min_p_lat)/bin_width_lat)+1
	      !! bin = floor((p-min_p)/bin_width)+1
	      weight = normalization*exp(-(r1_sq+r2_sq+r3_sq)/sigma_sq)
	      !! if (all(bin >= 1) .and. all (bin <= n_bins)) 
	      histo_count(bin_1,bin_2,bin_3) = histo_count(bin_1,bin_2,bin_3) + weight
	    end do
	  end do
	end do
      else
	! p_lat = at%g .mult. at%pos(:,i)
	! bin = floor((p_lat-min_p_lat)/bin_width_lat)+1
	bin = floor((p-min_p)/bin_width)+1
	if (all(bin >= 1) .and. all (bin <= n_bins)) histo_count(bin(1),bin(2),bin(3)) = histo_count(bin(1),bin(2),bin(3)) + 1.0_dp
      endif
    end do
  end subroutine accumulate_histo_count

  subroutine is_in_mask(mask_a, at, mask_str)
    type(Atoms), intent(in) :: at
    logical, intent(out) :: mask_a(at%N)
    character(len=*), optional, intent(in) :: mask_str

    integer :: i
    integer :: Zmask
    type(Table) :: atom_indices

    if (.not. present(mask_str)) then
      mask_a = .true.
      return
    endif

    if (len_trim(mask_str) == 0) then
      mask_a = .true.
      return
    endif

    mask_a = .false.
    if (mask_str(1:1)=='@') then
      call parse_atom_mask(mask_str,atom_indices)
      do i=1, atom_indices%N
	mask_a(atom_indices%int(1,i)) = .true.
      end do
    else if (scan(mask_str,'=')/=0) then
      call system_abort("property type mask not supported yet")
    else
      Zmask = Atomic_Number(mask_str)
      do i=1, at%N
	if (at%Z(i) == Zmask) mask_a(i) = .true.
      end do
    end if
  end subroutine is_in_mask

  subroutine reallocate_histos(histos, n, n_bins)
    real(dp), allocatable, intent(inout) :: histos(:,:,:,:)
    integer, intent(in) :: n, n_bins(3)

    integer :: new_size
    real(dp), allocatable :: t_histos(:,:,:,:)

    if (allocated(histos)) then
      if (n <= size(histos,4)) return
      allocate(t_histos(size(histos,1),size(histos,2),size(histos,3),size(histos,4)))
      t_histos = histos
      deallocate(histos)
      if (size(t_histos,4) <= 0) then
	new_size = 10
      else if (size(t_histos,4) < 1000) then
	new_size = 2*size(t_histos,4)
      else
	new_size = floor(1.25*size(t_histos,4))
      endif
      allocate(histos(size(t_histos,1), size(t_histos,2), size(t_histos,3), new_size))
      histos(1:size(t_histos,1),1:size(t_histos,2),1:size(t_histos,3),1:size(t_histos,4)) = &
	t_histos(1:size(t_histos,1),1:size(t_histos,2),1:size(t_histos,3),1:size(t_histos,4))
      histos(1:size(t_histos,1),1:size(t_histos,2),1:size(t_histos,3),size(t_histos,4)+1:size(histos,4)) = 0.0_dp
      deallocate(t_histos)
    else
      allocate(histos(n_bins(1), n_bins(2), n_bins(3), n))
    endif
  end subroutine reallocate_histos

  subroutine calc_mean_var_3array(a, mean, var)
    real(dp), intent(in) :: a(:,:,:,:)
    real(dp), intent(out) :: mean(:,:,:)
    real(dp), optional, intent(out) :: var(:,:,:)

    integer :: i1, i2, i3, n

    n = size(a,4)
    mean = sum(a, 4)/real(n,dp)
    if (present(var)) then
      do i1=1, size(a,1)
      do i2=1, size(a,2)
      do i3=1, size(a,3)
	var(i1,i2,i3) = sum((a(i1,i2,i3,:)-mean(i1,i2,i3))**2)/real(n,dp)
      end do
      end do
      end do
    end if
  end subroutine calc_mean_var_3array

end program density_new
