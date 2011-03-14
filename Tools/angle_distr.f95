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

! Calculates angle distribution in a radial mesh
! with error bars
! outputs: angle (1)  O-O distance (A)  distance from origin (A)  angle distribution (arb.units)

program angle_distr
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
  character(len=FIELD_LENGTH) :: mask_center
  integer :: decimation
  real(dp) :: min_time, max_time
  logical :: gaussian_smoothing
  real(dp) :: gaussian_sigma
  real(dp) :: sampling_sigma
  real(dp) :: gaussian_angle_sigma
  logical :: radial_histo, random_samples
  logical :: do_distance_dependence
  integer :: n_samples(2)
  integer :: nth_snapshots
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
  real(dp) :: m_bin_width

  real(dp) :: part_dens, cell_vol

  call system_initialise(PRINT_VERBOSE)

  call initialise(cli_params)
  call register_cli_params(cli_params,.true., infilename, infile_is_list, outfilename, commandfilename, &
    mask_str, mask_center, min_p, bin_width, n_bins, decimation, min_time, max_time, gaussian_smoothing, gaussian_sigma, gaussian_angle_sigma, sampling_sigma, &
    radial_histo, random_samples, n_samples, do_distance_dependence, nth_snapshots, &
    sort_Time, no_Time_dups, mean, mean_decorrelation_time, autocorrelation, autocorrelation_max_lag, quiet)
  if (.not. param_read_args(cli_params)) then
      call print_usage()
    call system_abort('could not parse argument line')
  end if

  call print("Reading configurations")
  no_compute_index=.false.
  if ((.not. sort_Time) .and. (.not. no_Time_dups)) no_compute_index=.true.
  call read_xyz(structure_ll, infilename, infile_is_list, decimation, min_time, max_time, sort_Time, no_Time_dups, quiet, no_compute_index)

  if (len_trim(commandfilename) == 0) then
    args_str = param_write_string(cli_params)
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
      mask_str, mask_center, min_p, bin_width, n_bins, decimation, min_time, max_time, gaussian_smoothing, gaussian_sigma, gaussian_angle_sigma, sampling_sigma, &
      radial_histo, random_samples, n_samples, do_distance_dependence, nth_snapshots, &
      sort_Time, no_Time_dups, mean, mean_decorrelation_time, autocorrelation, autocorrelation_max_lag, quiet)
    if (.not. param_read_line(cli_params, trim(args_str), ignore_unknown = .true.)) then
      call print_usage()
      call system_abort("could not parse argument line "//arg_line_no//" from file '"//trim(commandfilename//"'"))
    end if
    call finalise(cli_params)
  endif

  have_params = .true.
  do while (have_params)
    call print("infile " // trim(infilename) // " infile_is_list " // infile_is_list)
    call print("outfilename " // trim(outfilename))
    call print("decimation " // decimation // " min_time " // min_time // " max_time " // max_time)
    call print("AtomMask " // trim(mask_str))
    call print("CentMask " // trim(mask_center))
    call print("radial_histo " // radial_histo)
    if (radial_histo) then
      do i1 = 1, 3
       if (n_bins(i1) < 1) n_bins(i1) = 1
      end do
      call print("bin_width(1,2,3) "//bin_width(1)//" "//bin_width(2)//" "//bin_width(3)&
                 //" n_bins(1,2,3) " // n_bins(1)//" "//n_bins(2)//" "//n_bins(3))
      call print("random_samples " // random_samples)
      if (random_samples) then
	call print("n_samples " // n_samples(1))
      else
	call print("n_theta " // n_samples(1) // " n_phi " // n_samples(2))
      endif
    else
      call print("min_p " // min_p)
      call print("bin_width " // bin_width // " n_bins " // n_bins)
    endif
    call print("gaussian " // gaussian_smoothing // " sigma " // gaussian_sigma // " angle_sigma " // gaussian_angle_sigma // " sampling_sigma " // sampling_sigma)
    call print("do_distance_dependence " // do_distance_dependence)
    call print("nth_snapshots " // nth_snapshots)
    call print("sort_Time " // sort_Time // " no_Time_dups " // no_Time_dups)
    call print("mean " // mean // " mean_decorrelation_time " // mean_decorrelation_time)
    call print("autocorrelation " // autocorrelation // " autocorrelation_max_lag " // autocorrelation_max_lag)

    call print("Calculating angle distributions")
    call calc_histos(histo_raw, n_histos, min_p, bin_width, n_bins, structure_ll, mean_decorrelation_time, &
                     gaussian_smoothing, gaussian_sigma, gaussian_angle_sigma, sampling_sigma, &
                     radial_histo, random_samples, n_samples, mask_str, mask_center, do_distance_dependence, nth_snapshots)

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
       call print("# O--H-O angles (#/rad) ", file=outfile)
       call print("# r mean var n_samples", file=outfile)
       if(bin_width(1) == 0.0_dp) then
         if(n_bins(1) .ne. 1) then
           m_bin_width = PI/real(n_bins(1),dp)
         else
           m_bin_width = 0.0_dp
         end if
       else
         m_bin_width = bin_width(1)
       end if
       do i3=1, n_bins(3)
	do i1=1, n_bins(1)
          do i2=1, n_bins(2)
	    !call print(PI-(bin_width(1)*(real(i1,dp)-0.5)) // " " // (1.8_dp+(bin_width(2)*(real(i2,dp)-0.5))) // " " &
            !           //(bin_width(3)*(real(i3,dp)-0.5))// " " &
            !           // histo_raw(i1, i2, i3, 1) // " " //n_histos, file=outfile)
	    call print(PI-(m_bin_width*(real(i1,dp)-0.5)) // " " // (1.8_dp+(bin_width(2)*(real(i2,dp)-0.5))) // " " &
                       //(bin_width(3)*(real(i3,dp)-0.5))// " " &
                       // histo_mean(i1, i2, i3) // " " //  histo_var(i1, i2, i3)// " " //n_histos, file=outfile)
          end do
          call print('',file=outfile)
	end do
       end do
      else
	call print("# Density (#/vol)", file=outfile)
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
	  mask_str, mask_center, min_p, bin_width, n_bins, decimation, min_time, max_time, gaussian_smoothing, gaussian_sigma, gaussian_angle_sigma, sampling_sigma, &
	  radial_histo, random_samples, n_samples, do_distance_dependence, nth_snapshots, &
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

    call print("Usage: " // trim(my_exec_name)//" infile=filename [infile_is_list=logical(F)]", PRINT_ALWAYS)
    call print("       outfile=filename [commandfile=filename()] [AtomMask=species()] [min_p={x y z}]", PRINT_ALWAYS)
    call print("       bin_width={x y z}(y,z ignored for radial_histo) n_bins={nx ny nz}(y,z ignored for radial_histo)", PRINT_ALWAYS)
    call print("       [mean=logical(T)] [mean_decorrelation_time=t(0.0)]", PRINT_ALWAYS)
    call print("       [autocorrelation=logical(F)] [autocorrelation_max_lag=N(10000)]", PRINT_ALWAYS)
    call print("       [decimation=n(1)] [min_time=t(-1.0)] [max_time=t(-1.0)]", PRINT_ALWAYS)
    call print("       [gaussian=logical(F)(always true for radial_histo)] [sigma=s(1.0)] [radial_histo=logical(F)]", PRINT_ALWAYS)
    call print("       [random_samples=logical(F)(radial_histo only)] [n_samples={2 4})(component 2 ignored for random)]", PRINT_ALWAYS)
    call print("       [sort_Time(F)] [no_Time_dups(F)]", PRINT_ALWAYS)
    call print("       [quiet=logical(F)]", PRINT_ALWAYS)
  end subroutine print_usage

  subroutine register_cli_params(cli_params, initial, infilename, infile_is_list, outfilename, commandfilename, &
    mask_str, mask_center, min_p, bin_width, n_bins, decimation, min_time, max_time, gaussian_smoothing, gaussian_sigma, gaussian_angle_sigma, sampling_sigma, &
    radial_histo, random_samples, n_samples, do_distance_dependence, nth_snapshots, &
    sort_Time, no_Time_dups, mean, mean_decorrelation_time, autocorrelation, autocorrelation_max_lag, quiet)
    type(Dictionary), intent(inout) :: cli_params
    logical, intent(in) :: initial
    character(len=*), intent(inout) :: infilename
    logical, intent(inout) :: infile_is_list
    character(len=*), intent(inout) :: outfilename
    character(len=*), intent(inout) :: commandfilename
    character(len=*), intent(inout) :: mask_str
    character(len=*), intent(inout) :: mask_center
    real(dp), intent(inout) :: min_p(3), bin_width(3)
    integer, intent(inout) :: n_bins(3), decimation
    real(dp), intent(inout) :: min_time, max_time
    logical, intent(inout) :: gaussian_smoothing
    real(dp), intent(inout) :: gaussian_sigma
    real(dp), intent(inout) :: gaussian_angle_sigma
    real(dp), intent(inout) :: sampling_sigma
    logical, intent(inout) :: radial_histo
    logical, intent(inout) :: random_samples
    logical, intent(inout) :: do_distance_dependence
    integer, intent(inout) :: nth_snapshots
    integer, intent(inout) :: n_samples(2)
    logical, intent(inout) :: sort_Time, no_Time_dups
    logical, intent(inout) :: mean
    real(dp), intent(inout) :: mean_decorrelation_time
    logical, intent(inout) :: autocorrelation
    integer, intent(inout) :: autocorrelation_max_lag
    logical, intent(inout) :: quiet
    character(len=FIELD_LENGTH) :: t_field


    if (initial) then
      call param_register(cli_params, 'infile', param_mandatory, infilename, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(cli_params, 'infile_is_list', 'F', infile_is_list, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(cli_params, 'outfile', 'stdout', outfilename, help_string="No help yet.  This source file was $LastChangedBy$")
      mask_str = ""
      mask_center = ""
      commandfilename = ""
      call param_register(cli_params, 'commandfile', "", commandfilename, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(cli_params, 'AtomMask', "O", mask_str, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(cli_params, 'CentMask', "O", mask_center, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(cli_params, 'min_p', "0.0 0.0 0.0", min_p, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(cli_params, 'bin_width', PARAM_MANDATORY, bin_width, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(cli_params, 'n_bins', PARAM_MANDATORY, n_bins, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(cli_params, 'decimation', '1', decimation, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(cli_params, 'min_time', '-1.0', min_time, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(cli_params, 'max_time', '-1.0', max_time, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(cli_params, 'gaussian', 'F', gaussian_smoothing, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(cli_params, 'sigma', '0.0', gaussian_sigma, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(cli_params, 'do_distance_dependence', 'F', do_distance_dependence, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(cli_params, 'nth_snapshots', '1', nth_snapshots, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(cli_params, 'angle_sigma', '0.0', gaussian_angle_sigma, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(cli_params, 'sampling_sigma', '0.0', sampling_sigma, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(cli_params, 'radial_histo', 'F', radial_histo, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(cli_params, 'random_samples', 'F', random_samples, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(cli_params, 'n_samples', '2 4', n_samples, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(cli_params, 'sort_Time', 'F', sort_Time, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(cli_params, 'no_Time_dups', 'F', no_Time_dups, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(cli_params, 'mean', 'T', mean, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(cli_params, 'mean_decorrelation_time', '0.0', mean_decorrelation_time, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(cli_params, 'autocorrelation', 'F', autocorrelation, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(cli_params, 'autocorrelation_max_lag', '10000', autocorrelation_max_lag, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(cli_params, 'quiet', 'F', quiet, help_string="No help yet.  This source file was $LastChangedBy$")
    else
      t_field=infilename
      call param_register(cli_params, 'infile', trim(t_field), infilename, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(cli_params, 'infile_is_list', ""//infile_is_list, infile_is_list, help_string="No help yet.  This source file was $LastChangedBy$")
      t_field=outfilename
      call param_register(cli_params, 'outfile', trim(t_field), outfilename, help_string="No help yet.  This source file was $LastChangedBy$")
      t_field=commandfilename
      call param_register(cli_params, 'commandfile', trim(t_field), commandfilename, help_string="No help yet.  This source file was $LastChangedBy$")
      t_field=mask_str
      t_field=mask_center
      call param_register(cli_params, 'AtomMask', trim(t_field), mask_str, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(cli_params, 'CentMask', trim(t_field), mask_center, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(cli_params, 'min_p', ""//min_p, min_p, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(cli_params, 'bin_width', ""//bin_width, bin_width, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(cli_params, 'n_bins', ""//n_bins, n_bins, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(cli_params, 'min_time', ""//min_time, min_time, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(cli_params, 'max_time', ""//max_time, max_time, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(cli_params, 'gaussian', ""//gaussian_smoothing, gaussian_smoothing, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(cli_params, 'sigma', ""//gaussian_sigma, gaussian_sigma, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(cli_params, 'do_distance_dependence', ""//do_distance_dependence, do_distance_dependence, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(cli_params, 'nth_snapshots', ""//nth_snapshots, nth_snapshots, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(cli_params, 'angle_sigma', ""//gaussian_angle_sigma, gaussian_angle_sigma, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(cli_params, 'sampling_sigma', ""//sampling_sigma, sampling_sigma, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(cli_params, 'radial_histo', ""//radial_histo, radial_histo, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(cli_params, 'random_samples', ""//random_samples, random_samples, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(cli_params, 'n_samples', '2 4', n_samples, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(cli_params, 'mean', ""//mean, mean, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(cli_params, 'mean_decorrelation_time', ""//mean_decorrelation_time, mean_decorrelation_time, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(cli_params, 'autocorrelation', ""//autocorrelation, autocorrelation, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(cli_params, 'autocorrelation_max_lag', ""//autocorrelation_max_lag, autocorrelation_max_lag, help_string="No help yet.  This source file was $LastChangedBy$")
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

  subroutine calc_histos(histo_count, n_histos, min_p, bin_width, n_bins, structure_ll, interval, gaussian, &
                         gaussian_sigma, gaussian_angle_sigma, sampling_sigma, radial_histo, &
                         random_samples, n_samples, mask_str, mask_center, do_distance_dependence, nth_snapshots)
    real(dp), intent(inout), allocatable :: histo_count(:,:,:,:)
    integer, intent(out) :: n_histos
    real(dp), intent(in) :: min_p(3), bin_width(3)
    integer, intent(in) :: n_bins(3)
    type(atoms_ll), intent(in) :: structure_ll
    real(dp), intent(in) :: interval
    logical, intent(in) :: gaussian
    real(dp), intent(in) :: gaussian_sigma
    real(dp), intent(in) :: gaussian_angle_sigma
    real(dp), intent(in) :: sampling_sigma
    logical, intent(in) :: radial_histo
    logical :: random_samples
    integer :: n_samples(2)
    character(len=*), optional, intent(in) :: mask_str
    character(len=*), optional, intent(in) :: mask_center
    logical, intent(in) :: do_distance_dependence    
    integer, intent(in) :: nth_snapshots

    real(dp) :: last_time, cur_time
    type(atoms_ll_entry), pointer :: entry
    logical :: do_this_histo
    real(dp), allocatable :: theta(:), phi(:), w(:)

    n_histos = 0
    allocate(histo_count(n_bins(1),n_bins(2),n_bins(3),10))
    histo_count = 0.0_dp

    if (random_samples) then
      call calc_samples_random(n_samples(1), theta, phi, w)
    else
      call calc_samples_grid(n_samples(1), n_samples(2), theta, phi, w)
    endif
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
	if (mod(n_histos,10) == 0) then
          if(mod(n_histos,50) == 0) then
            write (mainlog%unit,'(i0)') n_histos
          else
            write (mainlog%unit,'(a,$)') '|'
          end if
        else
          write(mainlog%unit,'(a,$)') '.'
        end if
	if (mod(n_histos,1000) == 0) write (mainlog%unit,'(a)') " "
	call reallocate_histos(histo_count, n_histos, n_bins)
	if (radial_histo) then
          call calc_connect(entry%at)
          cell_vol = cell_volume(entry%at%lattice)
	  call accumulate_radial_histo_count(histo_count(:,:,:,n_histos), entry%at, bin_width(1), bin_width(2), bin_width(3), &
                                             n_bins(1), n_bins(2), n_bins(3), &
                                             gaussian, gaussian_sigma, gaussian_angle_sigma, sampling_sigma, &
                                             theta, phi, w, mask_str, mask_center, do_distance_dependence)
	endif
      endif
      entry => entry%next
    end do
  end subroutine calc_histos

  subroutine calc_samples_grid(n_t, n_p, theta, phi, w)
    integer, intent(in) :: n_t, n_p
    real(dp), allocatable, intent(inout) :: theta(:), phi(:), w(:)

    integer :: t_i, p_i, ii

    if (allocated(theta)) deallocate(theta)
    if (allocated(phi)) deallocate(phi)
    if (allocated(w)) deallocate(w)
    allocate(theta(n_t*n_p))
    allocate(phi(n_t*n_p))
    allocate(w(n_t*n_p))

    ii = 1
    do t_i=1, n_t
    do p_i=1, n_p
      phi(ii) = (p_i-1)*2.0_dp*PI/real(n_p,dp)
      if (mod(n_t,2) == 0) then
	theta(ii) = (t_i-1.5_dp-floor((n_t-1)/2.0_dp))*PI/real(n_t,dp)
      else
	theta(ii) = (t_i-1-floor((n_t-1)/2.0_dp))*PI/real(n_t,dp)
      endif
      ! not oversample top and bottom of spehre
      w(ii) = cos(theta(ii))
      ii = ii + 1
    end do
    end do
    w = w / ( (sampling_sigma*sqrt(2.0_dp*PI))**3 * sum(w) )
  end subroutine calc_samples_grid

  subroutine calc_samples_random(n, theta, phi, w)
    integer, intent(in) :: n
    real(dp), allocatable, intent(inout) :: theta(:), phi(:), w(:)

    integer :: ii
    real(dp) :: p(3), p_norm

    if (allocated(theta)) deallocate(theta)
    if (allocated(phi)) deallocate(phi)
    if (allocated(w)) deallocate(w)
    allocate(theta(n))
    allocate(phi(n))
    allocate(w(n))

    do ii = 1, n
      p_norm = 2.0_dp
      do while (p_norm > 1.0_dp)
	p = 0.0_dp
	call randomise(p, 2.0_dp)
	p_norm = norm(p)
      end do
      phi(ii) = atan2(p(2),p(1))
      theta(ii) = acos(p(3)/p_norm)
      w(ii) = 1.0_dp
    end do
    w = w / ( (sampling_sigma*sqrt(2.0_dp*PI))**3 * sum(w) )

  end subroutine calc_samples_random

  subroutine accumulate_radial_histo_count(histo_count, at, bin_width, dbin_width, sbin_width, n_bins, dn_bins, sn_bins, &
                                           gaussian, gaussian_sigma, gaussian_angle_sigma, sampling_sigma, &
                                           theta, phi, w, mask_str, mask_center, do_distance_dependence)
    real(dp), intent(inout) :: histo_count(:,:,:)
    type(Atoms), intent(in) :: at
    real(dp), intent(in)    :: bin_width, dbin_width, sbin_width
    integer, intent(in)     :: n_bins, dn_bins, sn_bins
    logical,intent(in)      :: gaussian
    real(dp), intent(in)    :: gaussian_sigma, gaussian_angle_sigma, sampling_sigma, theta(:), phi(:), w(:)
    character(len=*), optional, intent(in) :: mask_str
    character(len=*), optional, intent(in) :: mask_center
    logical,intent(in)      :: do_distance_dependence

    real(dp) :: bin_r, bin_a, bin_s, p(3), dist, r, d, sdist, distw, ang, v1(3), v2(3), m_bin_width
    logical, allocatable :: mask_a(:), mask_c(:)
    integer at_i, at_j, sample_i, bin_i, bin_j, bin_k, i, k, m
    real(dp) :: myNorm, fotPI, oostPI, sina, shell, sample_shell
    integer  :: myCount, sample_size
    real(dp) :: mycent(3)

    if(sbin_width == 0.0_dp) then
      sample_size = 1
    else
      sample_size = size(theta)
    end if
    
    allocate(mask_a(at%N))
    call is_in_mask(mask_a, at, mask_str)
    allocate(mask_c(at%N))
    call is_in_mask(mask_c, at, mask_center)

    m = 0
    do i = 1, at%N
     if (mask_a(i)) m=m+1
    end do
    part_dens=real(m,dp)/cell_vol

    m=0
    do i = 1, at%N
     if(mask_c(i)) then
       m=m+1
       k=i
     end if
    end do

    if(m==1) then
      mycent=at%pos(:,k)
    else
      mycent=(/0.0_dp,0.0_dp,0.0_dp/)
    end if

    myNorm = 1/(gaussian_angle_sigma*sqrt(2*PI))
    if(bin_width == 0.0_dp) then
      if(n_bins .ne. 1) then
        m_bin_width = PI/real(n_bins,dp)
      else
        m_bin_width = 0.0_dp
      end if
    else
      m_bin_width = bin_width
    end if
    fotPI = (4.0_dp * PI)/3.0_dp   
    oostPI = 1.0_dp/(sqrt(2.0_dp*PI))
! ratio of 20/4=5 is bad
! ratio of 20/3=6.66 is bad
! ratio of 20/2.5=8 is borderline (2e-4)
! ratio of 20/2.22=9 is fine  (error 2e-5)
! ratio of 20/2=10 is fine
    if ( (norm(at%lattice(:,1)) < 9.0_dp*sampling_sigma) .or. &
         (norm(at%lattice(:,2)) < 9.0_dp*sampling_sigma) .or. &
         (norm(at%lattice(:,3)) < 9.0_dp*sampling_sigma) ) &
      call print("WARNING: at%lattice may be too small for sigma, errors (noticeably too low a density) may result", PRINT_ALWAYS)
  if(.not. do_distance_dependence) then
   if(gaussian) then
     myCount=0
     do at_i=1, at%N
       if(.not. mask_c(at_i)) cycle
       myCount = myCount+1
       do bin_i = 1, n_bins !angles
       bin_a = PI - (real(bin_i,dp)-0.5)*m_bin_width
         do bin_j = 1, dn_bins !distances
           bin_r = 1.8_dp +(real(bin_j,dp)-0.5)*dbin_width
           shell = fotPI*( (bin_r+1.5_dp*gaussian_sigma)**3 - (max(0.0_dp,(bin_r-1.5_dp*gaussian_sigma)**3)) )*part_dens
           do at_j = 1, at%N
             if((at_j == at_i) .or. .not. mask_a(at_j)) cycle !(at%Z(at_j) .ne. 8)) cycle
             dist = distance_min_image(at,at%pos(:,at_i),at%pos(:,at_j))
              do i = 1, atoms_n_neighbours(at,at_j)
                k = atoms_neighbour(at,at_j,i)
                r = 0.0_dp
                if(is_nearest_neighbour(at,at_j,i) .and. at%Z(k) == 1) then
                  d = distance_min_image(at,at%pos(:,at_i),at%pos(:,k))
                  if(r == 0.0_dp .or. d < r) then
                    m = k !select the mth atom (hydrogen)
                    r = d
                  end if
                end if !nearest neighbour
              end do !neighbours of at_j
              v1 = at%pos(:,at_i) - at%pos(:,k)
              v2 = at%pos(:,at_j) - at%pos(:,k)
              ang = angle(v1,v2)
              sina = sin(PI-ang)
              if(sina .lt. 1.0E-10) sina=1.0E-10
              ang = myNorm*exp(-0.5_dp*((ang-bin_a)/gaussian_angle_sigma)**2)
              ang = ang/(2.0*PI*sina)
              histo_count(bin_i,bin_j,1) = histo_count(bin_i,bin_j,1) + &
                                           (ang*exp(-0.5_dp*((dist-bin_r)/gaussian_sigma)**2)*( (oostPI/gaussian_sigma) / shell ))
           end do !at_j
         end do !distances
       end do !angles
     end do ! at_i 
     do bin_i=1, n_bins
      do bin_j=1, dn_bins
       histo_count(bin_i,bin_j,1) = histo_count(bin_i,bin_j,1) / real(myCount,dp)
      end do
     end do 
   else if(.not. gaussian) then !use radial bins, weighted down by volume (assuming homogenous density)
    myCount = 0
    do at_i=1, at%N
      myCount = myCount+1
      if(.not. mask_c(at_i)) cycle
      do bin_i = 1, n_bins !angles
        bin_a = PI - (real(bin_i,dp)-0.5_dp)*m_bin_width
        do bin_j = 1, dn_bins !distances
          bin_r = 1.8_dp + ((real(bin_j,dp)-0.5_dp)*dbin_width)
          shell = fotPI*((bin_r+0.5_dp*dbin_width)**3 - (bin_r-0.5_dp*dbin_width)**3)*part_dens
          do at_j = 1, at%N
            if((at_j == at_i) .or. .not.mask_a(at_j)) cycle !(at%Z(at_j) .ne. 8)) cycle
            dist = distance_min_image(at,at%pos(:,at_i),at%pos(:,at_j))
            if( (dist > (bin_r-(dbin_width*0.5_dp))) .and. (dist < (bin_r+(dbin_width*0.5_dp))) ) then
              do i = 1, atoms_n_neighbours(at,at_j)
                k = atoms_neighbour(at,at_j,i)
                r = 0.0_dp
                if(is_nearest_neighbour(at,at_j,i) .and. at%Z(k) == 1) then
                  d = distance_min_image(at,at%pos(:,at_i),at%pos(:,k))
                  if(r == 0.0_dp .or. d < r) then
                    m = k !select the mth atom (hydrogen)
                    r = d
                  end if
                end if !nearest neighbour
              end do !neighbours of at_j
              v1 = at%pos(:,at_i) - at%pos(:,k)
              v2 = at%pos(:,at_j) - at%pos(:,k)
              ang = angle(v1,v2)
              sina=sin(PI-ang)
              if( (ang > (bin_a-(m_bin_width*0.5_dp))) .and. (ang < (bin_a+(m_bin_width*0.5_dp))) ) then
                histo_count(bin_i,bin_j,1) = histo_count(bin_i,bin_j,1) + (1.0_dp/(2.0*PI*sina))/shell
              end if !angle bin
            end if !distance bin
          end do !at_j
        end do !distances
      end do !angles
    end do ! at_i
    do bin_i=1, n_bins
     do bin_j = 1, dn_bins
      histo_count(bin_i,bin_j,1) = histo_count(bin_i,bin_j,1) / real(myCount,dp)
     end do
    end do 
   end if !no gaussian smoothing  
  else !distance dependence
   if(gaussian) then
    do bin_k = 1, sn_bins
      bin_s = (real(bin_k,dp)-0.5_dp)*sbin_width
     do at_i=1, at%N
       if(.not. mask_c(at_i)) cycle
         do sample_i = 1, sample_size 
           p(1) = mycent(1)+bin_s*cos(theta(sample_i))*cos(phi(sample_i))
           p(2) = mycent(2)+bin_s*cos(theta(sample_i))*sin(phi(sample_i))
           p(3) = mycent(3)+bin_s*sin(theta(sample_i))
           sdist = distance_min_image(at,p,at%pos(:,at_i))
           if( sdist < 3.0_dp*sampling_sigma) then
            distw = w(sample_i)*exp(-0.5_dp*(sdist/sampling_sigma)**2)
            do bin_i = 1, n_bins !angles
            bin_a = PI - (real(bin_i,dp)-0.5)*m_bin_width
              do bin_j = 1, dn_bins !distances
                bin_r = 1.8_dp +(real(bin_j,dp)-0.5_dp)*dbin_width
                shell = fotPI*( (bin_r+1.5_dp*gaussian_sigma)**3 - (max(0.0_dp,(bin_r-1.5_dp*gaussian_sigma)**3)) )*part_dens
                do at_j = 1, at%N
                  if((at_j == at_i) .or. .not.mask_a(at_j)) cycle !.or. (at%Z(at_j) .ne. 8)) cycle
                  dist = distance_min_image(at,at%pos(:,at_i),at%pos(:,at_j))
                   do i = 1, atoms_n_neighbours(at,at_j)
                     k = atoms_neighbour(at,at_j,i)
                     r = 0.0_dp
                     if(is_nearest_neighbour(at,at_j,i) .and. at%Z(k) == 1) then
                       d = distance_min_image(at,at%pos(:,at_i),at%pos(:,k))
                       if(r == 0.0_dp .or. d < r) then
                         m = k !select the mth atom (hydrogen)
                         r = d
                       end if
                     end if !nearest neighbour
                   end do !neighbours of at_j
                   v1 = at%pos(:,at_i) - at%pos(:,k)
                   v2 = at%pos(:,at_j) - at%pos(:,k)
                   ang = angle(v1,v2)
                   sina=sin(PI-ang)
                   if(sina .lt. 1.0E-10) sina=1.0E-10
                   ang = myNorm*exp(-0.5_dp*((ang-bin_a)/gaussian_angle_sigma)**2)
                   ang = ang/(2.0_dp*PI*sina)
                   histo_count(bin_i,bin_j,bin_k) = histo_count(bin_i,bin_j,bin_k) + &
                                                (ang*exp(-0.5_dp*((dist-bin_r)/gaussian_sigma)**2)*( (oostPI/gaussian_sigma) / shell )*distw)
                end do !at_j
              end do !distances
            end do !angles
          end if
         end do !samples
       end do ! at_i
     end do !sampling distances
   else !distance dependence, but not gaussian
    myCount = 0
    do bin_k = 1, sn_bins
      bin_s = (real(bin_k,dp)-0.5_dp)*sbin_width
      sample_shell = fotPI*((bin_s+0.5_dp*sbin_width)**3 - (bin_s-0.5_dp*sbin_width)**3)*part_dens
      do at_i=1, at%N
      !do at_i=991, 991
        if(.not. mask_c(at_i)) cycle
        !if(at%Z(at_i) .ne. 8) cycle !for oxygens only
        sdist = distance_min_image(at,(/0.0_dp,0.0_dp,0.0_dp/),at%pos(:,at_i))
        if( sdist > (bin_s-0.5_dp*sbin_width) .and. sdist < (bin_s+0.5_dp*sbin_width)) then
         do bin_i = 1, n_bins !angles
           bin_a = PI - (real(bin_i,dp)-0.5_dp)*m_bin_width
           do bin_j = 1, dn_bins !distances
             bin_r = 1.8_dp + ((real(bin_j,dp)-0.5_dp)*dbin_width)
             shell = fotPI*((bin_r+0.5_dp*dbin_width)**3 - (bin_r-0.5_dp*dbin_width)**3)*part_dens
             do at_j = 1, at%N
               if((at_j == at_i) .or. .not.mask_a(at_j)) cycle !(at%Z(at_j) .ne. 8)) cycle
               dist = distance_min_image(at,at%pos(:,at_i),at%pos(:,at_j))
               if( (dist > (bin_r-(dbin_width*0.5_dp))) .and. (dist < (bin_r+(dbin_width*0.5_dp))) ) then
                 do i = 1, atoms_n_neighbours(at,at_j)
                   k = atoms_neighbour(at,at_j,i)
                   r = 0.0_dp
                   if(is_nearest_neighbour(at,at_j,i) .and. at%Z(k) == 1) then
                     d = distance_min_image(at,at%pos(:,at_i),at%pos(:,k))
                     if(r == 0.0_dp .or. d < r) then
                       m = k !select the mth atom (hydrogen)
                       r = d
                     end if
                   end if !nearest neighbour
                 end do !neighbours of at_j
                 v1 = at%pos(:,at_i) - at%pos(:,k)
                 v2 = at%pos(:,at_j) - at%pos(:,k)
                 ang = angle(v1,v2)
                 sina=sin(PI-ang)
                 !cosa=cos(PI-ang)
                 if( (ang > (bin_a-(m_bin_width*0.5_dp))) .and. (ang < (bin_a+(m_bin_width*0.5_dp))) ) then
                   histo_count(bin_i,bin_j,bin_k) = histo_count(bin_i,bin_j,bin_k) + (1.0_dp/(2.0*PI*sina))/(shell*sample_shell)
                 end if !angle bin
               end if !distance bin
             end do !at_j
           end do !distances
         end do !angles 
       end if !if at_i is in sample_shell
      end do ! at_i
    end do !bin_k
   end if !gaussian in distance  dependence
  end if

    deallocate(mask_a)
  end subroutine accumulate_radial_histo_count

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

end program angle_distr
