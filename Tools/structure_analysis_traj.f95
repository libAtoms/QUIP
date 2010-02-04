! do various structural analysis to a trajcetory, save an intermediate file
! to be postprocessed (mean, variance, correlation, etc)
! #/vol density on a radial mesh
! #/vol density on a grid
! to be added: RDF, ADF

module structure_analysis_module
use libatoms_module
implicit none
private

type analysis
  logical :: density_radial, density_grid, rdfd, integrated_rdfd
  character(len=FIELD_LENGTH) :: outfilename
  character(len=FIELD_LENGTH) :: mask_str

  integer :: n_configs = 0

  real(dp) :: min_time, max_time
  integer :: min_frame, max_frame

  ! radial density stuff
  real(dp) :: radial_min_p, radial_bin_width
  integer :: radial_n_bins
  real(dp) :: radial_center(3)
  real(dp) :: radial_gaussian_sigma
  logical :: radial_random_angular_samples
  integer :: radial_n_angular_samples(2)
  real(dp), allocatable :: radial_histograms(:,:)
  real(dp), allocatable :: radial_pos(:)
  real(dp), allocatable :: radial_dr(:,:), radial_w(:)

  ! density grid stuff
  real(dp) :: grid_min_p(3), grid_bin_width(3)
  integer :: grid_n_bins(3)
  logical :: grid_gaussian_smoothing
  real(dp) :: grid_gaussian_sigma
  real(dp), allocatable :: grid_histograms(:,:,:,:)
  real(dp), allocatable :: grid_pos(:,:,:,:)

  ! rdfd stuff
  real(dp) :: rdfd_zone_center(3)
  character(FIELD_LENGTH) :: rdfd_center_mask_str, rdfd_neighbour_mask_str
  real(dp) :: rdfd_zone_width, rdfd_bin_width
  integer :: rdfd_n_zones, rdfd_n_bins
  logical :: rdfd_gaussian_smoothing
  real(dp) :: rdfd_gaussian_sigma
  logical :: rdfd_random_angular_samples
  integer :: rdfd_n_angular_samples(2)
  real(dp), allocatable :: rdfd_dr(:,:), rdfd_w(:)
  real(dp), allocatable :: rdfds(:,:,:)
  real(dp), allocatable :: rdfd_zone_pos(:), rdfd_bin_pos(:)

end type analysis

public :: analysis, analysis_read, check_analyses, do_analyses, print_analyses

interface reallocate_data
  module procedure reallocate_data_1d, reallocate_data_2d, reallocate_data_3d
end interface reallocate_data

contains

subroutine analysis_read(this, prev, args_str)
  type(analysis), intent(inout) :: this
  type(analysis), intent(in), optional :: prev
  character(len=*), optional, intent(in) :: args_str

  type(Dictionary) :: params
  integer :: dummy_i_1
  character(len=FIELD_LENGTH) :: dummy_c_1, dummy_c_2
  logical :: dummy_l_1, dummy_l_2

  call initialise(params)
  call param_register(params, 'infile', '', dummy_c_1)
  call param_register(params, 'commandfile', '', dummy_c_2)
  call param_register(params, 'decimation', '0', dummy_i_1)
  call param_register(params, 'infile_is_list', 'F', dummy_l_1)
  call param_register(params, 'quiet', 'F', dummy_l_2)

  if (.not. present(prev)) then
    ! general
    this%outfilename=''
    this%mask_str = ''
    call param_register(params, 'outfile', 'stdout', this%outfilename)
    call param_register(params, 'AtomMask', '', this%mask_str)
    call param_register(params, 'min_time', '-1.0', this%min_time)
    call param_register(params, 'max_time', '-1.0', this%max_time)
    call param_register(params, 'min_frame', '-1', this%min_frame)
    call param_register(params, 'max_frame', '-1', this%max_frame)
    call param_register(params, 'density_radial', 'F', this%density_radial)
    call param_register(params, 'density_grid', 'F', this%density_grid)
    call param_register(params, 'rdfd', 'F', this%rdfd)

    ! radial density
    call param_register(params, 'radial_min_p', "0.0", this%radial_min_p)
    call param_register(params, 'radial_bin_width', '-1.0', this%radial_bin_width)
    call param_register(params, 'radial_n_bins', '-1', this%radial_n_bins)
    call param_register(params, 'radial_center', '0.0 0.0 0.0', this%radial_center)
    call param_register(params, 'radial_sigma', '1.0', this%radial_gaussian_sigma)
    call param_register(params, 'radial_random_angular_samples', 'F', this%radial_random_angular_samples)
    call param_register(params, 'radial_n_angular_samples', '2 4', this%radial_n_angular_samples)

    ! grid density
    call param_register(params, 'grid_min_p', '0.0 0.0 0.0', this%grid_min_p)
    call param_register(params, 'grid_bin_width', '-1.0 -1.0 -1.0', this%grid_bin_width)
    call param_register(params, 'grid_n_bins', '-1 -1 -1', this%grid_n_bins)
    call param_register(params, 'grid_gaussian', 'F', this%grid_gaussian_smoothing)
    call param_register(params, 'grid_sigma', '1.0', this%grid_gaussian_sigma)

    ! rdfd
    call param_register(params, 'rdfd_zone_center', '0.0 0.0 0.0', this%rdfd_zone_center)
    call param_register(params, 'rdfd_bin_width', '-1', this%rdfd_bin_width)
    call param_register(params, 'rdfd_n_bins', '-1', this%rdfd_n_bins)
    call param_register(params, 'rdfd_zone_width', '-1.0', this%rdfd_zone_width)
    call param_register(params, 'rdfd_n_zones', '1', this%rdfd_n_zones)
    this%rdfd_center_mask_str=''
    call param_register(params, 'rdfd_center_mask', '', this%rdfd_center_mask_str)
    this%rdfd_neighbour_mask_str=''
    call param_register(params, 'rdfd_neighbour_mask', '', this%rdfd_neighbour_mask_str)
    call param_register(params, 'rdfd_gaussian', 'F', this%rdfd_gaussian_smoothing)
    call param_register(params, 'rdfd_sigma', '0.1', this%rdfd_gaussian_sigma)
    call param_register(params, 'rdfd_random_angular_samples', 'F', this%rdfd_random_angular_samples)
    call param_register(params, 'rdfd_n_angular_samples', '2 4', this%rdfd_n_angular_samples)

  else
    ! general
    call param_register(params, 'outfile', trim(prev%outfilename), this%outfilename)
    call param_register(params, 'AtomMask', trim(prev%mask_str), this%mask_str)
    call param_register(params, 'min_time', ''//prev%min_time, this%min_time)
    call param_register(params, 'max_time', ''//prev%max_time, this%max_time)
    call param_register(params, 'min_frame', ''//prev%min_frame, this%min_frame)
    call param_register(params, 'max_frame', ''//prev%max_frame, this%max_frame)
    call param_register(params, 'density_radial', ''//prev%density_radial, this%density_radial)
    call param_register(params, 'density_grid', ''//prev%density_grid, this%density_grid)
    call param_register(params, 'rdfd', ''//prev%rdfd, this%rdfd)

    ! radial density
    call param_register(params, 'radial_min_p', ''//this%radial_min_p, this%radial_min_p)
    call param_register(params, 'radial_bin_width', ''//this%radial_bin_width, this%radial_bin_width)
    call param_register(params, 'radial_n_bins', ''//this%radial_n_bins, this%radial_n_bins)
    call param_register(params, 'radial_center', ''//prev%radial_center, this%radial_center)
    call param_register(params, 'radial_sigma', ''//prev%radial_gaussian_sigma, this%radial_gaussian_sigma)
    call param_register(params, 'radial_random_angular_samples', ''//prev%radial_random_angular_samples, this%radial_random_angular_samples)
    call param_register(params, 'radial_n_angular_samples', ''//prev%radial_n_angular_samples, this%radial_n_angular_samples)

    ! grid density
    call param_register(params, 'grid_min_p', ''//prev%grid_min_p, this%grid_min_p)
    call param_register(params, 'grid_bin_width', ''//prev%grid_bin_width, this%grid_bin_width)
    call param_register(params, 'grid_n_bins', ''//prev%grid_n_bins, this%grid_n_bins)
    call param_register(params, 'grid_gaussian', ''//prev%grid_gaussian_smoothing, this%grid_gaussian_smoothing)
    call param_register(params, 'grid_sigma', ''//prev%grid_gaussian_sigma, this%grid_gaussian_sigma)

    ! rdfd
    call param_register(params, 'rdfd_zone_center', ''//prev%rdfd_zone_center, this%rdfd_zone_center)
    call param_register(params, 'rdfd_bin_width', ''//prev%rdfd_bin_width, this%rdfd_bin_width)
    call param_register(params, 'rdfd_n_bins', ''//prev%rdfd_n_bins, this%rdfd_n_bins)
    call param_register(params, 'rdfd_zone_width', ''//prev%rdfd_zone_width, this%rdfd_zone_width)
    call param_register(params, 'rdfd_n_zones', ''//prev%rdfd_n_zones, this%rdfd_n_zones)
    call param_register(params, 'rdfd_center_mask', trim(prev%rdfd_center_mask_str), this%rdfd_center_mask_str)
    call param_register(params, 'rdfd_neighbour_mask', trim(prev%rdfd_neighbour_mask_str), this%rdfd_neighbour_mask_str)
    call param_register(params, 'rdfd_gaussian', ''//prev%rdfd_gaussian_smoothing, this%rdfd_gaussian_smoothing)
    call param_register(params, 'rdfd_sigma', ''//prev%rdfd_gaussian_sigma, this%rdfd_gaussian_sigma)
    call param_register(params, 'rdfd_random_angular_samples', ''//prev%rdfd_random_angular_samples, this%rdfd_random_angular_samples)
    call param_register(params, 'rdfd_n_angular_samples', ''//prev%rdfd_n_angular_samples, this%rdfd_n_angular_samples)
  endif

  if (present(args_str)) then
    if (.not. param_read_line(params, trim(args_str), ignore_unknown=.false.)) &
      call system_abort("analysis_read failed to parse string '"//trim(args_str)//"'")
  else
    if (.not. param_read_args(params, do_check=.true.)) &
      call system_abort("analysis_read failed to parse command line arguments")
  endif

  if (count ( (/ this%density_radial, this%density_grid, this%rdfd/) ) /= 1) &
    call system_abort("Specified "//count( (/ this%density_radial, this%density_grid /) )//" types of analysis.  Possiblities: density_radial, density_grid, rdfd")

end subroutine analysis_read

subroutine check_analyses(a)
  type(analysis), intent(inout) :: a(:)

  integer :: i_a

  do i_a=1, size(a)
    if (a(i_a)%density_radial) then
      if (a(i_a)%radial_bin_width <= 0.0_dp) call system_abort("analysis " // i_a // " has radial_bin_width="//a(i_a)%radial_bin_width//" <= 0.0")
      if (a(i_a)%radial_n_bins <= 0) call system_abort("analysis " // i_a // " has radial_n_bins="//a(i_a)%radial_n_bins//" <= 0")
    else if (a(i_a)%density_grid) then
      if (any(a(i_a)%grid_bin_width <= 0.0_dp)) call system_abort("analysis " // i_a // " has grid_bin_width="//a(i_a)%grid_bin_width//" <= 0.0")
      if (any(a(i_a)%grid_n_bins <= 0)) call system_abort("analysis " // i_a // " has grid_n_bins="//a(i_a)%grid_n_bins//" <= 0")
    else if (a(i_a)%rdfd) then
      if (a(i_a)%rdfd_bin_width <= 0.0_dp) call system_abort("analysis " // i_a // " has rdfd_bin_width="//a(i_a)%rdfd_bin_width//" <= 0.0")
      if (a(i_a)%rdfd_n_bins <= 0) call system_abort("analysis " // i_a // " has rdfd_n_bins="//a(i_a)%rdfd_n_bins//" <= 0")
      if (a(i_a)%rdfd_n_zones <= 0) call system_abort("analysis " // i_a // " has rdfd_n_zones="//a(i_a)%rdfd_n_zones//" <= 0")
    else
      call system_abort("check_analyses: no type of analysis set for " // i_a)
    endif
  end do
end subroutine check_analyses

subroutine do_analyses(a, time, frame, at)
  type(analysis), intent(inout) :: a(:)
  real(dp), intent(in) :: time
  integer, intent(in) :: frame
  type(Atoms), intent(inout) :: at

  integer :: i_a

  call map_into_cell(at)
  at%travel = 0

  do i_a=1, size(a)

    if (do_this_analysis(a(i_a), time, frame)) then

      a(i_a)%n_configs = a(i_a)%n_configs + 1

      if (a(i_a)%density_radial) then
        call reallocate_data(a(i_a)%radial_histograms, a(i_a)%n_configs, a(i_a)%radial_n_bins)
        if (a(i_a)%radial_random_angular_samples) then
          call calc_angular_samples_random(a(i_a)%radial_n_angular_samples(1)*a(i_a)%radial_n_angular_samples(2), &
            a(i_a)%radial_dr, a(i_a)%radial_w, a(i_a)%radial_gaussian_sigma)
        endif
        if (a(i_a)%n_configs == 1) then
          if (.not. a(i_a)%radial_random_angular_samples) &
            call calc_angular_samples_grid(a(i_a)%radial_n_angular_samples(1), a(i_a)%radial_n_angular_samples(2), &
              a(i_a)%radial_dr, a(i_a)%radial_w, a(i_a)%radial_gaussian_sigma)
          allocate(a(i_a)%radial_pos(a(i_a)%radial_n_bins))
          call density_sample_radial_mesh_Gaussians(a(i_a)%radial_histograms(:,a(i_a)%n_configs), at, center_pos=a(i_a)%radial_center, &
            rad_bin_width=a(i_a)%radial_bin_width, n_rad_bins=a(i_a)%radial_n_bins, gaussian_sigma=a(i_a)%radial_gaussian_sigma, &
            dr=a(i_a)%radial_dr, w=a(i_a)%radial_w, mask_str=a(i_a)%mask_str, radial_pos=a(i_a)%radial_pos)
        else
          call density_sample_radial_mesh_Gaussians(a(i_a)%radial_histograms(:,a(i_a)%n_configs), at, center_pos=a(i_a)%radial_center, &
            rad_bin_width=a(i_a)%radial_bin_width, n_rad_bins=a(i_a)%radial_n_bins, gaussian_sigma= a(i_a)%radial_gaussian_sigma, &
            dr=a(i_a)%radial_dr, w=a(i_a)%radial_w, mask_str=a(i_a)%mask_str)
        endif
      else if (a(i_a)%density_grid) then
        call reallocate_data(a(i_a)%grid_histograms, a(i_a)%n_configs, a(i_a)%grid_n_bins)
        if (a(i_a)%grid_gaussian_smoothing) then
          if (a(i_a)%n_configs == 1) then
            allocate(a(i_a)%grid_pos(3,a(i_a)%grid_n_bins(1),a(i_a)%grid_n_bins(2),a(i_a)%grid_n_bins(3)))
            call density_sample_rectilinear_mesh_Gaussians(a(i_a)%grid_histograms(:,:,:,a(i_a)%n_configs), at, a(i_a)%grid_min_p, &
              a(i_a)%grid_bin_width, a(i_a)%grid_n_bins, a(i_a)%grid_gaussian_sigma, a(i_a)%mask_str, a(i_a)%grid_pos)
          else
            call density_sample_rectilinear_mesh_Gaussians(a(i_a)%grid_histograms(:,:,:,a(i_a)%n_configs), at, a(i_a)%grid_min_p, &
              a(i_a)%grid_bin_width, a(i_a)%grid_n_bins, a(i_a)%grid_gaussian_sigma, a(i_a)%mask_str)
          endif
        else
          if (a(i_a)%n_configs == 1) then
            allocate(a(i_a)%grid_pos(3,a(i_a)%grid_n_bins(1),a(i_a)%grid_n_bins(2),a(i_a)%grid_n_bins(3)))
            call density_bin_rectilinear_mesh(a(i_a)%grid_histograms(:,:,:,a(i_a)%n_configs), at, a(i_a)%grid_min_p, a(i_a)%grid_bin_width, &
              a(i_a)%grid_n_bins, a(i_a)%mask_str, a(i_a)%grid_pos)
          else
            call density_bin_rectilinear_mesh(a(i_a)%grid_histograms(:,:,:,a(i_a)%n_configs), at, a(i_a)%grid_min_p, a(i_a)%grid_bin_width, &
              a(i_a)%grid_n_bins, a(i_a)%mask_str)
          endif
        endif
      else if (a(i_a)%rdfd) then
        call reallocate_data(a(i_a)%rdfds, a(i_a)%n_configs, (/ a(i_a)%rdfd_n_bins, a(i_a)%rdfd_n_zones /) )
        if (a(i_a)%rdfd_random_angular_samples) then
          call calc_angular_samples_random(a(i_a)%rdfd_n_angular_samples(1)*a(i_a)%rdfd_n_angular_samples(2), &
            a(i_a)%rdfd_dr, a(i_a)%rdfd_w, a(i_a)%rdfd_gaussian_sigma)
        endif
        if (a(i_a)%n_configs == 1) then
          if (.not. a(i_a)%rdfd_random_angular_samples) &
            call calc_angular_samples_grid(a(i_a)%rdfd_n_angular_samples(1), a(i_a)%rdfd_n_angular_samples(2), &
              a(i_a)%rdfd_dr, a(i_a)%rdfd_w, a(i_a)%rdfd_gaussian_sigma)
          allocate(a(i_a)%rdfd_bin_pos(a(i_a)%rdfd_n_bins))
          allocate(a(i_a)%rdfd_zone_pos(a(i_a)%rdfd_n_zones))
          call rdfd_calc(a(i_a)%rdfds(:,:,a(i_a)%n_configs), at, a(i_a)%rdfd_zone_center, a(i_a)%rdfd_bin_width, a(i_a)%rdfd_n_bins, &
            a(i_a)%rdfd_zone_width, a(i_a)%rdfd_n_zones, a(i_a)%rdfd_gaussian_smoothing, a(i_a)%rdfd_gaussian_sigma, &
            a(i_a)%rdfd_dr, a(i_a)%rdfd_w, &
            a(i_a)%rdfd_center_mask_str, a(i_a)%rdfd_neighbour_mask_str, &
            a(i_a)%rdfd_bin_pos, a(i_a)%rdfd_zone_pos)
        else
          call rdfd_calc(a(i_a)%rdfds(:,:,a(i_a)%n_configs), at, a(i_a)%rdfd_zone_center, a(i_a)%rdfd_bin_width, a(i_a)%rdfd_n_bins, &
            a(i_a)%rdfd_zone_width, a(i_a)%rdfd_n_zones, a(i_a)%rdfd_gaussian_smoothing, a(i_a)%rdfd_gaussian_sigma, &
            a(i_a)%rdfd_dr, a(i_a)%rdfd_w, &
            a(i_a)%rdfd_center_mask_str, a(i_a)%rdfd_neighbour_mask_str)
        endif
      else 
        call system_abort("do_analyses: no type of analysis set for " // i_a)
      endif
    end if ! do this analysis
  end do
end subroutine do_analyses

function do_this_analysis(this, time, frame)
  type(analysis), intent(in) :: this
  real(dp), intent(in), optional :: time
  integer, intent(in), optional :: frame
  logical :: do_this_analysis

  do_this_analysis = .true.

if (present(time)) call print("do_this_analysis time " // time // " min " // this%min_time // " max " // this%max_time, ERROR)
  if (this%min_time > 0.0_dp) then
    if (.not. present(time)) call system_abort("analysis has non-zero min_time, but no time specified")
    if (time < 0.0_dp) call system_abort("analysis has non-zero min_time, but invalid time < 0")
    if (time < this%min_time) then
      do_this_analysis = .false.
    endif
  endif
  if (this%max_time > 0.0_dp) then
    if (.not. present(time)) call system_abort("analysis has non-zero max_time, but no time specified")
    if (time < 0.0_dp) call system_abort("analysis has non-zero max_time, but invalid time < 0")
    if (time > this%max_time) then
      do_this_analysis = .false.
    endif
  endif

  if (this%min_frame > 0) then
    if (.not. present(frame)) call system_abort("analysis has non-zero min_frame, but no frame specified")
    if (frame < 0) call system_abort("analysis has non-zero min_frame, but invalid frame < 0")
    if (frame < this%min_frame) then
      do_this_analysis = .false.
    endif
  endif
  if (this%max_frame > 0) then
    if (.not. present(frame)) call system_abort("analysis has non-zero max_frame, but no frame specified")
    if (frame < 0) call system_abort("analysis has non-zero max_frame, but invalid frame < 0")
    if (frame > this%max_frame) then
      do_this_analysis = .false.
    endif
  endif

end function do_this_analysis

subroutine print_analyses(a)
  type(analysis), intent(inout) :: a(:)

  type(inoutput) :: outfile
  integer :: i, i1, i2, i3, i_a
  real(dp), allocatable :: integrated_rdfds(:,:)

  do i_a=1, size(a)
    call initialise(outfile, a(i_a)%outfilename, OUTPUT)
    if (a(i_a)%outfilename == "stdout") then
      outfile%prefix="ANALYSIS_"//i_a
    endif

    if (a(i_a)%n_configs <= 0) then
      call print("# NO DATA", file=outfile)
    else
      if (a(i_a)%density_radial) then
        call print("# radial density histogram", file=outfile)
        call print("n_bins="//a(i_a)%radial_n_bins//" n_data="//a(i_a)%n_configs, file=outfile)
        do i=1, a(i_a)%radial_n_bins
          call print(a(i_a)%radial_pos(i), file=outfile)
        end do
        do i=1, a(i_a)%n_configs
          call print(a(i_a)%radial_histograms(:,i), file=outfile)
        end do
      else if (a(i_a)%density_grid) then
        call print("# grid density histogram", file=outfile)
        call print("n_bins="//a(i_a)%grid_n_bins(1)*a(i_a)%grid_n_bins(2)*a(i_a)%grid_n_bins(3)//" n_data="//a(i_a)%n_configs, file=outfile)
        do i1=1, a(i_a)%grid_n_bins(1)
        do i2=1, a(i_a)%grid_n_bins(2)
        do i3=1, a(i_a)%grid_n_bins(3)
          call print(""//a(i_a)%grid_pos(:,i1,i2,i3), file=outfile)
        end do
        end do
        end do
        do i=1, a(i_a)%n_configs
          call print(""//reshape(a(i_a)%grid_histograms(:,:,:,i), (/ a(i_a)%grid_n_bins(1)*a(i_a)%grid_n_bins(2)*a(i_a)%grid_n_bins(3) /) ), file=outfile)
        end do
      else if (a(i_a)%rdfd) then
        call print("# rdfd", file=outfile)
        call print("n_bins="//a(i_a)%rdfd_n_zones*a(i_a)%rdfd_n_bins//" n_data="//a(i_a)%n_configs, file=outfile)
        do i1=1, a(i_a)%rdfd_n_zones
        do i2=1, a(i_a)%rdfd_n_bins
          if (a(i_a)%rdfd_zone_width > 0.0_dp) then
            call print(""//a(i_a)%rdfd_zone_pos(i1)//" "//a(i_a)%rdfd_bin_pos(i2), file=outfile)
          else
            call print(""//a(i_a)%rdfd_bin_pos(i2), file=outfile)
          endif
        end do
        end do
        do i=1, a(i_a)%n_configs
          call print(""//reshape(a(i_a)%rdfds(:,:,i), (/ a(i_a)%rdfd_n_zones*a(i_a)%rdfd_n_bins /) ), file=outfile)
        end do

        call print("", file=outfile)
        call print("", file=outfile)
        call print("# integrated_rdfd", file=outfile)
        call print("n_bins="//a(i_a)%rdfd_n_zones*a(i_a)%rdfd_n_bins//" n_data="//a(i_a)%n_configs, file=outfile)
        do i1=1, a(i_a)%rdfd_n_zones
        do i2=1, a(i_a)%rdfd_n_bins
          if (a(i_a)%rdfd_zone_width > 0.0_dp) then
            call print(""//a(i_a)%rdfd_zone_pos(i1)//" "//a(i_a)%rdfd_bin_pos(i2), file=outfile)
          else
            call print(""//a(i_a)%rdfd_bin_pos(i2), file=outfile)
          endif
        end do
        end do
	allocate(integrated_rdfds(a(i_a)%rdfd_n_zones,a(i_a)%rdfd_n_bins))
        do i=1, a(i_a)%n_configs
	  integrated_rdfds = 0.0_dp
	  do i1=1, a(i_a)%rdfd_n_zones
	    do i2=2, a(i_a)%rdfd_n_bins
	      integrated_rdfds(i2,i1) = integrated_rdfds(i2-1,i1) + &
		(a(i_a)%rdfd_bin_pos(i2)-a(i_a)%rdfd_bin_pos(i2-1))* &
		4.0_dp*PI*((a(i_a)%rdfd_bin_pos(i2)**2)*a(i_a)%rdfds(i2,i1,i)+(a(i_a)%rdfd_bin_pos(i2-1)**2)*a(i_a)%rdfds(i2-1,i1,i))/2.0_dp
	    end do
	  end do
          call print(""//reshape(integrated_rdfds(:,:), (/ a(i_a)%rdfd_n_zones*a(i_a)%rdfd_n_bins /) ), file=outfile)
        end do
	deallocate(integrated_rdfds)

      else
        call system_abort("print_analyses: no type of analysis set for " // i_a)
      endif
    endif
    call finalise(outfile)
  end do
end subroutine print_analyses

subroutine calc_angular_samples_grid(n_t, n_p, dr, w, gaussian_sigma)
  integer, intent(in) :: n_t, n_p
  real(dp), allocatable, intent(inout) :: dr(:,:), w(:)
  real(dp), intent(in) :: gaussian_sigma

  integer :: t_i, p_i, ii
  real(dp) :: phi, theta

  if (allocated(dr)) deallocate(dr)
  if (allocated(w)) deallocate(w)
  allocate(dr(3,n_t*n_p))
  allocate(w(n_t*n_p))

  ii = 1
  do t_i=1, n_t
  do p_i=1, n_p
    phi = (p_i-1)*2.0_dp*PI/real(n_p,dp)
    if (mod(n_t,2) == 0) then
      theta = (t_i-1.5_dp-floor((n_t-1)/2.0_dp))*PI/real(n_t,dp)
    else
      theta = (t_i-1-floor((n_t-1)/2.0_dp))*PI/real(n_t,dp)
    endif
    ! not oversample top and bottom of spehre
    w(ii) = cos(theta)
    dr(1,ii) = cos(phi)*cos(theta)
    dr(2,ii) = sin(phi)*cos(theta)
    dr(3,ii) = sin(theta)
    ii = ii + 1
  end do
  end do
  w = w / ( (gaussian_sigma*sqrt(2.0_dp*PI))**3 * sum(w) )
end subroutine calc_angular_samples_grid

subroutine calc_angular_samples_random(n, dr, w, gaussian_sigma)
  integer, intent(in) :: n
  real(dp), allocatable, intent(inout) :: dr(:,:), w(:)
  real(dp), intent(in) :: gaussian_sigma

  integer :: ii
  real(dp) :: p(3), p_norm

  if (allocated(dr)) deallocate(dr)
  if (allocated(w)) deallocate(w)
  allocate(dr(3,n))
  allocate(w(n))

  do ii = 1, n
    p_norm = 2.0_dp
    do while (p_norm > 1.0_dp)
      p = 0.0_dp
      call randomise(p, 2.0_dp)
      p_norm = norm(p)
    end do
    p = p/p_norm
    dr(:,ii) = p
    w(ii) = 1.0_dp
  end do
  w = w / ( (gaussian_sigma*sqrt(2.0_dp*PI))**3 * sum(w) )

end subroutine calc_angular_samples_random

subroutine density_sample_radial_mesh_Gaussians(histogram, at, center_pos, center_i, rad_bin_width, n_rad_bins, gaussian_sigma, dr, w, mask_str, radial_pos, accumulate)
  real(dp), intent(inout) :: histogram(:)
  type(Atoms), intent(inout) :: at
  real(dp), intent(in), optional :: center_pos(3)
  integer, intent(in), optional :: center_i
  real(dp), intent(in) :: rad_bin_width
  integer, intent(in) :: n_rad_bins
  real(dp), intent(in) :: gaussian_sigma, dr(:,:), w(:)
  character(len=*), optional, intent(in) :: mask_str
  real(dp), intent(out), optional :: radial_pos(:)
  logical, optional, intent(in) :: accumulate

  logical :: my_accumulate
  real(dp) :: rad_sample_r, p(3), dist, r, t_center(3), exp_arg, diff(3)
  logical, allocatable :: mask_a(:)
  integer at_i, at_i_index, ang_sample_i, rad_sample_i
  real(dp), allocatable, save :: pos_mapped(:,:)

  if (present(center_pos) .and. present(center_i)) &
    call system_abort("density_sample_radial_mesh_Gaussians got both center_pos and center_i")
  if (.not. present(center_pos) .and. .not. present(center_i)) &
    call system_abort("density_sample_radial_mesh_Gaussians got neither center_pos nor center_i")

  my_accumulate = optional_default(.false., accumulate)
  if (.not. my_accumulate) histogram = 0.0_dp

    
  allocate(mask_a(at%N))
  call is_in_mask(mask_a, at, mask_str)

  if (present(radial_pos)) then
    do rad_sample_i=1, n_rad_bins
      rad_sample_r = (rad_sample_i-1)*rad_bin_width
      radial_pos(rad_sample_i) = rad_sample_r
    end do
  endif

  ! ratio of 20/4=5 is bad
  ! ratio of 20/3=6.66 is bad
  ! ratio of 20/2.5=8 is borderline (2e-4)
  ! ratio of 20/2.22=9 is fine  (error 2e-5)
  ! ratio of 20/2=10 is fine
  if ( (norm(at%lattice(:,1)) < 9.0_dp*gaussian_sigma) .or. &
       (norm(at%lattice(:,2)) < 9.0_dp*gaussian_sigma) .or. &
       (norm(at%lattice(:,3)) < 9.0_dp*gaussian_sigma) ) &
    call print("WARNING: at%lattice may be too small for sigma, errors (noticeably too low a density) may result", ERROR)

  if (present(center_pos)) then

    do at_i=1, at%N
      if (.not. mask_a(at_i)) cycle
      r = norm(at%pos(:,at_i))
      do rad_sample_i=1, n_rad_bins
        rad_sample_r = (rad_sample_i-1)*rad_bin_width
        do ang_sample_i=1,size(dr,2)
          p(:) = center_pos(:)+rad_sample_r*dr(:,ang_sample_i)
          dist = distance_min_image(at,p,at%pos(:,at_i))
!Include all the atoms, slow but minimises error
!	  if (dist > 4.0_dp*gaussian_sigma) cycle
          exp_arg = -0.5_dp*(dist/(gaussian_sigma))**2
          if (exp_arg > -20.0_dp) then ! good to about 1e-8
            histogram(rad_sample_i) = histogram(rad_sample_i) + exp(exp_arg)*w(ang_sample_i)
          endif
        end do ! ang_sample_i
      end do ! rad_sample_i
    end do ! at_i

  else ! no center_pos, must have center_i

    ! reallocate pos_mapped if necessary
    if (allocated(pos_mapped)) then
      if (size(pos_mapped,2) /= at%N) deallocate(pos_mapped)
    endif
    if (.not. allocated(pos_mapped)) allocate(pos_mapped(3,at%N))

    ! shift center_i to origin, and map into cell
    t_center=at%pos(:,center_i)
    do at_i=1, at%N
      pos_mapped(:,at_i) = at%pos(:,at_i) - t_center(:)
      call map_into_cell(pos_mapped(:,at_i), at%lattice, at%g)
    end do
!$OMP PARALLEL DO default(shared) private(at_i_index, at_i, rad_sample_i, rad_sample_r, ang_sample_i, p, dist,  exp_arg)
    do at_i_index=1, atoms_n_neighbours(at, center_i)
      at_i = atoms_neighbour(at, center_i, at_i_index, diff=diff)
      if (at_i == center_i) cycle
      if (.not. mask_a(at_i)) cycle
      do rad_sample_i=1, n_rad_bins
        rad_sample_r = (rad_sample_i-1)*rad_bin_width
        do ang_sample_i=1,size(dr,2)
          p(:) = pos_mapped(:,center_i)+rad_sample_r*dr(:,ang_sample_i)
          dist = norm(p-(pos_mapped(:,center_i)+diff))
          exp_arg = -0.5_dp*(dist/(gaussian_sigma))**2
          if (exp_arg > -20.0_dp) then ! good to about 1e-8
            histogram(rad_sample_i) = histogram(rad_sample_i) + exp(exp_arg)*w(ang_sample_i)
          endif
        end do ! ang_sample_i
      end do ! rad_sample_i
    end do
!$OMP END PARALLEL DO
  end if ! present(center_pos)

  deallocate(mask_a)
end subroutine density_sample_radial_mesh_Gaussians

subroutine rdfd_calc(rdfd, at, zone_center, bin_width, n_bins, zone_width, n_zones, gaussian_smoothing, gaussian_sigma, &
                     dr, w, center_mask_str, neighbour_mask_str, bin_pos, zone_pos)
  real(dp), intent(inout) :: rdfd(:,:)
  type(Atoms), intent(inout) :: at
  real(dp), intent(in) :: zone_center(3), bin_width, zone_width
  integer, intent(in) :: n_bins, n_zones
  logical, intent(in) :: gaussian_smoothing
  real(dp), intent(in) :: gaussian_sigma, dr(:,:), w(:)
  character(len=*), intent(in) :: center_mask_str, neighbour_mask_str
  real(dp), intent(inout), optional :: bin_pos(:), zone_pos(:)

  logical, allocatable :: center_mask_a(:), neighbour_mask_a(:)
  integer :: i_at, j_at, i_bin, i_zone
  integer, allocatable :: n_in_zone(:)
  real(dp) :: r, bin_inner_rad, bin_outer_rad

  allocate(center_mask_a(at%N))
  allocate(neighbour_mask_a(at%N))
  call is_in_mask(center_mask_a, at, center_mask_str)
  call is_in_mask(neighbour_mask_a, at, neighbour_mask_str)

  allocate(n_in_zone(n_zones))
  n_in_zone = 0

  if (present(zone_pos)) then
    if (zone_width > 0.0_dp) then
      do i_zone=1, n_zones
        zone_pos(i_zone) = (real(i_zone,dp)-0.5_dp)*zone_width
      end do
    else
      zone_pos(1) = -1.0_dp
    endif
  endif
  if (present(bin_pos)) then
    do i_bin=1, n_bins
      if (gaussian_smoothing) then
        bin_pos(i_bin) = (i_bin-1)*bin_width
      else
        bin_pos(i_bin) = (real(i_bin,dp)-0.5_dp)*bin_width
      endif
    end do
  endif

  if (gaussian_smoothing) then
    call set_cutoff(at, n_bins*bin_width+5.0_dp*gaussian_sigma)
    call calc_connect(at)
  endif

  rdfd = 0.0_dp
  do i_at=1, at%N ! loop over center atoms
    if (.not. center_mask_a(i_at)) cycle

    if (zone_width > 0.0_dp) then
      r = distance_min_image(at, zone_center, at%pos(:,i_at))
      i_zone = int(r/zone_width)+1
      if (i_zone > n_zones) cycle
    else
      i_zone = 1
    endif

    n_in_zone(i_zone) = n_in_zone(i_zone) + 1

    if (gaussian_smoothing) then
      call density_sample_radial_mesh_Gaussians(rdfd(:,i_zone), at, center_i=i_at, rad_bin_width=bin_width, n_rad_bins=n_bins, &
        gaussian_sigma=gaussian_sigma, dr=dr, w=w, mask_str=neighbour_mask_str, accumulate = .true.)
    else
      do j_at=1, at%N
        if (j_at == i_at) cycle
        if (.not. neighbour_mask_a(j_at)) cycle
        r = distance_min_image(at, i_at, j_at)
        i_bin = int(r/bin_width)+1
        if (i_bin <= n_bins) rdfd(i_bin,i_zone) = rdfd(i_bin,i_zone) + 1.0_dp
      end do ! j_at
    endif ! gaussian_smoothing
  end do ! i_at

  do i_zone=1, n_zones
    if (n_in_zone(i_zone) > 0) rdfd(:,i_zone) = rdfd(:,i_zone)/real(n_in_zone(i_zone),dp)
  end do

  if (.not. gaussian_smoothing) then
    do i_bin=1, n_bins
      bin_inner_rad = real(i_bin-1,dp)*bin_width
      bin_outer_rad = real(i_bin,dp)*bin_width
      rdfd(i_bin,:) = rdfd(i_bin,:)/(4.0_dp/3.0_dp*PI*bin_outer_rad**3 - 4.0_dp/3.0_dp*PI*bin_inner_rad**3)
    end do
  endif
  if (count(neighbour_mask_a) > 0) then
    rdfd = rdfd / (count(neighbour_mask_a)/cell_volume(at))
  endif

end subroutine rdfd_calc

subroutine density_sample_rectilinear_mesh_Gaussians(histogram, at, min_p, sample_dist, n_bins, gaussian_sigma, mask_str, grid_pos, accumulate)
  real(dp), intent(inout) :: histogram(:,:,:)
  type(Atoms), intent(in) :: at
  real(dp), intent(in) :: min_p(3), sample_dist(3)
  integer, intent(in) :: n_bins(3)
  real(dp), intent(in) :: gaussian_sigma
  character(len=*), optional, intent(in) :: mask_str
  real(dp), intent(out), optional :: grid_pos(:,:,:,:)
  logical, optional, intent(in) :: accumulate

  logical :: my_accumulate
  integer :: at_i, i1, i2, i3
  real(dp) :: p(3), w, dist
  logical, allocatable :: mask_a(:)

  my_accumulate = optional_default(.false., accumulate)
  if (.not. my_accumulate) histogram = 0.0_dp

  allocate(mask_a(at%N))
  call is_in_mask(mask_a, at, mask_str)

! ratio of 20/4=5 is bad
! ratio of 20/3=6.66 is bad
! ratio of 20/2.5=8 is borderline (2e-4)
! ratio of 20/2.22=9 is fine  (error 2e-5)
! ratio of 20/2=10 is fine
  if ( (norm(at%lattice(:,1)) < 9.0_dp*gaussian_sigma) .or. &
       (norm(at%lattice(:,2)) < 9.0_dp*gaussian_sigma) .or. &
       (norm(at%lattice(:,3)) < 9.0_dp*gaussian_sigma) ) &
    call print("WARNING: at%lattice may be too small for sigma, errors (noticeably too low a density) may result", ERROR)

  w = 1.0_dp / (gaussian_sigma*sqrt(2.0_dp*PI))**3

  if (present(grid_pos)) then
    do i1=1, n_bins(1)
      p(1) = min_p(1) + (i1-1)*sample_dist(1)
      do i2=1, n_bins(2)
        p(2) = min_p(2) + (i2-1)*sample_dist(2)
        do i3=1, n_bins(3)
          p(3) = min_p(3) + (i3-1)*sample_dist(3)
          grid_pos(:,i1,i2,i3) = p
        end do
      end do
    end do
  endif

  do at_i=1, at%N
    if (.not. mask_a(at_i)) cycle
    do i1=1, n_bins(1)
      p(1) = min_p(1) + (i1-1)*sample_dist(1)
      do i2=1, n_bins(2)
        p(2) = min_p(2) + (i2-1)*sample_dist(2)
        do i3=1, n_bins(3)
          p(3) = min_p(3) + (i3-1)*sample_dist(3)
          dist = distance_min_image(at,p,at%pos(:,at_i))
          histogram(i1,i2,i3) = histogram(i1,i2,i3) + exp(-0.5_dp*(dist/(gaussian_sigma))**2)*w
        end do ! i3
      end do ! i2
    end do ! i1
  end do ! at_i

  deallocate(mask_a)
end subroutine density_sample_rectilinear_mesh_Gaussians

subroutine density_bin_rectilinear_mesh(histogram, at, min_p, bin_width, n_bins, mask_str, grid_pos, accumulate)
  real(dp), intent(inout) :: histogram(:,:,:)
  type(Atoms), intent(in) :: at
  real(dp), intent(in) :: min_p(3), bin_width(3)
  integer, intent(in) :: n_bins(3)
  character(len=*), optional, intent(in) :: mask_str
  real(dp), intent(out), optional :: grid_pos(:,:,:,:)
  logical, optional, intent(in) :: accumulate

  logical :: my_accumulate
  integer :: i, i1, i2, i3, bin(3)
  logical, allocatable :: mask_a(:)
  real(dp) :: p(3)

  my_accumulate = optional_default(.false., accumulate)
  if (.not. my_accumulate) histogram = 0.0_dp

  if (present(grid_pos)) then
    do i1=1, n_bins(1)
      p(1) = min_p(1) + (real(i1,dp)-0.5_dp)*bin_width(1)
      do i2=1, n_bins(2)
        p(2) = min_p(2) + (real(i2,dp)-0.5_dp)*bin_width(2)
        do i3=1, n_bins(3)
          p(3) = min_p(3) + (real(i3,dp)-0.5_dp)*bin_width(3)
          grid_pos(:,i1,i2,i3) = p
        end do
      end do
    end do
  endif

  allocate(mask_a(at%N))
  call is_in_mask(mask_a, at, mask_str)

  do i=1, at%N
    if (.not. mask_a(i)) cycle
    bin = floor((at%pos(:,i)-min_p)/bin_width)+1
    if (all(bin >= 1) .and. all (bin <= n_bins)) histogram(bin(1),bin(2),bin(3)) = histogram(bin(1),bin(2),bin(3)) + 1.0_dp
  end do

  deallocate(mask_a)

end subroutine density_bin_rectilinear_mesh

subroutine is_in_mask(mask_a, at, mask_str)
  type(Atoms), intent(in) :: at
  logical, intent(out) :: mask_a(at%N)
  character(len=*), optional, intent(in) :: mask_str

  integer :: i_at, i_Z
  integer :: Zmask
  type(Table) :: atom_indices
  character(len=4) :: species(128)
  integer :: n_species

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
    do i_at=1, atom_indices%N
      mask_a(atom_indices%int(1,i_at)) = .true.
    end do
  else if (scan(mask_str,'=')/=0) then
    call system_abort("property type mask not supported yet")
  else
    call split_string(mask_str, ' ,', '""', species, n_species)
    do i_Z=1, n_species
      Zmask = Atomic_Number(species(i_Z))
      do i_at=1, at%N
        if (at%Z(i_at) == Zmask) mask_a(i_at) = .true.
      end do
    end do
  end if
end subroutine is_in_mask

subroutine reallocate_data_1d(data, n, n_bins)
  real(dp), allocatable, intent(inout) :: data(:,:)
  integer, intent(in) :: n, n_bins

  integer :: new_size
  real(dp), allocatable :: t_data(:,:)

  if (allocated(data)) then
    if (n <= size(data,2)) return
    allocate(t_data(size(data,1),size(data,2)))
    t_data = data
    deallocate(data)
    if (size(t_data,2) <= 0) then
      new_size = 10
    else if (size(t_data,2) < 1000) then
      new_size = 2*size(t_data,2)
    else
      new_size = floor(1.25*size(t_data,2))
    endif
    allocate(data(size(t_data,1), new_size))
    data(1:size(t_data,1),1:size(t_data,2)) = t_data(1:size(t_data,1),1:size(t_data,2))
    data(1:size(t_data,1),size(t_data,2)+1:size(data,2)) = 0.0_dp
    deallocate(t_data)
  else
    allocate(data(n_bins, n))
  endif
end subroutine reallocate_data_1d

subroutine reallocate_data_2d(data, n, n_bins)
  real(dp), allocatable, intent(inout) :: data(:,:,:)
  integer, intent(in) :: n, n_bins(2)
 
  integer :: new_size
  real(dp), allocatable :: t_data(:,:,:)

  if (allocated(data)) then
    if (n <= size(data,3)) return
    allocate(t_data(size(data,1),size(data,2),size(data,3)))
    t_data = data
    deallocate(data)
    if (size(t_data,3) <= 0) then
      new_size = 10
    else if (size(t_data,3) < 1000) then
      new_size = 2*size(t_data,3)
    else
      new_size = floor(1.25*size(t_data,3))
    endif
    allocate(data(size(t_data,1), size(t_data,2), new_size))
    data(1:size(t_data,1),1:size(t_data,2),1:size(t_data,3)) = &
      t_data(1:size(t_data,1),1:size(t_data,2),1:size(t_data,3))
    data(1:size(t_data,1),1:size(t_data,2),size(t_data,3)+1:size(data,3)) = 0.0_dp
    deallocate(t_data)
  else
    allocate(data(n_bins(1), n_bins(2), n))
  endif
end subroutine reallocate_data_2d

subroutine reallocate_data_3d(data, n, n_bins)
  real(dp), allocatable, intent(inout) :: data(:,:,:,:)
  integer, intent(in) :: n, n_bins(3)

  integer :: new_size
  real(dp), allocatable :: t_data(:,:,:,:)

  if (allocated(data)) then
    if (n <= size(data,4)) return
    allocate(t_data(size(data,1),size(data,2),size(data,3),size(data,4)))
    t_data = data
    deallocate(data)
    if (size(t_data,4) <= 0) then
      new_size = 10
    else if (size(t_data,4) < 1000) then
      new_size = 2*size(t_data,4)
    else
      new_size = floor(1.25*size(t_data,4))
    endif
    allocate(data(size(t_data,1), size(t_data,2), size(t_data,3), new_size))
    data(1:size(t_data,1),1:size(t_data,2),1:size(t_data,3),1:size(t_data,4)) = &
      t_data(1:size(t_data,1),1:size(t_data,2),1:size(t_data,3),1:size(t_data,4))
    data(1:size(t_data,1),1:size(t_data,2),1:size(t_data,3),size(t_data,4)+1:size(data,4)) = 0.0_dp
    deallocate(t_data)
  else
    allocate(data(n_bins(1), n_bins(2), n_bins(3), n))
  endif
end subroutine reallocate_data_3d

end module structure_analysis_module

program structure_analysis
use libatoms_module
use structure_analysis_module
implicit none

  type(Dictionary) :: cli_params

  character(len=FIELD_LENGTH) :: infilename
  integer :: decimation
  type(Inoutput) :: list_infile
  type (CInoutput) :: infile
  logical :: infile_is_list
  logical :: quiet

  character(len=FIELD_LENGTH) :: commandfilename
  type(Inoutput) :: commandfile

  character(len=10240) :: args_str, myline
  integer :: n_analysis_a
  type(analysis), allocatable :: analysis_a(:)

  logical :: more_files
  integer :: status, arg_line_no
  real(dp) :: time

  integer :: i_a, frame_count, raw_frame_count
  type(Atoms) :: structure
  logical :: do_verbose

  call system_initialise(NORMAL)

  call initialise(cli_params)
  call param_register(cli_params, "verbose", "F", do_verbose)
  if (.not. param_read_args(cli_params, ignore_unknown=.true.)) &
    call system_abort("Impossible failure to parse verbosity")
  call finalise(cli_params)
  if (do_verbose) then
    call system_initialise(verbosity=VERBOSE)
  endif

  call initialise(cli_params)
  commandfilename=''
  call param_register(cli_params, "commandfile", '', commandfilename)
  call param_register(cli_params, "infile", "stdin", infilename)
  call param_register(cli_params, "decimation", "1", decimation)
  call param_register(cli_params, "infile_is_list", "F", infile_is_list)
  call param_register(cli_params, "quiet", "F", quiet)
  if (.not. param_read_args(cli_params, ignore_unknown = .true., task="CLI")) then
    call system_abort("Failed to parse CLI")
  endif
  if (len_trim(commandfilename) == 0) then
    allocate(analysis_a(1))
    call analysis_read(analysis_a(1))
  else
    if (.not. param_read_args(cli_params, ignore_unknown = .false., task="CLI_again")) then
      call system_abort("Failed to parse CLI again after getting a commandfile, most likely passed in some analysis specific flags on command line")
    endif
    call initialise(commandfile, trim(commandfilename), INPUT)
    arg_line_no = 1
    myline = read_line(commandfile)
    read (unit=myline,fmt=*) n_analysis_a
    allocate(analysis_a(n_analysis_a))
    do i_a=1, n_analysis_a
      args_str = read_line(commandfile)
      if (i_a == 1) then
        call analysis_read(analysis_a(i_a), args_str=args_str)
      else
        call analysis_read(analysis_a(i_a), analysis_a(i_a-1), args_str)
      endif
    end do
  endif
  call finalise(cli_params)

  call check_analyses(analysis_a)

  more_files = .true.
  if (infile_is_list) then
    call initialise(list_infile, trim(infilename), INPUT)
    infilename = read_line(list_infile, status)
    more_files = .false.
    if (status == 0) more_files = .true.
  endif

  raw_frame_count = decimation-1
  frame_count = 0
  do while (more_files)

    call print(trim(infilename))
    call initialise(infile, infilename, INPUT)
    call read(structure, infile, status=status, frame=raw_frame_count)
    do while (status == 0)
      frame_count = frame_count + 1
      if (.not. quiet) then
        if (mod(frame_count,1000) == 1) write (mainlog%unit,'(I7,a,$)') frame_count," "
        if (mod(frame_count,10) == 0) write (mainlog%unit,'(I1,$)') mod(frame_count/10,10)
        if (mod(frame_count,1000) == 0) write (mainlog%unit,'(a)') " "
      endif

      if (.not. get_value(structure%params,"Time",time)) then
        time = -1.0_dp
      endif
      call do_analyses(analysis_a, time, frame_count, structure)
      call finalise(structure)

      raw_frame_count = raw_frame_count + decimation
      ! get ready for next structure
      call read(structure, infile, status=status, frame=raw_frame_count)
    end do
    raw_frame_count = raw_frame_count - decimation

    ! get ready for next file
    more_files = .false.
    call finalise(infile)
    if (infile_is_list) then
      infilename = read_line(list_infile, status)
      if (status == 0) then
        more_files = .true.
      endif
      if (infile%n_frame > 0) then
	raw_frame_count = (decimation-1)-(infile%n_frame-1-raw_frame_count)
      else
	raw_frame_count = decimation-1
      endif
    endif
    write (mainlog%unit,'(a)') " "
    frame_count = 0
  end do ! more_files
  if (infile_is_list) call finalise(list_infile)


  call print_analyses(analysis_a)

end program structure_analysis
