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

! do various structural analysis to a trajectory, save an intermediate file
! to be postprocessed (mean, variance, correlation, etc)
! #/vol density on a radial mesh
! #/vol density on a grid
! RDF, possibly as a function of distance of center from a fixed point
! to be added: ADF

module structure_analysis_module
use libatoms_module
use structure_analysis_routines_module
implicit none
private

type analysis
  character(len=STRING_LENGTH) :: type ! type of analysis, used to fill in list of logicals
  logical :: density_radial, density_grid, KE_density_radial, rdfd, integrated_rdfd, adfd, & !water
             KEdf_radial, propdf_radial, geometry, & !general
             density_axial_silica, num_hbond_silica, water_orientation_silica !silica-water interface
  character(len=STRING_LENGTH) :: outfilename
  character(len=STRING_LENGTH) :: mask_str

  integer :: n_configs = 0

  real(dp) :: min_time, max_time
  integer :: min_frame, max_frame

  ! radial density stuff
  real(dp) :: density_radial_bin_width
  integer :: density_radial_n_bins
  integer :: density_radial_center_at
  real(dp) :: density_radial_center(3)
  real(dp) :: density_radial_gaussian_sigma
  real(dp), allocatable :: density_radial_histograms(:,:)
  real(dp), allocatable :: density_radial_pos(:)

  ! density grid stuff
  real(dp) :: density_grid_min_p(3), density_grid_bin_width(3)
  integer :: density_grid_n_bins(3)
  logical :: density_grid_gaussian_smoothing
  real(dp) :: density_grid_gaussian_sigma
  real(dp), allocatable :: density_grid_histograms(:,:,:,:)
  real(dp), allocatable :: density_grid_pos(:,:,:,:)

  ! radial KE density stuff
  real(dp) :: KE_density_radial_bin_width
  integer :: KE_density_radial_n_bins
  integer :: KE_density_radial_center_at
  real(dp) :: KE_density_radial_center(3)
  real(dp) :: KE_density_radial_gaussian_sigma
  real(dp), allocatable :: KE_density_radial_histograms(:,:)
  real(dp), allocatable :: KE_density_radial_pos(:)

  ! rdfd stuff
  real(dp) :: rdfd_zone_center(3)
  integer :: rdfd_zone_atom_center
  character(STRING_LENGTH) :: rdfd_center_mask_str, rdfd_neighbour_mask_str
  real(dp) :: rdfd_zone_width, rdfd_bin_width
  integer :: rdfd_n_zones, rdfd_n_bins
  logical :: rdfd_gaussian_smoothing
  real(dp) :: rdfd_gaussian_sigma
  real(dp), allocatable :: rdfds(:,:,:)
  real(dp), allocatable :: rdfd_zone_pos(:), rdfd_bin_pos(:)

  ! adfd stuff
  real(dp) :: adfd_zone_center(3)
  character(STRING_LENGTH) :: adfd_center_mask_str, adfd_neighbour_1_mask_str, adfd_neighbour_2_mask_str
  logical :: adfd_dist_bin_rc2
  real(dp) :: adfd_neighbour_1_max_dist
  real(dp) :: adfd_zone_width, adfd_dist_bin_width
  integer :: adfd_n_zones, adfd_n_dist_bins, adfd_n_angle_bins
  real(dp), allocatable :: adfds(:,:,:,:)
  real(dp), allocatable :: adfd_zone_pos(:), adfd_dist_bin_pos(:), adfd_angle_bin_pos(:)

  ! r-dep KE distribution
  character(STRING_LENGTH) :: KEdf_radial_mask_str
  real(dp) :: KEdf_radial_gaussian_sigma
  real(dp) :: KEdf_radial_zone_center(3)
  integer :: KEdf_radial_zone_center_at
  real(dp) :: KEdf_radial_zone_width, KEdf_radial_bin_width
  integer :: KEdf_radial_n_zones, KEdf_radial_n_bins
  real(dp), allocatable :: KEdf_radial_bin_pos(:), KEdf_radial_zone_pos(:), KEdf_radial_histograms(:,:,:)

  ! r-dep |F| distribution
  character(STRING_LENGTH) :: propdf_radial_mask_str
  real(dp) :: propdf_radial_gaussian_sigma
  real(dp) :: propdf_radial_zone_center(3)
  integer :: propdf_radial_zone_center_at
  real(dp) :: propdf_radial_zone_width, propdf_radial_bin_width
  integer :: propdf_radial_n_zones, propdf_radial_n_bins
  real(dp), allocatable :: propdf_radial_bin_pos(:), propdf_radial_zone_pos(:), propdf_radial_histograms(:,:,:)
  character(len=STRING_LENGTH) :: propdf_radial_property

  !geometry
  character(STRING_LENGTH) :: geometry_filename
  type(Table) :: geometry_params
  integer :: geometry_central_atom
  real(dp), allocatable :: geometry_histograms(:,:)
  real(dp), allocatable :: geometry_pos(:)
  character(STRING_LENGTH), allocatable :: geometry_label(:)

  !uniaxial density - for silica-water interface
  integer :: density_axial_axis !x:1,y:2,z:3
  real(dp) :: density_axial_bin_width
  integer :: density_axial_n_bins
  integer :: density_axial_silica_atoms
  real(dp) :: density_axial_gaussian_sigma
  logical :: density_axial_gaussian_smoothing
  real(dp), allocatable :: density_axial_histograms(:,:)
  real(dp), allocatable :: density_axial_pos(:)

  !silica_numhb - hydrogen bond distribution for silica-water interface
  integer :: num_hbond_axis !x:1,y:2,z:3
  integer :: num_hbond_silica_atoms
  integer :: num_hbond_n_type=4
  integer :: num_hbond_n_bins
  logical :: num_hbond_gaussian_smoothing
  real(dp) :: num_hbond_gaussian_sigma
  real(dp), allocatable :: num_hbond_histograms(:,:,:)
  real(dp), allocatable :: integrated_num_hbond_histograms(:,:)
  integer, allocatable :: num_hbond_type_code(:)
  real(dp), allocatable :: num_hbond_bin_pos(:)
  character(STRING_LENGTH), allocatable :: num_hbond_type_label(:)

  !silica_water_orientation - water orientation distribution for silica-water interface
  integer :: water_orientation_axis !x:1,y:2,z:3
  integer :: water_orientation_silica_atoms
  integer :: water_orientation_n_angle_bins
  integer :: water_orientation_n_pos_bins
  logical :: water_orientation_gaussian_smoothing
  real(dp) :: water_orientation_pos_gaussian_sigma
  !real(dp) :: water_orientation_angle_gaussian_sigma
  logical :: water_orientation_use_dipole
  logical :: water_orientation_use_HOHangle_bisector
  real(dp), allocatable :: water_orientation_histograms(:,:,:)
  real(dp), allocatable :: integrated_water_orientation_histograms(:,:)
  real(dp), allocatable :: water_orientation_pos_bin(:)
  real(dp), allocatable :: water_orientation_angle_bin(:)
  real(dp), allocatable :: water_orientation_angle_bin_w(:)

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
  character(len=STRING_LENGTH) :: dummy_c_1, dummy_c_2
  logical :: dummy_l_1, dummy_l_2

  call initialise(params)
  ! dummy parameters required for some (undocumented) reason...
  call param_register(params, 'infile', '', dummy_c_1, help_string="(Dummy parameter)")
  call param_register(params, 'commandfile', '', dummy_c_2, help_string="(Dummy parameter)")
  call param_register(params, 'decimation', '0', dummy_i_1, help_string="(Dummy parameter)")
  call param_register(params, 'infile_is_list', 'F', dummy_l_1, help_string="(Dummy parameter)")
  call param_register(params, 'quiet', 'F', dummy_l_2, help_string="(Dummy parameter)")

  if (.not. present(prev)) then
    ! general
    call param_register(params, 'outfile', 'stdout', this%outfilename, help_string="Output file name (default: stdout)")
    call param_register(params, 'AtomMask', '', this%mask_str, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'min_time', '-1.0', this%min_time, help_string="Only include frames with Time >= min_time")
    call param_register(params, 'max_time', '-1.0', this%max_time, help_string="Only include frames with Time <= max_time")
    call param_register(params, 'min_frame', '-1', this%min_frame, help_string="Only include frames with number (not counting skipped) >= min_frame")
    call param_register(params, 'max_frame', '-1', this%max_frame, help_string="Only include frames with number (not counting skipped) <= max_frame")
    call param_register(params, 'type', '', this%type, help_string="[density_radial | density_grid | KE_density_radial | rdfd | adfd | KEdf_radial | propdf_radial | geometry | density_axial_silica | num_hbond_silica | water_orientation_silica ]")
  else
    ! general
    call param_register(params, 'outfile', trim(prev%outfilename), this%outfilename, help_string="")
    call param_register(params, 'AtomMask', trim(prev%mask_str), this%mask_str, help_string="")
    call param_register(params, 'min_time', ''//prev%min_time, this%min_time, help_string="")
    call param_register(params, 'max_time', ''//prev%max_time, this%max_time, help_string="")
    call param_register(params, 'min_frame', ''//prev%min_frame, this%min_frame, help_string="")
    call param_register(params, 'max_frame', ''//prev%max_frame, this%max_frame, help_string="")
    call param_register(params, 'type', ''//trim(prev%type), this%type, help_string="[density_radial | density_grid | KE_density_radial | rdfd | adfd | KEdf_radial | propdf_radial | geometry | density_axial_silica | num_hbond_silica | water_orientation_silica ]")
  endif

  if (.not. present(prev)) then
    ! radial density
    call param_register(params, 'density_radial_bin_width', '-1.0', this%density_radial_bin_width, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'density_radial_n_bins', '-1', this%density_radial_n_bins, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'density_radial_center', '0.0 0.0 0.0', this%density_radial_center, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'density_radial_center_at', '-1', this%density_radial_center_at, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'density_radial_sigma', '1.0', this%density_radial_gaussian_sigma, help_string="No help yet.  This source file was $LastChangedBy$")
  else
    ! radial density
    call param_register(params, 'density_radial_bin_width', ''//this%density_radial_bin_width, this%density_radial_bin_width, help_string="")
    call param_register(params, 'density_radial_n_bins', ''//this%density_radial_n_bins, this%density_radial_n_bins, help_string="")
    call param_register(params, 'density_radial_center', ''//prev%density_radial_center, this%density_radial_center, help_string="")
    call param_register(params, 'density_radial_center_at', ''//prev%density_radial_center_at, this%density_radial_center_at, help_string="")
    call param_register(params, 'density_radial_sigma', ''//prev%density_radial_gaussian_sigma, this%density_radial_gaussian_sigma, help_string="")
  endif

  if (.not. present(prev)) then
    ! grid density
    call param_register(params, 'density_grid_min_p', '0.0 0.0 0.0', this%density_grid_min_p, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'density_grid_bin_width', '-1.0 -1.0 -1.0', this%density_grid_bin_width, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'density_grid_n_bins', '-1 -1 -1', this%density_grid_n_bins, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'density_grid_gaussian', 'F', this%density_grid_gaussian_smoothing, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'density_grid_sigma', '1.0', this%density_grid_gaussian_sigma, help_string="No help yet.  This source file was $LastChangedBy$")
  else
    ! grid density
    call param_register(params, 'density_grid_min_p', ''//prev%density_grid_min_p, this%density_grid_min_p, help_string="")
    call param_register(params, 'density_grid_bin_width', ''//prev%density_grid_bin_width, this%density_grid_bin_width, help_string="")
    call param_register(params, 'density_grid_n_bins', ''//prev%density_grid_n_bins, this%density_grid_n_bins, help_string="")
    call param_register(params, 'density_grid_gaussian', ''//prev%density_grid_gaussian_smoothing, this%density_grid_gaussian_smoothing, help_string="")
    call param_register(params, 'density_grid_sigma', ''//prev%density_grid_gaussian_sigma, this%density_grid_gaussian_sigma, help_string="")
  endif

  if (.not. present(prev)) then
    ! radial KE density
    call param_register(params, 'KE_density_radial_bin_width', '-1.0', this%KE_density_radial_bin_width, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'KE_density_radial_n_bins', '-1', this%KE_density_radial_n_bins, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'KE_density_radial_center', '0.0 0.0 0.0', this%KE_density_radial_center, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'KE_density_radial_center_at', '-1', this%KE_density_radial_center_at, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'KE_density_radial_sigma', '1.0', this%KE_density_radial_gaussian_sigma, help_string="No help yet.  This source file was $LastChangedBy$")
  else
    ! radial KE_density
    !TODO why do these use 'this' for the default?
    call param_register(params, 'KE_density_radial_bin_width', ''//this%KE_density_radial_bin_width, this%KE_density_radial_bin_width, help_string="")
    call param_register(params, 'KE_density_radial_n_bins', ''//this%KE_density_radial_n_bins, this%KE_density_radial_n_bins, help_string="")
    !TODO why is KE_density_radial_center missing?
    call param_register(params, 'KE_density_radial_center_at', ''//prev%KE_density_radial_center_at, this%KE_density_radial_center_at, help_string="")
    call param_register(params, 'KE_density_radial_sigma', ''//prev%KE_density_radial_gaussian_sigma, this%KE_density_radial_gaussian_sigma, help_string="")
  endif

  if (.not. present(prev)) then
    ! rdfd
    call param_register(params, 'rdfd_zone_center', '0.0 0.0 0.0', this%rdfd_zone_center, help_string="Cartesian position to center zones for rdfd (ignored if rdfd_zone_atom_center is given)")
    call param_register(params, 'rdfd_zone_atom_center', '0', this%rdfd_zone_atom_center, help_string="If > 0, index of atom to center zones for rdfd")
    call param_register(params, 'rdfd_zone_width', '-1.0', this%rdfd_zone_width, help_string="Width of zones for rdfd")
    call param_register(params, 'rdfd_n_zones', '1', this%rdfd_n_zones, help_string="Number of zones for rdfd")
    call param_register(params, 'rdfd_bin_width', '-1', this%rdfd_bin_width, help_string="Width of bin (or distance between sample points) in each rdfd")
    call param_register(params, 'rdfd_n_bins', '-1', this%rdfd_n_bins, help_string="Number of bins (or sample points) in each rdfd")
    call param_register(params, 'rdfd_center_mask', '', this%rdfd_center_mask_str, help_string="Mask for atoms that are considered for rdfd centers")
    call param_register(params, 'rdfd_neighbour_mask', '', this%rdfd_neighbour_mask_str, help_string="Mask for atoms that are considered for rdfd neighbours")
    call param_register(params, 'rdfd_gaussian', 'F', this%rdfd_gaussian_smoothing, help_string="If true, use Gaussians to smear mass distributions being integrated for rdfd")
    call param_register(params, 'rdfd_sigma', '0.1', this%rdfd_gaussian_sigma, help_string="Width of Gaussians used to smear mass distributions for rdfd")
  else
    ! rdfd
    call param_register(params, 'rdfd_zone_center', ''//prev%rdfd_zone_center, this%rdfd_zone_center, help_string="")
    call param_register(params, 'rdfd_zone_atom_center', ''//prev%rdfd_zone_atom_center, this%rdfd_zone_atom_center, help_string="")
    call param_register(params, 'rdfd_zone_width', ''//prev%rdfd_zone_width, this%rdfd_zone_width, help_string="")
    call param_register(params, 'rdfd_n_zones', ''//prev%rdfd_n_zones, this%rdfd_n_zones, help_string="")
    call param_register(params, 'rdfd_bin_width', ''//prev%rdfd_bin_width, this%rdfd_bin_width, help_string="")
    call param_register(params, 'rdfd_n_bins', ''//prev%rdfd_n_bins, this%rdfd_n_bins, help_string="")
    call param_register(params, 'rdfd_center_mask', trim(prev%rdfd_center_mask_str), this%rdfd_center_mask_str, help_string="")
    call param_register(params, 'rdfd_neighbour_mask', trim(prev%rdfd_neighbour_mask_str), this%rdfd_neighbour_mask_str, help_string="")
    call param_register(params, 'rdfd_gaussian', ''//prev%rdfd_gaussian_smoothing, this%rdfd_gaussian_smoothing, help_string="")
    call param_register(params, 'rdfd_sigma', ''//prev%rdfd_gaussian_sigma, this%rdfd_gaussian_sigma, help_string="")
  endif

  if (.not. present(prev)) then
    ! adfd
    call param_register(params, 'adfd_zone_center', '0.0 0.0 0.0', this%adfd_zone_center, help_string="Cartesian position to center zones for adfd (ignored if rdfd_zone_atom_center is given)")
    call param_register(params, 'adfd_zone_width', '-1.0', this%adfd_zone_width, help_string="Width of zones for adfd")
    call param_register(params, 'adfd_n_zones', '1', this%adfd_n_zones, help_string="Number of zones for adfd")
    call param_register(params, 'adfd_dist_bin_width', '-1', this%adfd_dist_bin_width, help_string="Width of distance bin in adfd")
    call param_register(params, 'adfd_n_dist_bins', '-1', this%adfd_n_dist_bins, help_string="Number of distance bins in adfd")
    call param_register(params, 'adfd_n_angle_bins', '-1', this%adfd_n_angle_bins, help_string="Number of angular bins in adfd")
    call param_register(params, 'adfd_center_mask', '', this%adfd_center_mask_str, help_string="Mask for atoms that are considered for adfd centers")
    call param_register(params, 'adfd_neighbour_1_mask', '', this%adfd_neighbour_1_mask_str, help_string="Mask for atoms that are considered for the first neighbour in adfd")
    call param_register(params, 'adfd_neighbour_1_max_dist', '-1.0', this%adfd_neighbour_1_max_dist, help_string="Cutoff distance for atoms considered as the first neighbour")
    call param_register(params, 'adfd_neighbour_2_mask', '', this%adfd_neighbour_2_mask_str, help_string="Mask for atoms that are considered for the second neighbour in adfd")
    call param_register(params, 'adfd_dist_bin_rc2', '', this%adfd_dist_bin_rc2, help_string="If true, apply distance binning to r_center-neighbour2, otherwise to r_neighbour1-neighbour2")
  else
    ! adfd
    call param_register(params, 'adfd_zone_center', ''//prev%adfd_zone_center, this%adfd_zone_center, help_string="")
    call param_register(params, 'adfd_zone_width', ''//prev%adfd_zone_width, this%adfd_zone_width, help_string="")
    call param_register(params, 'adfd_n_zones', ''//prev%adfd_n_zones, this%adfd_n_zones, help_string="")
    call param_register(params, 'adfd_dist_bin_width', ''//prev%adfd_dist_bin_width, this%adfd_dist_bin_width, help_string="")
    call param_register(params, 'adfd_n_dist_bins', ''//prev%adfd_n_dist_bins, this%adfd_n_dist_bins, help_string="")
    call param_register(params, 'adfd_n_angle_bins', ''//prev%adfd_n_angle_bins, this%adfd_n_angle_bins, help_string="")
    call param_register(params, 'adfd_center_mask', trim(prev%adfd_center_mask_str), this%adfd_center_mask_str, help_string="")
    call param_register(params, 'adfd_neighbour_1_mask', trim(prev%adfd_neighbour_1_mask_str), this%adfd_neighbour_1_mask_str, help_string="")
    call param_register(params, 'adfd_neighbour_1_max_dist', ''//prev%adfd_neighbour_1_max_dist, this%adfd_neighbour_1_max_dist, help_string="")
    call param_register(params, 'adfd_neighbour_2_mask', trim(prev%adfd_neighbour_2_mask_str), this%adfd_neighbour_2_mask_str, help_string="")
    call param_register(params, 'adfd_dist_bin_rc2', ''//prev%adfd_dist_bin_rc2, this%adfd_dist_bin_rc2, help_string="If true, apply distance binning to r_center-neighbor2, otherwise to r_neighbor1-neighbor2")
  endif

  if (.not. present(prev)) then
    ! r-dep KE distribution
    call param_register(params, 'KEdf_radial_zone_center', '0.0 0.0 0.0', this%KEdf_radial_zone_center, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'KEdf_radial_zone_center_at', '-1', this%KEdf_radial_zone_center_at, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'KEdf_radial_zone_width', '-1.0', this%KEdf_radial_zone_width, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'KEdf_radial_bin_width', '0.002', this%KEdf_radial_bin_width, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'KEdf_radial_n_bins', '500', this%KEdf_radial_n_bins, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'KEdf_radial_n_zones', '1', this%KEdf_radial_n_zones, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'KEdf_radial_gaussian_sigma', '0.005', this%KEdf_radial_gaussian_sigma, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'KEdf_radial_mask', '', this%KEdf_radial_mask_str, help_string="No help yet.  This source file was $LastChangedBy$")
  else
    ! r-dep KE distribution
    call param_register(params, 'KEdf_radial_zone_center', ''//prev%KEdf_radial_zone_center, this%KEdf_radial_zone_center, help_string="")
    call param_register(params, 'KEdf_radial_zone_center_at', ''//prev%KEdf_radial_zone_center_at, this%KEdf_radial_zone_center_at, help_string="")
    call param_register(params, 'KEdf_radial_zone_width', ''//prev%KEdf_radial_zone_width, this%KEdf_radial_zone_width, help_string="")
    call param_register(params, 'KEdf_radial_bin_width', ''//prev%KEdf_radial_bin_width, this%KEdf_radial_bin_width, help_string="")
    call param_register(params, 'KEdf_radial_n_bins', ''//prev%KEdf_radial_n_bins, this%KEdf_radial_n_bins, help_string="")
    call param_register(params, 'KEdf_radial_n_zones', ''//prev%KEdf_radial_n_zones, this%KEdf_radial_n_zones, help_string="")
    call param_register(params, 'KEdf_radial_gaussian_sigma', ''//prev%KEdf_radial_gaussian_sigma, this%KEdf_radial_gaussian_sigma, help_string="")
    call param_register(params, 'KEdf_radial_mask', ''//trim(prev%KEdf_radial_mask_str), this%KEdf_radial_mask_str, help_string="")
  endif

  if (.not. present(prev)) then
    ! r-dep |F| distribution
    call param_register(params, 'propdf_radial_zone_center', '0.0 0.0 0.0', this%propdf_radial_zone_center, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'propdf_radial_zone_center_at', '-1', this%propdf_radial_zone_center_at, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'propdf_radial_zone_width', '-1.0', this%propdf_radial_zone_width, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'propdf_radial_n_zones', '1', this%propdf_radial_n_zones, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'propdf_radial_bin_width', '0.0', this%propdf_radial_bin_width, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'propdf_radial_n_bins', '0', this%propdf_radial_n_bins, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'propdf_radial_gaussian_sigma', '0.0', this%propdf_radial_gaussian_sigma, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'propdf_radial_mask', '', this%propdf_radial_mask_str, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'propdf_radial_property', '', this%propdf_radial_property, help_string="No help yet.  This source file was $LastChangedBy$")
  else
    ! r-dep |F| distribution
    call param_register(params, 'propdf_radial_zone_center', ''//prev%propdf_radial_zone_center, this%propdf_radial_zone_center, help_string="")
    call param_register(params, 'propdf_radial_zone_center_at', ''//prev%propdf_radial_zone_center_at, this%propdf_radial_zone_center_at, help_string="")
    call param_register(params, 'propdf_radial_zone_width', ''//prev%propdf_radial_zone_width, this%propdf_radial_zone_width, help_string="")
    call param_register(params, 'propdf_radial_n_zones', ''//prev%propdf_radial_n_zones, this%propdf_radial_n_zones, help_string="")
    call param_register(params, 'propdf_radial_bin_width', ''//prev%propdf_radial_bin_width, this%propdf_radial_bin_width, help_string="")
    call param_register(params, 'propdf_radial_n_bins', ''//prev%propdf_radial_n_bins, this%propdf_radial_n_bins, help_string="")
    call param_register(params, 'propdf_radial_gaussian_sigma', ''//prev%propdf_radial_gaussian_sigma, this%propdf_radial_gaussian_sigma, help_string="")
    call param_register(params, 'propdf_radial_mask', ''//trim(prev%propdf_radial_mask_str), this%propdf_radial_mask_str, help_string="")
    call param_register(params, 'propdf_radial_property', ''//trim(prev%propdf_radial_property), this%propdf_radial_property, help_string="")
  endif

  if (.not. present(prev)) then
    ! geometry
    call param_register(params, 'geometry_filename', '', this%geometry_filename, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'geometry_central_atom', '-1', this%geometry_central_atom, help_string="No help yet.  This source file was $LastChangedBy$")
  else
    ! geometry
    call param_register(params, 'geometry_filename', ''//trim(prev%geometry_filename), this%geometry_filename, help_string="")
    call param_register(params, 'geometry_central_atom', ''//prev%geometry_central_atom, this%geometry_central_atom, help_string="")
  endif

  if (.not. present(prev)) then
    ! uniaxial density silica
    call param_register(params, 'density_axial_n_bins', '-1', this%density_axial_n_bins, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'density_axial_axis', '-1', this%density_axial_axis, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'density_axial_silica_atoms', '-1', this%density_axial_silica_atoms, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'density_axial_gaussian', 'F', this%density_axial_gaussian_smoothing, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'density_axial_sigma', '1.0', this%density_axial_gaussian_sigma, help_string="No help yet.  This source file was $LastChangedBy$")
  else
    ! uniaxial density silica
    call param_register(params, 'density_axial_n_bins', ''//prev%density_axial_n_bins, this%density_axial_n_bins, help_string="")
    call param_register(params, 'density_axial_axis', ''//prev%density_axial_axis, this%density_axial_axis, help_string="")
    call param_register(params, 'density_axial_silica_atoms', ''//prev%density_axial_silica_atoms, this%density_axial_silica_atoms, help_string="")
    call param_register(params, 'density_axial_gaussian', ''//prev%density_axial_gaussian_smoothing, this%density_axial_gaussian_smoothing, help_string="")
    call param_register(params, 'density_axial_sigma', ''//prev%density_axial_gaussian_sigma, this%density_axial_gaussian_sigma, help_string="")
  endif

  if (.not. present(prev)) then
    ! num_hbond_silica
    call param_register(params, 'num_hbond_n_bins', '-1', this%num_hbond_n_bins, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'num_hbond_axis', '0', this%num_hbond_axis, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'num_hbond_silica_atoms', '0', this%num_hbond_silica_atoms, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'num_hbond_gaussian', 'F', this%num_hbond_gaussian_smoothing, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'num_hbond_sigma', '1.0', this%num_hbond_gaussian_sigma, help_string="No help yet.  This source file was $LastChangedBy$")
  else
    ! num_hbond_silica
    call param_register(params, 'num_hbond_n_bins', ''//prev%num_hbond_n_bins, this%num_hbond_n_bins, help_string="")
    call param_register(params, 'num_hbond_axis', ''//prev%num_hbond_axis, this%num_hbond_axis, help_string="")
    call param_register(params, 'num_hbond_silica_atoms', ''//prev%num_hbond_silica_atoms, this%num_hbond_silica_atoms, help_string="")
    call param_register(params, 'num_hbond_gaussian', ''//prev%num_hbond_gaussian_smoothing, this%num_hbond_gaussian_smoothing, help_string="")
    call param_register(params, 'num_hbond_sigma', ''//prev%num_hbond_gaussian_sigma, this%num_hbond_gaussian_sigma, help_string="")
  endif

  if (.not. present(prev)) then
    ! water_orientation_silica
    call param_register(params, 'water_orientation_n_pos_bins', '-1', this%water_orientation_n_pos_bins, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'water_orientation_n_angle_bins', '-1', this%water_orientation_n_angle_bins, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'water_orientation_axis', '0', this%water_orientation_axis, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'water_orientation_silica_atoms', '0', this%water_orientation_silica_atoms, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'water_orientation_gaussian', 'F', this%water_orientation_gaussian_smoothing, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'water_orientation_pos_sigma', '1.0', this%water_orientation_pos_gaussian_sigma, help_string="No help yet.  This source file was $LastChangedBy$")
    !call param_register(params, 'water_orientation_angle_sigma', '1.0', this%water_orientation_angle_gaussian_sigma, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'water_orientation_use_dipole', 'F', this%water_orientation_use_dipole, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'water_orientation_use_HOHangle_bisector', 'F', this%water_orientation_use_HOHangle_bisector, help_string="No help yet.  This source file was $LastChangedBy$")
  else
    ! water_orientation_silica
    call param_register(params, 'water_orientation_n_pos_bins', ''//prev%water_orientation_n_pos_bins, this%water_orientation_n_pos_bins, help_string="")
    call param_register(params, 'water_orientation_n_angle_bins', ''//prev%water_orientation_n_angle_bins, this%water_orientation_n_angle_bins, help_string="")
    call param_register(params, 'water_orientation_axis', ''//prev%water_orientation_axis, this%water_orientation_axis, help_string="")
    call param_register(params, 'water_orientation_silica_atoms', ''//prev%water_orientation_silica_atoms, this%water_orientation_silica_atoms, help_string="")
    call param_register(params, 'water_orientation_gaussian', ''//prev%water_orientation_gaussian_smoothing, this%water_orientation_gaussian_smoothing, help_string="")
    call param_register(params, 'water_orientation_pos_sigma', ''//prev%water_orientation_pos_gaussian_sigma, this%water_orientation_pos_gaussian_sigma, help_string="")
    !call param_register(params, 'water_orientation_angle_sigma', ''//prev%water_orientation_angle_gaussian_sigma, this%water_orientation_angle_gaussian_sigma, help_string="")
    call param_register(params, 'water_orientation_use_dipole', ''//prev%water_orientation_use_dipole, this%water_orientation_use_dipole, help_string="")
    call param_register(params, 'water_orientation_use_HOHangle_bisector', ''//prev%water_orientation_use_HOHangle_bisector, this%water_orientation_use_HOHangle_bisector, help_string="")
  endif

  if (present(args_str)) then
    if (.not. param_read_line(params, trim(args_str), ignore_unknown=.false.)) &
      call system_abort("analysis_read failed to parse string '"//trim(args_str)//"'")
  else
    if (.not. param_read_args(params)) &
      call system_abort("analysis_read failed to parse command line arguments")
  endif

  this%density_radial = .false.
  this%density_grid = .false.
  this%KE_density_radial = .false.
  this%rdfd = .false.
  this%adfd = .false.
  this%KEdf_radial = .false.
  this%propdf_radial = .false.
  this%geometry = .false.
  this%density_axial_silica = .false.
  this%num_hbond_silica = .false.
  this%water_orientation_silica = .false.
  select case(trim(this%type))
    case("density_radial")
      this%density_radial = .true.
    case("density_grid")
      this%density_grid = .true.
    case("KE_density_radial")
      this%KE_density_radial = .true.
    case("rdfd")
      this%rdfd = .true.
    case("adfd")
      this%adfd = .true.
    case("KEdf_radial")
      this%KEdf_radial = .true.
    case("propdf_radial")
      this%propdf_radial = .true.
    case("geometry")
      this%geometry = .true.
    case("density_axial_silica")
      this%density_axial_silica = .true.
    case("num_hbond_silica")
      this%num_hbond_silica = .true.
    case("water_orientation_silica")
      this%water_orientation_silica = .true.
    case default
      call system_abort("Unknown analysis type '"//trim(this%type)//"'")
  end select

  if (count ( (/ this%density_radial, this%density_grid, this%KE_density_radial, this%rdfd, this%adfd, &
		 this%KEdf_radial, this%propdf_radial, this%geometry, this%density_axial_silica, this%num_hbond_silica, &
		 this%water_orientation_silica /) ) /= 1) &
    call system_abort("Specified "//(/ this%density_radial, this%density_grid, this%KE_density_radial, &
      this%rdfd, this%adfd, this%KEdf_radial, this%propdf_radial, this%geometry, this%density_axial_silica, this%num_hbond_silica, &
      this%water_orientation_silica /)// &
      " types of analysis.  Possibilities: density_radial, density_grid, KE_density_radial, rdfd, adfd, KEdf_radial, propdf_radial, geometry, " // &
      " density_axial_silica, num_hbond_silica, water_orientation_silica.")

end subroutine analysis_read

subroutine check_analyses(a)
  type(analysis), intent(inout) :: a(:)

  integer :: i_a

  do i_a=1, size(a)
    if (a(i_a)%density_radial) then !density_radial
      if (a(i_a)%density_radial_bin_width <= 0.0_dp) call system_abort("analysis " // i_a // " has radial_bin_width="//a(i_a)%density_radial_bin_width//" <= 0.0")
      if (a(i_a)%density_radial_n_bins <= 0) call system_abort("analysis " // i_a // " has radial_n_bins="//a(i_a)%density_radial_n_bins//" <= 0")
    else if (a(i_a)%density_grid) then !density_grid
      if (any(a(i_a)%density_grid_bin_width <= 0.0_dp)) call system_abort("analysis " // i_a // " has grid_bin_width="//a(i_a)%density_grid_bin_width//" <= 0.0")
      if (any(a(i_a)%density_grid_n_bins <= 0)) call system_abort("analysis " // i_a // " has grid_n_bins="//a(i_a)%density_grid_n_bins//" <= 0")
    else if (a(i_a)%KE_density_radial) then !KE_density_radial
      if (a(i_a)%KE_density_radial_bin_width <= 0.0_dp) call system_abort("analysis " // i_a // " has radial_bin_width="//a(i_a)%KE_density_radial_bin_width//" <= 0.0")
      if (a(i_a)%KE_density_radial_n_bins <= 0) call system_abort("analysis " // i_a // " has radial_n_bins="//a(i_a)%KE_density_radial_n_bins//" <= 0")
    else if (a(i_a)%rdfd) then !rdfd
      if (a(i_a)%rdfd_bin_width <= 0.0_dp) call system_abort("analysis " // i_a // " has rdfd_bin_width="//a(i_a)%rdfd_bin_width//" <= 0.0")
      if (a(i_a)%rdfd_n_bins <= 0) call system_abort("analysis " // i_a // " has rdfd_n_bins="//a(i_a)%rdfd_n_bins//" <= 0")
      if (a(i_a)%rdfd_n_zones <= 0) call system_abort("analysis " // i_a // " has rdfd_n_zones="//a(i_a)%rdfd_n_zones//" <= 0")
    else if (a(i_a)%adfd) then !adfd
      if (a(i_a)%adfd_dist_bin_width <= 0.0_dp) call system_abort("analysis " // i_a // " has adfd_dist_bin_width="//a(i_a)%adfd_dist_bin_width//" <= 0.0")
      if (a(i_a)%adfd_n_dist_bins <= 0) call system_abort("analysis " // i_a // " has adfd_n_dist_bins="//a(i_a)%adfd_n_dist_bins//" <= 0")
      if (a(i_a)%adfd_n_angle_bins <= 0) call system_abort("analysis " // i_a // " has adfd_n_angle_bins="//a(i_a)%adfd_n_angle_bins//" <= 0")
      if (a(i_a)%adfd_n_zones <= 0) call system_abort("analysis " // i_a // " has adfd_n_zones="//a(i_a)%adfd_n_zones//" <= 0")
    else if (a(i_a)%KEdf_radial) then !KEdf_radial
      if (a(i_a)%KEdf_radial_bin_width <= 0.0_dp) call system_abort("analysis " // i_a // " has KEdf_radial_bin_width="//a(i_a)%KEdf_radial_bin_width//" <= 0.0")
      if (a(i_a)%KEdf_radial_n_bins <= 0) call system_abort("analysis " // i_a // " has KEdf_radial_n_bins="//a(i_a)%KEdf_radial_n_bins//" <= 0")
      if (a(i_a)%KEdf_radial_n_zones <= 0) call system_abort("analysis " // i_a // " has KEdf_radial_n_zones="//a(i_a)%KEdf_radial_n_zones//" <= 0")
    else if (a(i_a)%propdf_radial) then !propdf_radial
      if (a(i_a)%propdf_radial_bin_width <= 0.0_dp) call system_abort("analysis " // i_a // " has propdf_radial_bin_width="//a(i_a)%propdf_radial_bin_width//" <= 0.0")
      if (a(i_a)%propdf_radial_n_bins <= 0) call system_abort("analysis " // i_a // " has propdf_radial_n_bins="//a(i_a)%propdf_radial_n_bins//" <= 0")
      if (a(i_a)%propdf_radial_gaussian_sigma <= 0) call system_abort("analysis " // i_a // " has propdf_radial_gaussian_sigma="//a(i_a)%propdf_radial_gaussian_sigma//" <= 0")
      if (a(i_a)%propdf_radial_n_zones <= 0) call system_abort("analysis " // i_a // " has propdf_radial_n_zones="//a(i_a)%propdf_radial_n_zones//" <= 0")
    else if (a(i_a)%geometry) then !geometry
      if (trim(a(i_a)%geometry_filename)=="") call system_abort("analysis "//i_a//" has empty geometry_filename")
      !read geometry parameters to calculate from the file into a table
      call read_geometry_params(a(i_a),trim(a(i_a)%geometry_filename))
      if (a(i_a)%geometry_params%N==0) call system_abort("analysis "//i_a//" has no geometry parameters to calculate")
    else if (a(i_a)%density_axial_silica) then !density_axial_silica
      if (a(i_a)%density_axial_n_bins <= 0) call system_abort("analysis " // i_a // " has density_axial_n_bins="//a(i_a)%density_axial_n_bins//" <= 0")
      if (.not. any(a(i_a)%density_axial_axis == (/1,2,3/))) call system_abort("analysis " // i_a // " has density_axial_axis="//a(i_a)%density_axial_axis//" /= 1, 2 or 3")
    else if (a(i_a)%num_hbond_silica) then !num_hbond_silica
      if (a(i_a)%num_hbond_n_bins <= 0) call system_abort("analysis " // i_a // " has num_hbond_n_bins="//a(i_a)%num_hbond_n_bins//" <= 0")
      if (.not. any(a(i_a)%num_hbond_axis == (/1,2,3/))) call system_abort("analysis " // i_a // " has num_hbond_axis="//a(i_a)%num_hbond_axis//" /= 1, 2 or 3")
    else if (a(i_a)%water_orientation_silica) then !water_orientation_silica
      if (a(i_a)%water_orientation_n_pos_bins <= 0) call system_abort("analysis " // i_a // " has water_orientation_n_pos_bins="//a(i_a)%water_orientation_n_pos_bins//" <= 0")
      if (a(i_a)%water_orientation_n_angle_bins <= 0) call system_abort("analysis " // i_a // " has water_orientation_n_angle_bins="//a(i_a)%water_orientation_n_angle_bins//" <= 0")
      if (.not. any(a(i_a)%water_orientation_axis == (/1,2,3/))) call system_abort("analysis " // i_a // " has water_orientation_axis="//a(i_a)%water_orientation_axis//" /= 1, 2 or 3")
      if (.not. count((/a(i_a)%water_orientation_use_HOHangle_bisector,a(i_a)%water_orientation_use_dipole/)) == 1) call system_abort("Exactly one of water_orientation_use_HOHangle_bisector and water_orientation_use_dipole must be one.")
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
  real(dp) :: use_radial_center(3)

  integer :: i_a

  call map_into_cell(at)
!  at%t ravel = 0

  do i_a=1, size(a)

    if (do_this_analysis(a(i_a), time, frame)) then

      a(i_a)%n_configs = a(i_a)%n_configs + 1

      if (a(i_a)%density_radial) then !density_radial
        call reallocate_data(a(i_a)%density_radial_histograms, a(i_a)%n_configs, a(i_a)%density_radial_n_bins)
	if (a(i_a)%density_radial_center_at > 0) then
	  use_radial_center = at%pos(:,a(i_a)%density_radial_center_at)
	else
	  use_radial_center = a(i_a)%density_radial_center
	endif
        if (a(i_a)%n_configs == 1) then
          allocate(a(i_a)%density_radial_pos(a(i_a)%density_radial_n_bins))
          call density_sample_radial_mesh_Gaussians(a(i_a)%density_radial_histograms(:,a(i_a)%n_configs), at, center_pos=use_radial_center, &
            rad_bin_width=a(i_a)%density_radial_bin_width, n_rad_bins=a(i_a)%density_radial_n_bins, gaussian_sigma=a(i_a)%density_radial_gaussian_sigma, &
            mask_str=a(i_a)%mask_str, radial_pos=a(i_a)%density_radial_pos)
        else
          call density_sample_radial_mesh_Gaussians(a(i_a)%density_radial_histograms(:,a(i_a)%n_configs), at, center_pos=use_radial_center, &
            rad_bin_width=a(i_a)%density_radial_bin_width, n_rad_bins=a(i_a)%density_radial_n_bins, gaussian_sigma= a(i_a)%density_radial_gaussian_sigma, &
            mask_str=a(i_a)%mask_str)
        endif
      else if (a(i_a)%density_grid) then !density_grid
        call reallocate_data(a(i_a)%density_grid_histograms, a(i_a)%n_configs, a(i_a)%density_grid_n_bins)
        if (a(i_a)%density_grid_gaussian_smoothing) then
          if (a(i_a)%n_configs == 1) then
            allocate(a(i_a)%density_grid_pos(3,a(i_a)%density_grid_n_bins(1),a(i_a)%density_grid_n_bins(2),a(i_a)%density_grid_n_bins(3)))
            call density_sample_rectilinear_mesh_Gaussians(a(i_a)%density_grid_histograms(:,:,:,a(i_a)%n_configs), at, a(i_a)%density_grid_min_p, &
              a(i_a)%density_grid_bin_width, a(i_a)%density_grid_n_bins, a(i_a)%density_grid_gaussian_sigma, a(i_a)%mask_str, a(i_a)%density_grid_pos)
          else
            call density_sample_rectilinear_mesh_Gaussians(a(i_a)%density_grid_histograms(:,:,:,a(i_a)%n_configs), at, a(i_a)%density_grid_min_p, &
              a(i_a)%density_grid_bin_width, a(i_a)%density_grid_n_bins, a(i_a)%density_grid_gaussian_sigma, a(i_a)%mask_str)
          endif
        else
          if (a(i_a)%n_configs == 1) then
            allocate(a(i_a)%density_grid_pos(3,a(i_a)%density_grid_n_bins(1),a(i_a)%density_grid_n_bins(2),a(i_a)%density_grid_n_bins(3)))
            call density_bin_rectilinear_mesh(a(i_a)%density_grid_histograms(:,:,:,a(i_a)%n_configs), at, a(i_a)%density_grid_min_p, a(i_a)%density_grid_bin_width, &
              a(i_a)%density_grid_n_bins, a(i_a)%mask_str, a(i_a)%density_grid_pos)
          else
            call density_bin_rectilinear_mesh(a(i_a)%density_grid_histograms(:,:,:,a(i_a)%n_configs), at, a(i_a)%density_grid_min_p, a(i_a)%density_grid_bin_width, &
              a(i_a)%density_grid_n_bins, a(i_a)%mask_str)
          endif
        endif
      else if (a(i_a)%KE_density_radial) then ! KE_density_radial
        call reallocate_data(a(i_a)%KE_density_radial_histograms, a(i_a)%n_configs, a(i_a)%KE_density_radial_n_bins)
	if (a(i_a)%KE_density_radial_center_at > 0) then
	  use_radial_center = at%pos(:,a(i_a)%KE_density_radial_center_at)
	else
	  use_radial_center = a(i_a)%KE_density_radial_center
	endif
        if (a(i_a)%n_configs == 1) then
          allocate(a(i_a)%KE_density_radial_pos(a(i_a)%KE_density_radial_n_bins))
          call density_sample_radial_mesh_Gaussians(a(i_a)%KE_density_radial_histograms(:,a(i_a)%n_configs), at, center_pos=use_radial_center, &
            rad_bin_width=a(i_a)%KE_density_radial_bin_width, n_rad_bins=a(i_a)%KE_density_radial_n_bins, gaussian_sigma=a(i_a)%KE_density_radial_gaussian_sigma, &
            mask_str=a(i_a)%mask_str, radial_pos=a(i_a)%KE_density_radial_pos, quantity="KE")
        else
          call density_sample_radial_mesh_Gaussians(a(i_a)%KE_density_radial_histograms(:,a(i_a)%n_configs), at, center_pos=use_radial_center, &
            rad_bin_width=a(i_a)%KE_density_radial_bin_width, n_rad_bins=a(i_a)%KE_density_radial_n_bins, gaussian_sigma= a(i_a)%density_radial_gaussian_sigma, &
            mask_str=a(i_a)%mask_str, quantity="KE")
        endif
      else if (a(i_a)%rdfd) then !rdfd
        call reallocate_data(a(i_a)%rdfds, a(i_a)%n_configs, (/ a(i_a)%rdfd_n_bins, a(i_a)%rdfd_n_zones /) )
        if (a(i_a)%n_configs == 1) then
          allocate(a(i_a)%rdfd_bin_pos(a(i_a)%rdfd_n_bins))
          allocate(a(i_a)%rdfd_zone_pos(a(i_a)%rdfd_n_zones))
          call rdfd_calc(a(i_a)%rdfds(:,:,a(i_a)%n_configs), at, a(i_a)%rdfd_zone_center, a(i_a)%rdfd_zone_atom_center, a(i_a)%rdfd_bin_width, a(i_a)%rdfd_n_bins, &
            a(i_a)%rdfd_zone_width, a(i_a)%rdfd_n_zones, a(i_a)%rdfd_gaussian_smoothing, a(i_a)%rdfd_gaussian_sigma, &
            a(i_a)%rdfd_center_mask_str, a(i_a)%rdfd_neighbour_mask_str, &
            a(i_a)%rdfd_bin_pos, a(i_a)%rdfd_zone_pos)
        else
          call rdfd_calc(a(i_a)%rdfds(:,:,a(i_a)%n_configs), at, a(i_a)%rdfd_zone_center, a(i_a)%rdfd_zone_atom_center, a(i_a)%rdfd_bin_width, a(i_a)%rdfd_n_bins, &
            a(i_a)%rdfd_zone_width, a(i_a)%rdfd_n_zones, a(i_a)%rdfd_gaussian_smoothing, a(i_a)%rdfd_gaussian_sigma, &
            a(i_a)%rdfd_center_mask_str, a(i_a)%rdfd_neighbour_mask_str)
        endif
      else if (a(i_a)%adfd) then !adfd
        call reallocate_data(a(i_a)%adfds, a(i_a)%n_configs, (/ a(i_a)%adfd_n_angle_bins, a(i_a)%adfd_n_dist_bins, a(i_a)%adfd_n_zones /) )
        if (a(i_a)%n_configs == 1) then
          allocate(a(i_a)%adfd_dist_bin_pos(a(i_a)%adfd_n_dist_bins))
          allocate(a(i_a)%adfd_angle_bin_pos(a(i_a)%adfd_n_angle_bins))
          allocate(a(i_a)%adfd_zone_pos(a(i_a)%adfd_n_zones))
          call adfd_calc(a(i_a)%adfds(:,:,:,a(i_a)%n_configs), at, a(i_a)%adfd_zone_center, &
	    a(i_a)%adfd_n_angle_bins, a(i_a)%adfd_dist_bin_width, a(i_a)%adfd_n_dist_bins, &
            a(i_a)%adfd_zone_width, a(i_a)%adfd_n_zones, &
            a(i_a)%adfd_center_mask_str, a(i_a)%adfd_neighbour_1_mask_str, a(i_a)%adfd_neighbour_1_max_dist, a(i_a)%adfd_neighbour_2_mask_str, a(i_a)%adfd_dist_bin_rc2, &
            a(i_a)%adfd_angle_bin_pos, a(i_a)%adfd_dist_bin_pos, a(i_a)%adfd_zone_pos)
        else
          call adfd_calc(a(i_a)%adfds(:,:,:,a(i_a)%n_configs), at, a(i_a)%adfd_zone_center, &
	    a(i_a)%adfd_n_angle_bins, a(i_a)%adfd_dist_bin_width, a(i_a)%adfd_n_dist_bins, &
            a(i_a)%adfd_zone_width, a(i_a)%adfd_n_zones, &
            a(i_a)%adfd_center_mask_str, a(i_a)%adfd_neighbour_1_mask_str, a(i_a)%adfd_neighbour_1_max_dist, a(i_a)%adfd_neighbour_2_mask_str, a(i_a)%adfd_dist_bin_rc2)
        endif
      else if (a(i_a)%KEdf_radial) then
	call reallocate_data(a(i_a)%KEdf_radial_histograms, a(i_a)%n_configs, (/ a(i_a)%KEdf_radial_n_bins, a(i_a)%KEdf_radial_n_zones /) )
	if (a(i_a)%KEdf_radial_zone_center_at > 0) then
	  use_radial_center = at%pos(:,a(i_a)%KEdf_radial_zone_center_at)
	else
	  use_radial_center = a(i_a)%KEdf_radial_zone_center
	endif
	if (a(i_a)%n_configs == 1) then
	  allocate(a(i_a)%KEdf_radial_bin_pos(a(i_a)%KEdf_radial_n_bins))
	  allocate(a(i_a)%KEdf_radial_zone_pos(a(i_a)%KEdf_radial_n_zones))
	  call propdf_radial_calc(a(i_a)%KEdf_radial_histograms(:,:,a(i_a)%n_configs), at, a(i_a)%KEdf_radial_bin_width, a(i_a)%KEdf_radial_n_bins, &
	    use_radial_center, a(i_a)%KEdf_radial_zone_width, A(i_a)%KEdf_radial_n_zones, &
	    a(i_a)%KEdf_radial_gaussian_sigma, a(i_a)%KEdf_radial_mask_str, 'KE', a(i_a)%KEdf_radial_bin_pos, a(i_a)%KEdf_radial_zone_pos)
	else
	  call propdf_radial_calc(a(i_a)%KEdf_radial_histograms(:,:,a(i_a)%n_configs), at, a(i_a)%KEdf_radial_bin_width, a(i_a)%KEdf_radial_n_bins, &
	    use_radial_center, a(i_a)%KEdf_radial_zone_width, A(i_a)%KEdf_radial_n_zones, &
	    a(i_a)%KEdf_radial_gaussian_sigma, a(i_a)%KEdf_radial_mask_str, 'KE')
	endif
      else if (a(i_a)%propdf_radial) then
	call reallocate_data(a(i_a)%propdf_radial_histograms, a(i_a)%n_configs, (/ a(i_a)%propdf_radial_n_bins, a(i_a)%propdf_radial_n_zones /) )
	if (a(i_a)%propdf_radial_zone_center_at > 0) then
	  use_radial_center = at%pos(:,a(i_a)%propdf_radial_zone_center_at)
	else
	  use_radial_center = a(i_a)%propdf_radial_zone_center
	endif
	if (a(i_a)%n_configs == 1) then
	  allocate(a(i_a)%propdf_radial_bin_pos(a(i_a)%propdf_radial_n_bins))
	  allocate(a(i_a)%propdf_radial_zone_pos(a(i_a)%propdf_radial_n_zones))
	  call propdf_radial_calc(a(i_a)%propdf_radial_histograms(:,:,a(i_a)%n_configs), at, a(i_a)%propdf_radial_bin_width, a(i_a)%propdf_radial_n_bins, &
	    use_radial_center, a(i_a)%propdf_radial_zone_width, A(i_a)%propdf_radial_n_zones, &
	    a(i_a)%propdf_radial_gaussian_sigma, a(i_a)%propdf_radial_mask_str, A(i_a)%propdf_radial_property, &
	    a(i_a)%propdf_radial_bin_pos, a(i_a)%propdf_radial_zone_pos)
	else
	  call propdf_radial_calc(a(i_a)%propdf_radial_histograms(:,:,a(i_a)%n_configs), at, a(i_a)%propdf_radial_bin_width, a(i_a)%propdf_radial_n_bins, &
	    use_radial_center, a(i_a)%propdf_radial_zone_width, A(i_a)%propdf_radial_n_zones, &
	    a(i_a)%propdf_radial_gaussian_sigma, a(i_a)%propdf_radial_mask_str, a(i_a)%propdf_radial_property)
	endif
      else if (a(i_a)%geometry) then !geometry
        call reallocate_data(a(i_a)%geometry_histograms, a(i_a)%n_configs, a(i_a)%geometry_params%N)
        if (a(i_a)%n_configs == 1) then
          allocate(a(i_a)%geometry_pos(a(i_a)%geometry_params%N))
          allocate(a(i_a)%geometry_label(a(i_a)%geometry_params%N))
          call geometry_calc(a(i_a)%geometry_histograms(:,a(i_a)%n_configs), at, a(i_a)%geometry_params, a(i_a)%geometry_central_atom, &
               a(i_a)%geometry_pos(1:a(i_a)%geometry_params%N), a(i_a)%geometry_label(1:a(i_a)%geometry_params%N))
        else
          call geometry_calc(a(i_a)%geometry_histograms(:,a(i_a)%n_configs), at, a(i_a)%geometry_params, a(i_a)%geometry_central_atom)
        endif
      else if (a(i_a)%density_axial_silica) then !density_axial_silica
        call reallocate_data(a(i_a)%density_axial_histograms, a(i_a)%n_configs, a(i_a)%density_axial_n_bins)
        if (a(i_a)%n_configs == 1) then
           allocate(a(i_a)%density_axial_pos(a(i_a)%density_axial_n_bins))
           call density_axial_calc(a(i_a)%density_axial_histograms(:,a(i_a)%n_configs), at, &
             axis=a(i_a)%density_axial_axis, silica_center_i=a(i_a)%density_axial_silica_atoms, &
             n_bins=a(i_a)%density_axial_n_bins, &
             gaussian_smoothing=a(i_a)%density_axial_gaussian_smoothing, &
             gaussian_sigma=a(i_a)%density_axial_gaussian_sigma, &
             mask_str=a(i_a)%mask_str, axial_pos=a(i_a)%density_axial_pos)
        else
           call density_axial_calc(a(i_a)%density_axial_histograms(:,a(i_a)%n_configs), at, &
             axis=a(i_a)%density_axial_axis, silica_center_i=a(i_a)%density_axial_silica_atoms, &
             n_bins=a(i_a)%density_axial_n_bins, &
             gaussian_smoothing=a(i_a)%density_axial_gaussian_smoothing, &
             gaussian_sigma=a(i_a)%density_axial_gaussian_sigma, &
             mask_str=a(i_a)%mask_str)
        endif
      else if (a(i_a)%num_hbond_silica) then !num_hbond_silica
        a(i_a)%num_hbond_n_type = 4
        call reallocate_data(a(i_a)%num_hbond_histograms, a(i_a)%n_configs, (/a(i_a)%num_hbond_n_bins,a(i_a)%num_hbond_n_type/)) !4: ss, sw, ws, ww
        if (a(i_a)%n_configs == 1) then
           allocate(a(i_a)%num_hbond_bin_pos(a(i_a)%num_hbond_n_bins))
           allocate(a(i_a)%num_hbond_type_code(4))
           allocate(a(i_a)%num_hbond_type_label(4))
           call num_hbond_calc(a(i_a)%num_hbond_histograms(:,:,a(i_a)%n_configs), at, &
             axis=a(i_a)%num_hbond_axis, silica_center_i=a(i_a)%num_hbond_silica_atoms, &
             n_bins=a(i_a)%num_hbond_n_bins, &
             gaussian_smoothing=a(i_a)%num_hbond_gaussian_smoothing, &
             gaussian_sigma=a(i_a)%num_hbond_gaussian_sigma, &
             mask_str=a(i_a)%mask_str, num_hbond_pos=a(i_a)%num_hbond_bin_pos, &
             num_hbond_type_code=a(i_a)%num_hbond_type_code, &
             num_hbond_type_label=a(i_a)%num_hbond_type_label)
        else
           call num_hbond_calc(a(i_a)%num_hbond_histograms(:,:,a(i_a)%n_configs), at, &
             axis=a(i_a)%num_hbond_axis, silica_center_i=a(i_a)%num_hbond_silica_atoms, &
             n_bins=a(i_a)%num_hbond_n_bins, &
             gaussian_smoothing=a(i_a)%num_hbond_gaussian_smoothing, &
             gaussian_sigma=a(i_a)%num_hbond_gaussian_sigma, &
             mask_str=a(i_a)%mask_str)
        endif
      else if (a(i_a)%water_orientation_silica) then !water_orientation_silica
        call reallocate_data(a(i_a)%water_orientation_histograms, a(i_a)%n_configs, (/ a(i_a)%water_orientation_n_angle_bins, a(i_a)%water_orientation_n_pos_bins /) )
        if (a(i_a)%n_configs == 1) then
          allocate(a(i_a)%water_orientation_angle_bin(a(i_a)%water_orientation_n_angle_bins))
          allocate(a(i_a)%water_orientation_angle_bin_w(a(i_a)%water_orientation_n_angle_bins))
          allocate(a(i_a)%water_orientation_pos_bin(a(i_a)%water_orientation_n_pos_bins))
          call water_orientation_calc(a(i_a)%water_orientation_histograms(:,:,a(i_a)%n_configs), at, &
             axis=a(i_a)%water_orientation_axis, silica_center_i=a(i_a)%water_orientation_silica_atoms, &
            n_pos_bins=a(i_a)%water_orientation_n_pos_bins, n_angle_bins=a(i_a)%water_orientation_n_angle_bins, &
            gaussian_smoothing=a(i_a)%water_orientation_gaussian_smoothing, &
            pos_gaussian_sigma=a(i_a)%water_orientation_pos_gaussian_sigma, &
            !angle_gaussian_sigma=a(i_a)%water_orientation_angle_gaussian_sigma, &
            pos_bin=a(i_a)%water_orientation_pos_bin, angle_bin=a(i_a)%water_orientation_angle_bin, &
            angle_bin_w=a(i_a)%water_orientation_angle_bin_w, &
            use_dipole_rather_than_angle_bisector=a(i_a)%water_orientation_use_dipole)
        else
          call water_orientation_calc(a(i_a)%water_orientation_histograms(:,:,a(i_a)%n_configs), at, &
             axis=a(i_a)%water_orientation_axis, silica_center_i=a(i_a)%water_orientation_silica_atoms, &
            n_pos_bins=a(i_a)%water_orientation_n_pos_bins, n_angle_bins=a(i_a)%water_orientation_n_angle_bins, &
            gaussian_smoothing=a(i_a)%water_orientation_gaussian_smoothing, &
            pos_gaussian_sigma=a(i_a)%water_orientation_pos_gaussian_sigma, &
            !angle_gaussian_sigma=a(i_a)%water_orientation_angle_gaussian_sigma, &
            angle_bin_w=a(i_a)%water_orientation_angle_bin_w, &
            use_dipole_rather_than_angle_bisector=a(i_a)%water_orientation_use_dipole)
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
      if (a(i_a)%density_radial) then !density_radial
        call print("# radial density histogram", file=outfile)
        call print("n_bins="//a(i_a)%density_radial_n_bins//" n_data="//a(i_a)%n_configs, file=outfile)
        do i=1, a(i_a)%density_radial_n_bins
          call print(a(i_a)%density_radial_pos(i), file=outfile)
        end do
        do i=1, a(i_a)%n_configs
          call print(a(i_a)%density_radial_histograms(:,i), file=outfile)
        end do
      else if (a(i_a)%density_grid) then !density_grid
        call print("# grid density histogram", file=outfile)
        call print("n_bins="//a(i_a)%density_grid_n_bins(1)*a(i_a)%density_grid_n_bins(2)*a(i_a)%density_grid_n_bins(3)//" n_data="//a(i_a)%n_configs, file=outfile)
        do i1=1, a(i_a)%density_grid_n_bins(1)
        do i2=1, a(i_a)%density_grid_n_bins(2)
        do i3=1, a(i_a)%density_grid_n_bins(3)
          call print(""//a(i_a)%density_grid_pos(:,i1,i2,i3), file=outfile)
        end do
        end do
        end do
        do i=1, a(i_a)%n_configs
          call print(""//reshape(a(i_a)%density_grid_histograms(:,:,:,i), (/ a(i_a)%density_grid_n_bins(1)*a(i_a)%density_grid_n_bins(2)*a(i_a)%density_grid_n_bins(3) /) ), file=outfile)
        end do
      else if (a(i_a)%KE_density_radial) then !density_radial
        call print("# radial density histogram", file=outfile)
        call print("n_bins="//a(i_a)%KE_density_radial_n_bins//" n_data="//a(i_a)%n_configs, file=outfile)
        do i=1, a(i_a)%KE_density_radial_n_bins
          call print(a(i_a)%KE_density_radial_pos(i), file=outfile)
        end do
        do i=1, a(i_a)%n_configs
          call print(a(i_a)%KE_density_radial_histograms(:,i), file=outfile)
        end do
      else if (a(i_a)%rdfd) then !rdfd
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
	allocate(integrated_rdfds(a(i_a)%rdfd_n_bins,a(i_a)%rdfd_n_zones))
	integrated_rdfds = 0.0_dp
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

      else if (a(i_a)%adfd) then !adfd
        call print("# adfd", file=outfile)
        call print("n_bins="//a(i_a)%adfd_n_zones*a(i_a)%adfd_n_dist_bins*a(i_a)%adfd_n_angle_bins//" n_data="//a(i_a)%n_configs, file=outfile)
        do i1=1, a(i_a)%adfd_n_zones
        do i2=1, a(i_a)%adfd_n_dist_bins
        do i3=1, a(i_a)%adfd_n_angle_bins
          if (a(i_a)%adfd_zone_width > 0.0_dp) then
            call print(""//a(i_a)%adfd_zone_pos(i1)//" "//a(i_a)%adfd_dist_bin_pos(i2)//" "//a(i_a)%adfd_angle_bin_pos(i3), file=outfile)
          else
            call print(""//a(i_a)%adfd_dist_bin_pos(i2)//" "//a(i_a)%adfd_angle_bin_pos(i3), file=outfile)
          endif
        end do
        end do
        end do
        do i=1, a(i_a)%n_configs
          call print(""//reshape(a(i_a)%adfds(:,:,:,i), (/ a(i_a)%adfd_n_zones*a(i_a)%adfd_n_dist_bins*a(i_a)%adfd_n_angle_bins /) ), file=outfile)
        end do

      else if (a(i_a)%KEdf_radial) then !r-dep KE density
	call print("# r-dependent KE density", file=outfile)
	call print("n_bins="//a(i_a)%KEdf_radial_n_zones*a(i_a)%KEdf_radial_n_bins//" n_data="//a(i_a)%n_configs, file=outfile)

        do i1=1, a(i_a)%KEdf_radial_n_zones
        do i2=1, a(i_a)%KEdf_radial_n_bins
          if (a(i_a)%KEdf_radial_zone_width > 0.0_dp) then
            call print(""//a(i_a)%KEdf_radial_zone_pos(i1)//" "//a(i_a)%KEdf_radial_bin_pos(i2), file=outfile)
          else
            call print(""//a(i_a)%KEdf_radial_bin_pos(i2), file=outfile)
          endif
        end do
        end do
        do i=1, a(i_a)%n_configs
          call print(""//reshape(a(i_a)%KEdf_radial_histograms(:,:,i), (/ a(i_a)%KEdf_radial_n_zones*a(i_a)%KEdf_radial_n_bins /) ), file=outfile)
        end do

      else if (a(i_a)%propdf_radial) then !r-dep |F| density
	call print("# r-dependent prop density", file=outfile)
	call print("n_bins="//a(i_a)%propdf_radial_n_zones*a(i_a)%propdf_radial_n_bins//" n_data="//a(i_a)%n_configs, file=outfile)

        do i1=1, a(i_a)%propdf_radial_n_zones
        do i2=1, a(i_a)%propdf_radial_n_bins
          if (a(i_a)%propdf_radial_zone_width > 0.0_dp) then
            call print(""//a(i_a)%propdf_radial_zone_pos(i1)//" "//a(i_a)%propdf_radial_bin_pos(i2), file=outfile)
          else
            call print(""//a(i_a)%propdf_radial_bin_pos(i2), file=outfile)
          endif
        end do
        end do
        do i=1, a(i_a)%n_configs
          call print(""//reshape(a(i_a)%propdf_radial_histograms(:,:,i), (/ a(i_a)%propdf_radial_n_zones*a(i_a)%propdf_radial_n_bins /) ), file=outfile)
        end do

      else if (a(i_a)%geometry) then !geometry
        call print("# geometry histogram", file=outfile)
        call print("n_bins="//a(i_a)%geometry_params%N//" n_data="//a(i_a)%n_configs, file=outfile)
        do i=1, a(i_a)%geometry_params%N
!          call print(a(i_a)%geometry_pos(i), file=outfile)
          call print(trim(a(i_a)%geometry_label(i)), file=outfile)
        end do
        do i=1, a(i_a)%n_configs
          call print(a(i_a)%geometry_histograms(:,i), file=outfile)
        end do

      else if (a(i_a)%density_axial_silica) then !density_axial_silica
        call print("# uniaxial density histogram in direction "//a(i_a)%density_axial_axis, file=outfile)
        call print("n_bins="//a(i_a)%density_axial_n_bins//" n_data="//a(i_a)%n_configs, file=outfile)
        do i=1, a(i_a)%density_axial_n_bins
          call print(a(i_a)%density_axial_pos(i), file=outfile)
        end do
        do i=1, a(i_a)%n_configs
          call print(a(i_a)%density_axial_histograms(:,i), file=outfile)
        end do

      else if (a(i_a)%num_hbond_silica) then !num_hbond_silica
        !header
        call print("# num_hbond_silica in direction "//a(i_a)%num_hbond_axis, file=outfile)
        call print("n_bins="//a(i_a)%num_hbond_n_type*a(i_a)%num_hbond_n_bins//" n_data="//a(i_a)%n_configs, file=outfile)
        do i1=1, a(i_a)%num_hbond_n_type
          do i2=1, a(i_a)%num_hbond_n_bins
            !call print(""//a(i_a)%num_hbond_type_code(i1)//" "//a(i_a)%num_hbond_bin_pos(i2), file=outfile)
            call print(""//trim(a(i_a)%num_hbond_type_label(i1))//" "//a(i_a)%num_hbond_bin_pos(i2), file=outfile)
          end do
        end do
        !histograms
        do i=1, a(i_a)%n_configs
          call print(""//reshape(a(i_a)%num_hbond_histograms(:,:,i), (/ a(i_a)%num_hbond_n_type*a(i_a)%num_hbond_n_bins /) ), file=outfile)
        end do

        !integrated histograms header
        call print("", file=outfile)
        call print("", file=outfile)
        call print("# integrated_num_hbond_silica", file=outfile)
        call print("n_bins="//a(i_a)%num_hbond_n_type*a(i_a)%num_hbond_n_bins//" n_data="//a(i_a)%n_configs, file=outfile)
        do i1=1, a(i_a)%num_hbond_n_type
          do i2=1, a(i_a)%num_hbond_n_bins
            !call print(""//a(i_a)%num_hbond_type_code(i1)//" "//a(i_a)%num_hbond_bin_pos(i2), file=outfile)
            call print(""//trim(a(i_a)%num_hbond_type_label(i1))//" "//a(i_a)%num_hbond_bin_pos(i2), file=outfile)
          end do
        end do
        !integrated histograms
	allocate(a(i_a)%integrated_num_hbond_histograms(a(i_a)%num_hbond_n_bins,a(i_a)%num_hbond_n_type))
	a(i_a)%integrated_num_hbond_histograms = 0.0_dp
        do i=1, a(i_a)%n_configs
	  a(i_a)%integrated_num_hbond_histograms = 0.0_dp
	  do i1=1, a(i_a)%num_hbond_n_type
	    do i2=2, a(i_a)%num_hbond_n_bins
	      a(i_a)%integrated_num_hbond_histograms(i2,i1) = a(i_a)%integrated_num_hbond_histograms(i2-1,i1) + &
		(a(i_a)%num_hbond_bin_pos(i2)-a(i_a)%num_hbond_bin_pos(i2-1))* &
		4.0_dp*PI*((a(i_a)%num_hbond_bin_pos(i2)**2)*a(i_a)%num_hbond_histograms(i2,i1,i)+(a(i_a)%num_hbond_bin_pos(i2-1)**2)*a(i_a)%num_hbond_histograms(i2-1,i1,i))/2.0_dp
	    end do
	  end do
          call print(""//reshape(a(i_a)%integrated_num_hbond_histograms(:,:), (/ a(i_a)%num_hbond_n_type*a(i_a)%num_hbond_n_bins /) ), file=outfile)
        end do
	deallocate(a(i_a)%integrated_num_hbond_histograms)

      else if (a(i_a)%water_orientation_silica) then !water_orientation_silica
        !header
        call print("# water_orientation_silica in direction "//a(i_a)%water_orientation_axis, file=outfile)
        call print("n_bins="//a(i_a)%water_orientation_n_pos_bins*a(i_a)%water_orientation_n_angle_bins//" n_data="//a(i_a)%n_configs, file=outfile)
        do i1=1, a(i_a)%water_orientation_n_pos_bins
          do i2=1, a(i_a)%water_orientation_n_angle_bins
            call print(""//a(i_a)%water_orientation_pos_bin(i1)//" "//a(i_a)%water_orientation_angle_bin(i2), file=outfile)
          end do
        end do
        !histograms
        do i=1, a(i_a)%n_configs
          call print(""//reshape(a(i_a)%water_orientation_histograms(:,:,i), (/ a(i_a)%water_orientation_n_pos_bins*a(i_a)%water_orientation_n_angle_bins /) ), file=outfile)
        end do

      else
        call system_abort("print_analyses: no type of analysis set for " // i_a)
      endif
    endif
    call finalise(outfile)
  end do
end subroutine print_analyses

!read geometry parameters to calculate
!format:
!4               number of params
!anything        comment
!1 3             type #1: y-coord  of atom  #3
!2 4 5           type #2: distance of atoms #4--#5
!3 3 4 5         type #3: angle    of atoms #3--#4--#5
!4 3 1 2 4       type #4: dihedral of atoms #3--#1--#2--#4
subroutine read_geometry_params(this,filename)

  type(analysis), intent(inout) :: this
  character(*), intent(in) :: filename

  type(inoutput) :: geom_lib
  character(20), dimension(10) :: fields
  integer :: num_fields, status
  integer :: num_geom, i, geom_type
  character(STRING_LENGTH) :: comment

  call initialise(this%geometry_params,5,0,0,0,0) !type, atom1, atom2, atom3, atom4

  if (trim(filename)=="") then
     call print('WARNING! no file specified')
     return !with an empty table
  endif

  call initialise(geom_lib,trim(filename),action=INPUT)
  call parse_line(geom_lib,' ',fields,num_fields)
  if (num_fields < 1) then
     call print ('WARNING! empty file '//trim(filename))
     call finalise(geom_lib)
     return !with an empty table
  endif

  num_geom = string_to_int(fields(1))
  comment=""
  comment = read_line(geom_lib,status)
  call print(trim(comment),PRINT_VERBOSE)

  do i=1,num_geom
    call parse_line(geom_lib,' ',fields,num_fields)
    if (num_fields.gt.5 .or. num_fields.lt.2) call system_abort('read_geometry_params: 1 type and maximum 4 atoms must be in the geometry file')
    geom_type = string_to_int(fields(1))
    select case (geom_type)
      case (1) !y coord atom1
        if (num_fields.lt.2) call system_abort('type 1: coordinate, =1 atoms needed')
        call append(this%geometry_params,(/geom_type, &
                                           string_to_int(fields(2)), &
                                           0, &
                                           0, &
                                           0/) )
      case (2) !distance atom1-atom2
        if (num_fields.lt.3) call system_abort('type 2: bond length, =2 atoms needed')
        call append(this%geometry_params, (/geom_type, &
                                            string_to_int(fields(2)), &
                                            string_to_int(fields(3)), &
                                            0, &
                                            0/) )
      case (3) !angle atom1-atom2-atom3
        if (num_fields.lt.4) call system_abort('type 3: angle, =3 atoms needed')
        call append(this%geometry_params, (/geom_type, &
                                            string_to_int(fields(2)), &
                                            string_to_int(fields(3)), &
                                            string_to_int(fields(4)), &
                                            0/) )
      case (4) !dihedral atom1-atom2-atom3-atom4
        if (num_fields.lt.5) call system_abort('type 4: dihedral, =4 atoms needed')
        call append(this%geometry_params, (/geom_type, &
                                            string_to_int(fields(2)), &
                                            string_to_int(fields(3)), &
                                            string_to_int(fields(4)), &
                                            string_to_int(fields(5)) /) )
      case (5) !distance of any atom3(Z) from the atom1-atom2 bond
        if (num_fields.lt.4) call system_abort('type 5: Zatom3 distance from the atom1-atom2 bond, =3 atoms needed')
        call append(this%geometry_params, (/geom_type, &
                                            string_to_int(fields(2)), &
                                            string_to_int(fields(3)), &
                                            string_to_int(fields(4)), &
                                            0/) )
      case default
        call system_abort('unknown type '//geom_type//', must be one of 1(y coordinate), 2(bond length/distance), 3(angle), 4(dihedral).')
    end select
  enddo

  call finalise(geom_lib)

end subroutine read_geometry_params

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

  character(len=STRING_LENGTH) :: infilename
  integer :: decimation
  type(Inoutput) :: list_infile
  type (CInoutput) :: infile
  logical :: infile_is_list
  logical :: quiet

  character(len=STRING_LENGTH) :: commandfilename
  type(Inoutput) :: commandfile

  character(len=10240) :: args_str, myline
  integer :: n_analysis_a
  type(analysis), allocatable :: analysis_a(:)

  logical :: more_files
  integer :: status, arg_line_no
  integer :: error = ERROR_NONE
  real(dp) :: time

  integer :: i_a, frame_count, raw_frame_count
  type(Atoms) :: structure
  logical :: do_verbose

  call system_initialise(PRINT_NORMAL)

  call initialise(cli_params)
  call param_register(cli_params, "verbose", "F", do_verbose, help_string="Set verbosity level to VERBOSE.")
  if (.not. param_read_args(cli_params, ignore_unknown=.true.)) &
    call system_abort("Impossible failure to parse verbosity")
  call finalise(cli_params)
  if (do_verbose) then
    call system_initialise(verbosity=PRINT_VERBOSE)
  endif

  call initialise(cli_params)
  call param_register(cli_params, "commandfile", '', commandfilename, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params, "infile", "stdin", infilename, help_string="Input file name (default: stdin)")
  call param_register(cli_params, "decimation", "1", decimation, help_string="How many frames to advance in each step, i.e. every n-th frame will be analysed.")
  call param_register(cli_params, "infile_is_list", "F", infile_is_list, help_string="Treat infile as a list of filenames to process sequentially.")
  call param_register(cli_params, "quiet", "F", quiet, help_string="Suppress frame progress count.")
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
    call read(structure, infile, error=error, frame=raw_frame_count)
    do while (error == ERROR_NONE)
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
      call read(structure, infile, error=error, frame=raw_frame_count)
    end do
    ! We do not want to handle this error
    call clear_error(error)
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
