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

! do various structural analysis to a trajcetory, save an intermediate file
! to be postprocessed (mean, variance, correlation, etc)
! #/vol density on a radial mesh
! #/vol density on a grid
! RDF, possibly as a function of distance of center from a fixed point
! to be added: ADF

module structure_analysis_module
use libatoms_module
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
  call param_register(params, 'infile', '', dummy_c_1, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'commandfile', '', dummy_c_2, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'decimation', '0', dummy_i_1, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'infile_is_list', 'F', dummy_l_1, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'quiet', 'F', dummy_l_2, help_string="No help yet.  This source file was $LastChangedBy$")

  if (.not. present(prev)) then
    ! general
    call param_register(params, 'outfile', 'stdout', this%outfilename, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'AtomMask', '', this%mask_str, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'min_time', '-1.0', this%min_time, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'max_time', '-1.0', this%max_time, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'min_frame', '-1', this%min_frame, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'max_frame', '-1', this%max_frame, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'type', '', this%type, help_string="[density_radial | density_grid | KE_density_radial | rdfd | adfd | KEdf_radial | propdf_radial | geometry | density_axial_silica | num_hbond_silica | water_orientation_silica ]")
    ! radial density
    call param_register(params, 'density_radial_bin_width', '-1.0', this%density_radial_bin_width, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'density_radial_n_bins', '-1', this%density_radial_n_bins, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'density_radial_center', '0.0 0.0 0.0', this%density_radial_center, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'density_radial_center_at', '-1', this%density_radial_center_at, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'density_radial_sigma', '1.0', this%density_radial_gaussian_sigma, help_string="No help yet.  This source file was $LastChangedBy$")

    ! grid density
    call param_register(params, 'density_grid_min_p', '0.0 0.0 0.0', this%density_grid_min_p, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'density_grid_bin_width', '-1.0 -1.0 -1.0', this%density_grid_bin_width, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'density_grid_n_bins', '-1 -1 -1', this%density_grid_n_bins, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'density_grid_gaussian', 'F', this%density_grid_gaussian_smoothing, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'density_grid_sigma', '1.0', this%density_grid_gaussian_sigma, help_string="No help yet.  This source file was $LastChangedBy$")

    ! radial KE density
    call param_register(params, 'KE_density_radial_bin_width', '-1.0', this%KE_density_radial_bin_width, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'KE_density_radial_n_bins', '-1', this%KE_density_radial_n_bins, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'KE_density_radial_center', '0.0 0.0 0.0', this%KE_density_radial_center, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'KE_density_radial_center_at', '-1', this%KE_density_radial_center_at, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'KE_density_radial_sigma', '1.0', this%KE_density_radial_gaussian_sigma, help_string="No help yet.  This source file was $LastChangedBy$")

    ! rdfd
    call param_register(params, 'rdfd_zone_center', '0.0 0.0 0.0', this%rdfd_zone_center, help_string="Cartesian position to center zones for rdfd (ignored in rdfd_zone_center is usable)")
    call param_register(params, 'rdfd_zone_atom_center', '0', this%rdfd_zone_atom_center, help_string="If > 0, atom to center zones for rdfd")
    call param_register(params, 'rdfd_zone_width', '-1.0', this%rdfd_zone_width, help_string="Width of zones for rdfd")
    call param_register(params, 'rdfd_n_zones', '1', this%rdfd_n_zones, help_string="Number of zones for rdfd")
    call param_register(params, 'rdfd_bin_width', '-1', this%rdfd_bin_width, help_string="Width of bin (or distance between sample points) in each rdfd")
    call param_register(params, 'rdfd_n_bins', '-1', this%rdfd_n_bins, help_string="Number of bins (or sample points) in each rdfd")
    call param_register(params, 'rdfd_center_mask', '', this%rdfd_center_mask_str, help_string="Mask for atoms that are considered for rdfd centers")
    call param_register(params, 'rdfd_neighbour_mask', '', this%rdfd_neighbour_mask_str, help_string="Mask for atoms that are considered for rdfd neighbours")
    call param_register(params, 'rdfd_gaussian', 'F', this%rdfd_gaussian_smoothing, help_string="If true, use Gaussians to smear mass distributions being integrated for rdfd")
    call param_register(params, 'rdfd_sigma', '0.1', this%rdfd_gaussian_sigma, help_string="Width of Gaussians used to smear mass distributions for rdfd")

    ! adfd
    call param_register(params, 'adfd_zone_center', '0.0 0.0 0.0', this%adfd_zone_center, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'adfd_zone_width', '-1.0', this%adfd_zone_width, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'adfd_n_zones', '1', this%adfd_n_zones, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'adfd_dist_bin_width', '-1', this%adfd_dist_bin_width, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'adfd_n_dist_bins', '-1', this%adfd_n_dist_bins, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'adfd_n_angle_bins', '-1', this%adfd_n_angle_bins, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'adfd_center_mask', '', this%adfd_center_mask_str, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'adfd_neighbour_1_mask', '', this%adfd_neighbour_1_mask_str, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'adfd_neighbour_1_max_dist', '-1.0', this%adfd_neighbour_1_max_dist, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'adfd_neighbour_2_mask', '', this%adfd_neighbour_2_mask_str, help_string="No help yet.  This source file was $LastChangedBy$")

    ! r-dep KE distribution
    call param_register(params, 'KEdf_radial_zone_center', '0.0 0.0 0.0', this%KEdf_radial_zone_center, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'KEdf_radial_zone_center_at', '-1', this%KEdf_radial_zone_center_at, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'KEdf_radial_zone_width', '-1.0', this%KEdf_radial_zone_width, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'KEdf_radial_bin_width', '0.002', this%KEdf_radial_bin_width, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'KEdf_radial_n_bins', '500', this%KEdf_radial_n_bins, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'KEdf_radial_n_zones', '1', this%KEdf_radial_n_zones, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'KEdf_radial_gaussian_sigma', '0.005', this%KEdf_radial_gaussian_sigma, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'KEdf_radial_mask', '', this%KEdf_radial_mask_str, help_string="No help yet.  This source file was $LastChangedBy$")

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

    ! geometry
    call param_register(params, 'geometry_filename', '', this%geometry_filename, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'geometry_central_atom', '-1', this%geometry_central_atom, help_string="No help yet.  This source file was $LastChangedBy$")

    ! uniaxial density silica
    call param_register(params, 'density_axial_n_bins', '-1', this%density_axial_n_bins, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'density_axial_axis', '-1', this%density_axial_axis, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'density_axial_silica_atoms', '-1', this%density_axial_silica_atoms, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'density_axial_gaussian', 'F', this%density_axial_gaussian_smoothing, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'density_axial_sigma', '1.0', this%density_axial_gaussian_sigma, help_string="No help yet.  This source file was $LastChangedBy$")

    ! num_hbond_silica
    call param_register(params, 'num_hbond_n_bins', '-1', this%num_hbond_n_bins, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'num_hbond_axis', '0', this%num_hbond_axis, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'num_hbond_silica_atoms', '0', this%num_hbond_silica_atoms, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'num_hbond_gaussian', 'F', this%num_hbond_gaussian_smoothing, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'num_hbond_sigma', '1.0', this%num_hbond_gaussian_sigma, help_string="No help yet.  This source file was $LastChangedBy$")

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
    ! general
    call param_register(params, 'outfile', trim(prev%outfilename), this%outfilename, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'AtomMask', trim(prev%mask_str), this%mask_str, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'min_time', ''//prev%min_time, this%min_time, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'max_time', ''//prev%max_time, this%max_time, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'min_frame', ''//prev%min_frame, this%min_frame, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'max_frame', ''//prev%max_frame, this%max_frame, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'type', ''//trim(prev%type), this%type, help_string="[density_radial | density_grid | KE_density_radial | rdfd | adfd | KEdf_radial | propdf_radial | geometry | density_axial_silica | num_hbond_silica | water_orientation_silica ]")

    ! radial density
    call param_register(params, 'density_radial_bin_width', ''//this%density_radial_bin_width, this%density_radial_bin_width, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'density_radial_n_bins', ''//this%density_radial_n_bins, this%density_radial_n_bins, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'density_radial_center', ''//prev%density_radial_center, this%density_radial_center, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'density_radial_center_at', ''//prev%density_radial_center_at, this%density_radial_center_at, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'density_radial_sigma', ''//prev%density_radial_gaussian_sigma, this%density_radial_gaussian_sigma, help_string="No help yet.  This source file was $LastChangedBy$")

    ! grid density
    call param_register(params, 'density_grid_min_p', ''//prev%density_grid_min_p, this%density_grid_min_p, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'density_grid_bin_width', ''//prev%density_grid_bin_width, this%density_grid_bin_width, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'density_grid_n_bins', ''//prev%density_grid_n_bins, this%density_grid_n_bins, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'density_grid_gaussian', ''//prev%density_grid_gaussian_smoothing, this%density_grid_gaussian_smoothing, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'density_grid_sigma', ''//prev%density_grid_gaussian_sigma, this%density_grid_gaussian_sigma, help_string="No help yet.  This source file was $LastChangedBy$")

    ! radial KE_density
    call param_register(params, 'KE_density_radial_bin_width', ''//this%KE_density_radial_bin_width, this%KE_density_radial_bin_width, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'KE_density_radial_n_bins', ''//this%KE_density_radial_n_bins, this%KE_density_radial_n_bins, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'KE_density_radial_center_at', ''//prev%KE_density_radial_center_at, this%KE_density_radial_center_at, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'KE_density_radial_sigma', ''//prev%KE_density_radial_gaussian_sigma, this%KE_density_radial_gaussian_sigma, help_string="No help yet.  This source file was $LastChangedBy$")

    ! rdfd
    call param_register(params, 'rdfd_zone_center', ''//prev%rdfd_zone_center, this%rdfd_zone_center, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'rdfd_zone_atom_center', ''//prev%rdfd_zone_atom_center, this%rdfd_zone_atom_center, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'rdfd_zone_width', ''//prev%rdfd_zone_width, this%rdfd_zone_width, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'rdfd_n_zones', ''//prev%rdfd_n_zones, this%rdfd_n_zones, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'rdfd_bin_width', ''//prev%rdfd_bin_width, this%rdfd_bin_width, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'rdfd_n_bins', ''//prev%rdfd_n_bins, this%rdfd_n_bins, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'rdfd_center_mask', trim(prev%rdfd_center_mask_str), this%rdfd_center_mask_str, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'rdfd_neighbour_mask', trim(prev%rdfd_neighbour_mask_str), this%rdfd_neighbour_mask_str, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'rdfd_gaussian', ''//prev%rdfd_gaussian_smoothing, this%rdfd_gaussian_smoothing, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'rdfd_sigma', ''//prev%rdfd_gaussian_sigma, this%rdfd_gaussian_sigma, help_string="No help yet.  This source file was $LastChangedBy$")

    ! r-dep KE distribution
    call param_register(params, 'KEdf_radial_zone_center', ''//prev%KEdf_radial_zone_center, this%KEdf_radial_zone_center, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'KEdf_radial_zone_center_at', ''//prev%KEdf_radial_zone_center_at, this%KEdf_radial_zone_center_at, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'KEdf_radial_zone_width', ''//prev%KEdf_radial_zone_width, this%KEdf_radial_zone_width, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'KEdf_radial_bin_width', ''//prev%KEdf_radial_bin_width, this%KEdf_radial_bin_width, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'KEdf_radial_n_bins', ''//prev%KEdf_radial_n_bins, this%KEdf_radial_n_bins, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'KEdf_radial_n_zones', ''//prev%KEdf_radial_n_zones, this%KEdf_radial_n_zones, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'KEdf_radial_gaussian_sigma', ''//prev%KEdf_radial_gaussian_sigma, this%KEdf_radial_gaussian_sigma, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'KEdf_radial_mask', ''//trim(prev%KEdf_radial_mask_str), this%KEdf_radial_mask_str, help_string="No help yet.  This source file was $LastChangedBy$")

    ! r-dep |F| distribution
    call param_register(params, 'propdf_radial_zone_center', ''//prev%propdf_radial_zone_center, this%propdf_radial_zone_center, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'propdf_radial_zone_center_at', ''//prev%propdf_radial_zone_center_at, this%propdf_radial_zone_center_at, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'propdf_radial_zone_width', ''//prev%propdf_radial_zone_width, this%propdf_radial_zone_width, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'propdf_radial_n_zones', ''//prev%propdf_radial_n_zones, this%propdf_radial_n_zones, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'propdf_radial_bin_width', ''//prev%propdf_radial_bin_width, this%propdf_radial_bin_width, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'propdf_radial_n_bins', ''//prev%propdf_radial_n_bins, this%propdf_radial_n_bins, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'propdf_radial_gaussian_sigma', ''//prev%propdf_radial_gaussian_sigma, this%propdf_radial_gaussian_sigma, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'propdf_radial_mask', ''//trim(prev%propdf_radial_mask_str), this%propdf_radial_mask_str, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'propdf_radial_property', ''//trim(prev%propdf_radial_property), this%propdf_radial_property, help_string="No help yet.  This source file was $LastChangedBy$")

    ! adfd
    call param_register(params, 'adfd_zone_center', ''//prev%adfd_zone_center, this%adfd_zone_center, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'adfd_zone_width', ''//prev%adfd_zone_width, this%adfd_zone_width, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'adfd_n_zones', ''//prev%adfd_n_zones, this%adfd_n_zones, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'adfd_dist_bin_width', ''//prev%adfd_dist_bin_width, this%adfd_dist_bin_width, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'adfd_n_dist_bins', ''//prev%adfd_n_dist_bins, this%adfd_n_dist_bins, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'adfd_n_angle_bins', ''//prev%adfd_n_angle_bins, this%adfd_n_angle_bins, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'adfd_center_mask', trim(prev%adfd_center_mask_str), this%adfd_center_mask_str, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'adfd_neighbour_1_mask', trim(prev%adfd_neighbour_1_mask_str), this%adfd_neighbour_1_mask_str, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'adfd_neighbour_1_max_dist', ''//prev%adfd_neighbour_1_max_dist, this%adfd_neighbour_1_max_dist, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'adfd_neighbour_2_mask', trim(prev%adfd_neighbour_2_mask_str), this%adfd_neighbour_2_mask_str, help_string="No help yet.  This source file was $LastChangedBy$")

    ! geometry
    call param_register(params, 'geometry_filename', ''//trim(prev%geometry_filename), this%geometry_filename, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'geometry_central_atom', ''//prev%geometry_central_atom, this%geometry_central_atom, help_string="No help yet.  This source file was $LastChangedBy$")

    ! uniaxial density silica
    call param_register(params, 'density_axial_n_bins', ''//prev%density_axial_n_bins, this%density_axial_n_bins, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'density_axial_axis', ''//prev%density_axial_axis, this%density_axial_axis, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'density_axial_silica_atoms', ''//prev%density_axial_silica_atoms, this%density_axial_silica_atoms, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'density_axial_gaussian', ''//prev%density_axial_gaussian_smoothing, this%density_axial_gaussian_smoothing, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'density_axial_sigma', ''//prev%density_axial_gaussian_sigma, this%density_axial_gaussian_sigma, help_string="No help yet.  This source file was $LastChangedBy$")

    ! num_hbond_silica
    call param_register(params, 'num_hbond_n_bins', ''//prev%num_hbond_n_bins, this%num_hbond_n_bins, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'num_hbond_axis', ''//prev%num_hbond_axis, this%num_hbond_axis, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'num_hbond_silica_atoms', ''//prev%num_hbond_silica_atoms, this%num_hbond_silica_atoms, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'num_hbond_gaussian', ''//prev%num_hbond_gaussian_smoothing, this%num_hbond_gaussian_smoothing, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'num_hbond_sigma', ''//prev%num_hbond_gaussian_sigma, this%num_hbond_gaussian_sigma, help_string="No help yet.  This source file was $LastChangedBy$")

    ! water_orientation_silica
    call param_register(params, 'water_orientation_n_pos_bins', ''//prev%water_orientation_n_pos_bins, this%water_orientation_n_pos_bins, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'water_orientation_n_angle_bins', ''//prev%water_orientation_n_angle_bins, this%water_orientation_n_angle_bins, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'water_orientation_axis', ''//prev%water_orientation_axis, this%water_orientation_axis, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'water_orientation_silica_atoms', ''//prev%water_orientation_silica_atoms, this%water_orientation_silica_atoms, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'water_orientation_gaussian', ''//prev%water_orientation_gaussian_smoothing, this%water_orientation_gaussian_smoothing, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'water_orientation_pos_sigma', ''//prev%water_orientation_pos_gaussian_sigma, this%water_orientation_pos_gaussian_sigma, help_string="No help yet.  This source file was $LastChangedBy$")
    !call param_register(params, 'water_orientation_angle_sigma', ''//prev%water_orientation_angle_gaussian_sigma, this%water_orientation_angle_gaussian_sigma, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'water_orientation_use_dipole', ''//prev%water_orientation_use_dipole, this%water_orientation_use_dipole, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'water_orientation_use_HOHangle_bisector', ''//prev%water_orientation_use_HOHangle_bisector, this%water_orientation_use_HOHangle_bisector, help_string="No help yet.  This source file was $LastChangedBy$")

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
      " types of analysis.  Possiblities: density_radial, density_grid, KE_density_radial, rdfd, adfd, KEdf_radial, propdf_radial, geometry, " // &
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
            a(i_a)%adfd_center_mask_str, a(i_a)%adfd_neighbour_1_mask_str, a(i_a)%adfd_neighbour_1_max_dist, a(i_a)%adfd_neighbour_2_mask_str, &
            a(i_a)%adfd_angle_bin_pos, a(i_a)%adfd_dist_bin_pos, a(i_a)%adfd_zone_pos)
        else
          call adfd_calc(a(i_a)%adfds(:,:,:,a(i_a)%n_configs), at, a(i_a)%adfd_zone_center, &
	    a(i_a)%adfd_n_angle_bins, a(i_a)%adfd_dist_bin_width, a(i_a)%adfd_n_dist_bins, &
            a(i_a)%adfd_zone_width, a(i_a)%adfd_n_zones, &
            a(i_a)%adfd_center_mask_str, a(i_a)%adfd_neighbour_1_mask_str, a(i_a)%adfd_neighbour_1_max_dist, a(i_a)%adfd_neighbour_2_mask_str)
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

subroutine density_sample_radial_mesh_Gaussians(histogram, at, center_pos, center_i, rad_bin_width, n_rad_bins, gaussian_sigma, mask_str, radial_pos, accumulate, quantity)
  real(dp), intent(inout) :: histogram(:)
  type(Atoms), intent(inout) :: at
  real(dp), intent(in), optional :: center_pos(3)
  integer, intent(in), optional :: center_i
  real(dp), intent(in) :: rad_bin_width
  integer, intent(in) :: n_rad_bins
  real(dp), intent(in) :: gaussian_sigma
  character(len=*), optional, intent(in) :: mask_str
  real(dp), intent(out), optional :: radial_pos(:)
  logical, optional, intent(in) :: accumulate
  character(len=*), optional, intent(in) :: quantity

  logical :: my_accumulate

  real(dp) :: use_center_pos(3), d, r0, s_sq, s_cu, h_val
  real(dp), parameter :: SQROOT_PI = sqrt(PI), PI_THREE_HALVES = PI**(3.0_dp/2.0_dp)
  logical, allocatable :: mask_a(:)
  integer :: at_i, rad_sample_i
  real(dp) :: rad_sample_r, exp_arg, ep, em
  logical :: quantity_1, quantity_KE
  real(dp) :: sample_weight

  if (present(center_pos) .and. present(center_i)) then
    call system_abort("density_sample_radial_mesh_Gaussians received both center_pos and center_i")
  else if (.not. present(center_pos) .and. .not. present(center_i)) then
    call system_abort("density_sample_radial_mesh_Gaussians received neither center_pos nor center_i")
  endif

  quantity_1 = .false.
  quantity_KE = .false.
  if (present(quantity)) then
    select case (quantity)
      case("1")
	quantity_1 = .true.
      case("KE")
	quantity_KE = .true.
      case default
	call system_abort("density_sample_radial_mesh_Gaussians called with unknown quantity='"//trim(quantity)//"'")
    end select
  else
    quantity_1 = .true.
  endif
  if (quantity_1) sample_weight = 1.0_dp

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
  ! ratio of 20/2=10 is definitely fine
  if (min(norm(at%lattice(:,1)), norm(at%lattice(:,2)), norm(at%lattice(:,3))) < 9.0_dp*gaussian_sigma) &
    call print("WARNING: at%lattice may be too small for sigma, errors (noticeably too low a density) may result", PRINT_ALWAYS)

  if (present(center_i)) then
    use_center_pos = at%pos(:,center_i)
  else
    use_center_pos = center_pos
  end if

  s_sq = gaussian_sigma**2
  s_cu = s_sq*gaussian_sigma
  do at_i=1, at%N
    if (.not. mask_a(at_i)) cycle
    if (present(center_i)) then
      if (at_i == center_i) cycle
    endif
    if (quantity_KE) then
      if (associated(at%velo)) then
	sample_weight = kinetic_energy(ElementMass(at%Z(at_i)), at%velo(1:3,at_i))
      else
	sample_weight = 0.0_dp
      endif
    endif
    d = distance_min_image(at,use_center_pos,at%pos(:,at_i))
    ! skip if atom is outside range
    if (d > (n_rad_bins-1)*rad_bin_width+6.0_dp*gaussian_sigma) cycle
    do rad_sample_i=1, n_rad_bins
      rad_sample_r = (rad_sample_i-1)*rad_bin_width
      r0 = rad_sample_r
      ep = 0.0_dp
      exp_arg = -(r0+d)**2/s_sq
      if (exp_arg > -20.0_dp) ep = exp(exp_arg)
      em = 0.0_dp
      exp_arg = -(r0-d)**2/s_sq
      if (exp_arg > -20.0_dp) em = exp(exp_arg)
      ! should really fix d->0 limit
      if (d .feq. 0.0_dp) then
	h_val = 0.0_dp
	exp_arg = -r0**2/s_sq
	if (exp_arg > -20.0_dp) &
	  h_val = exp(exp_arg) / (PI_THREE_HALVES * s_cu)
      else if (r0 .feq. 0.0_dp) then
	h_val = 0.0_dp
	exp_arg = -d**2/s_sq
	if (exp_arg > -20.0_dp) &
	  h_val = exp(exp_arg) / (PI_THREE_HALVES * s_cu)
      else
	h_val = (r0/(SQROOT_PI * gaussian_sigma * d) * (em - ep)) /  (4.0_dp * PI * r0**2)
      endif
      histogram(rad_sample_i) = histogram(rad_sample_i) +  sample_weight*h_val
    end do ! rad_sample_i
  end do ! at_i

end subroutine density_sample_radial_mesh_Gaussians

subroutine rdfd_calc(rdfd, at, zone_center, zone_atom_center, bin_width, n_bins, zone_width, n_zones, gaussian_smoothing, gaussian_sigma, &
                     center_mask_str, neighbour_mask_str, bin_pos, zone_pos)
  real(dp), intent(inout) :: rdfd(:,:)
  type(Atoms), intent(inout) :: at
  real(dp), intent(in) :: zone_center(3), bin_width, zone_width
  integer :: zone_atom_center
  integer, intent(in) :: n_bins, n_zones
  logical, intent(in) :: gaussian_smoothing
  real(dp), intent(in) :: gaussian_sigma
  character(len=*), intent(in) :: center_mask_str, neighbour_mask_str
  real(dp), intent(inout), optional :: bin_pos(:), zone_pos(:)

  logical, allocatable :: center_mask_a(:), neighbour_mask_a(:)
  integer :: i_at, j_at, i_bin, i_zone
  integer, allocatable :: n_in_zone(:)
  real(dp) :: r, bin_inner_rad, bin_outer_rad, my_zone_center(3)

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

!  if (gaussian_smoothing) then
!    call set_cutoff(at, n_bins*bin_width+5.0_dp*gaussian_sigma)
!    call calc_connect(at)
!  endif

  rdfd = 0.0_dp
  do i_at=1, at%N ! loop over center atoms
    if (.not. center_mask_a(i_at)) cycle

    !calc which zone the atom is in
    if (zone_width > 0.0_dp) then
      if (zone_atom_center > 0) then
	 if (zone_atom_center <= at%N) then
	    my_zone_center = at%pos(:,zone_atom_center)
	 else
	    call system_abort("rdfd_calc got zone_atom_center="//zone_atom_center//" out of range (1.."//at%N//")")
	 endif
      else
	 my_zone_center = zone_center
      endif
      r = distance_min_image(at, my_zone_center, at%pos(:,i_at))
      i_zone = int(r/zone_width)+1
      if (i_zone > n_zones) cycle
    else
      i_zone = 1
    endif

    !count the number of atoms in that zone
    n_in_zone(i_zone) = n_in_zone(i_zone) + 1

    !calc rdfd in each bin for this zone
    if (gaussian_smoothing) then
      call density_sample_radial_mesh_Gaussians(rdfd(:,i_zone), at, center_i=i_at, rad_bin_width=bin_width, n_rad_bins=n_bins, &
        gaussian_sigma=gaussian_sigma, mask_str=neighbour_mask_str, accumulate = .true.)
    else
      !loop over atoms and advance the bins
      do j_at=1, at%N
        if (j_at == i_at) cycle
        if (.not. neighbour_mask_a(j_at)) cycle
        r = distance_min_image(at, i_at, j_at)
        i_bin = int(r/bin_width)+1
        if (i_bin <= n_bins) rdfd(i_bin,i_zone) = rdfd(i_bin,i_zone) + 1.0_dp
      end do ! j_at
    endif ! gaussian_smoothing
  end do ! i_at

  if (.not. gaussian_smoothing) then
    !calculate local density by dividing bins with their volumes
    do i_bin=1, n_bins
      bin_inner_rad = real(i_bin-1,dp)*bin_width
      bin_outer_rad = real(i_bin,dp)*bin_width
      rdfd(i_bin,:) = rdfd(i_bin,:)/(4.0_dp/3.0_dp*PI*bin_outer_rad**3 - 4.0_dp/3.0_dp*PI*bin_inner_rad**3)
    end do
  end if

  ! normalise zones by the number of atoms in that zone
  do i_zone=1, n_zones
    if (n_in_zone(i_zone) > 0) rdfd(:,i_zone) = rdfd(:,i_zone)/real(n_in_zone(i_zone),dp)
  end do
  ! normalise with the global density
  if (count(neighbour_mask_a) > 0) then
    rdfd = rdfd / (count(neighbour_mask_a)/cell_volume(at))
  endif

end subroutine rdfd_calc

subroutine propdf_radial_calc(histograms, at, bin_width, n_bins, &
  zone_center, zone_width, n_zones, gaussian_sigma, mask_str, property, bin_pos, zone_pos)
  real(dp), intent(inout) :: histograms(:,:)
  type(Atoms), intent(in) :: at
  real(dp), intent(in) :: bin_width, zone_width
  integer, intent(in) :: n_bins, n_zones
  real(dp), intent(in) :: gaussian_sigma
  real(dp), intent(in) :: zone_center(3)
  character(len=*), intent(in) :: mask_str, property
  real(dp), intent(inout), optional :: bin_pos(:), zone_pos(:)

  logical, allocatable :: mask_a(:)
  integer :: i_at, i_bin_ctr, i_bin, i_zone
  integer, allocatable :: n_in_zone(:)
  real(dp) :: n, quant, r, bin_ctr
  logical :: has_mass
  logical :: doing_KE
  character(len=100) :: use_property
  real(dp), pointer :: prop_a(:), prop_3a(:,:)
  logical :: prop_is_scalar

  allocate(mask_a(at%N))
  call is_in_mask(mask_a, at, mask_str)

  doing_KE = (trim(property) == 'KE') 

  if (len_trim(property) == 0) &
    call system_abort("propdf_radial_calc has no property specified.  Are you sure you passed a propdf_radial_property value?")

  if (doing_KE) then ! KE is a special case
    has_mass = has_property(at, 'mass')
    if (.not. has_property(at, 'velo')) &
      call system_abort("propdf_radial_calc has no 'velo' property in Atoms structure")
  else ! assign pointers for property
    if (.not. has_property(at, property)) &
       call system_abort("propdf_radial_calc has no '"//trim(property)//"' property in Atoms structure")

    prop_is_scalar = .true.
    if (.not. assign_pointer(at, trim(property), prop_a)) then
       prop_is_scalar = .false.
       if (.not. assign_pointer(at, trim(property), prop_3a)) then
	  call system_abort("propdf_radial_calc failed to assign 1- or 3-array real(dp) pointer for property '"//trim(property))
       endif
    endif
  endif

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
      bin_pos(i_bin) = (i_bin-1)*bin_width
    end do
  endif

  histograms = 0.0_dp
  do i_at=1, at%N ! loop over atoms
    if (.not. mask_a(i_at)) cycle

    !calc which zone the atom is in
    if (zone_width > 0.0_dp) then
      r = distance_min_image(at, zone_center, at%pos(:,i_at))
      i_zone = int(r/zone_width)+1
      if (i_zone > n_zones) cycle
    else
      i_zone = 1
    endif

    !count the number of atoms in that zone
    n_in_zone(i_zone) = n_in_zone(i_zone) + 1

    if (doing_KE) then
      !calc KE for this atom
      if (has_mass) then
	quant = kinetic_energy(at%mass(i_at), at%velo(:,i_at))
      else
	quant = kinetic_energy(ElementMass(at%Z(i_at)), at%velo(:,i_at))
      endif
    else
      if (prop_is_scalar) then
	quant = prop_a(i_at)
      else
	quant = norm(prop_3a(:,i_at))
      endif
    endif
    i_bin_ctr = quant/bin_width
    do i_bin=max(1,i_bin_ctr-floor(6*gaussian_sigma/bin_width)), min(n_bins,i_bin_ctr+floor(6*gaussian_sigma/bin_width))
      bin_ctr = (i_bin-1)*bin_width
      histograms(i_bin,i_zone) = histograms(i_bin,i_zone) + exp(-0.5*(quant-bin_ctr)**2/gaussian_sigma**2)
    end do
  end do ! i_at

  ! normalise zones by the number of atoms in that zone
  n = gaussian_sigma/sqrt(2.0_dp*PI)
  do i_zone=1, n_zones
    if (n_in_zone(i_zone) > 0) histograms(:,i_zone) = n*histograms(:,i_zone)/real(n_in_zone(i_zone),dp)
  end do

end subroutine propdf_radial_calc

subroutine adfd_calc(adfd, at, zone_center, n_angle_bins, dist_bin_width, n_dist_bins, zone_width, n_zones, &
		     center_mask_str, neighbour_1_mask_str, neighbour_1_max_dist, neighbour_2_mask_str, & 
		     angle_bin_pos, dist_bin_pos, zone_pos)
  real(dp), intent(inout) :: adfd(:,:,:)
  type(Atoms), intent(inout) :: at
  real(dp), intent(in) :: zone_center(3), dist_bin_width, zone_width
  integer, intent(in) :: n_angle_bins, n_dist_bins, n_zones
  character(len=*), intent(in) :: center_mask_str, neighbour_1_mask_str, neighbour_2_mask_str
  real(dp), intent(in) :: neighbour_1_max_dist
  real(dp), intent(inout), optional :: angle_bin_pos(:), dist_bin_pos(:), zone_pos(:)

  logical, allocatable :: center_mask_a(:), neighbour_1_mask_a(:), neighbour_2_mask_a(:)
  integer :: i_at, j_at, k_at, i_dist_bin, i_angle_bin, i_zone, ji, ki
  integer, allocatable :: n_in_zone(:)
  real(dp) :: dist_bin_inner_rad, dist_bin_outer_rad, angle_bin_width, angle_bin_min, angle_bin_max
  real(dp) :: r, r_ij, r_jk, jik_angle

  allocate(center_mask_a(at%N), neighbour_1_mask_a(at%N), neighbour_2_mask_a(at%N))
  call is_in_mask(center_mask_a, at, center_mask_str)
  call is_in_mask(neighbour_1_mask_a, at, neighbour_1_mask_str)
  call is_in_mask(neighbour_2_mask_a, at, neighbour_2_mask_str)

  allocate(n_in_zone(n_zones))
  n_in_zone = 0

  if (present(zone_pos)) then
    if (zone_width > 0.0_dp) then
      do i_zone= 1, n_zones
	zone_pos(i_zone) = (real(i_zone,dp)-0.5_dp)*zone_width
      end do
    else
      zone_pos(1) = -1.0_dp
    end if
  end if
  if (present(dist_bin_pos)) then
    do i_dist_bin=1, n_dist_bins
      dist_bin_pos(i_dist_bin) = (real(i_dist_bin,dp)-0.5_dp)*dist_bin_width
    end do
  end if
  angle_bin_width = PI/real(n_angle_bins,dp)
  if (present(angle_bin_pos)) then
    do i_angle_bin=1, n_angle_bins
      angle_bin_pos(i_angle_bin) = (real(i_angle_bin,dp)-0.5_dp)*angle_bin_width
    end do
  end if

  adfd = 0.0_dp
  call set_cutoff(at, n_dist_bins*dist_bin_width)
  call calc_connect(at)
  do i_at=1, at%N ! center atom
    if (.not. center_mask_a(i_at)) cycle
    if (zone_width > 0.0_dp) then
      r = distance_min_image(at, zone_center, at%pos(:,i_at))
      i_zone = int(r/zone_width)+1
      if (i_zone > n_zones) cycle
    else
      i_zone = 1
    endif
    n_in_zone(i_zone) = n_in_zone(i_zone) + 1
    do ji=1, atoms_n_neighbours(at, i_at)
      j_at = atoms_neighbour(at, i_at, ji, distance=r_ij)
      if (j_at == i_at) cycle
      if (.not. neighbour_1_mask_a(j_at)) cycle
      if (neighbour_1_max_dist > 0.0_dp) then
	if (r_ij > neighbour_1_max_dist) cycle
      else
	if (r_ij > bond_length(at%Z(i_at),at%Z(j_at))*at%nneightol) cycle
      endif
      do ki=1, atoms_n_neighbours(at, i_at)
	k_at = atoms_neighbour(at, i_at, ki)
	if (k_at == i_at .or. k_at == j_at) cycle
	if (.not. neighbour_2_mask_a(k_at)) cycle
	r_jk = distance_min_image(at, j_at, k_at)
	i_dist_bin = int(r_jk/dist_bin_width)+1
	if (i_dist_bin > n_dist_bins) cycle
! call print("doing triplet ijk " // i_at // " " // j_at // " "// k_at, PRINT_ALWAYS)
! call print("  Zijk " // at%Z(i_at) // " " // at%Z(j_at) // " " // at%Z(k_at), PRINT_ALWAYS)
! call print("  pi " // at%pos(:,i_at), PRINT_ALWAYS)
! call print("  pj " // at%pos(:,j_at), PRINT_ALWAYS)
! call print("  pk " // at%pos(:,k_at), PRINT_ALWAYS)
! call print("  r_ij " // diff_min_image(at,i_at,j_at) // "     " // r_ij, PRINT_ALWAYS)
! call print("  r_jk " // diff_min_image(at,j_at,k_at) // "     " // distance_min_image(at, j_at, k_at), PRINT_ALWAYS)
! call print("  r_ik " // diff_min_image(at,i_at,k_at) // "     " // distance_min_image(at, i_at, k_at), PRINT_ALWAYS)
        jik_angle = angle(diff_min_image(at,i_at,j_at), diff_min_image(at,i_at,k_at))
! call print("  r_ij " // r_ij // " r_jk " // r_jk // " jik_angle " // (jik_angle*180.0/PI), PRINT_ALWAYS)
	i_angle_bin = int(jik_angle/angle_bin_width)+1
	if (i_angle_bin > n_angle_bins) i_angle_bin = n_angle_bins
	adfd(i_angle_bin,i_dist_bin,i_zone) = adfd(i_angle_bin,i_dist_bin,i_zone) + 1.0_dp ! /(2.0_dp*PI*sin(jik_angle))
      end do ! k_at
    end do ! j_at
  end do ! i_at

  do i_angle_bin=1, n_angle_bins
    angle_bin_min = real(i_angle_bin-1,dp)*angle_bin_width
    angle_bin_max = real(i_angle_bin,dp)*angle_bin_width
    do i_dist_bin=1, n_dist_bins
      dist_bin_inner_rad = real(i_dist_bin-1,dp)*dist_bin_width
      dist_bin_outer_rad = real(i_dist_bin,dp)*dist_bin_width
      adfd(i_angle_bin,i_dist_bin,:) = adfd(i_angle_bin,i_dist_bin,:) / (2.0_dp*PI * ((dist_bin_outer_rad**3 - dist_bin_inner_rad**3)/3.0_dp) * (cos(angle_bin_min)-cos(angle_bin_max)))
      ! adfd(i_angle_bin,i_dist_bin,:) = adfd(i_angle_bin,i_dist_bin,:) / (4.0_dp/3.0_dp*PI * (dist_bin_outer_rad**3 - dist_bin_inner_rad**3))
    end do
  end do

  ! normalise zones by the number of atoms in that zone
  do i_zone=1, n_zones
    if (n_in_zone(i_zone) > 0) adfd(:,:,i_zone) = adfd(:,:,i_zone)/real(n_in_zone(i_zone),dp)
  end do
  ! normalise with the global density
  if (count(neighbour_2_mask_a) > 0) then
    adfd = adfd / (count(neighbour_2_mask_a)/cell_volume(at))
  endif

end subroutine adfd_calc

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
    call print("WARNING: at%lattice may be too small for sigma, errors (noticeably too low a density) may result", PRINT_ALWAYS)

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

!Calculates distances, angles, dihedrals of the given atoms
subroutine geometry_calc(histogram, at, geometry_params, central_atom, geometry_pos, geometry_label)

  real(dp), intent(inout) :: histogram(:)
  type(Atoms), intent(inout) :: at
  type(Table), intent(in) :: geometry_params
  integer, intent(in) :: central_atom
  real(dp), intent(out), optional :: geometry_pos(:)
  character(STRING_LENGTH), intent(out), optional :: geometry_label(:)

  integer :: i, j, geom_type
  integer :: atom1, atom2, atom3, atom4
  real(dp) :: shift(3), bond12(3),bond23(3),bond34(3)
  real(dp) :: a, b, c, dist, min_dist
  logical :: found

  !center around central_atom if requested
  if (central_atom.gt.at%N) call system_abort('central atom is greater than atom number '//at%N)
  if (central_atom.gt.0) then !center around central_atom
     shift = at%pos(1:3,central_atom)
     do j=1,at%N
        at%pos(1:3,j) = at%pos(1:3,j) - shift(1:3)
     enddo
     call map_into_cell(at) !only in this case, otherwise it has been mapped
  endif

  !loop over the parameters to calculate
  do i=1, geometry_params%N
     geom_type=geometry_params%int(1,i)
     atom1 = geometry_params%int(2,i)
     atom2 = geometry_params%int(3,i)
     atom3 = geometry_params%int(4,i)
     atom4 = geometry_params%int(5,i)
     select case (geom_type)
       case(1) !y coord atom1
         if (atom1<1.or.atom1>at%N) call system_abort('atom1 must be >0 and < '//at%N)
         histogram(i) = at%pos(2,atom1)
       case(2) !distance atom1-atom2
         if (atom1<1.or.atom1>at%N) call system_abort('atom1 must be >0 and < '//at%N)
         if (atom2<1.or.atom2>at%N) call system_abort('atom2 must be >0 and < '//at%N)
         !histogram(i) = norm(at%pos(1:3,atom1)-at%pos(1:3,atom2))
         histogram(i) = distance_min_image(at,atom1,atom2)
       case(3) !angle atom1-atom2-atom3
         if (atom1<1.or.atom1>at%N) call system_abort('atom1 must be >0 and < '//at%N)
         if (atom2<1.or.atom2>at%N) call system_abort('atom2 must be >0 and < '//at%N)
         if (atom3<1.or.atom2>at%N) call system_abort('atom3 must be >0 and < '//at%N)
         !histogram(i) = angle(at%pos(1:3,atom1)-at%pos(1:3,atom2), &
         !                     at%pos(1:3,atom3)-at%pos(1:3,atom2))
         histogram(i) = angle(diff_min_image(at,atom2,atom1), &
                              diff_min_image(at,atom2,atom3))
       case(4) !dihedral atom1-(bond12)->atom2-(bond23)->atom3-(bond34)->atom4
         if (atom1<1.or.atom1>at%N) call system_abort('atom1 must be >0 and < '//at%N)
         if (atom2<1.or.atom2>at%N) call system_abort('atom2 must be >0 and < '//at%N)
         if (atom3<1.or.atom2>at%N) call system_abort('atom3 must be >0 and < '//at%N)
         if (atom4<1.or.atom2>at%N) call system_abort('atom4 must be >0 and < '//at%N)
         !bond12(1:3) = at%pos(1:3,atom2)-at%pos(1:3,atom1)
         bond12(1:3) = diff_min_image(at,atom1,atom2)
         !bond23(1:3) = at%pos(1:3,atom3)-at%pos(1:3,atom2)
         bond23(1:3) = diff_min_image(at,atom2,atom3)
         !bond34(1:3) = at%pos(1:3,atom4)-at%pos(1:3,atom3)
         bond34(1:3) = diff_min_image(at,atom3,atom4)
         histogram(i) = atan2(norm(bond23(1:3)) * bond12(1:3).dot.(bond23(1:3).cross.bond34(1:3)), &
                              (bond12(1:3).cross.bond23(1:3)) .dot. (bond23(1:3).cross.bond34(1:3)))
       case(5) !distance of atom3-Z from atom1-atom2 bond
         min_dist=HUGE(1._dp)
         found=.false.
         do j=1, at%N
            if (j==atom1 .or. j==atom2) cycle
            if (at%Z(j)/=atom3) cycle
            if (j<=6) cycle !CH3Cl2- specific!!!
            !only close ones
            a = distance_min_image(at,atom1,atom2)
            b = distance_min_image(at,atom1,j)
            c = distance_min_image(at,atom2,j)
!            if (b>4.0_dp .or. c>4.0_dp) cycle
            if (b>a .or. c>a) cycle
            dist = 0.25_dp*sqrt((a+b+c)*(b+c-a)*(c+a-b)*(a+b-c)) / (0.5_dp*a)
            if (dist<min_dist) then
               found=.true.
               min_dist=dist
            endif
            if (found) then
               histogram(i) = min_dist
            else
               histogram(i) = 0._dp
            endif
         enddo
       case default
         call system_abort("geometry_calc: unknown geometry type "//geom_type)
     end select
     if (present(geometry_pos)) geometry_pos(i) = real(i,dp)
     if (present(geometry_label)) geometry_label(i) = geom_type//'=='//atom1//'--'//atom2//'--'//atom3//'--'//atom4
  enddo

end subroutine geometry_calc


!silica-water inerface: put half-half of the silica slab to the 2 edges of the cell along the axis  normal to the surface
subroutine shift_silica_to_edges(at, axis, silica_center_i, mask_str)
  type(Atoms), intent(inout) :: at
  integer, intent(in) :: axis
  integer, intent(in) :: silica_center_i
  character(len=*), optional, intent(in) :: mask_str

  logical, allocatable :: mask_a(:), mask_silica(:)
  integer :: i
  integer :: counter
  integer, allocatable :: Si_atoms(:)
  real(dp) :: com(3), shift(3)
  real(dp), pointer :: mass_p(:)

  if (.not.silica_center_i>0) return

  !find the com of the silica (using the Si atoms) and move that to the edges in this direction
  shift(1:3) = at%pos(1:3,1)
  do i=1,at%N
     at%pos(1:3,i) = at%pos(1:3,i) - shift(1:3)
  enddo
  call map_into_cell(at)
  allocate(mask_silica(at%N))
  call is_in_mask(mask_silica, at, "Si")
  allocate(Si_atoms(count(mask_silica)))
  counter = 0
  do i=1,at%N
     if (mask_silica(i)) then
        counter = counter + 1
        Si_atoms(counter) = i
     endif
  enddo
  !allocate(at%mass(at%N))
  !at%mass = ElementMass(at%Z)
  call add_property(at,"mass",0._dp)
  if (.not.(assign_pointer(at, "mass", mass_p))) call system_abort('??')
  mass_p = ElementMass(at%Z)
  com = centre_of_mass(at,index_list=Si_atoms(1:size(Si_atoms)),origin=1)
  !shift axis to the edge (-0.5*edge_length)
  at%pos(axis,1:at%N) = at%pos(axis,1:at%N) - 0.5_dp * at%lattice(axis,axis) - com(axis)
  call map_into_cell(at) !everyone becomes -0.5b<y<0.5b
  !NO !shift everyone to positive coordinate along axis
  !at%pos(axis,1:at%N) = at%pos(axis,1:at%N) + 0.5_dp * at%lattice(axis,axis)
  deallocate(Si_atoms)
  deallocate(mask_silica)

end subroutine shift_silica_to_edges

subroutine density_axial_calc(histogram, at, axis, silica_center_i,n_bins, gaussian_smoothing, gaussian_sigma, mask_str, axial_pos, accumulate)
  real(dp), intent(inout) :: histogram(:)
  type(Atoms), intent(inout) :: at
  integer, intent(in) :: axis
  integer, intent(in) :: silica_center_i
  integer, intent(in) :: n_bins
  logical, intent(in) :: gaussian_smoothing
  real(dp), intent(in) :: gaussian_sigma
  character(len=*), optional, intent(in) :: mask_str
  real(dp), intent(out), optional :: axial_pos(:)
  logical, optional, intent(in) :: accumulate

  logical :: my_accumulate
  real(dp) :: ax_sample_r, dist, r, exp_arg
  logical, allocatable :: mask_a(:)
  integer at_i, ax_sample_i, i
  real(dp) :: bin_width

  my_accumulate = optional_default(.false., accumulate)
  if (.not. my_accumulate) histogram = 0.0_dp

  !atommask
  allocate(mask_a(at%N))
  call is_in_mask(mask_a, at, mask_str)

  !bin labels
  bin_width = at%lattice(axis,axis) / n_bins
  if (present(axial_pos)) then
    do ax_sample_i=1, n_bins
      ax_sample_r = (real(ax_sample_i,dp)-0.5_dp)*bin_width !the middle of the bin
      axial_pos(ax_sample_i) = ax_sample_r
    end do
  endif

  if (silica_center_i>0) then ! silica
     call shift_silica_to_edges(at, axis, silica_center_i, mask_str)
  endif

  !now bins are from -axis/2 to axis/2
  !shift everyone to positive coordinate along axis
  at%pos(axis,1:at%N) = at%pos(axis,1:at%N) + 0.5_dp * at%lattice(axis,axis)


  !simply check distances and bin them


  ! ratio of 20/4=5 is bad
  ! ratio of 20/3=6.66 is bad
  ! ratio of 20/2.5=8 is borderline (2e-4)
  ! ratio of 20/2.22=9 is fine  (error 2e-5)
  ! ratio of 20/2=10 is fine
  if ( gaussian_smoothing .and. (at%lattice(axis,axis) < 9.0_dp*gaussian_sigma) ) &
    call print("WARNING: at%lattice may be too small for sigma, errors (noticeably too low a density) may result", PRINT_ALWAYS)

  if (gaussian_smoothing) then

    do at_i=1, at%N
      if (silica_center_i>0 .and. at_i<=silica_center_i) cycle !ignore silica atoms
      if (.not. mask_a(at_i)) cycle
      r = at%pos(axis,at_i)
      do ax_sample_i=1, n_bins

        ax_sample_r = (real(ax_sample_i,dp)-0.5_dp)*bin_width
        dist = abs(r - ax_sample_r)
!Include all the atoms, slow but minimises error
!	  if (dist > 4.0_dp*gaussian_sigma) cycle
          exp_arg = -0.5_dp*(dist/(gaussian_sigma))**2
          if (exp_arg > -20.0_dp) then ! good to about 1e-8
            histogram(ax_sample_i) = histogram(ax_sample_i) + exp(exp_arg)/(gaussian_sigma*sqrt(2.0_dp*PI)) !Gaussian in 1 dimension
          endif

      end do ! ax_sample_i
    end do ! at_i

  else !no gaussian_smoothing

    do at_i=1, at%N
      if (silica_center_i>0 .and. at_i<=silica_center_i) cycle !ignore silica atoms
      if (.not. mask_a(at_i)) cycle
      r = at%pos(axis,at_i)

        histogram(int(r/bin_width)+1) = histogram(int(r/bin_width)+1) + 1

    end do ! at_i

  endif

  deallocate(mask_a)
end subroutine density_axial_calc

!
! Calculates the number of H-bonds along an axis
!  -- 1st center around atom 1
!  -- then calculate COM of Si atoms
!  -- then shift the centre of mass along the y axis to y=0
!  -- calculate H bonds for
!        water - water
!        water - silica
!        silica - water
!        silica - silica
!     interactions
!  -- the definition of a H-bond (O1-H1 - - O2):
!        d(O1,O2) < 3.5 A
!        d(O2,H1) < 2.45 A
!        angle(H1,O1,O2) < 30 degrees
!     source: P. Jedlovszky, J.P. Brodholdt, F. Bruni, M.A. Ricci and R. Vallauri, J. Chem. Phys. 108, 8525 (1998)
!
subroutine num_hbond_calc(histogram, at, axis, silica_center_i,n_bins, gaussian_smoothing, gaussian_sigma, mask_str, num_hbond_pos, num_hbond_type_code, num_hbond_type_label,accumulate)
  real(dp),                          intent(inout) :: histogram(:,:)
  type(Atoms),                       intent(inout) :: at
  integer,                           intent(in)    :: axis
  integer,                           intent(in)    :: silica_center_i
  integer,                           intent(in)    :: n_bins
  logical,                           intent(in)    :: gaussian_smoothing
  real(dp),                          intent(in)    :: gaussian_sigma
  character(len=*),        optional, intent(in)    :: mask_str
  real(dp),                optional, intent(out)   :: num_hbond_pos(:)
  integer,                 optional, intent(out)   :: num_hbond_type_code(:)
  character(STRING_LENGTH), optional, intent(out)   :: num_hbond_type_label(:)
  logical,                 optional, intent(in)    :: accumulate

  logical :: my_accumulate
  real(dp) :: num_hbond_sample_r, dist, r, exp_arg
  logical, allocatable :: mask_a(:)
  integer :: num_hbond_sample_i
  real(dp) :: bin_width
  real(dp), parameter                   :: dist_O2_H1 = 2.45_dp
  real(dp), parameter                   :: dist_O1_O2 = 3.5_dp
  real(dp), parameter                   :: angle_H1_O1_O2 = 30._dp
real(dp) :: min_distance, distance, HOO_angle
integer :: H1, O1, i, j, k, O2, num_atoms, hbond_type

  my_accumulate = optional_default(.false., accumulate)
  if (.not. my_accumulate) histogram = 0.0_dp

  !atommask
  allocate(mask_a(at%N))
  call is_in_mask(mask_a, at, mask_str)

  !bin labels
  bin_width = at%lattice(axis,axis) / n_bins
  if (present(num_hbond_pos)) then
    do num_hbond_sample_i=1, n_bins
      num_hbond_sample_r = (real(num_hbond_sample_i,dp)-0.5_dp)*bin_width !the middle of the bin
      num_hbond_pos(num_hbond_sample_i) = num_hbond_sample_r
    end do
  endif
  if (present(num_hbond_type_label)) then
    num_hbond_type_label(1)="water-water"
    num_hbond_type_label(2)="water-silica"
    num_hbond_type_label(3)="silica-water"
    num_hbond_type_label(4)="silica-silica"
  endif
  if (present(num_hbond_type_code)) then
    num_hbond_type_code(1)=11
    num_hbond_type_code(2)=10
    num_hbond_type_code(3)=01
    num_hbond_type_code(4)=00
  endif

  if (silica_center_i>0) then ! silica
     call shift_silica_to_edges(at, axis, silica_center_i, mask_str)
  endif

  !!calc_connect now, before shifting positions to positive, because it would remap the positions!!
  !call calc_connect including the H-bonds
  call set_cutoff(at,dist_O2_H1)
  call calc_connect(at)

  !now bins are from -axis/2 to axis/2
  !shift everyone to positive coordinate along axis
  at%pos(axis,1:at%N) = at%pos(axis,1:at%N) + 0.5_dp * at%lattice(axis,axis)

  !simply check hbonds and bin them


  ! ratio of 20/4=5 is bad
  ! ratio of 20/3=6.66 is bad
  ! ratio of 20/2.5=8 is borderline (2e-4)
  ! ratio of 20/2.22=9 is fine  (error 2e-5)
  ! ratio of 20/2=10 is fine
  if ( gaussian_smoothing .and. (at%lattice(axis,axis) < 9.0_dp*gaussian_sigma) ) &
    call print("WARNING: at%lattice may be too small for sigma, errors (noticeably too low a density) may result", PRINT_ALWAYS)

!  call set_cutoff(at,dist_O2_H1)
!  call calc_connect(at)

  num_atoms = 0
  do H1=1, at%N
     if(at%Z(H1)/=1) cycle !find H: H1
     !Count the atoms
     call print('Found H'//H1//ElementName(at%Z(H1)),PRINT_ANAL)
     num_atoms = num_atoms + 1

     !find closest O: O1
     min_distance = huge(1._dp)
     O1 = 0
     k = 0
     do i = 1, atoms_n_neighbours(at,H1)
        j = atoms_neighbour(at,H1,i,distance)
        if (distance<min_distance) then
           min_distance = distance
           k = i  !the closest neighbour is the k-th one
           O1 = j
        endif
     enddo
     if (O1==0) call system_abort('H has no neighbours.')
     !if (at%Z(O1).ne.8) call system_abort('H'//H1//' has not O closest neighbour '//ElementName(at%Z(O1))//O1//'.')
     if (.not.mask_a(O1)) call system_abort('H'//H1//' has not O closest neighbour '//ElementName(at%Z(O1))//O1//'.')

     !loop over all other Os: O2
     do i = 1, atoms_n_neighbours(at,H1)
        if (i.eq.k) cycle
        O2 = atoms_neighbour(at,H1,i)
        !if (at%Z(O2).ne.8) cycle !only keep O
        if (.not. mask_a(O2)) cycle
        !check O1-O2 distance for definition
        if (distance_min_image(at,O1,O2).gt.dist_O1_O2) cycle
        !check H1-O1-O2 angle for definition
        HOO_angle = angle(diff_min_image(at,O1,H1), &
                          diff_min_image(at,O1,O2)) *180._dp/PI
        call print('HOO_ANGLE '//ElementName(at%Z(H1))//H1//' '//ElementName(at%Z(O1))//O1//' '//ElementName(at%Z(O2))//O2//' '//HOO_angle,PRINT_ANAL)
        if (HOO_angle.gt.angle_H1_O1_O2) cycle

        !We've found a H-bond.

        !Find out the type (what-to-what)
        if (O1>silica_center_i .and. O2>silica_center_i) then  ! water - water
           call print('Found water-water H-bond.',PRINT_ANAL)
           hbond_type = 1
        elseif (O1>silica_center_i .and. O2<=silica_center_i) then  ! water - silica
           call print('Found water-silica H-bond.',PRINT_ANAL)
           hbond_type = 2
        elseif (O1<=silica_center_i .and. O2>silica_center_i) then  ! silica - water
           call print('Found silica-water H-bond.',PRINT_ANAL)
           hbond_type = 3
        elseif (O1<=silica_center_i .and. O2<=silica_center_i) then  ! silica - silica
           call print('Found silica-silica H-bond.',PRINT_ANAL)
           hbond_type = 4
        endif

        !Build histogram
        r = at%pos(axis,H1) !the position of H1

        if (gaussian_smoothing) then !smear the position along axis
           do num_hbond_sample_i=1, n_bins
             num_hbond_sample_r = (real(num_hbond_sample_i,dp)-0.5_dp)*bin_width
             dist = abs(r - num_hbond_sample_r)
!!!!!!Include all the atoms, slow but minimises error
!!!!!!	  if (dist > 4.0_dp*gaussian_sigma) cycle
               exp_arg = -0.5_dp*(dist/(gaussian_sigma))**2
               if (exp_arg > -20.0_dp) then ! good to about 1e-8
                 histogram(num_hbond_sample_i,hbond_type) = histogram(num_hbond_sample_i,hbond_type) + exp(exp_arg)/(gaussian_sigma*sqrt(2.0_dp*PI)) !Gaussian in 1 dimension
               endif
           end do ! num_hbond_sample_i
        else !no gaussian_smoothing
           histogram(int(r/bin_width)+1,hbond_type) = histogram(int(r/bin_width)+1,hbond_type) + 1
        endif
     end do ! i, atom_neighbours

  end do ! H1

  deallocate(mask_a)

end subroutine num_hbond_calc

!
! Calculates the orientation of water molecules along the y axis
!  -- 1st center around atom 1
!  -- then calculate COM of Si atoms
!  -- then shift the centre of mass along the y axis to y=0
!  -- calculate the orientation of the {dipole moment} / {angle half line} of the water: if skip_atoms is set to the last atom of the silica
!
subroutine water_orientation_calc(histogram, at, axis, silica_center_i,n_pos_bins, n_angle_bins, gaussian_smoothing, pos_gaussian_sigma, pos_bin, angle_bin, angle_bin_w, use_dipole_rather_than_angle_bisector, accumulate)
!subroutine water_orientation_calc(histogram, at, axis, silica_center_i,n_pos_bins, n_angle_bins, gaussian_smoothing, pos_gaussian_sigma, angle_gaussian_sigma, pos_bin, angle_bin, angle_bin_w, use_dipole_rather_than_angle_bisector, accumulate)
  real(dp),                          intent(inout) :: histogram(:,:)
  type(Atoms),                       intent(inout) :: at
  integer,                           intent(in)    :: axis
  integer,                           intent(in)    :: silica_center_i
  integer,                           intent(in)    :: n_pos_bins, n_angle_bins
  logical,                           intent(in)    :: gaussian_smoothing
  real(dp),                          intent(in)    :: pos_gaussian_sigma !, angle_gaussian_sigma
  real(dp),                optional, intent(out)   :: pos_bin(:), angle_bin(:)
  real(dp),                optional, intent(inout) :: angle_bin_w(:)
  logical,                 optional, intent(in)    :: use_dipole_rather_than_angle_bisector
  logical,                 optional, intent(in)    :: accumulate

  logical :: my_accumulate
  real(dp) :: sample_r, sample_angle, r
  logical, allocatable :: mask_a(:)
  integer :: sample_i
  real(dp) :: pos_bin_width, angle_bin_width
  real(dp) :: sum_w
integer :: n, num_atoms
integer :: O, H1, H2
  real(dp) :: surface_normal(3)
  logical :: use_dipole
real(dp) :: vector_OH1(3), vector_OH2(3)
real(dp) :: bisector_vector(3)
real(dp) :: dipole(3)
real(dp) :: orientation_angle
    real(dp), parameter                   :: charge_O = -0.834_dp
    real(dp), parameter                   :: charge_H = 0.417_dp
real(dp) :: sum_counts
real(dp) :: dist, exp_arg

  my_accumulate = optional_default(.false., accumulate)
  if (.not. my_accumulate) histogram = 0.0_dp

  use_dipole = optional_default(.true.,use_dipole_rather_than_angle_bisector)
  if (use_dipole) then
     call print("Using dipole to calculate angle with the surface normal.",PRINT_VERBOSE)
  else
     call print("Using HOH angle bisector to calculate angle with the surface normal.",PRINT_VERBOSE)
  endif

  !pos bin labels along axis
  pos_bin_width = at%lattice(axis,axis) / real(n_pos_bins,dp)
  if (present(pos_bin)) then
    do sample_i=1, n_pos_bins
      sample_r = (real(sample_i,dp)-0.5_dp)*pos_bin_width !the middle of the bin
      pos_bin(sample_i) = sample_r
    end do
  endif

 !angle bin labels
  angle_bin_width = PI/n_angle_bins
  if (present(pos_bin)) then
    sum_w = 0._dp
    do sample_i=1, n_angle_bins
      sample_angle = (real(sample_i,dp)-0.5_dp)*angle_bin_width !the middle of the bin
      angle_bin(sample_i) = sample_angle
      !the normalised solid angle is (1/4pi) * 2pi * sin((fi_north)-sin(fi_south)) where fi e [-pi/2,pi/2]
      angle_bin_w(sample_i) = 0.5_dp * ( sin(0.5_dp*PI - real(sample_i-1,dp)*angle_bin_width) - &
                                         sin(0.5_dp*PI - real(sample_i, dp)*angle_bin_width) )
      sum_w = sum_w + angle_bin_w(sample_i)
    end do
    angle_bin_w(1:n_angle_bins) = angle_bin_w(1:n_angle_bins) / sum_w
  endif

  !shift silica slab to the edges, water in the middle    || . . . ||
  if (silica_center_i>0) then ! silica
     call shift_silica_to_edges(at, axis, silica_center_i)
  endif

  !!calc_connect now, before shifting positions to positive, because it would remap the positions!!
  !call calc_connect including the H-bonds
  call set_cutoff(at,0._dp)
  call calc_connect(at)

  !now bins are from -axis/2 to axis/2
  !shift everyone to positive coordinate along axis
  at%pos(axis,1:at%N) = at%pos(axis,1:at%N) + 0.5_dp * at%lattice(axis,axis)

  !simply check water orientations and bin them


!  if (gaussian_smoothing) call system_abort('not implemented.')

  surface_normal(1:3) = 0._dp
  surface_normal(axis) = 1._dp

  num_atoms = 0
  do O=silica_center_i+1, at%N !only check water molecules
     if(at%Z(O)==1) cycle !find O
     !Count the atoms
     call print('Found O'//O//ElementName(at%Z(O)),PRINT_ANAL)
     num_atoms = num_atoms + 1

     !find H neighbours
     n = atoms_n_neighbours(at,O)
     if (n.ne.2) then ! O with =2 nearest neighbours
        call print("WARNING! water(?) oxygen with "//n//"/=2 neighbours will be skipped!")
        cycle
     endif
     H1 = atoms_neighbour(at,O,1)
     H2 = atoms_neighbour(at,O,2)
     if ((at%Z(H1).ne.1).or.(at%Z(H2).ne.1)) then !2 H neighbours
        call print("WARNING! water(?) oxygen with non H neighbour will be skipped!")
        cycle
     endif
 
     !We've found a water molecule.

     !Build histogram
     r = at%pos(axis,H1) !the position of O
     !HOH_angle = angle(diff_min_image(at,O,H1), &
     !                  diff_min_image(at,O,H2)) !the H-O-H angle

     if (.not. use_dipole) then
     !VERSION 1.
     !the direction of the HH->O vector, the bisector of the HOH angle
     !vector that is compared to the surface normal:
     !  point from the bisector of the 2 Hs (scaled to have the same bond length)
     !        to the O
         vector_OH1(1:3) = diff_min_image(at, O, H1)
         if (norm(vector_OH1) > 1.2_dp) &
            call system_abort('too long OH bond? '//O//' '//H1//' '//norm(vector_OH1))
         vector_OH2(1:3) = diff_min_image(at, O, H2)
         if (norm(vector_OH2) > 1.2_dp) &
            call system_abort('too long OH bond? '//O//' '//H2//' '//norm(vector_OH1))
         bisector_vector(1:3) = vector_OH1(1:3) / norm(vector_OH1) * norm(vector_OH2)

         ! a.dot.b = |a|*|b|*cos(angle)
         orientation_angle = dot_product((bisector_vector(1:3)),surface_normal(1:3)) / &
                             sqrt(dot_product(bisector_vector(1:3),bisector_vector(1:3))) / &
                             sqrt(dot_product(surface_normal(1:3),surface_normal(1:3)))
     else ! use_dipole

     !VERSION 2.
     !the dipole of the water molecule = sum(q_i*r_i)
     !Calculate the dipole and its angle compared to the surface normal

         !
         dipole(1:3) = ( diff_min_image(at,O,H1)*charge_H + &
                         diff_min_image(at,O,H2)*charge_H )
         if (norm(diff_min_image(at,O,H1)).gt.1.2_dp) call system_abort('too long O-H1 bond (atoms '//O//'-'//H1//'): '//norm(diff_min_image(at,O,H1)))
         if (norm(diff_min_image(at,O,H2)).gt.1.2_dp) call system_abort('too long O-H2 bond (atoms '//O//'-'//H2//'): '//norm(diff_min_image(at,O,H2)))
!call print ('dipole '//dipole(1:3))
    
         ! a.dot.b = |a|*|b|*cos(angle)
         orientation_angle = dot_product((dipole(1:3)),surface_normal(1:3)) / &
                             sqrt(dot_product(dipole(1:3),dipole(1:3))) / &
                             sqrt(dot_product(surface_normal(1:3),surface_normal(1:3)))
     endif

     if (orientation_angle.gt.1._dp) then
        call print('WARNING | correcting cos(angle) to 1.0 = '//orientation_angle)
        orientation_angle = 1._dp
     else if (orientation_angle.lt.-1._dp) then
        call print('WARNING | correcting cos(angle) to -1.0 = '//orientation_angle)
        orientation_angle = -1._dp
     endif
     orientation_angle = acos(orientation_angle)
     if (orientation_angle.lt.0._dp) then
        call print('WARNING | correcting angle to 0.0: '//orientation_angle)
        orientation_angle = 0._dp
     endif
     if (orientation_angle.gt.PI) then
        call print('WARNING | correcting angle to pi : '//orientation_angle)
        orientation_angle = PI
     endif
!call print ('angle '//(orientation_angle*180._dp/pi))

     call print('Storing angle for water '//O//'--'//H1//'--'//H2//' with reference = '//round(orientation_angle,5)//'degrees',PRINT_ANAL)
     call print('   with distance -1/2 b -- '//O//' = '//round(r,5)//'A',PRINT_ANAL)

     if (gaussian_smoothing) then !smear the position along axis
        !call system_abort('not implemented.')
        do sample_i=1, n_pos_bins
          sample_r = (real(sample_i,dp)-0.5_dp)*pos_bin_width
          dist = abs(r - sample_r)
          !Include all the atoms, slow but minimises error
          !	  if (dist > 4.0_dp*gaussian_sigma) cycle
            exp_arg = -0.5_dp*(dist/(pos_gaussian_sigma))**2
            if (exp_arg > -20.0_dp) then ! good to about 1e-8
              histogram(int(orientation_angle/angle_bin_width)+1,sample_i) = histogram(int(orientation_angle/angle_bin_width)+1,sample_i) + exp(exp_arg)/(pos_gaussian_sigma*sqrt(2.0_dp*PI)) !Gaussian in 1 dimension
            endif
        end do ! sample_i

     else !no gaussian_smoothing
        histogram(int(orientation_angle/angle_bin_width)+1,int(r/pos_bin_width)+1) = histogram(int(orientation_angle/angle_bin_width)+1,int(r/pos_bin_width)+1) + 1._dp
     endif

  end do ! O

  !normalise for the number of molecules in each pos_bin
  do sample_i=1,n_pos_bins
     sum_counts = sum(histogram(1:n_angle_bins,sample_i))
     if (sum_counts /= 0) histogram(1:n_angle_bins,sample_i) = histogram(1:n_angle_bins,sample_i) / sum_counts
  enddo

  !normalise for different solid angles of each angle_bin
  do sample_i=1, n_angle_bins
     histogram(sample_i,1:n_pos_bins) = histogram(sample_i,1:n_pos_bins) / angle_bin_w(sample_i)
  enddo 

end subroutine water_orientation_calc

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
  call param_register(cli_params, "verbose", "F", do_verbose, help_string="No help yet.  This source file was $LastChangedBy$")
  if (.not. param_read_args(cli_params, ignore_unknown=.true.)) &
    call system_abort("Impossible failure to parse verbosity")
  call finalise(cli_params)
  if (do_verbose) then
    call system_initialise(verbosity=PRINT_VERBOSE)
  endif

  call initialise(cli_params)
  call param_register(cli_params, "commandfile", '', commandfilename, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params, "infile", "stdin", infilename, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params, "decimation", "1", decimation, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params, "infile_is_list", "F", infile_is_list, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params, "quiet", "F", quiet, help_string="No help yet.  This source file was $LastChangedBy$")
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
