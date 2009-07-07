! move displacement field along a shift vector (defined in reference config) to
! figure out which shifted atom's displacement corresponds to which original atom
program move_displacement_field
use libatoms_module
implicit none
  real(dp) :: shift_vec(3)
  real(dp) :: cutoff
  real(dp), allocatable :: new_pos(:,:)
  type(Atoms) :: config, ref_config
  character(len=FIELD_LENGTH) :: config_filename, ref_config_filename
  type(Inoutput) :: config_io, ref_config_io
  type(Dictionary) :: cli_params
  real(dp) :: shifted_displacement(3)
  real(dp) :: dist
  integer :: i, shifted_i
  integer :: cell_image_Na, cell_image_Nb, cell_image_Nc

  call system_initialise()
  call initialise(cli_params)
  call param_register(cli_params, "config_filename", PARAM_MANDATORY, config_filename)
  call param_register(cli_params, "ref_config_filename", PARAM_MANDATORY, ref_config_filename)
  call param_register(cli_params, "shift_vec", PARAM_MANDATORY, shift_vec)
  call param_register(cli_params, "cutoff", "5.0_dp", cutoff)
  if (.not. param_read_args(cli_params, do_check = .true.)) then
    call print("Usage: move_displacement_field config_filename=file ref_config_filename=file shift_vec={x y z}", ERROR)
    call system_abort("confused by command line arguments")
  endif
  call finalise(cli_params)

  call initialise(config_io, config_filename, action=INPUT)
  call read_xyz(config, config_io)
  call finalise(config_io)

  call initialise(ref_config_io, ref_config_filename, action=INPUT)
  call read_xyz(ref_config, ref_config_io)
  call finalise(ref_config_io)

  if (config%N /= ref_config%N) call system_abort("number of atoms mismatch " // config%N // " " // ref_config%N)

  allocate(new_pos(3,config%N))

  call set_cutoff(ref_config, cutoff)
  call calc_connect(ref_config)

  call fit_box_in_cell(cutoff, cutoff, cutoff, ref_config%lattice, cell_image_Na, cell_image_Nb, cell_image_Nc)
  cell_image_Na = max(1,(cell_image_Na+1)/2)
  cell_image_Nb = max(1,(cell_image_Nb+1)/2)
  cell_image_Nc = max(1,(cell_image_Nc+1)/2)

  do i=1, config%N
    ! atom closest to displaced position
    shifted_i = closest_atom(ref_config, ref_config%pos(:,i)-shift_vec, cell_image_Na, cell_image_Nb, cell_image_Nc, dist=dist)

    if (shifted_i > 0 .and. dist < 0.1_dp) then ! found a close match
      ! current displacement of atom shifted_i
      shifted_displacement(1:3) = config%pos(:,shifted_i) - ref_config%pos(:,shifted_i)
      new_pos(:,i) = ref_config%pos(:,i) + shifted_displacement
    else
      new_pos(:,i) = config%pos(:,shifted_i)
    endif

  end do

  config%pos(1:3,1:config%N) = new_pos(1:3,1:config%N)

  mainlog%prefix="SHIFTED_CONFIG"
  call print_xyz(config, mainlog, all_properties=.true.)
  mainlog%prefix=""

  call system_finalise()
end program
