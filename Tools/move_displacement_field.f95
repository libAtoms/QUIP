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

! move displacement field along a shift vector (defined in reference config) to
! figure out which shifted atom's displacement corresponds to which original atom
program move_displacement_field
use libatoms_module
implicit none
  real(dp) :: shift_vec(3)
  real(dp) :: cutoff, dist_tol
  real(dp), allocatable :: new_pos(:,:)
  type(Atoms) :: config, ref_config
  character(len=FIELD_LENGTH) :: config_filename, ref_config_filename
  type(CInoutput) :: config_io, ref_config_io
  type(Dictionary) :: cli_params
  real(dp) :: cur_displacement(3)
  real(dp) :: dist
  integer :: i, shifted_i
  integer :: cell_image_Na, cell_image_Nb, cell_image_Nc

  call system_initialise()
  call initialise(cli_params)
  call param_register(cli_params, "config_filename", PARAM_MANDATORY, config_filename, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params, "ref_config_filename", PARAM_MANDATORY, ref_config_filename, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params, "shift_vec", PARAM_MANDATORY, shift_vec, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params, "cutoff", "5.0", cutoff, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params, "dist_tol", "0.3", dist_tol, help_string="No help yet.  This source file was $LastChangedBy$")
  if (.not. param_read_args(cli_params)) then
    call print("Usage: move_displacement_field config_filename=file ref_config_filename=file shift_vec={x y z}", PRINT_ALWAYS)
    call print("       cutoff=r(5.0) dist_tol=t(0.3)", PRINT_ALWAYS)
    call system_abort("confused by command line arguments")
  endif
  call finalise(cli_params)

  call print("parameters:")
  call print("config_filename " // trim(config_filename))
  call print("ref_config_filename " // trim(ref_config_filename))
  call print("shift_vec " // shift_vec)
  call print("cutoff " // cutoff)
  call print("dist_tol " // dist_tol)
  call print("")

  call initialise(config_io, config_filename, action=INPUT)
  call read(config, config_io)
  call finalise(config_io)

  call initialise(ref_config_io, ref_config_filename, action=INPUT)
  call read(ref_config, ref_config_io)
  call finalise(ref_config_io)

  if (config%N /= ref_config%N) call system_abort("number of atoms mismatch " // config%N // " " // ref_config%N)

  allocate(new_pos(3,config%N))

  call set_cutoff(ref_config, cutoff)
  call calc_connect(ref_config)

  call fit_box_in_cell(cutoff, cutoff, cutoff, ref_config%lattice, cell_image_Na, cell_image_Nb, cell_image_Nc)
  cell_image_Na = max(1,(cell_image_Na+1)/2)
  cell_image_Nb = max(1,(cell_image_Nb+1)/2)
  cell_image_Nc = max(1,(cell_image_Nc+1)/2)

  new_pos(1:3,1:config%N) = config%pos(1:3,1:config%N)

  do i=1, config%N
    if (mod(i,1000) == 1) write (*,'(a,$)') "."
    ! atom closest to displaced position
    shifted_i = closest_atom(ref_config, ref_config%pos(:,i)+shift_vec, cell_image_Na, cell_image_Nb, cell_image_Nc, dist=dist)

    if (current_verbosity() >= PRINT_VERBOSE) then
      if (shifted_i > 0) then
	call print("mapping displacement cur " // i // " " // config%pos(:,i) // " to " // shifted_i // " " // config%pos(:,shifted_i) // " dist " // dist)
      endif
    endif

    if (shifted_i > 0 .and. dist < dist_tol) then ! found a close match
      ! current displacement of atom i
      cur_displacement(1:3) = config%pos(:,i) - ref_config%pos(:,i)
      ! new position of atom shifted_i from current displacement of atom i
      new_pos(:,shifted_i) = ref_config%pos(:,shifted_i) + cur_displacement
    else
      if (current_verbosity() >= PRINT_VERBOSE) then
	call print("no mapping for " // i // " " // config%pos(:,i))
      endif
    endif

  end do
  call print("")

  config%pos(1:3,1:config%N) = new_pos(1:3,1:config%N)

  call write(config, "stdout", prefix="SHIFTED_CONFIG")

  call system_finalise()
end program
