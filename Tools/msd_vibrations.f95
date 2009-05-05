program msd
use libatoms_module
  character(len=FIELD_LENGTH) :: infile_name
  type(inoutput) :: infile
  type(Atoms) :: current, mean, prev
  type(Dictionary) :: cli_params
  integer :: stat
  integer :: i, n_config
  integer :: shift(3)
  real(dp) :: dummy
  real(dp), pointer :: dr(:)

  call system_initialise()

  call initialise(cli_params)
  call param_register(cli_params, 'infile', PARAM_MANDATORY, infile_name)

  call print("n_args " // cmd_arg_count())

  if (.not. param_read_args(cli_params, do_check = .true.)) then
    call print('Usage: msd infile=filename.xyz')
    call system_abort("Bad CLI parameter")
  endif
  call finalise(cli_params)

  call initialise(infile, infile_name)

  call read_xyz(current, infile, status=stat)

  n_config = 0
  do while (stat == 0)
    n_config = n_config + 1
    if (mean%N == 0) then ! first config
      mean = current
      prev = current
! call print("mean_calc initial current")
! call print_xyz(current, mainlog)
    else
      do i=1, current%N
	dummy = distance_min_image(current, current%pos(:,i), prev%pos(:,i), shift=shift)
	if (any(shift /= 0)) then
	  current%pos(:,i) = current%pos(:,i) - (current%lattice .mult. shift)
	endif
      end do
      mean%pos = mean%pos + current%pos
      prev = current
! call print("mean_calc wrapped current")
! call print_xyz(current, mainlog)
    endif
    call read_xyz(current, infile, status=stat)
  end do

  mean%pos = mean%pos/real(n_config,dp)

! call print("mean")
! call print_xyz(mean, mainlog)

  call finalise(infile)

  call add_property(mean, "dr", 0.0_dp)
  if (.not. assign_pointer(mean, "dr", dr)) &
    call system_abort("Impossible failure to assign pointer for dr")
  dr(:) = 0.0_dp

  call initialise(infile, infile_name)

  call read_xyz(current, infile, status=stat)
  prev = current
  do while (stat == 0)
! call print("msd_calc pre-wrap prev")
! call print_xyz(prev, mainlog)
! call print("msd_calc pre-wrap current")
! call print_xyz(current, mainlog)
    do i=1, current%N
      dummy = distance_min_image(current, current%pos(:,i), prev%pos(:,i), shift=shift)
! call print("dummy " // dummy // " shift " // shift)
      if (any(shift /= 0)) then
	current%pos(:,i) = current%pos(:,i) - (current%lattice .mult. shift)
      endif
    end do
    prev = current

! call print("msd_calc wrapped current")
! call print_xyz(current, mainlog)
    do i=1, current%N
      dr(i) = dr(i) + sum((current%pos(:,i)-mean%pos(:,i))**2)
    end do

    call read_xyz(current, infile, status=stat)
  end do

  dr = dr/real(n_config,dp)

  call finalise(infile)

  mainlog%prefix="MSD"
  call print_xyz(mean, mainlog, all_properties=.true.)
  mainlog%prefix=""

  call system_finalise()
end program msd
