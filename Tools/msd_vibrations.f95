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
  call param_register(cli_params, 'infile', PARAM_MANDATORY, infile_name, help_string="No help yet.  This source file was $LastChangedBy$")

  call print("n_args " // cmd_arg_count())

  if (.not. param_read_args(cli_params)) then
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
