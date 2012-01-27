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

#include "error.inc"

program rings
  use libatoms_module

  implicit none

  real(DP), parameter :: CC_cutoff = 1.80_DP
  integer, parameter  :: MAX_RING_LEN     = 20
  integer, parameter  :: DIST_CUTOFF      = MAX_RING_LEN+1
  integer, parameter  :: MAX_GROUPS       = 20

  type(Dictionary) :: cli_params
  character(len=STRING_LENGTH) :: infilename
  type(CInOutput) :: infile
  type(Atoms) :: at
  integer i
  integer :: error = ERROR_NONE

  integer                :: diameter, o

  integer, allocatable   :: dist(:, :)
  integer                :: ring_stat(MAX_RING_LEN)

  logical, allocatable   :: mask(:)

  call system_initialise(PRINT_NORMAL)

  call initialise(cli_params)
  call param_register(cli_params, 'infile', 'stdin', infilename, help_string="No help yet.  This source file was $LastChangedBy$")
  if (.not. param_read_args(cli_params)) then
    call system_abort('could not parse argument line')
  end if
  call finalise(cli_params)

  call print("infile " // trim(infilename))

  call initialise(infile, infilename, action=INPUT, no_compute_index=.true.)

  i = 1
  call read(infile, at, error=error)
  HANDLE_ERROR(error)

  call set_cutoff(at, CC_cutoff)
  call calc_connect(at)

  allocate(dist(at%N, at%N))

  !
  ! Create the graph of connected C-C
  !

  ring_stat = 0

  allocate(mask(at%N))

  !mask = (at%Z == 6) .and. (at%pos(3, :) > 11.0_DP)
  mask = (at%Z == 6)

  call distance_map(at, dist, mask=mask, diameter=diameter)
  call print("# Graph diameter: " // diameter)
  call count_sp_rings(at, CC_cutoff, dist, MAX_RING_LEN, &
       ring_stat, mask=mask)

  o = sum(ring_stat)
  do i = 1, MAX_RING_LEN
     write (*, '(I5,I10,F10.5)')  i, ring_stat(i), 1.0_DP*ring_stat(i)/o
  enddo

  deallocate(mask)

  deallocate(dist)

  call finalise(infile)
  call system_finalise()

end program rings
