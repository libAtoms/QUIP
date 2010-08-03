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

program randomise_calc

use libAtoms_module
use Potential_module
use Potential_module

implicit none

  type(Potential) pot
  type(Atoms):: at, at2
  real(dp)::e
  real(dp):: lat(3,3), v(3,3)
  real(dp), allocatable::f(:,:)
  real(dp):: pos_delta, lattice_delta
  logical:: lattice_pull, lattice_pull1, lattice_pull2, lattice_pull3, scale_positions, dryrun
  integer :: i
  type(Dictionary) :: cli
  character(len=FIELD_LENGTH) :: infile, outfile, pot_init_args
  integer :: rng_seed, n_configs

  call system_initialise(enable_timing=.true.)

  call initialise(cli)
  call param_register(cli, "pot_init_args", PARAM_MANDATORY, pot_init_args)
  call param_register(cli, "infile", "stdin", infile)
  call param_register(cli, "outfile", "stdout", outfile)
  call param_register(cli, "rng_seed", "-1", rng_seed)
  call param_register(cli, "n_configs", "10", n_configs)
  call param_register(cli, "positions_delta", "0.2", pos_delta)
  call param_register(cli, "lattice_delta", "0.2", lattice_delta)
  call param_register(cli, "lattice_pull1", "F", lattice_pull1)
  call param_register(cli, "lattice_pull2", "F", lattice_pull2)
  call param_register(cli, "lattice_pull3", "F", lattice_pull3)
  call param_register(cli, "scale_positions", "T", scale_positions)
  call param_register(cli, "dryrun", "F", dryrun)

  if (.not. param_read_args(cli, do_check=.true., ignore_unknown=.false.)) &
    call system_abort ("Failed to parse command line arguments")
  call finalise(cli)

  if (rng_seed >= 0) call system_set_random_seeds(rng_seed)

  call Initialise(pot, trim(pot_init_args))

  call read(at2, trim(infile))
  allocate(f(3,at2%N))

  lattice_pull = .false.
  if(lattice_pull1 .or. lattice_pull2 .or. lattice_pull3) lattice_pull = .true.

  at=at2
  do i=1,n_configs
     if(.not. lattice_pull) at=at2

     if(pos_delta /= 0.0_dp) call randomise(at%pos, pos_delta)
     if(lattice_delta /= 0.0_dp) then
        if(lattice_pull) then
           lat = at%lattice
           if(lattice_pull1) lat(:,1) = lat(:,1)*(1.0_dp+lattice_delta)
           if(lattice_pull2) lat(:,2) = lat(:,2)*(1.0_dp+lattice_delta)
           if(lattice_pull3) lat(:,3) = lat(:,3)*(1.0_dp+lattice_delta)
        else
           lat = at%lattice
           call randomise(lat, lattice_delta)
        end if

        call set_lattice(at, lat, scale_positions=scale_positions)
     end if
     if(.not. dryrun) then
        call system_timer("calc")
        call calc(pot, at, e=e, f=f, virial=v)
        call system_timer("calc")
     end if
     call add_property(at, 'f', f)
     call set_value(at%properties,'Energy', e)
     call set_value(at%properties,'virial', reshape(v, (/9 /)))
     call write(at, outfile, properties='pos:f', append=.true.)
  end do
  deallocate(f)

  call system_finalise()

end program randomise_calc
