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

implicit none

  type(Potential) pot
  type(Atoms):: at, at2
  real(dp)::e
  real(dp):: lat(3,3), v(3,3), strain(3,3)
  real(dp), allocatable::f(:,:)
  real(dp):: pos_delta, lattice_delta, normal_strain_delta, shear_strain_delta, strain_vector(6)
  logical:: lattice_pull, lattice_pull1, lattice_pull2, lattice_pull3, scale_positions, dryrun, no_11, no_12, no_13, no_22, no_23, no_33
  integer :: i
  type(Dictionary) :: cli
  character(len=FIELD_LENGTH) :: infile, outfile, pot_init_args
  integer :: rng_seed, n_configs

  call system_initialise(enable_timing=.true.)

  call initialise(cli)
  call param_register(cli, "pot_init_args", PARAM_MANDATORY, pot_init_args, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli, "infile", "stdin", infile, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli, "outfile", "stdout", outfile, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli, "rng_seed", "-1", rng_seed, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli, "n_configs", "10", n_configs, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli, "positions_delta", "0.2", pos_delta, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli, "lattice_delta", "0.2", lattice_delta, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli, "lattice_pull1", "F", lattice_pull1, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli, "lattice_pull2", "F", lattice_pull2, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli, "lattice_pull3", "F", lattice_pull3, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli, "scale_positions", "T", scale_positions, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli, "dryrun", "F", dryrun, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli, "normal_strain_delta", "0.0", normal_strain_delta, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli, "shear_strain_delta", "0.0", shear_strain_delta, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli, "strain_vector", "0.0 0.0 0.0 0.0 0.0 0.0", strain_vector, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli, "no_11", "F", no_11, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli, "no_12", "F", no_12, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli, "no_13", "F", no_13, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli, "no_22", "F", no_22, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli, "no_23", "F", no_23, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli, "no_33", "F", no_33, help_string="No help yet.  This source file was $LastChangedBy$")

  if (.not. param_read_args(cli, do_check=.true., ignore_unknown=.false.)) &
    call system_abort ("Failed to parse command line arguments")
  call finalise(cli)

  if (rng_seed >= 0) call system_set_random_seeds(rng_seed)

  if(.not. dryrun) call Initialise(pot, trim(pot_init_args))

  call read(at2, trim(infile))
  if(.not. dryrun) allocate(f(3,at2%N))

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

     if ((normal_strain_delta .fne. 0.0_dp) .or. (shear_strain_delta .fne. 0.0_dp) .or. (strain_vector .fne. (/0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp/))) then
        strain = 0.0_dp
        call add_identity(strain)

        if (normal_strain_delta .fne. 0.0_dp) then
           if (.not. no_11) strain(1,1) = strain(1,1) + (ran_uniform()-0.5_dp)*normal_strain_delta
           if (.not. no_22) strain(2,2) = strain(2,2) + (ran_uniform()-0.5_dp)*normal_strain_delta
           if (.not. no_33) strain(3,3) = strain(3,3) + (ran_uniform()-0.5_dp)*normal_strain_delta
        end if

        if (shear_strain_delta .fne. 0.0_dp) then
           if (.not. no_12) strain(1,2) = strain(1,2) + (ran_uniform()-0.5_dp)*shear_strain_delta
           if (.not. no_13) strain(1,3) = strain(1,3) + (ran_uniform()-0.5_dp)*shear_strain_delta
           if (.not. no_23) strain(2,3) = strain(2,3) + (ran_uniform()-0.5_dp)*shear_strain_delta
        end if

        if (strain_vector .fne. (/0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp/)) then
           if (.not. no_11) strain(1,1) = strain(1,1) + strain_vector(1)
           if (.not. no_22) strain(2,2) = strain(2,2) + strain_vector(2)
           if (.not. no_33) strain(3,3) = strain(3,3) + strain_vector(3)
           if (.not. no_12) strain(1,2) = strain(1,2) + strain_vector(4)
           if (.not. no_13) strain(1,3) = strain(1,3) + strain_vector(5)
           if (.not. no_23) strain(2,3) = strain(2,3) + strain_vector(6)
        end if

        call print("configuration: " // i)
        call print("strain matrix:")
        call print(strain)

        call set_lattice(at, strain .mult. at%lattice, scale_positions = .true., remap = .true.)
     end if

     if(.not. dryrun) then
        call system_timer("calc")
        call calc(pot, at, energy=e, force=f, virial=v)
        call system_timer("calc")
        call add_property(at, 'f', f)
        call set_value(at%properties,'Energy', e)
        call set_value(at%properties,'virial', reshape(v, (/9 /)))
        call write(at, outfile, properties='species:pos:f', append=.true.)
     else
        call write(at, outfile, properties='species:pos', append=.true.)
     end if
  end do
  if(.not. dryrun) deallocate(f)

  call finalise(at)
  call finalise(at2)
  call finalise(pot)

  call system_finalise()

end program randomise_calc
