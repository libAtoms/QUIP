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
module analyze_md_phonons_mod
use libatoms_module
implicit none

contains

subroutine phonon_evec_to_displacements(phonon, Z)
  real(dp) :: phonon(:,:)
  integer :: Z(:)

  integer i

  do i=1, size(phonon,2)
    phonon(:,i) = phonon(:,i) / sqrt(ElementMass(Z(i)))
  end do
end subroutine
end module analyze_md_phonons_mod

program analyze_md_phonons
use libatoms_module
use analyze_md_phonons_mod
implicit none
  type(atoms) :: config
  integer n_phonons
  real(dp), allocatable :: phonons(:,:,:)
  type(CInOutput) :: phonons_io, md_io
  real(dp), pointer :: phonon_displ(:,:)
  integer i, config_i, j
  logical done
  integer error
  real(dp), allocatable :: ke_proj(:)
  real(dp) :: vel_mag, ke_tot, p(3), L(3), mode_vel
  real(dp) :: MoI(3,3), MoI_evecs(3,3), MoI_evals(3), MoI_evecs_orig(3,3), R(3,3)
  character(len=FIELD_LENGTH) phonons_file
  logical :: fix_rotation, regular_eigenproblem
  type(dictionary) :: cli_params


  call system_initialise()

  call initialise(cli_params)
    call param_register(cli_params, "phonons_file", "phonons.xyz", phonons_file, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(cli_params, "fix_rotation", "F", fix_rotation, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(cli_params, "regular_eigenproblem", "F", regular_eigenproblem, help_string="No help yet.  This source file was $LastChangedBy$")
    if (.not. param_read_args(cli_params)) then
      call print("Usage: analyze_md_phonons [phonons_file=file(phonons.xyz)] [fix_rotation=T/F(F)]", PRINT_ALWAYS)
      call print("       [regular_eigenproblem=T/F(F)]", PRINT_ALWAYS)
      call system_abort("Confused by parameters")
    endif
  call finalise(cli_params)

  call print("phonons_file '" // trim(phonons_file) // "'")
  call print("fix_rotation " // fix_rotation)
  call initialise(phonons_io, phonons_file)

  call read(config, phonons_io)
  n_phonons = 3*config%N
  allocate(phonons(3, config%N, n_phonons))
  if (.not. assign_pointer(config, "phonon", phonon_displ)) &
    call system_abort("Couldn't find field phonon in phonon config " // 1)
  if (regular_eigenproblem) call phonon_evec_to_displacements(phonon_displ, config%Z)
  phonons(:,:,1) = phonon_displ
  call print("phonon norm 1 " // sum(ElementMass(config%Z)*sum(phonons(:,:,1)**2,1)))

  do i=2, n_phonons
    call read(config, phonons_io)
    if (.not. assign_pointer(config, "phonon", phonon_displ)) &
      call system_abort("Couldn't find field phonon in phonon config " // i)
    if (regular_eigenproblem) call phonon_evec_to_displacements(phonon_displ, config%Z)
    phonons(:,:,i) = phonon_displ
    call print("phonon norm " // i // " " // sum(ElementMass(config%Z)*sum(phonons(:,:,i)**2,1)))
  end do

  call finalise(phonons_io)

  call initialise(md_io, "stdin")

  allocate(ke_proj(n_phonons))

  call print("Masses " // ElementMass(config%Z))

call add_property(config, "mass", 0.0_dp, 1)
call atoms_repoint(config)
config%mass = ElementMass(config%Z)
MoI = moment_of_inertia_tensor(config, centre_of_mass(config))
call diagonalise(MoI, MoI_evals, MoI_evecs_orig)
call print(MoI_evecs_orig, PRINT_VERBOSE)

  done = .false.
  config_i = 0
  do while (.not. done)
     config_i = config_i + 1
     call read(config, md_io, error=error)
     if (error /= 0) then
        if (error == ERROR_IO_EOF) then
	   done = .true.
	else
	   HANDLE_ERROR(error)
	endif
     else
      call add_property(config, "mass", 0.0_dp, 1)
      call atoms_repoint(config)
      config%mass = ElementMass(config%Z)
      if (.not. associated(config%velo)) &
	call system_abort("config " // config_i // " has no velo pointer")
      ke_tot = kinetic_energy(config)
      p = momentum(config)
      L = angular_momentum(config, centre_of_mass(config))
      MoI = moment_of_inertia_tensor(config, centre_of_mass(config))

      if (fix_rotation) then
	call diagonalise(MoI, MoI_evals, MoI_evecs)
	do i=1,3
	if ((MoI_evecs(:,i) .dot. MoI_evecs_orig(:,i)) < 0) MoI_evecs(:,i) = - MoI_evecs(:,i)
	end do
	call print(MoI_evecs, PRINT_VERBOSE)
	R = transpose(MoI_evecs_orig) .mult. MoI_evecs
	config%velo = R .mult. config%velo
      endif

      vel_mag = sqrt(sum(config%velo**2))
      if (vel_mag .feq. 0.0_dp) then
	ke_proj = 0.0_dp
      else
	do i=1, n_phonons
	  mode_vel = 0.0_dp
	  do j=1, config%N
	    mode_vel = mode_vel + ElementMass(config%Z(j))*sum(phonons(:,j,i)*config%velo(:,j))
	  end do
	  ke_proj(i) = mode_vel**2/2.0_dp
	end do
      endif
      call print("check " // config_i // " p " // p // " L ", PRINT_VERBOSE)
      call print("check ke_tot " // ke_tot // " sum(ke_proj) " // sum(ke_proj) // &
	" sum(ke_proj(1:6)) " // sum(ke_proj(1:6)) // &
	" sum(ke_proj(7:)) " // sum(ke_proj(7:)) // " ratio " // &
	(sum(ke_proj(7:))/sum(ke_proj)), PRINT_VERBOSE)
      call print("config " // config_i // " ke_proj " // ke_proj)
    endif
  end do

  call system_finalise()

end program analyze_md_phonons

