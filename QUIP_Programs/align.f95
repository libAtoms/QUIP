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

! Align a configuration (presumably a molecule) so that moment of inertia
! axes are aligned with cartesian axes
program align_prog
use libatoms_module
implicit none
  type(Atoms) :: at
  type(Dictionary) :: cli_params
  real(dp) :: cutoff_factor

  real(dp) :: CoM(3), MoI(3,3), MoI_evecs(3,3), MoI_evals(3)
  real(dp) :: rot_mat(3,3)
  integer :: i, ii(1)
  real(dp), allocatable :: orig_mass(:)
  real(dp), pointer :: mass(:)
  character(len=2048) :: props

  call system_initialise()

  call initialise(cli_params)
  call param_register(cli_params, 'cutoff_factor', '1.0', cutoff_factor, help_string="cutoff factor for nearest-neighbors to figure out what's a single molecule")
  if (.not. param_read_args(cli_params)) then
    call print("Usage: align [cutoff_factor=1.0]", PRINT_ALWAYS)
    call system_abort("Confused by CLI parameters")
  endif
  call finalise(cli_params)

  call read(at, "stdin")
  props = prop_names_string(at)

  if (.not.(assign_pointer(at,'mass', mass))) then
    call add_property(at,'mass',0.0_dp)
    if (.not.(assign_pointer(at,'mass', mass))) &
      call system_abort("ERROR: Impossible failure to add mass property to atoms")
    mass = 1.0_dp
  else
    allocate(orig_mass(at%N))
    orig_mass = mass
    mass = 1.0_dp
  endif

  call atoms_repoint(at)

  call set_cutoff_factor(at, cutoff_factor)
  call calc_connect(at)

  call coalesce_in_one_periodic_image(at)

  CoM = centre_of_mass(at)
  do i=1, at%N
    at%pos(:,i) = at%pos(:,i) - CoM(:)
  end do

  MoI = moment_of_inertia_tensor(at)
  call diagonalise(MoI, MoI_evals, MoI_evecs)

  do i=1, 3
    ii = maxloc(MoI_evals)
    rot_mat(i,:) = MoI_evecs(:,ii(1))
    MoI_evals(ii(1)) = -1.0e38_dp
  end do

  do i=1, at%N
    at%pos(:,i) = matmul(rot_mat,at%pos(:,i))
  end do

  if (allocated(orig_mass)) then
    mass = orig_mass
    deallocate(orig_mass)
  endif

  at%lattice = 0.0_dp
  at%lattice(1,1) = (maxval(at%pos(1,:)) - minval(at%pos(1,:)))*2.0_dp+5.0_dp
  at%lattice(2,2) = (maxval(at%pos(2,:)) - minval(at%pos(2,:)))*2.0_dp+5.0_dp
  at%lattice(3,3) = (maxval(at%pos(3,:)) - minval(at%pos(3,:)))*2.0_dp+5.0_dp

  ! call print("props" // trim(props))
  call write(at, "stdout", properties=trim(props), prefix="ALIGNED")

  call system_finalise()
end program
