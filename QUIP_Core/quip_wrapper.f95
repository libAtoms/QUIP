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

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X quip_wrapper subroutine
!X
!% wrapper to make it nicer for non-QUIP programs to use a QUIP potential
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

subroutine quip_wrapper(N,lattice,symbol,coord,args_str,args_str_len,energy,force,virial,do_energy,do_force,do_virial)

  use system_module, only : dp, print, system_initialise, system_abort, PRINT_NORMAL, PRINT_SILENT, verbosity_push, verbosity_pop
  use dictionary_module
  use periodictable_module
  use mpi_context_module
  use atoms_types_module, only : assign_pointer
  use atoms_module

  use potential_module

  implicit none

  integer, intent(in) :: N
  real(dp), dimension(3,3), intent(inout) :: lattice
  character(len=3), dimension(N), intent(in) :: symbol
  integer, intent(in) :: args_str_len
  character(len=args_str_len) :: args_str
  real(dp), dimension(3,N), intent(in) :: coord
  real(dp), intent(out) :: energy
  real(dp), dimension(3,N), intent(out) :: force
  real(dp), dimension(3,3), intent(out) :: virial
  logical, intent(in) :: do_energy, do_force, do_virial
  
  type(atoms), save       :: at
  type(Potential), save   :: pot
  type(MPI_context), save :: mpi_glob

  integer :: i
  real(dp), dimension(:,:), pointer :: quip_wrapper_force
  real(dp) :: cut, clustlenx, clustleny, clustlenz, boxlenx, boxleny, boxlenz
  real(dp), dimension(3) :: upper, lower
  character(len=STRING_LENGTH) :: calc_args

  logical, save :: first_run = .true.

  if( first_run ) then
     call system_initialise(verbosity=PRINT_SILENT)
     call Initialise(mpi_glob)
     call Potential_Filename_Initialise(pot, args_str=trim(args_str), param_filename="quip_params.xml",mpi_obj=mpi_glob)
     call verbosity_push(PRINT_NORMAL)
     call Print(pot)
     call verbosity_pop()
     call initialise(at,N,lattice)
  endif
  
  if( .not. first_run .and. (N /= at%N) ) then
     call finalise(at)
     call initialise(at,N,lattice)
  endif

  if(lattice(1,1) == 0.0_dp) then
     lower = minval(coord, dim=2)
     upper = maxval(coord, dim=2)
     cut = cutoff(pot)

     clustlenx = abs(upper(1) - lower(1))
     clustleny = abs(upper(2) - lower(2))
     clustlenz = abs(upper(3) - lower(3))

     boxlenx = clustlenx + 2*cut + 1.0_dp
     boxleny = clustleny + 2*cut + 1.0_dp
     boxlenz = clustlenz + 2*cut + 1.0_dp

     lattice = reshape((/ boxlenx, 0.0_dp, 0.0_dp, 0.0_dp, boxleny, 0.0_dp, 0.0_dp, 0.0_dp, boxlenz /), shape(lattice))
  endif

  call set_lattice(at,lattice, scale_positions=.false.)
  
  do i = 1, at%N
     at%Z(i) = atomic_number_from_symbol(symbol(i))
  enddo 
  at%pos = coord

  call set_cutoff(at,cutoff(pot)+0.5_dp)
  call calc_connect(at)

  calc_args = ""
  if(do_energy) calc_args = trim(calc_args)//" energy=quip_wrapper_energy "
  if(do_force) calc_args = trim(calc_args)//" force=quip_wrapper_force "
  if(do_virial) calc_args = trim(calc_args)//" virial=quip_wrapper_virial "

  call calc(pot,at,args_str=trim(calc_args))

  if(do_energy) call get_param_value(at, "quip_wrapper_energy", energy)
  if(do_virial) call get_param_value(at, "quip_wrapper_virial", virial)
  if(do_force) then
     if(.not. assign_pointer(at,"quip_wrapper_force",quip_wrapper_force) ) call system_abort("Could not calculate forces")
     force = quip_wrapper_force
  endif

  first_run = .false.

endsubroutine quip_wrapper
