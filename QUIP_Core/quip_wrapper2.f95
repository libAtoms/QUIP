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

! frac_coord: fractional_coordinates
! lattice: LAttice vectors: L(:,1), L(:,2), L(:,3)
! energy: eV
! forces: eV/A
! stress: kbar
subroutine quip_wrapper2(N,Z, frac_coord,lattice,args_str,mpi_communicator,energy,force,stress)

  use libatoms_module
  use potential_module

  implicit none

  integer, intent(in) :: N
  real(dp), dimension(3,3), intent(in) :: lattice
  integer, dimension(N), intent(in) :: Z
  character(len=STRING_LENGTH) :: args_str
  real(dp), dimension(3,N), intent(in) :: frac_coord
  real(dp), intent(out) :: energy
  real(dp), dimension(3,N), intent(out) :: force
  real(dp), dimension(3,3), intent(out) :: stress
  integer, intent(in) :: mpi_communicator 


  type(atoms), save       :: at
  type(Potential), save   :: pot
  type(MPI_context), save :: mpi_glob

  integer :: i

  logical, save :: first_run = .true.

  if( first_run ) then
     call system_initialise(verbosity=PRINT_VERBOSE)
     call Initialise(mpi_glob, communicator=mpi_communicator)
     call Potential_Filename_Initialise(pot, args_str=trim(args_str), param_filename="quip_params.xml",mpi_obj=mpi_glob)
!     call verbosity_push(PRINT_NORMAL)
!     call Print(pot)
!     call verbosity_pop()
     call print(N)
     call print(lattice)
     call initialise(at,N,lattice)
     first_run = .false.
  endif
  
  if( .not. first_run .and. (N /= at%N) ) then
     call finalise(at)
     call initialise(at,N,lattice)
  endif

  call set_lattice(at,lattice, scale_positions=.false.)
  
  do i = 1, at%N
     at%Z(i) = Z(i)
  enddo 
  at%pos = matmul(lattice,frac_coord)

  call set_cutoff(at,cutoff(pot)+0.5_dp)
  call calc_connect(at)

  call calc(pot,at,energy=energy,force=force,virial=stress)
  stress = stress/cell_volume(at)*GPA*10 ! in kbar

endsubroutine quip_wrapper2
