subroutine gap_wrapper(N,lattice,symbol,coord,energy,force,stress)
  use libatoms_module
  use quip_module

  implicit none

  integer, intent(in) :: N
  real(dp), dimension(3,3), intent(in) :: lattice
  character(len=8), dimension(N), intent(in) :: symbol
  real(dp), dimension(3,N), intent(in) :: coord
  real(dp), intent(out) :: energy
  real(dp), dimension(3,N), intent(out) :: force
  real(dp), dimension(3,3), intent(out) :: stress
  
  type(atoms), save     :: at
  type(Potential), save :: pot

  integer :: i
  real(dp), dimension(:), pointer :: charge

  logical, save :: first_run = .true.

  call system_initialise(verbosity=SILENT)

  if( first_run ) then
     call Initialise(pot, "IP GAP", "" )
     call initialise(at,N,lattice*BOHR)
  endif
  
  if( .not. first_run .and. (N /= at%N) ) then
     call finalise(at)
     call initialise(at,N,lattice*BOHR)
  endif

  call set_lattice(at,lattice*BOHR)
  
  do i = 1, at%N
     at%Z(i) = atomic_number_from_symbol(symbol(i))
  enddo 
  at%pos = coord*BOHR

  call atoms_set_cutoff(at,cutoff(pot)+0.5_dp)
  call calc_connect(at)

  call calc(pot,at,e=energy,f=force,virial=stress)

  energy = energy / HARTREE
  force = force / HARTREE * BOHR
  stress = stress / HARTREE * (BOHR**3)

  first_run = .false.
  call system_finalise()

endsubroutine gap_wrapper
