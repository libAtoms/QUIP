program test_DFTB

use System_module
use MPI_context_module
use Atoms_module
use TBCalculate_module

implicit none

  type(TBCalculate) tbc
  type(Atoms) at
  real(dp) :: energy
  integer i

  type(MPI_context) mpi_glob

  type (Inoutput) params

  call system_initialise()

  call Initialise(params, "quip_params.xml")

  call Initialise(mpi_glob)

  call Print ("******************************************** DFTB **************************")
  call Read_xyz(at, "atoms_dftb.xyz")
  call Print ("Atoms", verbosity=VERBOSE)
  call Print (at)

  call Print ("Initialise TBCalculate DFTB")
  call Initialise(tbc, "DFTB", params, mpi_obj=mpi_glob)
  call Print ("Post initialise")
  ! call Print(tbc)

  at%cutoff = tbc%tbsys%tbmodel%cutoff
  at%use_uniform_cutoff = .true.
  call calcConnectFast(at)

  call Print ("call setup_atoms")
  call setup_atoms(tbc, at)

  call Print ("call solve_diag")
  call solve_diag(tbc, at)

  call Print ("Post solve_diag evals")
  call Print(tbc%evals)

  call Print ("call calc_diag")
  energy = calc_diag(tbc, at)
  call Print ("post calc_diag")
  write (line, '("got energy ",f20.10)') energy
  call Print (line)

  call grad_test(tbc, at)

  call Print ("Finalise")
  call Finalise(tbc)
  call Finalise(at)

  call system_finalise()

end program

subroutine grad_test(tbc, at)
use System_module
use Atoms_module
use TBCalculate_module
implicit none
  type(TBCalculate) :: tbc
  type(Atoms) :: at

  real(dp), allocatable :: pert_pos(:,:), forces0(:,:)
  real(dp) p0, dpos
  real(dp) energy0, energyp, energym, fd_f

  type (MatrixD) :: dH0, Hp, Hm, fd_dH

  integer :: ml(2)

  allocate(pert_pos(3,at%N))
  allocate(forces0(3,at%N))

  call random_number(pert_pos)
  pert_pos = 0.1*(pert_pos-0.5D0)
  at%pos = at%pos + pert_pos

  p0 = at%pos(1,3)
  dpos = 0.0001_dp

  call calcDists(at)
  call unfreshen_everything(tbc)
  energy0 = calc_diag(tbc, at, forces=forces0)
  write (line, '("energy 0 ",2f30.20)') dpos, energy0
  call Print (line)

  call Initialise(dH0, tbc%tbsys%N)
  call copy (dH0, tbc%tbsys%dH(1)%sdata_d(1))

  at%pos(1,3) = p0+dpos
  call calcDists(at)
  call unfreshen_everything(tbc)
  energyp = calc_diag(tbc, at)
  write (line, '("energy + ",2f30.20)') dpos, energyp
  call Print (line)

  call Initialise(Hp, tbc%tbsys%H%data_d(1)%N)
  Hp%data = tbc%tbsys%H%data_d(1)%data

  at%pos(1,3) = p0-dpos
  call calcDists(at)
  call unfreshen_everything(tbc)
  energym = calc_diag(tbc, at)
  write (line, '("energy - ",2f30.20)') dpos, energym
  call Print (line)

  call Initialise(Hm, tbc%tbsys%H%data_d(1)%N)
  Hm%data = tbc%tbsys%H%data_d(1)%data

  fd_f = -(energyp-energym)/(2*dpos)
  write (line, '("force analytical ",f30.20," fd ",f30.20)') forces0(1,3), fd_f
  call Print(line)

  call Initialise(fd_dH, tbc%tbsys%H%data_d(1)%N)
  fd_dH%data = (Hp%data-Hm%data)/(2*dpos)

  call Print (dH0, name="dH0")
  call Print (fd_dH, name="fd_dH")

  write (line, '("maxval(dH0-fd_dH) ", f30.20)') maxval(abs(fd_dH%data-dH0%data))
  call Print(line)
  write (line, '("maxloc(dH0-fd_dH) ", 2(i0,x))') maxloc(abs(fd_dH%data-dH0%data))
  call Print(line)

  ml = maxloc(abs(fd_dH%data-dH0%data))

  print '("dH0 maxloc ", 4f20.15)', dH0%data(ml(1), ml(2)), dH0%data(ml(2), ml(1))
  print '("fd_dH maxloc ", 4f20.15)', fd_dH%data(ml(1), ml(2)), fd_dH%data(ml(2), ml(1))

end subroutine
