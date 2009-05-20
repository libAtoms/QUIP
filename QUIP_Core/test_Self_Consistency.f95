program test_Self_Consistency

use System_module
use Atoms_module

use QUIP_module

implicit none

  type(TBCalculate) tbc
  type(Atoms) at
  real(dp) :: energy
  integer i
  real(dp), allocatable :: w_n(:)

  type(MPI_context) mpi_glob

  type (Inoutput) params

  call system_initialise()

  call Initialise(params, "quip_params.xml")

  call Initialise(mpi_glob)

  call Print ("******************************************** NRL-TB **************************")
  call Read_xyz(at, "atoms_self_cons.xyz")
  call Print ("Atoms", verbosity=VERBOSE)
  call Print (at)

  call Print ("Initialise TBCalculate NRL-TB")
  call Initialise(tbc, "NRL-TB Silicon SCF_NONLOCAL_U_DFTB", params, mpi_obj=mpi_glob)
  call Print ("Post initialise")
  call Print(tbc)

  call set_cutoff(at, tbc%tbsys%tbmodel%cutoff)
  call calc_connect(at)

  call Print ("call setup_atoms")
  call setup_atoms(tbc, at)

  call Print ("call calc_diag")
  energy = calc_diag(tbc, at)
  call Print ("post calc_diag")
  write (line, '("got energy ",f20.10)') energy
  call Print (line)

  if (tbc%tbsys%scf%type == SCF_GLOBAL_U .or. &
      tbc%tbsys%scf%type == SCF_GCN) then
    allocate(w_n(at%N))
    w_n = 0.0_dp
    w_n(1) = 1.0_dp
    w_n(2) = 1.0_dp
    energy = calc_diag(tbc, at, w_n = w_n)
    write (line, '("energy 0 ",f30.20)') energy
    call Print(line)
    deallocate(w_n)
  else
    call grad_test(tbc, at, 1d-1)
    call grad_test(tbc, at, 1d-2)
    call grad_test(tbc, at, 1d-3)
    call grad_test(tbc, at, 1d-4)
  endif

  call Print ("Finalise")
  call Finalise(tbc)
  call Finalise(at)

  call Print ("******************************************** Bowler **************************")
  call Read_xyz(at, "atoms_self_cons.xyz")
  call Print ("Atoms", verbosity=VERBOSE)
  call Print (at)

  call rewind(params)
  call Print ("Initialise TBCalculate Bowler")
  call Initialise(tbc, "Bowler SCF_NONLOCAL_U_DFTB", params, mpi_obj=mpi_glob)
  call Print ("Post initialise")
  call Print(tbc)

  call set_cutoff(at, tbc%tbsys%tbmodel%cutoff)
  call calc_connect(at)

  call Print ("call setup_atoms")
  call setup_atoms(tbc, at)

  call Print ("call calc_diag")
  energy = calc_diag(tbc, at)
  call Print ("post calc_diag")
  write (line, '("got energy ",f20.10)') energy
  call Print (line)

  if (tbc%tbsys%scf%type == SCF_GLOBAL_U .or. &
      tbc%tbsys%scf%type == SCF_GCN) then
    allocate(w_n(at%N))
    w_n = 0.0_dp
    w_n(1) = 1.0_dp
    w_n(2) = 1.0_dp
    energy = calc_diag(tbc, at, w_n = w_n)
    write (line, '("energy 0 ",f30.20)') energy
    call Print(line)
    deallocate(w_n)
  else
    call grad_test(tbc, at, 1d-1)
    call grad_test(tbc, at, 1d-2)
    call grad_test(tbc, at, 1d-3)
    call grad_test(tbc, at, 1d-4)
  endif

  call Print ("Finalise")
  call Finalise(tbc)
  call Finalise(at)

  call system_finalise()

end program

subroutine grad_test(tbc, at, dpos)
use System_module
use Atoms_module
use QUIP_module
implicit none
  type(TBCalculate) :: tbc
  type(Atoms) :: at
  real(dp) :: dpos

  real(dp), allocatable :: pert_pos(:,:), forces0(:,:)
  real(dp) p0
  real(dp) energy0, energyp, energym, fd_f

  type (MatrixD) :: dH0, Hp, Hm, fd_dH
  integer fd_ind

integer i, j
  integer :: ml(2)

  allocate(pert_pos(3,at%N))
  allocate(forces0(3,at%N))

  fd_ind = 2

  call random_number(pert_pos)
  pert_pos = 0.1*(pert_pos-0.5D0)
!NB  at%pos = at%pos + pert_pos

  p0 = at%pos(1,fd_ind)

  call calc_dists(at)
  energy0 = calc_diag(tbc, at, forces=forces0)
  write (line, '("energy 0 ",2f30.20)') dpos, energy0
  call Print (line)

  print '("forces0(:,1) ", 3F20.10)', forces0(1:3,1)
!call system_abort ("end")

  call Initialise(dH0, tbc%tbsys%N)
  call copy (dH0, tbc%tbsys%dH(1)%sdata_d(1))

  at%pos(1,fd_ind) = p0+dpos
  call calc_dists(at)
  energyp = calc_diag(tbc, at)
  write (line, '("energy + ",2f30.20)') dpos, energyp
  call Print (line)

  call Initialise(Hp, tbc%tbsys%H%data_d(1)%N)
  Hp%data = tbc%tbsys%H%data_d(1)%data

  at%pos(1,fd_ind) = p0-dpos
  call calc_dists(at)
  energym = calc_diag(tbc, at)
  write (line, '("energy - ",2f30.20)') dpos, energym
  call Print (line)

  call Initialise(Hm, tbc%tbsys%H%data_d(1)%N)
  Hm%data = tbc%tbsys%H%data_d(1)%data

  fd_f = -(energyp-energym)/(2*dpos)
  write (line, '("force analytical ",f30.20," fd ",f30.20, " diff ", f30.20)') forces0(1,fd_ind), fd_f, forces0(1,fd_ind)-fd_f
  call Print(line)

  at%pos(1,fd_ind) = p0

  ! call Initialise(fd_dH, tbc%tbsys%H%data_d(1)%N)
  ! fd_dH%data = (Hp%data-Hm%data)/(2*dpos)

  ! call Print (dH0, name="dH0")
  ! call Print (fd_dH, name="fd_dH")

  ! write (line, '("maxval(dH0-fd_dH) ", f30.20)') maxval(abs(fd_dH%data-dH0%data))
  ! call Print(line)
  ! write (line, '("maxloc(dH0-fd_dH) ", 2(i0,x))') maxloc(abs(fd_dH%data-dH0%data))
  ! call Print(line)

  ! ml = maxloc(abs(fd_dH%data-dH0%data))

  ! print '("dH0 @ maxloc ", 4f20.15)', dH0%data(ml(1), ml(2)), dH0%data(ml(2), ml(1))
  ! print '("fd_dH @ maxloc ", 4f20.15)', fd_dH%data(ml(1), ml(2)), fd_dH%data(ml(2), ml(1))

end subroutine
