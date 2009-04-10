program test_TBModel

use System_module
use MPI_context_module
use TBSystem_module

implicit none

  type(TBSystem) tbs
  type(Atoms) at
  real(dp) :: local_E
  integer i

  type(MPI_context) mpi_glob

  type (Inoutput) params

  call system_initialise()

  call Initialise(params, "quip_params.xml")

  call Initialise(mpi_glob)

  call Print ("******************************************** NRL-TB **************************")
  call Read_xyz(at, "atoms_nrl_tb.xyz")
  call Print ("Atoms", verbosity=VERBOSE)
  call Print (at)

  call Print ("Initialize TBSystem NRL-TB")
  call Initialise(tbs, "NRL-TB", params, mpi_obj=mpi_glob)
  call Print ("Post initialize")
  call Print(tbs)

  at%cutoff = tbs%tbmodel%cutoff
  at%use_uniform_cutoff = .true.
  write (line, '(a,f10.5)') "at%cutoff ", at%cutoff
  call Print(line)
  call calc_connect(at)

  call Print ("call Setup_atoms")
  call Setup_atoms(tbs, at)
  call Print ("Post setup_atoms")
  call Print(tbs)

  call Print ("call fill_matrices")
  call fill_matrices(tbs, at)
  call Print ("Post fill_matrices")
  call Print(tbs)

  call Print ("local_E")
  do i=1, at%N
    local_E = get_local_rep_E(tbs%tbmodel, at, i)
    call Print ('local_E ' // local_E)
  end do

  call Print ("Finalise")
  call Finalise(tbs)
  call Finalise(at)

  call rewind(params)

  call Print ("******************************************** Bowler **************************")
  call Read_xyz(at, "atoms_bowler.xyz")
  call Print ("Atoms", verbosity=VERBOSE)
  call Print (at)

  call Print ("Initialize TBSystem Bowler")
  call Initialise(tbs, "Bowler", params, mpi_obj=mpi_glob)
  call Print ("Post initialize")
  call Print(tbs)

  at%cutoff = tbs%tbmodel%cutoff
  at%use_uniform_cutoff = .true.
  write (line, '(a,f10.5)') "at%cutoff ", at%cutoff
  call Print(line)
  call calc_connect(at)

  call Print ("call Setup_atoms")
  call Setup_atoms(tbs, at)
  call Print ("Post setup_atoms")
  call Print(tbs)

  call Print ("call fill_matrices")
  call fill_matrices(tbs, at)
  call Print ("Post fill_matrices")
  call Print(tbs)

  call Print ("local_E")
  do i=1, at%N
    local_E = get_local_rep_E(tbs%tbmodel, at, i)
    call Print ('local_E ' // local_E)
  end do

  call Print ("Finalise")
  call Finalise(tbs)
  call Finalise(at)

  call rewind(params)

  call Print ("******************************************** DFTB **************************")
  call Read_xyz(at, "atoms_dftb.xyz")
  call Print ("Atoms", verbosity=VERBOSE)
  call Print (at)

  call Print ("Initialize TBSystem DFTB")
  call Initialise(tbs, "DFTB", params, mpi_obj=mpi_glob)
  call Print ("Post initialize")
  call Print(tbs)

  at%cutoff = tbs%tbmodel%cutoff
  at%use_uniform_cutoff = .true.
  write (line, '(a,f10.5)') "at%cutoff ", at%cutoff
  call Print(line)
  call calc_connect(at)

  call Print ("call Setup_atoms")
  call Setup_atoms(tbs, at)
  call Print ("Post setup_atoms")
  call Print(tbs)

  call Print ("call fill_matrices")
  call fill_matrices(tbs, at)
  call Print ("Post fill_matrices")
  call Print(tbs)

  call Print ("local_E")
  do i=1, at%N
    local_E = get_local_rep_E(tbs%tbmodel, at, i)
    call Print ('local_E ' // local_E)
  end do

  call Print ("Finalise")
  call Finalise(tbs)
  call Finalise(at)


  call system_finalise()

end program
