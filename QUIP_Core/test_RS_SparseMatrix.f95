program test_RS_SparseMAtrix

use System_module
use Matrix_module
use RS_SparseMatrix_module

implicit none

  type(RS_SparseMatrixD) sm, smsm, sm2
  type(MatrixD) dm, dmdm, dmsm, dm2
  type(Atoms) at

  integer i
  integer, allocatable :: first_orb_of_atom(:)

  type (Inoutput) params

  call system_initialise()

  call Initialise(params, "quip_params.xml")

  call Read_xyz(at, "atoms_rs_sparsematrix.xyz")
  call Print ("Atoms", verbosity=VERBOSE)
  call Print (at)

  call set_cutoff(at,3.0_dp)
  call calcConnectFast(at, rebuild=.true.)
  call Print(at)


  allocate(first_orb_of_atom(at%N+1))
  first_orb_of_atom(1) = 1
  do i=2, at%N+1
    if (i .eq. at%N+1) then
      first_orb_of_atom(i) = first_orb_of_atom(i-1) + 1
    else
      first_orb_of_atom(i) = first_orb_of_atom(i-1) + 4
    endif
  end do

  call Print ("call initialise")
  call Initialise(sm, at, first_orb_of_atom)
  sm%data = 1.0_dp
  call Print(sm, name="sm")

  print *, "calling Initialise(dm)", first_orb_of_atom(at%N+1)-1
  call Initialise(dm, first_orb_of_atom(at%N+1)-1)
  call Print(dm)
  call copy(dm, sm)
  call Print(dm, name="dm=sm")

  call Initialise(dmdm, dm%N)
  call Initialise(dmsm, dm%N)
  call matrix_product(dmdm, dm, dm)
  call matrix_product(dmsm, dm, sm)

  write (line, '("max dm*dm - dm*sm ",f20.10)'), maxval(dmdm%data-dmsm%data)
  call Print(line)

  call Print(dmdm, name="dm*dm")

  call Initialise(smsm, at, first_orb_of_atom)
  call matrix_product(smsm, sm, sm)
  call Print(smsm, name="sm*sm with structure of sm")

  call Initialise(sm2, sm%l, sm%l)
  call matrix_product(sm2, sm, sm)
  call Print(sm2, name="sm*sm with proper structure")

  call Initialise(dm2, dm%N)
  call copy(dm2, sm2)
  write (line, '("max dm*dm - sm*sm ",f20.10)'), maxval(dmdm%data-dmsm%data)
  call Print(line)


end program
