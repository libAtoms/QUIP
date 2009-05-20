program test_GreensFunc
use System_module
use Atoms_module
use QUIP_module
implicit none

  type(TBCalculate) tbc
  type(ApproxFermi) AF
  type(GreensFunc) GF
  type(Atoms) at

  integer i

  type(MPI_context) mpi_glob

  real(dp) energy_exact, energy_exact_no_loc, energy_approx, energy_GF
  real(dp) Fermi_E, band_width
  real(dp) E, emin, emax

  real(dp) energy_approx_plus, energy_approx_minus
  real(dp) energy_GF_plus, energy_GF_minus
  real(dp) f0
  real(dp), allocatable :: p0(:)

  real(dp), allocatable :: local_e_exact(:), local_e_approx(:), local_e_GF(:)
  real(dp), allocatable :: local_n_exact(:), local_n_approx(:), local_n_GF(:)
  real(dp), allocatable :: GF_forces(:,:)
  real(dp), allocatable :: w_e(:), w_n(:)

  ! type(MatrixZ) zS_H, inv_test
  ! complex(dp) :: mone

  type (Inoutput) params

  call system_initialise()

  call Initialise(params, "quip_params.xml")

  call Initialise(mpi_glob)

  call Print ("******************************************** DFTB **************************")
  call Read_xyz(at, "atoms_GF.xyz")
  call Print ("Atoms", verbosity=VERBOSE)
  call Print (at)


  call Print ("Initialize TBSystem DFTB")
  call Initialise(tbc, "DFTB SCF_NONLOCAL_U_DFTB", params, mpi_obj=mpi_glob)
!NB  call Print ("Post initialize")
!NB  call Print(tbc)

  call atoms_set_cutoff(at, tbc%tbsys%tbmodel%cutoff)
  call calc_connect(at)

  call Setup_atoms(tbc, at)

  allocate (w_n(at%N))
  w_n = 0.0_dp
  w_n(1) = 1.0_dp
  w_n(2) = 1.0_dp

  allocate(local_e_exact(at%N))
  allocate(local_n_exact(at%N))
  energy_exact_no_loc = calc_diag(tbc, at, fermi_T = 0.1_dp, w_n = w_n)
  energy_exact = calc_diag(tbc, at, fermi_T = 0.1_dp, local_e = local_e_exact, w_n = w_n, local_n = local_n_exact)

!NB  call Print(tbc%evals, name="tbc%evals")

  Fermi_E = tbc%Fermi_E
  band_width = (Fermi_E-min(tbc%tbsys%kpoints, tbc%evals%data_d(1,:)))*1.2_dp
print '("got Fermi_E, band_width ", 2F20.10)', Fermi_E, band_width
Fermi_E = -4.5_dp
band_width = 12.0_dp

  call Print("Fermi_E " // fermi_E // " band_width " // band_width)

  call Initialise(AF, Fermi_E, 0.1_dp, band_width)
  call Print("got AF%n_poles " // AF%n_poles)

  call Print ("Initialize GF")
  call Initialise(GF, AF%z, AF%a, tbc%tbsys)
  call Print ("Post initialize")
  ! call Print(GF)

  emin = (Fermi_E-band_width)*1.5_dp
  emax = Fermi_E+band_width
!NB  do i=0, 100
!NB    E = emin + dble(i)/100_dp * (emax-emin)
!NB    call Print ("approx_f_fermi " // E // " " // f_fermi(Fermi_E, 0.1_dp, E) // &
!NB      " " // approx_f_fermi(AF, E))
!NB  end do

  ! call Print ("GF%tbsys")
  ! call Print(GF%tbsys)

!NB  call calc_Gs(GF, at)
!NB  call Print ("Post calc_Gs")
  ! call Print(GF)

  ! call Initialise(inv_test, GF%G(1)%N)
  ! call Initialise(zS_H, GF%G(1)%N)
  ! mone = -1.0_dp
  ! call scaled_sum(zS_H, GF%z(1), GF%tbsys%S%data_z(1), mone, GF%tbsys%H%data_z(1))
  ! inv_test%data = matmul(GF%G(1)%data_z(1)%data, zS_H%data)
  ! print *, "GF * (zS-H)"
  ! call Print(inv_test)


  allocate(local_e_approx(at%N))
  allocate(local_n_approx(at%N))

  energy_approx = calc_diag(tbc, at, use_fermi_E = .true., local_e = local_e_approx, w_n = w_n, local_n = local_n_approx, AF = AF)

  allocate (w_e(at%N))
  w_e = 0.0_dp
  w_e(1) = 1.0_dp

  allocate(local_e_GF(at%N))
  allocate(local_n_GF(at%N))
  allocate(GF_forces(3,at%N))
  energy_GF = calc_GF(tbc, at, .true., AF%fermi_E, 0.1_dp, band_width, w_e = w_e, local_e = local_e_GF, &
    w_n = w_n, local_n = local_n_GF, forces = GF_forces)

  call Print ("exact energy no local " // energy_exact_no_loc)
  call Print ("exact energy " // energy_exact)
  call Print (local_e_exact)
  call Print (local_n_exact)

  call Print ("approx energy " // energy_approx)
  call Print(local_e_approx)
  call Print(local_n_approx)

  call print ("GF energy " // energy_GF)
  call print (local_e_GF)
  call print (local_n_GF)

  do i=1, at%N
    call Print("GF force " // GF_forces(:,i))
  end do

  f0 = GF_forces(1,1)
  if (at%N > 6) then
    f0 = f0 + GF_forces(1,7) + GF_forces(1,13) + GF_forces(1,19)
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  allocate(p0(at%N))

  p0 = at%pos(1,:)

  at%pos(1,:) = p0
  at%pos(1,1) = at%pos(1,1) + 1.0e-3_dp
  if (at%N > 6) then
    at%pos(1,7) = at%pos(1,7) + 1.0e-3_dp
    at%pos(1,13) = at%pos(1,13) + 1.0e-3_dp
    at%pos(1,19) = at%pos(1,19) + 1.0e-3_dp
  endif
  call calc_dists(at)
  energy_approx_plus = calc_diag(tbc, at, use_fermi_E = .true., w_e = w_e, w_n = w_n, AF = AF)
  energy_GF_plus = calc_GF(tbc, at, .true., AF%fermi_E, 0.1_dp, band_width, w_e = w_e, w_n = w_n)

  at%pos(1,:) = p0
  at%pos(1,1) = at%pos(1,1) - 1.0e-3_dp
  if (at%N > 6) then
    at%pos(1,7) = at%pos(1,7) - 1.0e-3_dp
    at%pos(1,13) = at%pos(1,13) - 1.0e-3_dp
    at%pos(1,19) = at%pos(1,19) - 1.0e-3_dp
  endif
  call calc_dists(at)
  energy_approx_minus = calc_diag(tbc, at, use_fermi_E = .true., w_e = w_e, w_n = w_n, AF = AF)
  energy_GF_minus = calc_GF(tbc, at, .true., AF%fermi_E, 0.1_dp, band_width, w_e = w_e, w_n = w_n)

  call Print ("approx Ep " // energy_approx_plus // " Em " // energy_approx_minus)
  call Print ("GF Ep " // energy_GF_plus // " Em " // energy_GF_minus)
  call Print ("f0 " // f0 // " FD " // -(energy_GF_plus-energy_GF_minus)/2.0e-3_dp)

  at%pos(1,:) = p0
  at%pos(1,1) = at%pos(1,1) + 1.0e-4_dp
  if (at%N > 6) then
    at%pos(1,7) = at%pos(1,7) + 1.0e-4_dp
    at%pos(1,13) = at%pos(1,13) + 1.0e-4_dp
    at%pos(1,19) = at%pos(1,19) + 1.0e-4_dp
  endif
  call calc_dists(at)
  energy_approx_plus = calc_diag(tbc, at, use_fermi_E = .true., w_e = w_e, w_n = w_n, AF = AF)
  energy_GF_plus = calc_GF(tbc, at, .true., AF%fermi_E, 0.1_dp, band_width, w_e = w_e, w_n = w_n)

  at%pos(1,:) = p0
  at%pos(1,1) = at%pos(1,1) - 1.0e-4_dp
  if (at%N > 6) then
    at%pos(1,7) = at%pos(1,7) - 1.0e-4_dp
    at%pos(1,13) = at%pos(1,13) - 1.0e-4_dp
    at%pos(1,19) = at%pos(1,19) - 1.0e-4_dp
  endif
  call calc_dists(at)
  energy_approx_minus = calc_diag(tbc, at, use_fermi_E = .true., w_e = w_e, w_n = w_n, AF = AF)
  energy_GF_minus = calc_GF(tbc, at, .true., AF%fermi_E, 0.1_dp, band_width, w_e = w_e, w_n = w_n)

  call Print ("approx Ep " // energy_approx_plus // " Em " // energy_approx_minus)
  call Print ("GF Ep " // energy_GF_plus // " Em " // energy_GF_minus)
  call Print ("f0 " // f0 // " FD " // -(energy_GF_plus-energy_GF_minus)/2.0e-4_dp)

end program
