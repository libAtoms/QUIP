program test_potential
use System_module
use Atoms_module
use Extendable_Str_module
use QUIP_module

implicit none

  type(Potential) pot
  type(MPI_context) mpi_glob
  type(inoutput) params
  type(extendable_str) params_str
  type(Atoms) at

  call system_initialise()

  call Initialise(params, "quip_params.xml")
  call Initialise(params_str)
  call read(params_str, params%unit)

  call Initialise(mpi_glob)

  call read_xyz(at, "atoms_potential_cuau.xyz")
  at%use_uniform_cutoff = .true.

  call print("Atoms")
  call print(at)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call Initialise(pot, "IP LJ", string(params_str), mpi_obj=mpi_glob)
  call Print("post initialise pot LJ")
  call Print(pot)

  at%cutoff = cutoff(pot)
  call calc_connect(at)

  call Print (pot)
  call test_local_e(pot, at)
!NB  call test_grad(pot, at, .true.)
  call test_grad(pot, at, .false.)

  call Finalise(pot)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call finalise(at)
  call read_xyz(at, "atoms_potential_sw.xyz")
  at%use_uniform_cutoff = .true.

  call Initialise(pot, "IP SW", string(params_str), mpi_obj=mpi_glob)
  call Print("post initialise pot SW")
  call Print(pot)

  at%cutoff = cutoff(pot)
  call calc_connect(at)

  call Print (pot)
  call test_local_e(pot, at)
!  call test_grad(pot, at, .true.)
  call test_grad(pot, at, .false.)

  call Finalise(pot)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call finalise(at)
  call read_xyz(at, "atoms_potential_tersoff.xyz")
  at%use_uniform_cutoff = .true.

  call Initialise(pot, "IP Tersoff", string(params_str), mpi_obj=mpi_glob)
  call Print("post initialise pot Tersoff")
  call Print(pot)

  at%cutoff = cutoff(pot)
  call calc_connect(at)

  call Print (pot)
  call test_local_e(pot, at)
!NB  call test_grad(pot, at, .true.)
  call test_grad(pot, at, .false.)

  call Finalise(pot)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call finalise(at)
  call read_xyz(at, "atoms_potential_sw.xyz")
  at%use_uniform_cutoff = .true.

  call Initialise(pot, "TB Bowler", string(params_str), mpi_obj=mpi_glob)
  call Print("post initialise pot Bowler Si")

  at%cutoff = cutoff(pot)
  call calc_connect(at)

!  call Print (pot)
  call test_local_e(pot, at)
  call test_grad(pot, at, .false.)
  call Print(pot%tb%evals)

  call Finalise(pot)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call finalise(at)
  call read_xyz(at, "atoms_potential_sw.xyz")
  at%use_uniform_cutoff = .true.

  call Initialise(pot, "TB DFTB", string(params_str), mpi_obj=mpi_glob)
  call Print("post initialise pot DFTB Si")

  at%cutoff = cutoff(pot)
  call calc_connect(at)

!  call Print (pot)
  call test_local_e(pot, at)
  call test_grad(pot, at, .false.)
  call Print(pot%tb%evals)

  call Finalise(pot)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call finalise(at)
  call read_xyz(at, "atoms_potential_sw.xyz")
  at%use_uniform_cutoff = .true.

  call Initialise(pot, "TB NRL-TB Silicon", string(params_str), mpi_obj=mpi_glob)
  call Print("post initialise pot NRL-TB Si")

  at%cutoff = cutoff(pot)
  call calc_connect(at)

  call Print (pot)
  call test_local_e(pot, at)
  call test_grad(pot, at, .false.)
  call Print(pot%tb%evals)

  call Finalise(pot)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call finalise(at)
  call read_xyz(at, "atoms_potential_sw.xyz")
  at%use_uniform_cutoff = .true.

  call Initialise(pot, "TB DFTB SCF_NONLOCAL_U_DFTB", string(params_str), mpi_obj=mpi_glob)
  call Print("post initialise pot DFTB Si SCF_NONLOCAL_U_DFTB")

  at%cutoff = cutoff(pot)
  call calc_connect(at)

!  call Print (pot)
  call test_local_e(pot, at)
  call test_grad(pot, at, .false.)
  call Print(pot%tb%evals)

  call Finalise(pot)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call finalise(at)
  call read_xyz(at, "atoms_potential_sw.xyz")
  at%use_uniform_cutoff = .true.

  call Initialise(pot, "TB NRL-TB Silicon SCF_NONLOCAL_U_DFTB", string(params_str), mpi_obj=mpi_glob)
  call Print("post initialise pot NRL-TB Si SCF_NONLOCAL_U_DFTB")

  at%cutoff = cutoff(pot)
  call calc_connect(at)

  call Print (pot)
  call test_local_e(pot, at)
  call test_grad(pot, at, .false.)
  call Print(pot%tb%evals)

  call Finalise(pot)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call finalise(at)
  call read_xyz(at, "atoms_potential_cuau.xyz")
  at%use_uniform_cutoff = .true.

  call Initialise(pot, "TB NRL-TB Cu-Au--PW91", string(params_str), mpi_obj=mpi_glob)
  call Print("post initialise pot NRL-TB Ag")

  at%cutoff = cutoff(pot)
  call calc_connect(at)

  call Print (pot)
  call test_local_e(pot, at)
  call test_grad(pot, at, .false.)

  call Finalise(pot)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  call Finalise(at)
  call system_finalise()

end program test_potential

subroutine test_local_e(pot, at)
use System_module
use Atoms_module
use Potential_module
implicit none
  type(Potential), intent(inout) :: pot
  type(Atoms), intent(inout) :: at

  integer i
  real(dp) :: e0, e1
  real(dp), allocatable :: local_e(:)
  real(dp), pointer :: w_e(:)
  logical dummy

  allocate(local_e(at%N))
  allocate(w_e(at%N))

  if (.not. assign_pointer(at, "weight", w_e)) then
    call add_property(at, "weight", 0.0_dp)
    dummy = assign_pointer(at, "weight", w_e) 
  endif

  w_e = 0.0
  w_e(1) = 1.0

  call calc_connect(at)

  call calc(pot, at, e = E0)
  call calc(pot, at, e = E1, local_e = local_e)

  do i=1, at%N
    print '("LOCE local_e ", (i0," ",f20.10,x))', i, local_e(i)
  end do
  print '("LOCE e0 e1 sum(local_e) ", 3(f20.10,x))', e0, e1, sum(local_e)

  call calc(pot, at, e = E0)
  call calc(pot, at, e = E1, local_e = local_e)

  do i=1, at%N
    print '("LOCE w/w_e local_e ", (f20.10,x))', local_e(i)
  end do
  print '("LOCE w/w_e e0 e1 sum(local_e*w_e) ", 3(f20.10,x))', e0, e1, sum(local_e*w_e)

end subroutine

subroutine test_grad(pot, at, do_w_e)
use System_module
use Atoms_module
use Potential_module
implicit none
  type(Potential), intent(inout) :: pot
  type(Atoms), intent(inout) :: at
  logical :: do_w_e

  real(dp), allocatable :: p0(:,:), pert_pos(:,:), dir(:,:)
  real(dp) :: E0, Ep, Em
  real(dp), allocatable :: F0(:,:)

  real(dp) :: lat0(3,3)
  real(dp) :: vir0(3,3)
  real(dp) :: deform_grad(3,3)

  integer ii, jj
  integer fdi, i
  real(dp) fd

  logical dummy
  real(dp), pointer :: w_e(:)

  call calc_connect(at)

  allocate(p0(3,at%N))
  allocate(pert_pos(3,at%N))
  allocate(F0(3,at%N))
  allocate(dir(3,at%N))

  if (do_w_e) then
    if (.not. assign_pointer(at,"weight",w_e)) then
      call add_property(at,"weight",0.0_dp)
      dummy = assign_pointer(at,"weight",w_e)
    endif

    w_e = 0.0_dp
    w_e(1) = 1.0_dp
    w_e(2) = 0.5_dp
  else
    call remove_property(at,"weight")
  endif

  call random_number(pert_pos)
  pert_pos = 0.1_dp*(pert_pos-0.5_dp)

  p0 = at%pos
  lat0 = at%lattice

  at%pos = at%pos + pert_pos
  call calc_dists(at)

  call random_number(dir)
  dir = dir-0.5_dp
  dir = dir/sqrt(sum(dir*dir))


  call calc(pot, at, e = E0, f = F0, virial = vir0)
  print '("GT E0 ", F30.20)', E0

  do i=1, at%N
    print '("F0 ",I3," ",3F20.10)', i, F0(:,i)
  end do

  do fdi=1,8
    fd = 10.0_dp**(-fdi)

    at%pos = p0 + pert_pos - fd*dir
    call calc_dists(at)
    call calc(pot, at, e = Em)

    at%pos = p0 + pert_pos + fd*dir
    call calc_dists(at)
    call calc(pot, at, e = Ep)

    print '("GT dr ",F20.10," fd ",F30.20," analyt ",F30.20, " diff ", F30.20)', fd, (Em-Ep)/(2.0_dp*fd), sum(F0*dir), &
      (Em-Ep)/(2.0_dp*fd)-sum(F0*dir)
  end do

  do i=1, 3
    print '("Virial0 ",3F20.10)', vir0(i,:)
  end do

  ii = 1; jj = 1
  do fdi=1,8
    fd = 10.0_dp**(-fdi)

    deform_grad = 0.0_dp
    call add_identity(deform_grad)
    deform_grad(ii,jj) = deform_grad(ii,jj) - fd
    deform_grad(jj,ii) = deform_grad(jj,ii) - fd
    at%lattice = deform_grad .mult. lat0
    at%pos = deform_grad .mult. (p0 + pert_pos)
    call calc_dists(at)
    call calc(pot, at, e = Em)

    deform_grad = 0.0_dp
    call add_identity(deform_grad)
    deform_grad(ii,jj) = deform_grad(ii,jj) + fd
    deform_grad(jj,ii) = deform_grad(jj,ii) + fd
    at%lattice = deform_grad .mult. lat0
    at%pos = deform_grad .mult. (p0 + pert_pos)
    call calc_dists(at)
    call calc(pot, at, e = Ep)

    print '("VT("I0,",",I0,") dr ",F20.10," fd ",F30.20," analyt ",F30.20," diff ",F30.20)', ii, jj, fd, (Em-Ep)/(4.0_dp*fd), vir0(ii,jj), &
      (Em-Ep)/(4.0_dp*fd) - vir0(ii,jj)

  end do

  ii = 1; jj = 2
  do fdi=1,8
    fd = 10.0_dp**(-fdi)

    deform_grad = 0.0_dp
    call add_identity(deform_grad)
    deform_grad(ii,jj) = deform_grad(ii,jj) - fd
    deform_grad(jj,ii) = deform_grad(jj,ii) - fd
    at%lattice = deform_grad .mult. lat0
    at%pos = deform_grad .mult. (p0 + pert_pos)
    call calc_dists(at)
    call calc(pot, at, e = Em)

    deform_grad = 0.0_dp
    call add_identity(deform_grad)
    deform_grad(ii,jj) = deform_grad(ii,jj) + fd
    deform_grad(jj,ii) = deform_grad(jj,ii) + fd
    at%lattice = deform_grad .mult. lat0
    at%pos = deform_grad .mult. (p0 + pert_pos)
    call calc_dists(at)
    call calc(pot, at, e = Ep)

    print '("VT("I0,",",I0,") dr ",F20.10," fd ",F30.20," analyt ",F30.20," diff ",F30.20)', ii, jj, fd, (Em-Ep)/(4.0_dp*fd), vir0(ii,jj), &
      (Em-Ep)/(4.0_dp*fd) - vir0(ii,jj)

  end do

  at%pos = p0
  at%lattice = lat0
  call calc_dists(at)
end subroutine
