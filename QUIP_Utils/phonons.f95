! IR intensities from K. Jackson, M. R. Pederson, and D. Porezag,
! Z. Hajnal, and T. Frauenheim, Phys. Rev. B v. 55, 2549 (1997).
module phonons_module
use libatoms_module
use quip_module
use libatoms_misc_utils_module
implicit none
private

public phonons, eval_frozen_phonon

contains

function eval_frozen_phonon(metapot, at, dx, evec, calc_args)
  type(MetaPotential), intent(inout) :: metapot
  type(Atoms), intent(inout) :: at
  real(dp), intent(in) :: dx
  real(dp), intent(in) :: evec(:)
  character(len=*), intent(in), optional :: calc_args
  real(dp) :: eval_frozen_phonon ! result

  real(dp) :: Ep, E0, Em
  real(dp), allocatable :: pos0(:,:), dpos(:,:)

  allocate(pos0(3,at%N))
  allocate(dpos(3,at%N))

  pos0 = at%pos

  dpos = reshape(evec, (/ 3, at%N /) )
  dpos = dpos / sqrt(sum(dpos**2))

  call calc_dists(at)
  call calc(metapot, at, E=E0, args_str=calc_args)
  mainlog%prefix="FROZ_E0"
  call print_xyz(at, mainlog, real_format='f14.10', properties="pos:phonon")
  mainlog%prefix=""

  at%pos = pos0 + dx*dpos

  call calc_dists(at)
  call calc(metapot, at, E=Ep, args_str=calc_args)
  mainlog%prefix="FROZ_EP"
  call print_xyz(at, mainlog, real_format='f14.10', properties="pos:phonon")
  mainlog%prefix=""

  at%pos = pos0 - dx*dpos

  call calc_dists(at)
  call calc(metapot, at, E=Em, args_str=calc_args)
  mainlog%prefix="FROZ_EM"
  call print_xyz(at, mainlog, real_format='f14.10', properties="pos:phonon")
  mainlog%prefix=""

  call print("Em " // Em // " E0 " // E0 // " Ep " // Ep, VERBOSE)

  eval_frozen_phonon = ((Ep-E0)/dx - (E0-Em)/dx)/dx

  at%pos = pos0
  call calc_dists(at)

  deallocate(dpos)
  deallocate(pos0)

end function eval_frozen_phonon

subroutine phonons(metapot, at, dx, evals, evecs, effective_masses, calc_args, IR_intensities, do_parallel, &
		   zero_translation, zero_rotation, force_const_mat)
  type(MetaPotential), intent(inout) :: metapot
  type(Atoms), intent(inout) :: at
  real(dp), intent(in) :: dx
  real(dp), intent(out) :: evals(at%N*3)
  real(dp), intent(out), optional :: evecs(at%N*3,at%N*3)
  real(dp), intent(out), optional :: effective_masses(at%N*3)
  character(len=*), intent(in), optional :: calc_args
  real(dp), intent(out), optional :: IR_intensities(:)
  logical, intent(in), optional :: do_parallel
  logical, intent(in), optional :: zero_translation, zero_rotation
  real(dp), intent(out), optional :: force_const_mat(:,:)

  integer i, j, alpha, beta
  integer err
  real(dp), allocatable :: pos0(:,:), fp(:,:), fm(:,:)
  real(dp), allocatable :: dm(:,:)

  logical :: override_zero_freq_phonons = .true.
  real(dp) :: phonon_norm, CoM(3), axis(3), dr_proj(3)
  real(dp), allocatable :: zero_overlap_inv(:,:)
  real(dp), allocatable :: phonon(:,:), zero_phonon(:,:), zero_phonon_p(:,:), P(:,:), dm_t(:,:)
  real(dp) :: mu_m(3), mu_p(3), dmu_dq(3)
  real(dp), allocatable :: dmu_dr(:,:,:)
  real(dp), pointer :: local_dn(:), mass(:)

  integer :: n_zero
  logical :: do_zero_translation, do_zero_rotation

  integer :: ind
  logical :: my_do_parallel
  type(MPI_context) :: mpi_glob

  my_do_parallel = optional_default(.false., do_parallel)

  if (my_do_parallel) then
    call initialise(mpi_glob)
  endif

  do_zero_translation = optional_default(zero_translation, .true.)
  do_zero_rotation = optional_default(zero_rotation, .false.)

  allocate(pos0(3,at%N))
  allocate(fp(3,at%N))
  allocate(fm(3,at%N))
  allocate(dm(at%N*3,at%N*3))
  allocate(dmu_dr(3, at%N, 3))

  if (present(force_const_mat)) then
    if (size(force_const_mat,1) /= size(dm,1) .or. size(force_const_mat,2) /= size(dm,2)) &
      call system_abort("phonons received force_const_mat, shape="//shape(force_const_mat) // &
			" which doesn't match shape(dm)="//shape(dm))
  endif

  if (present(IR_intensities)) then
    if (.not. assign_pointer(at, 'local_dn', local_dn)) then
      call add_property(at, 'local_dn', 0.0_dp, 1)
    endif
  endif

  pos0 = at%pos

  if (dx == 0.0_dp) &
    call system_abort("phonons called with dx == 0.0")

  if (my_do_parallel) then
    dm = 0.0_dp
    dmu_dr = 0.0_dp
  endif

  call set_cutoff(at, cutoff(metapot))
  call calc_connect(at)
  ! calculate dynamical matrix with finite differences
  ind = -1
  do i=1, at%N
    do alpha=1,3
      ind = ind + 1
      if (my_do_parallel) then
	if (mod(ind, mpi_glob%n_procs) /= mpi_glob%my_proc) cycle
      endif

      at%pos = pos0
      at%pos(alpha,i) = at%pos(alpha,i) + dx
      call calc_dists(at)
      call calc(metapot, at, f=fp, args_str=calc_args)
      if (present(IR_intensities)) then
	if (.not. assign_pointer(at, 'local_dn', local_dn)) &
	  call system_abort("phonons impossible failure to assign pointer for local_dn")
	mu_p = dipole_moment(at%pos, local_dn)
      endif

      at%pos = pos0
      at%pos(alpha,i) = at%pos(alpha,i) - dx
      call calc_dists(at)
      call calc(metapot, at, f=fm, args_str=calc_args)
      if (present(IR_intensities)) then
	if (.not. assign_pointer(at, 'local_dn', local_dn)) &
	  call system_abort("phonons impossible failure to assign pointer for local_dn")
	mu_m = dipole_moment(at%pos, local_dn)
      endif

      dmu_dr(alpha, i, :) = (mu_p-mu_m)/(2.0_dp*dx)

      do j=1, at%N
	do beta=1,3
	  dm((i-1)*3+alpha,(j-1)*3+beta) = -((fp(beta,j)-fm(beta,j))/(2.0_dp*dx))
	end do
      end do

    end do
  end do

  at%pos = pos0
  call calc_dists(at)

  if (my_do_parallel) then
    call sum_in_place(mpi_glob, dm)
    call sum_in_place(mpi_glob, dmu_dr)
  endif

  if (.not. assign_pointer(at, 'mass', mass)) &
    call add_property(at, 'mass', 0.0_dp, 1)
  if (.not. assign_pointer(at, 'mass', mass)) &
    call system_abort("impossible failure to assign pointer for mass")
  do i=1, at%N
    mass(i) = ElementMass(at%Z(i))
  end do

  if (present(force_const_mat)) then
    force_const_mat = dm
  endif

  ! transform from generalized eigenproblem to regular eigenproblem
  do i=1, at%N
    do j=1, at%N
      dm((i-1)*3+1:(i-1)*3+3,(j-1)*3+1:(j-1)*3+3) = dm((i-1)*3+1:(i-1)*3+3,(j-1)*3+1:(j-1)*3+3) / &
					sqrt(ElementMass(at%Z(i))*ElementMass(at%Z(j)))
    end do
  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  n_zero = 0
  if (do_zero_translation) n_zero = n_zero + 3
  if (do_zero_rotation) n_zero = n_zero + 3

  if (n_zero > 0) then

    allocate(zero_phonon(at%N*3,n_zero))
    allocate(zero_phonon_p(at%N*3,n_zero))
    allocate(phonon(3,at%N))
    do i=1, n_zero
      if (do_zero_translation .and. i <= 3) then ! translation
	beta = i
	phonon = 0.0_dp
	phonon(beta,:) = 1.0_dp
      else ! rotation
	beta = i-3
	CoM = centre_of_mass(at)
	axis = 0.0_dp; axis(beta) = 1.0_dp
	do j=1, at%N
	  dr_proj = at%pos(:,j)-CoM
	  dr_proj(beta) = 0.0_dp
	  phonon(:,j) = dr_proj .cross. axis
	end do
      endif
      phonon_norm=sqrt(sum(ElementMass(at%Z)*sum(phonon**2,1)))
      zero_phonon(:,i) = reshape(phonon/phonon_norm, (/ 3*at%N /) )
    end do ! i
    deallocate(phonon)

    ! transform from generalized eigenproblem to regular eigenproblem
    do i=1, at%N
      zero_phonon((i-1)*3+1:(i-1)*3+3,:) = zero_phonon((i-1)*3+1:(i-1)*3+3,:)*sqrt(ElementMass(at%Z(i)))
    end do

    allocate(zero_overlap_inv(n_zero,n_zero))
    ! project out zero frequency modes
    do i=1, n_zero
      do j=1, n_zero
	zero_overlap_inv(i,j) = sum(zero_phonon(:,i)*zero_phonon(:,j))
      end do
    end do
    call inverse(zero_overlap_inv)

    zero_phonon_p = 0.0_dp; call matrix_product_sub(zero_phonon_p, zero_phonon, zero_overlap_inv)
    deallocate(zero_overlap_inv)

    allocate(dm_t(at%N*3,at%N*3))
    allocate(P(at%N*3,at%N*3))
    P = 0.0_dp; call matrix_product_sub(P, zero_phonon_p, zero_phonon, .false., .true.)
    deallocate(zero_phonon_p)
    P = -P
    call add_identity(P)

    dm_t = 0.0_dp; call matrix_product_sub(dm_t, dm, P)
    dm = 0.0_dp; call matrix_product_sub(dm, P, dm_t)
    deallocate(dm_t)
    deallocate(P)
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! symmetrize dynamical matrix exactly
  do i=1, 3*at%N
    do j=i+1, 3*at%N
      dm(i,j) = dm(j,i)
    end do
  end do

  call print("dm", NERD)
  call print(dm, NERD)

  ! diagonalise dynamical matrix
  call diagonalise(dm, evals, evecs, err)
  if (err /= 0) then
    call system_abort("calc_phonons got error " // err // " in diagonalise")
  endif

  if (override_zero_freq_phonons .and. do_zero_rotation) then
    zero_phonon(:,n_zero-1) = zero_phonon(:,n_zero-1)-zero_phonon(:,n_zero-2)*sum(zero_phonon(:,n_zero-1)*zero_phonon(:,n_zero-2))
    zero_phonon(:,n_zero-1) = zero_phonon(:,n_zero-1)/sqrt(sum(zero_phonon(:,n_zero-1)**2))
    zero_phonon(:,n_zero) = zero_phonon(:,n_zero)-zero_phonon(:,n_zero-2)*sum(zero_phonon(:,n_zero)*zero_phonon(:,n_zero-2))
    zero_phonon(:,n_zero) = zero_phonon(:,n_zero)-zero_phonon(:,n_zero-1)*sum(zero_phonon(:,n_zero)*zero_phonon(:,n_zero-1))
    zero_phonon(:,n_zero) = zero_phonon(:,n_zero)/sqrt(sum(zero_phonon(:,n_zero)**2))
    evecs(:,1:n_zero) = zero_phonon
  endif
  deallocate(zero_phonon)

  ! transform from evecs of regular eigenproblem to evecs of original generalized eigenproblem
  if (present(evecs)) then
    do i=1, at%N
      evecs((i-1)*3+1:(i-1)*3+3,:) = evecs((i-1)*3+1:(i-1)*3+3,:) / sqrt(ElementMass(at%Z(i))) 
    end do
  endif

  ! calculate effective masses
  if (present(effective_masses)) then
    do i=1, 3*at%N
      effective_masses(i) = 0.0_dp
      do j=1, at%N
	effective_masses(i) = 1.0_dp/sum(evecs(:,i)**2)
      end do
    end do
  endif

  if (present(IR_intensities)) then
    do i=1, at%N*3
      dmu_dq(1) = sum(dmu_dr(:,:,1)*reshape(evecs(:,i), (/ 3, at%N /) ) )
      dmu_dq(2) = sum(dmu_dr(:,:,2)*reshape(evecs(:,i), (/ 3, at%N /) ) )
      dmu_dq(3) = sum(dmu_dr(:,:,3)*reshape(evecs(:,i), (/ 3, at%N /) ) )
      IR_intensities(i) = 3.0_dp/(PI*3.0*10**3)*sum(dmu_dq**2)
    end do
  endif

  call print("evals", VERBOSE)
  call print(evals, VERBOSE)
  if (present(evecs)) then
    call print("evecs", NERD)
    call print(evecs, NERD)
  endif
  if (present(effective_masses)) then
    call print("effective masses", VERBOSE)
    call print(effective_masses, VERBOSE)
  endif
  if (present(IR_intensities)) then
    call print("IR intensities", VERBOSE)
    call print(IR_intensities, VERBOSE)
  endif

  deallocate(dmu_dr)
  deallocate(pos0)
  deallocate(fp)
  deallocate(fm)
  deallocate(dm)

end subroutine phonons

end module phonons_module
