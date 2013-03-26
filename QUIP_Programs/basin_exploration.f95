program basin_exploration
use libatoms_module
use gp_basic_module
implicit none

type (gp_basic) :: gp
integer :: i_sample, n_samples
real(dp), allocatable :: sample_pos(:,:), sample_grad(:,:), sample_noise(:, :)
integer :: pos_init(2)

integer :: i_cycle, n_cycles, i_candidate, n_candidates
real(dp) :: rv(2), dpos(2),cur_pos(2)
real(dp) :: simplex_step_length, shooting_step_length, gaussian_width
integer, allocatable :: candidate_i_sample(:)
real(dp), allocatable :: candidate_pos(:,:), candidate_E(:), candidate_prob(:)
real(dp) :: sampling_T

call system_initialise(verbosity=PRINT_SILENT)

n_cycles = 500
n_samples = 500
n_candidates = 10
sampling_T = 0.1_dp

pos_init=(/ 0.0_dp, 0.0_dp /)

allocate(sample_pos(2,3+n_cycles), sample_grad(2,3+n_cycles), sample_noise(2, 3+n_cycles))

simplex_step_length=0.25_dp
shooting_step_length=0.10_dp
gaussian_width=0.5_dp

sample_pos(:,1) = pos_init
call eval_grad(sample_pos(:,1), sample_grad(:,1), sample_noise(:,1))
rv = (/ ran_normal(), ran_normal () /)
rv = simplex_step_length*rv/norm(rv)
sample_pos(:,2) = rv
call eval_grad(sample_pos(:,2), sample_grad(:,2), sample_noise(:,2))
rv = (/ ran_normal(), ran_normal () /)
dpos = sample_pos(:,1) - sample_pos(:,2)
rv = rv - (rv .dot. dpos)*dpos/(norm(dpos)**2)
sample_pos(:,3) = 0.5_dp*(sample_pos(:,1)+sample_pos(:,2)) + simplex_step_length*sqrt(3.0_dp)/2.0_dp*rv/norm(rv)
call eval_grad(sample_pos(:,3), sample_grad(:,3), sample_noise(:,3))
n_samples = 3

! print *, "sample_pos ", sample_pos(:,1)
! print *, "sample_pos ", sample_pos(:,2)
! print *, "sample_pos ", sample_pos(:,3)

call initialise(gp, (/ gaussian_width, gaussian_width /), (/ 0.0_dp, 0.0_dp /), 0.25_dp, SE_kernel_r_rr, &
  g_r=sample_pos(:,1:n_samples), g_v=sample_grad(:,1:n_samples), g_n=sample_noise(:,1:n_samples))
call print_gp(gp)
do i_sample=1, n_samples
  print *, sample_pos(1:2,i_sample)
end do
print *, ""
print *, ""

allocate(candidate_i_sample(n_candidates), candidate_pos(2,n_candidates), candidate_E(n_candidates))
allocate(candidate_prob(0:n_candidates))

do i_cycle=1, n_cycles
  ! generate candidates
  do i_candidate=1, n_candidates
    candidate_i_sample(i_candidate) = int(ran_uniform()*n_samples)+1
    call generate_sample_pos(candidate_pos(:,i_candidate), sample_pos(:, candidate_i_sample(i_candidate)), gp)
    candidate_E(i_candidate) = f_predict(gp, candidate_pos(:, i_candidate), SE_kernel_r_rr)
  end do

  ! select candidate
  candidate_prob(0) = 0.0_dp
  candidate_prob(1:n_candidates) = exp(-candidate_E/sampling_T)
  candidate_prob = candidate_prob / sum(candidate_prob)
  do i_candidate=1, n_candidates
    candidate_prob(i_candidate) = candidate_prob(i_candidate) + candidate_prob(i_candidate-1)
  end do
  rv(1) = ran_uniform()
  do i_candidate=1, n_candidates
    if (rv(1) >= candidate_prob(i_candidate-1) .and. rv(1) < candidate_prob(i_candidate)) then
      cur_pos = candidate_pos(:,i_candidate)
      exit
    end if
  end do

  ! evaluate at cur_pos
  n_samples = n_samples + 1
  sample_pos(:, n_samples) = cur_pos
  call eval_grad(sample_pos(:,n_samples), sample_grad(:,n_samples), sample_noise(:,n_samples))
  call initialise(gp, (/ gaussian_width, gaussian_width /), (/ 0.0_dp, 0.0_dp /), 0.25_dp, SE_kernel_r_rr, &
    g_r=sample_pos(:,1:n_samples), g_v=sample_grad(:,1:n_samples), g_n=sample_noise(:,1:n_samples))
  call print_gp(gp)
  do i_sample=1, n_samples
    print *, sample_pos(1:2,i_sample)
  end do
  print *, ""
  print *, ""
end do

call system_finalise()

contains

subroutine generate_sample_pos(new_pos, init_pos, gp)
  real(dp), intent(out) :: new_pos(2)
  real(dp), intent(in) :: init_pos(2)
  type(gp_basic), intent(inout) :: gp

  rv = (/ ran_normal(), ran_normal() /)
  rv = rv / norm(rv)
  new_pos = sample_pos(:, i_sample)
  new_pos = new_pos + shooting_step_length*rv
  do while (f_predict_var(gp, new_pos, SE_kernel_r_rr) < 0.10_dp)
    new_pos = new_pos + shooting_step_length*rv
  end do
end subroutine generate_sample_pos

subroutine print_gp(gp)
  type(gp_basic), intent(inout) :: gp

  integer:: i, j
  real(dp) :: pos(2)
  real(dp) :: fv

  integer :: grid_size = 40

  do i=0, grid_size
  do j=0, grid_size
    pos = (/ 5.0*(2.0*(real(i, dp)/real(grid_size, dp)-0.5_dp)), 5.0*(2.0*(real(j, dp)/real(grid_size, dp)-0.5_dp)) /)
    call eval_func(pos, fv)
    write (unit=*, fmt=*) pos, f_predict(gp, pos, SE_kernel_r_rr), f_predict_var(gp, pos, SE_kernel_r_rr), fv
  end do
  end do
  write (unit=*, fmt='(A)') ""
  write (unit=*, fmt='(A)') ""
end subroutine

subroutine eval_func(pos, func)
  real(dp), intent(in) :: pos(2)
  real(dp), intent(out) :: func

  real(dp) :: r0(2) = (/ 0.0_dp, 0.0_dp /)
  real(dp) :: r1(2) = (/ 3.0_dp, 0.0_dp /)

  func = -exp(-normsq(pos-r0)) - exp(-normsq(pos-r1))

end subroutine eval_func

subroutine eval_grad(pos, grad, noise)
  real(dp), intent(in) :: pos(2)
  real(dp), intent(out) :: grad(2), noise(2)

  real(dp) :: r0(2) = (/ 0.0_dp, 0.0_dp /)
  real(dp) :: r1(2) = (/ 3.0_dp, 0.0_dp /)

  ! E = -exp(-|pos-0|^2) - exp(-|pos-(2,0)|^2)
  grad = exp(-normsq(pos-r0))*2.0_dp*(pos-r0) + exp(-normsq(pos-r1))*2.0_dp*(pos-r1) + 0.01_dp*(/ ran_normal(), ran_normal() /)
  noise = 0.01_dp**2

end subroutine eval_grad

end program
