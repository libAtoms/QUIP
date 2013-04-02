program basin_exploration
use libatoms_module
use k_means_clustering_module
use gp_basic_module
implicit none

integer :: n_dim

type (gp_basic), target :: gp
type (gp_basic), pointer :: gp_sparse
integer :: i_sample, n_samples
real(dp), allocatable :: sample_pos(:,:), sample_grad(:,:), sample_noise(:, :)

integer :: i_cycle, n_cycles, i_candidate, n_candidates
real(dp), allocatable :: rv(:), cur_pos(:)
real(dp) :: simplex_step_length, shooting_step_length, gaussian_width, step_var_target
integer, allocatable :: candidate_i_sample(:)
real(dp), allocatable :: candidate_pos(:,:), candidate_E(:), candidate_prob(:)
real(dp) :: sampling_T
real(dp), allocatable :: len_scale(:), periodicity(:), sparse_len_scale(:)

integer :: n_sparse_samples
integer, allocatable :: g_sparse_set(:)

logical :: doubling_interp, doubling_interp_repeat
integer :: gp_guess_sparsify

call system_initialise(verbosity=PRINT_SILENT, enable_timing=.true.)
call verbosity_push(PRINT_NORMAL)

read *, n_dim
read *, n_cycles, n_candidates
read *, sampling_T
read *, simplex_step_length, shooting_step_length, step_var_target
read *, gaussian_width
read *, doubling_interp, doubling_interp_repeat, gp_guess_sparsify

print *, "n_dim ", n_dim
print *, "n_cycles n_candidates ", n_cycles, n_candidates
print *, "sampling_T ", sampling_T
print *, "simplex_step_length shooting_step_length ",simplex_step_length, shooting_step_length
print *, "gaussian_width ", gaussian_width
print *, "doubling_interp, doubling_interp_repeat, gp_guess_sparsify ", &
          doubling_interp, doubling_interp_repeat, gp_guess_sparsify

allocate(rv(n_dim), cur_pos(n_dim))
allocate(len_scale(n_dim), periodicity(n_dim), sparse_len_scale(n_dim))

! n_cycles = 500
! n_candidates = 1000
! sampling_T = 0.01_dp
! 
! simplex_step_length=0.25_dp
! shooting_step_length=0.10_dp
! gaussian_width=0.5_dp

len_scale = gaussian_width
periodicity = 0.0_dp

! sparsify_factor^(1/dim) = len_scale_factor
sparse_len_scale = len_scale * (gp_guess_sparsify**(1.0_dp/real(n_dim, dp)))
! sparse_len_scale = 2.0_dp*len_scale

print *, "len_scale ", len_scale
print *, "sparse_len_scale ", sparse_len_scale

allocate(sample_pos(n_dim,n_dim+1+n_cycles), sample_grad(n_dim,n_dim+1+n_cycles), sample_noise(n_dim, n_dim+1+n_cycles))

call create_simplex(sample_pos(:,1:n_dim+1))
sample_pos(:,1:n_dim+1) = simplex_step_length*sample_pos(:,1:n_dim+1) 
n_samples = n_dim+1
do i_sample=1, n_samples
   call eval_grad(sample_pos(:,i_sample), sample_grad(:,i_sample), sample_noise(:,i_sample))
end do

if (gp_guess_sparsify > 1) then
   allocate(gp_sparse)
else
   gp_sparse => gp
endif

call initialise(gp, len_scale, periodicity, 0.25_dp, SE_kernel_r_rr, &
  g_r=sample_pos(:,1:n_samples), g_v=sample_grad(:,1:n_samples), g_n=sample_noise(:,1:n_samples), jitter=1.0e-10_dp)
call print_gp(gp, gp_sparse, sample_pos(:, 1:n_samples), 0)
if (gp_guess_sparsify > 1) then
   call initialise(gp_sparse, sparse_len_scale,  periodicity, 0.25_dp, SE_kernel_r_rr, &
     g_r=sample_pos(:,1:n_samples), g_v=sample_grad(:,1:n_samples), g_n=sample_noise(:,1:n_samples), jitter=1.0e-10_dp)
endif

allocate(candidate_i_sample(n_candidates), candidate_pos(n_dim,n_candidates), candidate_E(n_candidates))
allocate(candidate_prob(0:n_candidates))

do i_cycle=1, n_cycles
call system_timer("generate")
  ! generate candidates
  do i_candidate=1, n_candidates
    candidate_i_sample(i_candidate) = int(ran_uniform()*n_samples)+1
    call generate_sample_pos(candidate_pos(:,i_candidate), sample_pos(:, candidate_i_sample(i_candidate)), gp, gp_sparse, step_var_target, &
      doubling_interp, doubling_interp_repeat)
    candidate_E(i_candidate) = f_predict(gp, candidate_pos(:, i_candidate), SE_kernel_r_rr)
  end do
call system_timer("generate")

call system_timer("select")
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
      print *, i_cycle, " selected sample pos, E ",cur_pos, candidate_E(i_candidate)
      exit
    end if
  end do
call system_timer("select")

call system_timer("evaluate_teach")
  ! evaluate at cur_pos
  n_samples = n_samples + 1
  sample_pos(:, n_samples) = cur_pos
  call eval_grad(sample_pos(:,n_samples), sample_grad(:,n_samples), sample_noise(:,n_samples))
  call initialise(gp, len_scale, periodicity, 0.25_dp, SE_kernel_r_rr, &
    g_r=sample_pos(:,1:n_samples), g_v=sample_grad(:,1:n_samples), g_n=sample_noise(:,1:n_samples), jitter=1.0e-10_dp)
call system_timer("evaluate_teach")
call system_timer("print")
  call print_gp(gp, gp_sparse, sample_pos(:, 1:n_samples), i_cycle)
call system_timer("print")
  if (gp_guess_sparsify > 1) then
call system_timer("evaluate_teach_sparse")
     if (n_samples > gp_guess_sparsify*10) then
	n_sparse_samples = n_samples/gp_guess_sparsify
	allocate(g_sparse_set(n_sparse_samples))
	call k_means_clustering_pick(sample_pos(:,1:n_samples), periodicity, g_sparse_set)
	call initialise(gp_sparse, sparse_len_scale,  periodicity, 0.25_dp, SE_kernel_r_rr, &
	  g_r=sample_pos(:,1:n_samples), g_v=sample_grad(:,1:n_samples), g_n=sample_noise(:,1:n_samples), &
	  g_sparse_set=g_sparse_set, jitter=1.0e-10_dp)
	deallocate(g_sparse_set)
     else
	call initialise(gp_sparse, sparse_len_scale, periodicity, 0.25_dp, SE_kernel_r_rr, &
	  g_r=sample_pos(:,1:n_samples), g_v=sample_grad(:,1:n_samples), g_n=sample_noise(:,1:n_samples), jitter=1.0e-10_dp)
     endif
call system_timer("evaluate_teach_sparse")
  endif

end do

call system_finalise()

contains

subroutine generate_sample_pos(new_pos, init_pos, gp, gp_for_var, step_var_target, doubling_interp, doubling_interp_repeat)
  real(dp), intent(out) :: new_pos(:)
  real(dp), intent(in) :: init_pos(:)
  type(gp_basic), intent(inout) :: gp, gp_for_var
  real(dp), intent(in) :: step_var_target
  logical, intent(in) :: doubling_interp, doubling_interp_repeat

  real(dp) :: rv(size(new_pos)), new_interpolated_pos(size(new_pos))
  real(dp) :: prev_var, new_var
  real(dp) :: step_size, step_frac

  integer :: n_dim, i_dim

  n_dim = size(new_pos)

  step_size = shooting_step_length

  do i_dim=1, n_dim
    rv(i_dim) = ran_normal()
  end do
  rv = rv / norm(rv)
  new_pos = sample_pos(:, i_sample) + step_size*rv
  prev_var = 0.0_dp
  new_var = f_predict_var(gp_for_var, new_pos, SE_kernel_r_rr)
  print *, "predict_var initial ", new_var
  do while (new_var < step_var_target)
    if (doubling_interp) step_size = 2.0_dp * step_size
    new_pos = new_pos + step_size*rv
    prev_var = new_var
    new_var = f_predict_var(gp_for_var, new_pos, SE_kernel_r_rr)
    print *, "predict_var shooting ", new_var
  end do

  if (doubling_interp) then
     print *, "final vars ", prev_var, new_var
     ! prev_var + step_frac*(new_var-prev_var) = step_var_target
     step_frac = (step_var_target-prev_var)/(new_var-prev_var)
     new_interpolated_pos = new_pos - (1.0_dp-step_frac)*step_size*rv
     if (doubling_interp_repeat) then
	prev_var = f_predict_var(gp_for_var, new_interpolated_pos, SE_kernel_r_rr)
	print *, "predict_var initial interpolation ", prev_var
	do while (prev_var < step_var_target)
	   step_size = (1.0_dp-step_frac)*step_size
	   step_frac = 0.5_dp
	   new_interpolated_pos = new_pos - (1.0_dp-step_frac)*step_size*rv
	   prev_var = f_predict_var(gp_for_var, new_interpolated_pos, SE_kernel_r_rr)
	   print *, "predict_var interpolation ", prev_var
	end do
     endif ! doubling_interp_repeat
     new_pos = new_interpolated_pos
  endif ! doubling_interp
  !
  new_var = f_predict_var(gp_for_var, new_pos, SE_kernel_r_rr)
  print *, "candidate final pos var ", new_pos, new_var, f_predict_var(gp, new_pos, SE_kernel_r_rr)
  !

end subroutine generate_sample_pos

subroutine print_gp(gp, gp_sparse, sample_pos, i_label)
  type(gp_basic), intent(inout) :: gp, gp_sparse
  real(dp), intent(in) :: sample_pos(:,:)
  integer, intent(in) :: i_label

  integer:: i, j
  real(dp) :: pos(size(sample_pos,1)), origin(size(sample_pos,1))
  real(dp) :: true_fval, err_offset, plot_offset

  integer :: grid_size = 40
  integer :: err_count
  real(dp) :: err


  ! align origin (minimum) for error offset
  origin = 0.0_dp ! ; origin(1) = -5.0_dp
  call eval_func(origin, err_offset)
  err_offset = err_offset - f_predict(gp, origin, SE_kernel_r_rr)

  ! align far field for plot offset
  origin = 0.0_dp; origin(1) = -5.0_dp
  call eval_func(origin, plot_offset)
  plot_offset = plot_offset - f_predict(gp, origin, SE_kernel_r_rr)

  open (unit=100, file="gp."//i_label, status="unknown")
  err = 0.0_dp
  err_count = 0
  do i=0, grid_size
  do j=0, grid_size
    pos = 0.0_dp
    pos(1:2) = (/ 5.0*(2.0*(real(i, dp)/real(grid_size, dp)-0.5_dp)), 5.0*(2.0*(real(j, dp)/real(grid_size, dp)-0.5_dp)) /)
    call eval_func(pos, true_fval)
    write (unit=100, fmt=*) f_predict(gp, pos, SE_kernel_r_rr)+plot_offset, f_predict_var(gp, pos, SE_kernel_r_rr), &
       f_predict_var(gp_sparse, pos, SE_kernel_r_rr), true_fval, pos(1:2)
    if (f_predict_var(gp, pos, SE_kernel_r_rr) < 0.1_dp) then
       err = err + (f_predict(gp, pos, SE_kernel_r_rr)+err_offset-true_fval)**2
       err_count = err_count + 1
    endif
  end do
  end do
  write (unit=100, fmt='(A)') ""
  write (unit=100, fmt='(A)') ""
  do i=1, size(sample_pos, 2)
    write(unit=100, fmt=*) sample_pos(:, i)
  end do
  write (unit=100, fmt='(A)') ""
  write (unit=100, fmt='(A)') ""
  err = sqrt(err/real(err_count,dp))
  write (unit=100, fmt='("ERR ",I6,F10.5)') err_count, err

  close (unit=100)

end subroutine print_gp

subroutine eval_func(pos, func)
  real(dp), intent(in) :: pos(:)
  real(dp), intent(out) :: func

  real(dp) :: r0(size(pos)), r1(size(pos))

  r0 = 0.0_dp
  r1 = 0.0_dp; r1(1) = 3.0_dp

  func = -exp(-normsq(pos-r0)) - exp(-normsq(pos-r1))

end subroutine eval_func

subroutine eval_grad(pos, grad, noise)
  real(dp), intent(in) :: pos(:)
  real(dp), intent(out) :: grad(:), noise(:)

  integer  :: i_dim
  real(dp) :: r0(size(pos)), r1(size(pos))

  n_dim = size(pos)

  if (size(grad) /= n_dim ) then
     call system_abort("eval_grad got grad with bad dimensions "//shape(grad))
  endif
  if (size(noise) /= n_dim ) then
     call system_abort("eval_grad got noise with bad dimensions "//shape(noise))
  endif

  r0 = 0.0_dp
  r1 = 0.0_dp; r1(1) = 3.0_dp

  grad = exp(-normsq(pos-r0))*2.0_dp*(pos-r0) + exp(-normsq(pos-r1))*2.0_dp*(pos-r1)
  do i_dim=1, n_dim
    grad(i_dim) = grad(i_dim) + 0.01_dp*ran_normal()
  end do
  noise = 0.01_dp**2

end subroutine eval_grad

subroutine create_simplex(pos)
   real(dp) :: pos(:, :)

   integer :: n_dim, i_dim
   real(dp) :: target_dot, cur_dot

   n_dim = size(pos, 1)
   if (size(pos, 2) /= n_dim + 1) then
      call system_abort("create_simple got pos with bad dimensions "//shape(pos))
   endif

   target_dot = -1.0_dp / real(n_dim, dp)

   pos = 0.0_dp
   ! first vector
   pos(:,1) = 0.0_dp; pos(1,1) = 1.0_dp
   ! all others
   do i_dim=2, n_dim+1
      if (i_dim == 2) then
	 cur_dot = 0.0_dp
      else
	 cur_dot = pos(1:i_dim-2, i_dim) .dot. pos(1:i_dim-2, i_dim-1)
      endif
      pos(i_dim-1, i_dim:n_dim+1) = (target_dot - cur_dot)/pos(i_dim-1,i_dim-1)
      pos(i_dim, i_dim) = sqrt(1.0_dp-normsq(pos(1:i_dim-1,i_dim)))
   end do

   pos = pos / norm(pos(:,1)-pos(:,2))
end subroutine create_simplex

end program
