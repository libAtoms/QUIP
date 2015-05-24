! H0 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! H0 X
! H0 X   libAtoms+QUIP: atomistic simulation library
! H0 X
! H0 X   Portions of this code were written by
! H0 X     Albert Bartok-Partay, Silvia Cereda, Gabor Csanyi, James Kermode,
! H0 X     Ivan Solt, Wojciech Szlachta, Csilla Varnai, Steven Winfield.
! H0 X
! H0 X   Copyright 2006-2010.
! H0 X
! H0 X   These portions of the source code are released under the GNU General
! H0 X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html
! H0 X
! H0 X   If you would like to license the source code under different terms,
! H0 X   please contact Gabor Csanyi, gabor@csanyi.net
! H0 X
! H0 X   Portions of this code were written by Noam Bernstein as part of
! H0 X   his employment for the U.S. Government, and are not subject
! H0 X   to copyright in the USA.
! H0 X
! H0 X
! H0 X   When using this software, please cite the following reference:
! H0 X
! H0 X   http://www.libatoms.org
! H0 X
! H0 X  Additional contributions by
! H0 X    Alessio Comisso, Chiara Gattinoni, and Gianpietro Moras
! H0 X
! H0 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine matrix_general_diagonlise_lobpcg_sandia(A, M, evals, evecs, block_size, &
        				      precond, no_guess, err)
    real(dp),intent(in), dimension(:,:) :: A
    real(dp),intent(in), dimension(:,:) :: M
    real(dp),intent(inout), dimension(:) :: evals
    real(dp),intent(inout), dimension(:,:) :: evecs
    integer, intent(in), optional :: block_size
    real(dp),intent(in), dimension(:,:), optional :: precond(:,:)
    logical, intent(in), optional :: no_guess
    integer, intent(out), optional :: err

    real(dp), allocatable :: X_twid(:,:), X(:,:), RI(:,:), R(:), PI(:,:), &
      HI(:,:), S(:,:), Y(:,:), YI(:,:), A_X(:,:), M_X(:,:), rr_M(:,:), M_HI(:,:), M_YI(:,:)
    real(dp), allocatable :: evals_block(:), Y_evals(:)

    logical :: done, failed
    logical :: my_no_guess
    integer :: n_evecs, n_full_blocks, extra_block_size
    integer :: N, i, iter, max_n_iter, bs, cbs, i_block, ii_block, nrr, prev_cbs
    integer :: first_eval, last_eval, max_bs, max_n_new_block, got_n_evals
    real(dp) :: max_r, max_r_hist(10)
    integer :: n_needed

    call print("matrix_general_diagonlise_lobpcg_sandia", PRINT_VERBOSE)

    N = size(A, 1)

    my_no_guess = optional_default(.true., no_guess)
    max_bs = optional_default(1, block_size)

    n_evecs = size(evals)
    if (n_evecs /= size(evecs,2)) &
      call system_abort("called matrix_general_digonalise_lobpcg_sandia with mismatched evals and evecs size")

    call print("matrix_general_diagonalise_lpbpcg_sandia n_evecs " // &
      n_evecs // " block_size " // max_bs, PRINT_VERBOSE)

    call print("n_full_blocks " // n_full_blocks // " extra_block_size " // extra_block_size, PRINT_VERBOSE)

    ! randomise evecs
    if (my_no_guess) then
      evecs = 0.0_dp
      do i=1, n_evecs
        call randomise(evecs(:,i), 1.0_dp)
      end do
    endif

    max_n_iter = 200
    max_n_new_block = 30

    allocate(X_twid(N,max_bs))
    allocate(X(N,max_bs))
    allocate(A_X(N,max_bs))
    allocate(M_X(N,max_bs))
    allocate(evals_block(max_bs))
    allocate(R(N))
    allocate(RI(N,max_bs))
    allocate(PI(N,max_bs))
    allocate(HI(N,max_bs))
    allocate(M_HI(N,max_bs))
    allocate(Y(3*max_bs,max_bs))
    allocate(YI(3*max_bs,max_bs))
    allocate(M_YI(3*max_bs,max_bs))
    allocate(Y_evals(max_bs))
    allocate(S(N,3*max_bs))
    allocate(rr_M(3*max_bs,3*max_bs))

    got_n_evals = 0
    i_block = 1
    failed = .false.
    do while (got_n_evals < n_evecs .and. .not. failed)

      first_eval = got_n_evals + 1
      last_eval = min(first_eval+max_bs-1, N)
      bs = last_eval-first_eval+1

      call print("doing block " // first_eval // " - " // last_eval, PRINT_VERBOSE)

      n_needed = n_evecs-first_eval
      ! copy block of eigenvectors
      X_twid(1:N,1:min(bs,n_needed)) = evecs(1:N,first_eval:min(last_eval,n_evecs))
      if (n_needed < bs) then
        do i=n_needed+1, bs
          call randomise(X_twid(1:N,i), 1.0_dp)
        end do
      endif

      if (first_eval > 1) then
	call matrix_product_sub(X(1:N,1:bs), M(1:N,1:N), X_twid(1:N,1:bs))
        call orthogonalise_col_vectors(X_twid(1:N,1:bs), evecs(1:N,1:first_eval-1), X(1:N,1:bs))
      endif

      call rr_solve(X_twid(1:N,1:bs), A(1:N,1:N), M(1:N,1:N), Y(1:bs,1:bs), Y_evals(1:bs))
      call matrix_product_sub(X(1:N,1:bs),X_twid(1:N,1:bs),Y(1:bs,1:bs))
      call matrix_product_sub(A_X(1:N,1:bs), A(1:N,1:N), X(1:N,1:bs))
      call matrix_product_sub(M_X(1:N,1:bs), M(1:N,1:N), X(1:N,1:bs))
      do i=1, bs
        RI(1:N,i) = A_X(1:N,i) - Y_evals(i)*M_X(1:N,i)
      end do
      PI(1:N,1:bs) = 0.0_dp
      cbs = bs
      iter = 1
      done = .false.
      max_r_hist = 1e38_dp
      do while (iter < max_n_iter .and. .not. done .and. .not. failed)
        call print("doing iter "//iter, PRINT_NERD)
        if (present(precond)) then
          call matrix_product_sub(HI(1:N,1:cbs), precond(1:N,1:N), RI(1:N,1:cbs))
        else
          HI(1:N,1:cbs) = RI(1:N,1:cbs)
        endif

        S(1:N,1:bs) = X(1:N,1:bs)
        S(1:N,bs+1:bs+cbs) = PI(1:N,1:cbs)

	call matrix_product_sub(M_HI(1:N,1:cbs), M(1:N,1:N), HI(1:N,1:cbs))
        call orthogonalise_col_vectors(HI(1:N,1:cbs), S(1:N,1:bs+cbs), M_HI(1:N,1:cbs))

        if (first_eval > 1) then
          call orthogonalise_col_vectors(HI(1:N,1:cbs), evecs(1:N,1:first_eval-1), M_HI(1:N,1:cbs))
        end if

        S(1:N,1:bs) = X(1:N,1:bs)
        S(1:N,bs+1:bs+cbs) = HI(1:N,1:cbs)
        S(1:N,bs+cbs+1:bs+2*cbs) = PI(1:N,1:cbs)
        if (iter == 1) then
          nrr = bs+cbs
        else
          nrr = bs+2*cbs
        endif
        call rr_solve(S(1:N,1:nrr), A(1:N,1:N), M(1:N,1:N), Y(1:nrr,1:bs), Y_evals(1:bs), &
                      rr_M=rr_M(1:nrr,1:nrr))

        call matrix_product_sub(X(1:N,1:bs), S(1:N,1:nrr), Y(1:nrr,1:bs))

        call matrix_product_sub(A_X(1:N,1:bs), A(1:N,1:N), X(1:N,1:bs))
        call matrix_product_sub(M_X(1:N,1:bs), M(1:N,1:N), X(1:N,1:bs))
        if (iter == 1) then
          YI = 0.0_dp
        endif
	prev_cbs = cbs
        cbs = 0
        max_r = 0.0_dp
        do i=1, bs
          R(1:N) = A_X(1:N,i) - Y_evals(i)*M_X(1:N,i)
          max_r = max(max_R, maxval(abs(R)))
          if (maxval(abs(R(1:N))) > 1.0e-7_dp) then ! column not converged
            cbs = cbs + 1
            RI(1:N,cbs) = R(1:N)
            YI(1:nrr,cbs) = Y(1:nrr,i)
          endif
        end do
        call print("residual max " // max_r // " cbs " // cbs, PRINT_VERBOSE)
        max_r_hist(1:9) = max_r_hist(2:10)
        max_r_hist(10) = max_r

	if (prev_cbs < cbs) then
          call print("diverging, but block partially converged " // (bs-cbs), PRINT_VERBOSE)
          done = .true.
          n_needed = n_evecs-first_eval+1
          got_n_evals = got_n_evals + (bs-cbs)
          ! copy all evecs and converged evals
          evecs(1:N,first_eval:min(first_eval+(bs)-1,n_evecs)) = X(1:N,1:min(bs,n_needed))
          evals(first_eval:min(first_eval+(bs-cbs)-1,n_evecs)) = Y_evals(1:min(bs-cbs,n_needed))
	else if (cbs == 0) then ! all converged
          done = .true.
          n_needed = n_evecs-first_eval+1
          got_n_evals = got_n_evals + bs
          ! copy converged evecs/evals
          evecs(1:N,first_eval:min(last_eval,n_evecs)) = X(1:N,1:min(bs,n_needed))
          evals(first_eval:min(last_eval,n_evecs)) = Y_evals(1:min(bs,n_needed))
	else if (got_n_evals+(bs-cbs) >= n_evecs) then
          call print("enough evecs done, but block partially converged " // (bs-cbs), PRINT_VERBOSE)
          done = .true.
          n_needed = n_evecs-first_eval+1
          got_n_evals = got_n_evals + (bs-cbs)
          ! copy all evecs and converged evals
          evecs(1:N,first_eval:min(first_eval+(bs)-1,n_evecs)) = X(1:N,1:min(bs,n_needed))
          evals(first_eval:min(first_eval+(bs-cbs)-1,n_evecs)) = Y_evals(1:min(bs-cbs,n_needed))
        else if (cbs /= bs .and. (iter > max_n_new_block .or. max_r_hist(10) > max_r_hist(5)*0.75_dp)) then
          call print("block partially converged " // (bs-cbs), PRINT_VERBOSE)
          done = .true.
          n_needed = n_evecs-first_eval+1
          got_n_evals = got_n_evals + (bs-cbs)
          ! copy all evecs and converged evals
          evecs(1:N,first_eval:min(first_eval+(bs)-1,n_evecs)) = X(1:N,1:min(bs,n_needed))
          evals(first_eval:min(first_eval+(bs-cbs)-1,n_evecs)) = Y_evals(1:min(bs-cbs,n_needed))
        else
          done = .false.
          YI(1:bs,1:cbs) = 0.0_dp
	  call matrix_product_sub(M_YI(1:nrr,1:bs), rr_M(1:nrr,1:nrr), YI(1:nrr, 1:bs))
          call orthogonalise_col_vectors(YI(1:nrr,1:cbs), Y(1:nrr,1:bs), M_YI(1:nrr,1:cbs))

          call matrix_product_sub(PI(1:N,1:cbs), S(1:N,1:nrr), YI(1:nrr,1:cbs))
        endif

        iter = iter + 1

      end do ! iter
      call print("block " // i_block // " took " // iter // " iterations")
      if (iter >= max_n_iter) then
        failed = .true.
      endif
    end do ! i_block

    deallocate(X_twid)
    deallocate(X)
    deallocate(A_X)
    deallocate(M_X)
    deallocate(evals_block)
    deallocate(R)
    deallocate(RI)
    deallocate(PI)
    deallocate(HI)
    deallocate(M_HI)
    deallocate(Y)
    deallocate(YI)
    deallocate(M_YI)
    deallocate(Y_evals)
    deallocate(S)
    deallocate(rr_M)

    if (present(err)) err = 0

    if (failed) then
      if (present(err)) then
        err = 1
        call print("matrix_general_diagonalise_lobpcg_sandia failed")
        return
      else
        call system_abort("matrix_general_diagonalise_lobpcg_sandia failed")
      endif
    endif

  end subroutine matrix_general_diagonlise_lobpcg_sandia

  subroutine rr_solve(X, A, M, Y, evals, rr_A, rr_M)
    real(dp), intent(in) :: X(:,:), A(:,:), M(:,:)
    real(dp), intent(out) :: Y(:,:)
    real(dp), intent(out) :: evals(:)
    real(dp), intent(out), target, optional :: rr_A(:,:), rr_M(:,:)

    integer :: N, rr_size
    real(dp), allocatable :: YY(:,:), ev(:)
    real(dp), allocatable :: A_X(:,:), M_X(:,:)
    real(dp), pointer :: l_rr_A(:,:), l_rr_M(:,:)

    if (size(A,1) /= size(A,2)) call system_abort("called rr_solve with non-square A")
    if (size(M,1) /= size(M,2)) call system_abort("called rr_solve with non-square M")
    if (size(X,1) /= size(A,1)) call system_abort("called rr_solve with X size not matching A/M size")
    if (size(X,2) /= size(Y,1)) call system_abort("called rr_solve with X size not matching Y size")
    if (size(X,2) < size(Y,2)) call system_abort("called rr_solve with X size too small for  Y size")
    if (size(X,2) < size(evals)) call system_abort("called rr_solve with Y size too small for evals size")

    N = size(X,1)
    rr_size = size(X,2)

    if (present(rr_A)) then
      if (size(rr_A,1) /= size(rr_A,2)) call system_abort("called rr_solve with non-square rr_A")
      if (size(rr_A,1) /= rr_size) call system_abort("called rr_solve with rr_A size not matching X size")
      l_rr_A => rr_A
    else
      allocate(l_rr_A(rr_size,rr_size))
    endif
    if (present(rr_M)) then
      if (size(rr_M,1) /= size(rr_M,2)) call system_abort("called rr_solve with non-square rr_M")
      if (size(rr_M,1) /= rr_size) call system_abort("called rr_solve with rr_M size not matching X size")
      l_rr_M => rr_M
    else
      allocate(l_rr_M(rr_size,rr_size))
    endif

    allocate(YY(rr_size,rr_size))
    allocate(ev(rr_size))
    allocate(A_X(N,rr_size))
    allocate(M_X(N,rr_size))

    call matrix_product_sub(A_X(1:N,1:rr_size), A(1:N,1:N), X(1:N,1:rr_size))
    call matrix_product_sub(M_X(1:N,1:rr_size), M(1:N,1:N), X(1:N,1:rr_size))
    call matrix_product_sub(l_rr_A(1:rr_size,1:rr_size), X(1:N,1:rr_size), A_X(1:N,1:rr_size), &
      m1_transpose=.true., m2_transpose=.false.)
    call matrix_product_sub(l_rr_M(1:rr_size,1:rr_size), X(1:N,1:rr_size), M_X(1:N,1:rr_size), &
      m1_transpose=.true., m2_transpose=.false.)

    call diagonalise(l_rr_A, l_rr_M, ev, YY)

    evals(:) = ev(1:size(evals))
    Y(1:rr_size,:) = YY(1:rr_size,1:size(Y,2))

    deallocate(YY)
    deallocate(ev)
    deallocate(A_X)
    deallocate(M_X)
    if (.not. present(rr_A)) deallocate(l_rr_A)
    if (.not. present(rr_M)) deallocate(l_rr_M)
  end subroutine rr_solve

!  subroutine matrix_general_diagonlise_lobpcg(B, A, evals, evecs, n_evecs, block_size, &
!					      precond, no_guess, largest, err)
!    real(dp),intent(in), dimension(:,:) :: A
!    real(dp),intent(in), dimension(:,:) :: B
!    real(dp),intent(inout), dimension(:) :: evals
!    real(dp),intent(inout), dimension(:,:) :: evecs
!    integer, intent(in), optional :: n_evecs
!    integer, intent(in), optional :: block_size
!    real(dp),intent(in), dimension(:,:), optional :: precond(:,:)
!    logical, intent(in), optional :: no_guess
!    logical, intent(in), optional :: largest
!    integer, intent(out), optional :: err
!
!    logical :: my_no_guess, my_largest
!    integer :: my_n_evecs, my_block_size
!    integer :: i_block, ii_block, n_full_blocks, extra_block_size
!    integer :: i, j
!    integer :: N
!    integer :: iter, max_n_iter
!    real(dp), allocatable :: evals_block(:)
!    real(dp), allocatable :: RR_evals(:)
!    real(dp), allocatable :: tt(:,:), R(:,:), W(:,:), P(:,:), RR_basis(:,:), RR_basis_AorB(:,:), &
!      RR_mat_A(:,:), RR_mat_B(:,:), evecs_block(:,:), A_evecs_block(:,:), B_evecs_block(:,:)
!    real(dp), allocatable :: RR_evecs(:,:)
!    logical :: done
!
!    integer :: i_max, rr_min, rr_max, rr_step
!    real(dp) :: mag
!
!    call print("matrix_general_diagonlise_lobpcg", PRINT_VERBOSE)
!
!    N = size(A, 1)
!    call print("N " // N, PRINT_VERBOSE)
!
!    my_no_guess = optional_default(.true., no_guess)
!    my_largest = optional_default(.false., largest)
!    my_n_evecs = optional_default(N, n_evecs)
!    my_block_size = optional_default(1, block_size)
!    n_full_blocks = my_n_evecs/my_block_size
!    extra_block_size = my_n_evecs - n_full_blocks*my_block_size
!    if (extra_block_size /= 0) &
!      call system_abort("lobpcg can't handle non integer number of blocks yet")
!
!    call print("matrix_general_diagonalise_lobpcg n_evecs " // &
!      my_n_evecs // " block_size " // my_block_size, PRINT_VERBOSE)
!
!    call print("n_full_blocks " // n_full_blocks // " extra_block_size " // extra_block_size, PRINT_VERBOSE)
!
!    allocate(evecs_block(N, my_block_size))
!    allocate(A_evecs_block(N, my_block_size))
!    allocate(B_evecs_block(N, my_block_size))
!    allocate(evals_block(my_block_size))
!    allocate(tt(N, my_block_size))
!    allocate(R(N, my_block_size))
!    allocate(W(N, my_block_size))
!    allocate(P(N, my_block_size))
!    allocate(RR_basis(N, 3*my_block_size))
!    allocate(RR_basis_AorB(N, 3*my_block_size))
!    allocate(RR_mat_A(3*my_block_size, 3*my_block_size))
!    allocate(RR_mat_B(3*my_block_size, 3*my_block_size))
!    allocate(RR_evecs(3*my_block_size, 3*my_block_size))
!    allocate(RR_evals(3*my_block_size))
!
!    ! randomise evecs
!    if (my_no_guess) then
!      evecs = 0.0_dp
!      do i=1, n_evecs
!	call randomise(evecs(:,i), 1.0_dp)
!      end do
!    endif
!
!    max_n_iter = 100
!
!    do i_block=1, n_full_blocks
!      call print("doing block " // i_block, PRINT_VERBOSE)
!
!      ! copy block of eigenvectors
!      evecs_block = evecs(:,(i_block-1)*my_block_size+1:i_block*my_block_size)
!      call print("raw evecs_block", PRINT_ANAL)
!      call print(evecs_block, PRINT_ANAL)
!
!      ! A orthogonalize evecs_block w.r.t. previous blocks
!      do ii_block=1, i_block-1
!	call orthogonalise_col_vectors(evecs_block, &
!	    evecs(:,(ii_block-1)*my_block_size+1:ii_block*my_block_size), A)
!      end do
!
!      call print("orthogonalized evecs_block", PRINT_ANAL)
!      call print(evecs_block, PRINT_ANAL)
!
!      P = 0.0_dp
!      iter = 1
!      done = .false.
!      do while (iter <= max_n_iter .and. .not. done)
!        call print("doing iter " // iter, PRINT_NERD)
!
!	! evaluate eigenvalues
!	call matrix_product_sub(A_evecs_block, A, evecs_block)
!	call matrix_product_sub(B_evecs_block, B, evecs_block)
!	if (iter == 1) &
!	  evals_block(:) = sum(evecs_block*B_evecs_block,1)/ &
!			   sum(evecs_block*A_evecs_block,1)
!
!	call print("evals_block", PRINT_NERD)
!	call print(evals_block, PRINT_NERD)
!	call print("evecs_block", PRINT_ANAL)
!	call print(evecs_block, PRINT_ANAL)
!
!	! evaluate residual
!	do i=1, block_size
!	  R(:,i) = B_evecs_block(:,i) - evals_block(i)*A_evecs_block(:,i)
!	end do
!
!	call print("R", PRINT_ANAL)
!	call print(R, PRINT_ANAL)
!
!	mainlog%default_real_precision=16
!	call print("residual max " // maxval(abs(R)), PRINT_VERBOSE)
!	mainlog%default_real_precision=8
!	done = (maxval(abs(R)) < 1.0e-8_dp)
!
!	if (.not. done) then
!
!	  ! precondition
!	  if (present(precond)) then
!	    call matrix_product_sub(W, precond, R)
!	  else
!	    W = R
!	  endif
!
!          call print("W", PRINT_ANAL)
!          call print(W, PRINT_ANAL)
!
!	  ! A orthogonalize W w.r.t. previous blocks
!	  do ii_block=1, i_block-1
!	    call orthogonalise_col_vectors(W, &
!		evecs(:,(ii_block-1)*my_block_size+1:ii_block*my_block_size), A)
!	  end do
!
!	  call print("orthogonalized W", PRINT_ANAL)
!	  call print(W, PRINT_ANAL)
!
!	  call print("P", PRINT_ANAL)
!	  call print(P, PRINT_ANAL)
!
!	  ! evaluate Rayleigh-Ritz basis
!	  RR_basis(:,1:my_block_size) = evecs_block
!	  RR_basis(:,my_block_size+1:2*my_block_size) = W
!	  RR_basis(:,2*my_block_size+1:3*my_block_size) = P
!
!
!	  call print("RR_basis", PRINT_ANAL)
!	  call print(RR_basis, PRINT_ANAL)
!
!	  ! evaluate Rayleigh-Ritz matrices
!	  call matrix_product_sub(RR_basis_AorB, A, RR_basis)
!
!	  call matrix_product_sub(RR_mat_A, RR_basis_AorB, RR_basis, &
!	    m1_transpose=.true., m2_transpose=.false.)
!	  call print("RR_mat_A", PRINT_ANAL)
!	  call print(RR_mat_A, PRINT_ANAL)
!	  call matrix_product_sub(RR_basis_AorB, B, RR_basis)
!	  call matrix_product_sub(RR_mat_B, RR_basis_AorB, RR_basis, &
!	    m1_transpose=.true., m2_transpose=.false.)
!	  call print("RR_mat_B", PRINT_ANAL)
!	  call print(RR_mat_B, PRINT_ANAL)
!
!	  ! diagonalise Rayleigh-Ritz matrix
!	  if (iter == 1) then
!	    call diagonalise(RR_mat_B(1:2*block_size,1:2*block_size), &
!	      RR_mat_A(1:2*block_size,1:2*block_size), &
!	      RR_evals(1:2*block_size), RR_evecs(1:2*block_size,1:2*block_size))
!	      call print("RR_evals", PRINT_NERD)
!	      call print(RR_evals(1:2*block_size), PRINT_NERD)
!	      call print("RR_evecs", PRINT_ANAL)
!	      call print(RR_evecs(1:2*block_size,1:2*block_size), PRINT_ANAL)
!	  else
!	    call diagonalise(RR_mat_B, RR_mat_A, RR_evals, RR_evecs)
!	    call print("RR_evals", PRINT_NERD)
!	    call print(RR_evals, PRINT_NERD)
!	    call print("RR_evecs", PRINT_ANAL)
!	    call print(RR_evecs, PRINT_ANAL)
!	  endif
!
!	  if (iter == 1) then
!	    i_max = 2*my_block_size
!	    if (my_largest) then
!	      rr_min = 2*my_block_size
!	      rr_max = my_block_size+1
!	      rr_step = -1
!	    else
!	      rr_min = my_block_size+1
!	      rr_max = 2*my_block_size
!	      rr_step = 1
!	    endif
!	  else
!	    i_max = 3*my_block_size
!	    if (my_largest) then
!	      rr_min = 3*my_block_size
!	      rr_max = 2*my_block_size+1
!	      rr_step = -1
!	    else
!	      rr_min = 2*my_block_size+1
!	      rr_max = 3*my_block_size
!	      rr_step = 1
!	    endif
!	  endif
!
!	  ! Create new evecs from Rayleigh-Ritz vectors
!	  call matrix_product_sub(evecs_block, RR_basis(:,1:i_max), &
!					       RR_evecs(1:i_max,rr_min:rr_max:rr_step))
!
!	  ! copy Rayleigh-Ritz evals (ignored, because of orthogonalization)
!	  evals_block = RR_evals(rr_min:rr_max:rr_step)
!
!	  ! Create new Ps from Rayleigh-Ritz vectors
!	  RR_evecs(1:my_block_size,:) = 0.0_dp
!	  call matrix_product_sub(P, RR_basis(:,1:i_max), RR_evecs(1:i_max,rr_min:rr_max:rr_step))
!
!	  iter = iter + 1
!	end if
!      end do ! while not converged
!
!      if (.not. done) then
!	if (present(err)) then
!	  err = 1
!	  deallocate(evecs_block)
!	  deallocate(A_evecs_block)
!	  deallocate(B_evecs_block)
!	  deallocate(evals_block)
!	  deallocate(tt)
!	  deallocate(R)
!	  deallocate(W)
!	  deallocate(P)
!	  deallocate(RR_basis)
!	  deallocate(RR_basis_AorB)
!	  deallocate(RR_mat_A)
!	  deallocate(RR_mat_B)
!	  deallocate(RR_evecs)
!	  deallocate(RR_evals)
!	  return
!	else
!	  call system_abort("matrix_general_diagonalise_lobpcg failed to converge in " // &
!	    iter // " iterations")
!	endif
!      endif
!
!      ! copy blocks back into main arrays
!      evecs(:,(i_block-1)*my_block_size+1:i_block*my_block_size) = evecs_block
!      evals((i_block-1)*my_block_size+1:i_block*my_block_size) = evals_block
!      call print("block " // i_block // " took " // iter // " iterations")
!      
!    end do ! i_block
!
!    if (present(err)) err = 0
!
!    deallocate(evecs_block)
!    deallocate(A_evecs_block)
!    deallocate(B_evecs_block)
!    deallocate(evals_block)
!    deallocate(tt)
!    deallocate(R)
!    deallocate(W)
!    deallocate(P)
!    deallocate(RR_basis)
!    deallocate(RR_basis_AorB)
!    deallocate(RR_mat_A)
!    deallocate(RR_mat_B)
!    deallocate(RR_evecs)
!    deallocate(RR_evals)
!
!  end subroutine matrix_general_diagonlise_lobpcg
!
!
!  ! assume A is symmetric
!  subroutine orthogonalise_col_vectors(V1, V2, A)
!    real(dp), intent(inout) :: V1(:,:)
!    real(dp), intent(in) :: V2(:,:)
!    real(dp), intent(in), optional :: A(:,:)
!
!    real(dp), allocatable :: A_V1(:,:), V2_V2T(:,:)
!    integer :: N, n_v1, n_v2
!integer i
!
!    if (present(A)) then
!      if (size(A,2) /= size(V1,1)) &
!	call system_abort("called orthogonalise_col_vectors with V1 size not matching A size")
!      if (size(A,1) /= size(V2,1)) &
!	call system_abort("called orthogonalise_col_vectors with V2 size not matching A size")
!    else
!      if (size(V1,1) /= size(V2,1)) &
!	call system_abort("called orthogonalise_col_vectors with V1 size not matching V2 size")
!    endif
!
!    if (present(A)) then
!      allocate(A_V1(size(A,1),size(V1,2)))
!      call matrix_product_sub(A_V1, A, V1)
!    endif
!
!    allocate(V2_V2T(size(V2,1), size(V2,1)))
!    call matrix_product_sub(V2_V2T, V2, V2, m1_transpose=.false., m2_transpose=.true.)
!
!    if (present(A)) then
!      call matrix_product_sub(V1, V2_V2T, A_V1, lhs_factor=1.0_dp, rhs_factor=-1.0_dp)
!call matrix_product_sub(A_V1, A, V1)
!do i=1, size(V1,2)
!  V1(:,i) = V1(:,i)/sqrt(sum(V1(:,i)*A_V1(:,i)))
!end do
!      deallocate(A_V1)
!    else
!      call matrix_product_sub(V1, V2_V2T, V1, lhs_factor=1.0_dp, rhs_factor=-1.0_dp)
!do i=1, size(V1,2)
!  V1(:,i) = V1(:,i)/sqrt(sum(V1(:,i)*V1(:,i)))
!end do
!    endif
!    deallocate(V2_V2T)
!
!  end subroutine orthogonalise_col_vectors

  ! assume A is symmetric
  subroutine orthogonalise_col_vectors(V1, V2, A_V1)
    real(dp), intent(inout) :: V1(:,:)
    real(dp), intent(in) :: V2(:,:)
    real(dp), intent(inout), optional :: A_V1(:,:)

    real(dp), allocatable :: V2_V2T(:,:)
    real(dp) :: mag
    integer :: i

    if (present(A_V1)) then
      if (size(A_V1,1) /= size(V1,1)) &
	call system_abort("called orthogonalise_col_vectors with V1 size not matching A_V1 size")
      if (size(A_V1,2) /= size(V1,2)) &
	call system_abort("called orthogonalise_col_vectors with V1 size not matching A_V1 size")
    else
      if (size(V1,1) /= size(V2,1)) &
	call system_abort("called orthogonalise_col_vectors with V1 size not matching V2 size")
    endif

    allocate(V2_V2T(size(V2,1), size(V2,1)))
    call matrix_product_sub(V2_V2T, V2, V2, m1_transpose=.false., m2_transpose=.true.)

    if (present(A_V1)) then
      call matrix_product_sub(V1, V2_V2T, A_V1, lhs_factor=1.0_dp, rhs_factor=-1.0_dp)
      do i=1, size(V1,2)
	mag = sqrt(sum(V1(:,i)*A_V1(:,i)))
	V1(:,i) = V1(:,i)/mag
	A_V1(:,i) = A_V1(:,i)/mag
      end do
    else
      call matrix_product_sub(V1, V2_V2T, V1, lhs_factor=1.0_dp, rhs_factor=-1.0_dp)
      do i=1, size(V1,2)
	V1(:,i) = V1(:,i)/sqrt(sum(V1(:,i)*V1(:,i)))
      end do
    endif

    deallocate(V2_V2T)

  end subroutine orthogonalise_col_vectors
