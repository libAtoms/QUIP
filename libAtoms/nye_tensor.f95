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

! Nye Tensor calculations from C. S. Hartley and Y. Mishin, Acta Mat. v. 53, p. 1313 (2005)

module nye_tensor_module
  use linearalgebra_module
  use atoms_module
implicit none

contains
  subroutine calc_nye_tensor(at, ref_lat, alpha, G)
    type(Atoms), intent(inout) :: at
    type(Atoms), intent(inout) :: ref_lat
    real(dp), intent(out) :: alpha(:,:,:)
    real(dp), intent(out), target, optional :: G(:,:,:)

    real(dp), pointer :: G_p(:,:,:)
    real(dp), allocatable :: Qplus(:,:,:)
    integer, allocatable :: match_list(:,:)
    integer, allocatable :: n_matches(:)

    integer, allocatable :: ref_nn_list(:)
    real(dp), allocatable :: ref_nn_vec(:,:), ref_nn_vec_normalised(:,:)
    integer, allocatable :: nn_list(:)
    real(dp), allocatable :: nn_vec(:,:), nn_vec_normalised(:,:)

    integer :: ref_n_neighbours, n_neighbours, max_n_neighbours
    integer :: i_at

    integer i_gamma, II, MM, i, j, k, m
    real(dp), allocatable :: Delta_G(:,:,:)
    real(dp), allocatable :: Delta_G_vec(:)
    real(dp) :: A(3,3,3)
    real(dp) :: eps(3,3,3)

    call calc_connect(ref_lat)
    call calc_connect(at)

    max_n_neighbours = 0
    do i_at = 1, at%N
      n_neighbours = atoms_n_neighbours(at, i_at)
      if (n_neighbours > max_n_neighbours) max_n_neighbours = n_neighbours
    end do

    if (.not.present(G)) then
      allocate(G_p(3,3,at%N))
    else
      G_p => G
    endif
    allocate(Qplus(3,max_n_neighbours,at%N))
    allocate(match_list(max_n_neighbours,at%N))
    allocate(n_matches(at%N))
    allocate(Delta_G(3,3,max_n_neighbours))
    allocate(Delta_G_vec(max_n_neighbours))
    
    call get_nn_list(ref_lat, 1, ref_n_neighbours, ref_nn_list, ref_nn_vec, ref_nn_vec_normalised)

    do i_at=1, at%N
      call get_nn_list(at, i_at, n_neighbours, nn_list, nn_vec, nn_vec_normalised)
      call find_lattice_correspondence(n_neighbours, nn_list, nn_vec, nn_vec_normalised, &
				       ref_n_neighbours, ref_nn_list, ref_nn_vec, ref_nn_vec_normalised, &
				       n_matches(i_at), match_list(:,i_at), Qplus(:,:,i_at), G_p(:,:,i_at))
      call print("atom i_at " // i_at, PRINT_NERD)
      call print("n_matches " // n_matches(i_at) // " match_list " // match_list(1:n_matches(i_at),i_at), PRINT_NERD)
      call print("Qplus", PRINT_NERD)
      call print(Qplus(1:3,1:n_matches(i_at),i_at), PRINT_NERD)
      call print("G", PRINT_NERD)
      call print(G_p(1:3,1:3,i_at), PRINT_NERD)
    end do

    eps = permutation_symbol()

    do i_at = 1, at%N
      do i_gamma = 1, n_matches(i_at)
	Delta_G(:,:,i_gamma) = G_p(:,:,match_list(i_gamma,i_at)) - G_p(:,:,i_at)
	call print("Delta_G(:,:,"//i_gamma//")", PRINT_NERD)
	call print(Delta_G(:,:,i_gamma), PRINT_NERD)
      end do
      do II = 1, 3
      do MM = 1, 3
	Delta_G_vec(1:n_matches(i_at)) = Delta_G(II,MM,1:n_matches(i_at))
	call print("II " // II // " MM " // MM, PRINT_NERD)
	call print("Delta_G_vec " // Delta_G_vec, PRINT_NERD)
	A(II,MM,1:3) = matmul(Qplus(1:3,1:n_matches(i_at),i_at), Delta_G_vec(1:n_matches(i_at)))
	call print("A " // A(II, MM, 1:3), PRINT_NERD)
      end do
      end do

      alpha(:,:,i_at) = 0.0_dp
      do j=1, 3
      do k=1, 3
      do m=1, 3
      do i=1, 3
	alpha(j,k,i_at) = alpha(j,k,i_at) - eps(j,i,m)*A(i,m,k)
      end do
      end do
      end do
      end do

      call print("alpha(:,:,"//i_at//")", PRINT_VERBOSE)
      call print(alpha(1,:,i_at), PRINT_VERBOSE)
      call print(alpha(2,:,i_at), PRINT_VERBOSE)
      call print(alpha(3,:,i_at), PRINT_VERBOSE)

    end do

    if (.not. present(G)) then
      deallocate(G_p)
    endif
    deallocate(match_list, Qplus)
    deallocate(n_matches)
    deallocate(Delta_G, Delta_G_vec)
  end subroutine calc_nye_tensor

  subroutine find_lattice_correspondence(n_neighbours, nn_list, nn_vec, nn_vec_normalised, &
      ref_n_neighbours, ref_nn_list, ref_nn_vec, ref_nn_vec_normalised, &
      n_matches, match_list, Qplus, G)
    integer, intent(in) :: n_neighbours
    integer, intent(in) :: nn_list(:)
    real(dp), intent(in) :: nn_vec(:,:), nn_vec_normalised(:,:)
    integer, intent(in) :: ref_n_neighbours
    integer, intent(in) :: ref_nn_list(:)
    real(dp), intent(in) :: ref_nn_vec(:,:), ref_nn_vec_normalised(:,:)
    integer, intent(out) :: n_matches
    integer, intent(out) :: match_list(:)
    real(dp), intent(out) :: Qplus(:,:), G(:,:)

    real(dp), allocatable :: P(:,:), Q(:,:)
    real(dp) :: QTQ_inv(3,3)

    integer :: i_beta, i_gamma, i_gammap
    integer :: match_i, t_n_matches
    integer, allocatable :: t_match_list(:)
    real(dp), allocatable :: t_match_length(:)
    real(dp) :: smallest_angle
    real(dp), allocatable :: angles(:,:)
    integer :: save_match, i_gamma_best_match(1), i_match

    !real(dp) :: max_angle_dev = 27.0_dp*PI/180.0_dp

    allocate(angles(n_neighbours, ref_n_neighbours))
    allocate(t_match_list(ref_n_neighbours))
    allocate(t_match_length(ref_n_neighbours))


    ! calculate all angles between lat and ref neighbor vectors
    do i_gamma = 1, n_neighbours
      do i_beta = 1, ref_n_neighbours
	angles(i_gamma, i_beta) = acos(sum(nn_vec_normalised(1:3,i_gamma)*ref_nn_vec_normalised(1:3,i_beta)))
      end do
    end do


    ! for each lat vector i_gamma, find matching ref vector and save in match_list(i_gamma), save 0 if more than 1 ref. vector matches
    do i_gamma = 1, n_neighbours
      smallest_angle = 1.0e38_dp
      n_matches = 0 ! ref. vector matches to this real neighbor vector
      do i_beta=1, ref_n_neighbours
	if (angles(i_gamma, i_beta) .feq. smallest_angle) then
	  n_matches = n_matches + 1
	else if (angles(i_gamma, i_beta) < smallest_angle) then
	  n_matches = 1
	  match_i = i_beta
	  smallest_angle = angles(i_gamma, i_beta)
	endif
      end do
      if (n_matches /= 1 .or. smallest_angle > 27.0_dp*PI/180.0_dp) then
	match_list(i_gamma) = 0
      else
	match_list(i_gamma) = match_i
      endif
    end do

    ! look for multiple lat vectors that match same ref vector, and pick best length match - zero other entries of match_list
    do i_gamma = 1, n_neighbours
      if (match_list(i_gamma) == 0) cycle
      t_n_matches = 0
      do i_gammap = 1, n_neighbours
	if (match_list(i_gammap) == match_list(i_gamma)) then
	  t_n_matches = t_n_matches + 1
	  t_match_list(t_n_matches) = i_gammap
	  t_match_length(t_n_matches) = norm(nn_vec(:,i_gammap))
	endif
      end do
      if (t_n_matches > 1) then
	i_gamma_best_match = minloc(abs(t_match_length-norm(ref_nn_vec(:,match_list(i_gamma)))))
	save_match = match_list(i_gamma)
	match_list(t_match_list(1:t_n_matches)) = 0
	match_list(i_gamma_best_match(1)) = save_match
      endif
    end do

    ! fill in P and Q vectors for matched pairs
    n_matches = count(match_list /= 0)
    allocate(P(n_matches,3))
    allocate(Q(n_matches,3))
    i_match = 0
    do i_gamma=1, n_neighbours
      if (match_list(i_gamma) /= 0) then
	i_match = i_match + 1
	Q(i_match,1:3) = nn_vec(1:3,i_gamma)
	P(i_match,1:3) = ref_nn_vec(1:3,match_list(i_gamma))
      endif
    end do

    ! compute Qplus and G
    QTQ_inv(1:3,1:3) = matmul(transpose(Q(1:n_matches,1:3)),Q(1:n_matches,1:3))
    call inverse(QTQ_inv)
    Qplus(1:3,1:n_matches) = matmul(QTQ_inv(1:3,1:3),transpose(Q(1:n_matches,1:3)))
    G(1:3,1:3) = matmul(Qplus(1:3,1:n_matches), P(1:n_matches,1:3))

    call print("P", PRINT_NERD)
    call print(P, PRINT_NERD)
    call print("Q", PRINT_NERD)
    call print(Q, PRINT_NERD)
    call print("Q - P", PRINT_NERD)
    call print((Q - P), PRINT_NERD)
    call print("Q.G - P", PRINT_NERD)
    call print(((Q.mult.G) - P), PRINT_NERD)
    call print("rms Q-P " // sum((Q-P)**2), PRINT_NERD)
    call print("rms Q.G-P " // sum(((Q.mult.G)-P)**2), PRINT_NERD)

    deallocate(P,Q)

    deallocate(angles)
    deallocate(t_match_list)
    deallocate(t_match_length)

    ! make output match_list(), which contains indices of atoms for each matched bond
    allocate(t_match_list(n_neighbours))
    i_gammap = 0
    do i_gamma = 1, n_neighbours
       if (match_list(i_gamma) /= 0) then
	  i_gammap = i_gammap + 1
	  t_match_list(i_gammap) = nn_list(i_gamma)
       end if
    end do
    match_list = 0
    match_list(1:n_matches) = t_match_list(1:n_matches)
    deallocate(t_match_list)
  end subroutine find_lattice_correspondence

  subroutine get_nn_list(at, i_at, n_neighbours, nn_list, nn_vec, nn_vec_normalised)
    type(Atoms), intent(inout) :: at
    integer, intent(in) :: i_at
    integer, intent(out) :: n_neighbours
    integer, allocatable, intent(inout) :: nn_list(:)
    real(dp), allocatable, intent(inout) :: nn_vec(:,:), nn_vec_normalised(:,:)

    real(dp) :: nn_v(3), nn_v_normalised(3)
    integer i_neigh

    n_neighbours = atoms_n_neighbours(at, i_at)

    if (n_neighbours > size(nn_list)) then
      if (allocated(nn_list)) deallocate(nn_list)
      allocate(nn_list(n_neighbours))
      if (allocated(nn_vec)) deallocate(nn_vec)
      allocate(nn_vec(3,n_neighbours))
      if (allocated(nn_vec_normalised)) deallocate(nn_vec_normalised)
      allocate(nn_vec_normalised(3,n_neighbours))
    endif

    do i_neigh = 1, n_neighbours
      nn_list(i_neigh) = atoms_neighbour(at, i_at, i_neigh, diff = nn_v, cosines = nn_v_normalised)
      nn_vec(:, i_neigh) = nn_v
      nn_vec_normalised(:, i_neigh) = nn_v_normalised
    end do

  end subroutine get_nn_list
      
end module nye_tensor_module
