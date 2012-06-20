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

!X
!X Ewald module 
!X
!% This module computes the Ewald sum.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
module Ewald_module

use system_module
use units_module

use functions_module

implicit none

private

integer max_n_r, max_n_q
integer n_r, n_q
integer, allocatable:: r_list(:,:), q_list(:,:)

integer:: last_N_A = 0
integer:: n_pairs
integer, allocatable:: pair_list(:,:)

public add_madelung_matrix, add_dmadelung_matrix, add_dmadelung_matrix_dr

contains

!| subroutine add_madelung_matrix
!| calculate Madelung matrix (interactions of pt charges in 3-D periodic supercell)
!|
!|I integer N_A: number of atoms
!|I real(dp) a1(3), a2(3), a3(3): periodic supercell vectors
!|I real(dp) A(3,N_A): atomic positions
!|O real(dp) madelung_matrix(N_A,N_A): Madelung matrix
!|I logical redo_lattice: if true, redo whole calculation, rather than
!|    just summing over periodic images that had significant contributions 
!|    last time
!|
!| Noam Bernstein 9/12/2001
!| From Fundamental Formulas of Physics, p. 599
!% Calculate Madelung matrix (interactions of pt charges in 3-D periodic supercell).
!% It is possible to choose between redoing the whole calculation and summing over those periodic images 
!% that had significant contributions last time by setting the logical variable redo_lattice.
subroutine add_madelung_matrix(N_A, a1_in, a2_in, a3_in, A_in, madelung_matrix, redo_lattice)
    integer, intent(in) :: N_A                           !% number of atoms (input) 
    real(dp), intent(in) :: a1_in(3), a2_in(3), a3_in(3) !% periodic supercell vectors (input)
    real(dp), intent(in) :: A_in(3,N_A)                  !% atomic positions (input)   
    real(dp), intent(inout) :: madelung_matrix(N_A,N_A)  !% Madelung matrix  (output)
    logical, intent(in), optional :: redo_lattice        !% if true, redo whole calculation, rather than just summing over periodic images that had significant contributions last time (input)

    integer i, j, range_a(3)
    integer s_a(3), s2, s3
    real(dp) b1(3), b2(3), b3(3), h(3,3)
    real(dp) vol
    real(dp) K(3), ri(3), rj(3), t, tv(3), Kmagsq, Rmag

    real(dp) :: a1(3), a2(3), a3(3)
    real(dp), allocatable :: A(:,:)

    integer rr, qq

    real(dp) epsilon, fourepsilonsq
    real(dp) recip_space_prefac, K_dep_factor
    real(dp) phase
    real(dp) conv_crit, prev_conv_crit

    integer ni, nip1, nip2
    real(dp) b1mag, b2mag, b3mag
    real(dp) a1mag, a2mag, a3mag

    logical use_this_one
    real(dp) convf

    logical u_redo_lattice

    integer i_pair

    integer nv

    allocate(A(3,N_A))
    a1 = a1_in / Bohr
    a2 = a2_in / Bohr
    a3 = a3_in / Bohr
    A = A_in / Bohr

    u_redo_lattice = .true.
    if (present(redo_lattice)) u_redo_lattice = redo_lattice

    !convf = 1.0D-10
    convf = 1.0D-14

    h(1:3,1) = a1(1:3)
    h(1:3,2) = a2(1:3)
    h(1:3,3) = a3(1:3)
    vol = calc_volume(h)
    call cross_3 (a1, a2, b3)
    call cross_3 (a2, a3, b1)
    call cross_3 (a3, a1, b2)
    b1 = 2.0_dp*PI* b1 / vol
    b2 = 2.0_dp*PI* b2 / vol
    b3 = 2.0_dp*PI* b3 / vol

    !epsilon = 0.7_dp/nn_dist

    epsilon = 1.5_dp / ( 1.0_dp/sqrt(PI) * 2.0_dp * ((3.0_dp*vol)/(4.0_dp*PI))**(1.0_dp/3.0_dp) )
    !epsilon = 2.25_dp / ( 1.0_dp/sqrt(PI) * 2.0_dp * ((3.0_dp*vol)/(4.0_dp*PI))**(1.0_dp/3.0_dp) )

    fourepsilonsq = 4.0_dp*epsilon**2

    recip_space_prefac = 4.0_dp*PI/vol

    prev_conv_crit = 0.0_dp
    conv_crit = prev_conv_crit+1.0_dp

    b1mag = sqrt(sum(b1**2))
    b2mag = sqrt(sum(b2**2))
    b3mag = sqrt(sum(b3**2))

    nv = 0
    range_a = 0

    if (u_redo_lattice) then

	if (N_A .ne. last_N_A) then
	    if (allocated(pair_list)) deallocate(pair_list)
	    call assign_atom_pairs(N_A, n_pairs, pair_list)
	    last_N_A = N_A
	end if

	max_n_r = 2000
	max_n_q = 2000
	if (allocated(r_list)) deallocate(r_list)
	if (allocated(q_list)) deallocate(q_list)
	allocate(r_list(3,max_n_r))
	allocate(q_list(3,max_n_q))
	n_r = 0
	n_q = 0


	! reciprocal space
	do while (abs(conv_crit-prev_conv_crit) .gt. 1.0e-14_dp) 

	    !!!!!
	    do ni=1,3
	    range_a(ni) = range_a(ni) + 1
	    nip1 = mod(ni,3)+1
	    nip2 = mod(ni+1,3)+1

	    s_a(ni) = -range_a(ni)
	    do s2 = -range_a(nip1), range_a(nip1)
		s_a(nip1) = s2
		do s3 = -range_a(nip2), range_a(nip2)
		    s_a(nip2) = s3

		    use_this_one = .false.

		    K = s_a(1)*b1 + s_a(2)*b2 + s_a(3)*b3
		    Kmagsq = sum(K**2)
		    K_dep_factor = recip_space_prefac*exp(-Kmagsq/fourepsilonsq)/Kmagsq

		    do i_pair = 1, n_pairs

			i = pair_list(1,i_pair)
			j = pair_list(2,i_pair)
			ri = A(1:3,i)
			rj = A(1:3,j)

			phase = sum(K*(ri-rj))
			t = Hartree*K_dep_factor*cos(phase)
			madelung_matrix(i,j) = madelung_matrix(i,j) + 2*t

			if (i /= j) madelung_matrix(j,i) = madelung_matrix(j,i) + 2*t

			if (t /= 0.0_dp .and. abs((2.0_dp*t)/madelung_matrix(i,j)) > convf) then
! workaround for buggy Cray MTA compiler
#ifndef cray_mta
			    use_this_one = .true.
#endif
			end if

		    end do

		    if (use_this_one) then
			n_q = n_q + 1
			q_list(1:3,n_q) = s_a
		    end if

		    nv = nv + 1

		end do ! s3
	    end do ! s2
	    end do ! ni

	    prev_conv_crit = conv_crit
	    !NO_PARALLEL t_conv_crit = sum(abs(madelung_matrix))
	    !NO_PARALLEL call global_sum(t_conv_crit, conv_crit)
	    conv_crit = sum(abs(madelung_matrix))
	    ! print *, "G conv ", range_a, conv_crit, prev_conv_crit

	end do ! while not conv
	! print *, "nv ", nv


	!NO_PARALLEL t_conv_crit = sum(abs(madelung_matrix))
	!NO_PARALLEL call global_sum(t_conv_crit,prev_conv_crit)
	prev_conv_crit = sum(abs(madelung_matrix))
	conv_crit = prev_conv_crit+1.0_dp

	! real space shift=0
	a1mag = sqrt(sum(a1**2))
	a2mag = sqrt(sum(a2**2))
	a3mag = sqrt(sum(a3**2))

	tv = 0.0_dp

	do i_pair = 1, n_pairs
	    i = pair_list(1,i_pair)
	    j = pair_list(2,i_pair)
	    if (i .eq. j) cycle

	    ri = A(1:3,i)
	    rj = A(1:3,j)

	    Rmag = sqrt(sum((rj+tv-ri)**2))
	    t = Hartree * erfc(epsilon*Rmag)/Rmag
	    madelung_matrix(i,j) = madelung_matrix(i,j) + t
	    if (i .ne. j) madelung_matrix(j,i) = madelung_matrix(j,i) + t
	end do


	prev_conv_crit = conv_crit
	!NO_PARALLEL t_conv_crit = sum(dabs(madelung_matrix))
	!NO_PARALLEL call global_sum(t_conv_crit, conv_crit)
	conv_crit = sum(abs(madelung_matrix))
	! print *, "R conv ", range, conv_crit, prev_conv_crit

	nv = 0
	range_a = 0

	! real space shift /= 0
	do while (abs(conv_crit-prev_conv_crit) .gt. 1.0e-14_dp) 

	    ! amag_a(1) = a1mag*(range_a(1)+1)
	    ! amag_a(2) = a2mag*(range_a(2)+1)
	    ! amag_a(3) = a3mag*(range_a(3)+1)
	    ! next_increment = minloc(amag_a)
	    ! ni = next_increment(1)

	    !!!!!
	    do ni=1,3

	    range_a(ni) = range_a(ni) + 1
	    nip1 = mod(ni,3)+1
	    nip2 = mod(ni+1,3)+1

	    s_a(ni) = -range_a(ni)
	    do s2 = -range_a(nip1), range_a(nip1)
		s_a(nip1) = s2
		do s3 = -range_a(nip2), range_a(nip2)
		    s_a(nip2) = s3

		    use_this_one = .false.

		    tv = s_a(1)*a1 + s_a(2)*a2 + s_a(3)*a3

		    do i_pair = 1, n_pairs
			i = pair_list(1,i_pair)
			j = pair_list(2,i_pair)
			ri = A(1:3,i)
			rj = A(1:3,j)

			Rmag = sqrt(sum((rj+tv-ri)**2))
			t = Hartree * erfc(epsilon*Rmag)/Rmag
			madelung_matrix(i,j) = madelung_matrix(i,j) + t
			if (i .ne. j) madelung_matrix(j,i) = madelung_matrix(j,i) + t
			if (t .ne. 0 .and. abs(t/madelung_matrix(i,j)) .gt. convf) then
			    use_this_one = .true.
			end if

		    end do

		    !NO_PARALLEL call global_or(use_this_one, g_use_this_one)
		    !NO_PARALLEL use_this_one = g_use_this_one

		    if (use_this_one) then
			n_r = n_r + 1
			r_list(1:3,n_r) = s_a
		    end if

		    nv = nv + 1
		end do
	    end do

	    s_a(ni) = range_a(ni)
	    do s2 = -range_a(nip1), range_a(nip1)
		s_a(nip1) = s2
		do s3 = -range_a(nip2), range_a(nip2)
		    s_a(nip2) = s3

		    use_this_one = .false.

		    tv = s_a(1)*a1 + s_a(2)*a2 + s_a(3)*a3
		    do i_pair = 1, n_pairs
			i = pair_list(1,i_pair)
			j = pair_list(2,i_pair)
			ri = A(1:3,i)
			rj = A(1:3,j)

			Rmag = sqrt(sum((rj+tv-ri)**2))
			t = Hartree * erfc(epsilon*Rmag)/Rmag
			madelung_matrix(i,j) = madelung_matrix(i,j) + t
			if (i .ne. j) madelung_matrix(j,i) = madelung_matrix(j,i) + t

			if (t .ne. 0 .and. abs(t/madelung_matrix(i,j)) .gt. convf) then
			    use_this_one = .true.
			end if

		    end do

		    !NO_PARALLEL call global_or(use_this_one, g_use_this_one)
		    !NO_PARALLEL use_this_one = g_use_this_one

		    if (use_this_one) then
			n_r = n_r + 1
			r_list(1:3,n_r) = s_a
		    end if

		    nv = nv + 1
		end do
	    end do

	    !!!!!
	    end do ! ni


	    prev_conv_crit = conv_crit
	    !NO_PARALLEL t_conv_crit = sum(abs(madelung_matrix))
	    !NO_PARALLEL call global_sum(t_conv_crit, conv_crit)
	    conv_crit = sum(abs(madelung_matrix))
	    ! print *, "R conv ", range, conv_crit, prev_conv_crit

	end do ! while conv_crit

	! print *, "nv ", nv
	! print *, "n_r, n_q ", n_r, n_q
    else ! don't redo lattice

	do qq=1, n_q

	    s_a = q_list(1:3,qq)

	    K = s_a(1)*b1 + s_a(2)*b2 + s_a(3)*b3
	    Kmagsq = sum(K**2)
	    K_dep_factor = recip_space_prefac*exp(-Kmagsq/fourepsilonsq)/Kmagsq
	    do i_pair = 1, n_pairs
		i = pair_list(1,i_pair)
		j = pair_list(2,i_pair)
		ri = A(1:3,i)
		rj = A(1:3,j)

		phase = sum(K*(ri-rj))
		t = Hartree * K_dep_factor*cos(phase)
		madelung_matrix(i,j) = madelung_matrix(i,j) + 2*t
		if (i .ne. j) madelung_matrix(j,i) = madelung_matrix(j,i) + 2*t
	    end do

	end do

	tv = 0.0_dp
	do i_pair = 1, n_pairs
	    i = pair_list(1,i_pair)
	    j = pair_list(2,i_pair)
	    if (i .eq. j) cycle
	    ri = A(1:3,i)
	    rj = A(1:3,j)

	    Rmag = sqrt(sum((rj+tv-ri)**2))
	    t = Hartree * erfc(epsilon*Rmag)/Rmag
	    madelung_matrix(i,j) = madelung_matrix(i,j) + t
	    if (i .ne. j) madelung_matrix(j,i) = madelung_matrix(j,i) + t
	end do

	do rr=1, n_r
	    s_a = r_list(1:3,rr)

	    tv = s_a(1)*a1 + s_a(2)*a2 + s_a(3)*a3
	    do i_pair = 1, n_pairs
		i = pair_list(1,i_pair)
		j = pair_list(2,i_pair)
		ri = A(1:3,i)
		rj = A(1:3,j)

		Rmag = sqrt(sum((rj+tv-ri)**2))
		t = Hartree * erfc(epsilon*Rmag)/Rmag
		madelung_matrix(i,j) = madelung_matrix(i,j) + t
		if (i .ne. j) madelung_matrix(j,i) = madelung_matrix(j,i) + t
	    end do
	end do

    end if

    conv_crit = sum(abs(madelung_matrix))
    ! print *, "overall conv ", conv_crit

    do i=1, N_A
	madelung_matrix(i,i) = madelung_matrix(i,i) - Hartree * 2.0_dp*epsilon/sqrt(PI)
    end do

    madelung_matrix = madelung_matrix - Hartree * PI/(vol*epsilon**2)

    deallocate(A)

end subroutine add_madelung_matrix

!| subroutine add_dmadelung_matrix
!| calculate derivative of Madelung matrix (interactions of pt charges in 
!|    3-D periodic supercell)
!|
!|I integer N_A: number of atoms
!|I real(dp) a1(3), a2(3), a3(3): periodic supercell vectors
!|I real(dp) A(3,N_A): atomic positions
!|O real(dp) dmadelung_matrix(N_A,N_A): derivative of Madelung matrix
!|
!| Noam Bernstein 9/12/2001
!% This subroutine calculate derivative of Madelung matrix (interactions of pt charges in 
!%  3-D periodic supercell).
subroutine add_dmadelung_matrix(N_A, a1_in, a2_in, a3_in, A_in, dmadelung_matrix)
    integer  :: N_A                             !% number of atoms                                   
    real(dp) :: a1_in(3), a2_in(3), a3_in(3)   !% periodic supercell vectors  
    real(dp) :: A_in(3,N_A)                    !% atomic positions                               
    real(dp) :: dmadelung_matrix(N_A,N_A,3)    !% derivative of Madelung matrix 

    !NO PARALLEL real(dp), allocatable:: t_dmadelung_matrix(:,:,:)

    integer i, j
    integer s_a(3)
    real(dp) b1(3), b2(3), b3(3), h(3,3)
    real(dp) vol
    real(dp) K(3), ri(3), rj(3), dt, tv(3), Kmagsq, Rmag, dr(3)

    real(dp) :: a1(3), a2(3), a3(3)
    real(dp), allocatable :: A(:,:)

    integer i_pair
    integer rr, qq

    real(dp) epsilon, fourepsilonsq
    real(dp) recip_space_prefac, K_dep_factor
    real(dp) phase

    allocate(A(3,N_A))
    a1 = a1_in / Bohr
    a2 = a2_in / Bohr
    a3 = a3_in / Bohr
    A = A_in / Bohr

    h(1:3,1) = a1(1:3)
    h(1:3,2) = a2(1:3)
    h(1:3,3) = a3(1:3)
    vol = calc_volume(h)
    call cross_3 (a1, a2, b3)
    call cross_3 (a2, a3, b1)
    call cross_3 (a3, a1, b2)
    b1 = 2.0_dp*PI* b1 / vol
    b2 = 2.0_dp*PI* b2 / vol
    b3 = 2.0_dp*PI* b3 / vol

    !epsilon = 0.7_dp/nn_dist

    epsilon = 1.5_dp / ( 1.0_dp/sqrt(PI) * 2.0_dp * ((3.0_dp*vol)/(4.0_dp*PI))**(1.0_dp/3.0_dp) )
    !epsilon = 2.25_dp / ( 1.0_dp/sqrt(PI) * 2.0_dp * ((3.0_dp*vol)/(4.0_dp*PI))**(1.0_dp/3.0_dp) )

    fourepsilonsq = 4.0_dp*epsilon**2

    recip_space_prefac = 4.0_dp*PI/vol

    do qq=1, n_q

	s_a = q_list(1:3,qq)

	K = s_a(1)*b1 + s_a(2)*b2 + s_a(3)*b3
	Kmagsq = sum(K**2)
	K_dep_factor = recip_space_prefac*exp(-Kmagsq/fourepsilonsq)/Kmagsq
	do i_pair=1, n_pairs
	    i = pair_list(1,i_pair)
	    j = pair_list(2,i_pair)
	    ri = A(1:3,i)
	    rj = A(1:3,j)

	    phase = sum(K*(ri-rj))
	    dt = Hartree/Bohr * K_dep_factor*(-sin(phase))
	    dmadelung_matrix(i,j,1:3) = dmadelung_matrix(i,j,1:3) + 2.0_dp*dt*K(1:3)
	    if (i .ne. j) dmadelung_matrix(j,i,1:3) = dmadelung_matrix(j,i,1:3) - 2.0_dp*dt*K(1:3)
	end do

    end do

    tv = 0.0_dp
    do i_pair=1, n_pairs
	i = pair_list(1,i_pair)
	j = pair_list(2,i_pair)
	if (i .eq. j) cycle
	ri = A(1:3,i)
	rj = A(1:3,j)

	dr = ri-tv-rj
	Rmag = sqrt(sum(dr**2))
	dt = Hartree/Bohr * ( (-2.0_dp/sqrt(PI))*Rmag*epsilon*exp(-(epsilon*Rmag)**2) - erfc(epsilon*Rmag) ) / Rmag**2
	dmadelung_matrix(i,j,1:3) = dmadelung_matrix(i,j,1:3) + dt*dr/Rmag
	if (i .ne. j) dmadelung_matrix(j,i,1:3) = dmadelung_matrix(j,i,1:3) - dt*dr/Rmag

    end do

    do rr=1, n_r
	s_a = r_list(1:3,rr)

	tv = s_a(1)*a1 + s_a(2)*a2 + s_a(3)*a3
	do i_pair=1, n_pairs
	    i = pair_list(1,i_pair)
	    j = pair_list(2,i_pair)
	    ri = A(1:3,i)
	    rj = A(1:3,j)

	    dr = ri-tv-rj
	    Rmag = sqrt(sum(dr**2))
	    dt = Hartree/Bohr * ( (-2.0_dp/sqrt(PI))*Rmag*epsilon*exp(-(epsilon*Rmag)**2) - erfc(epsilon*Rmag) ) / Rmag**2
	    dmadelung_matrix(i,j,1:3) = dmadelung_matrix(i,j,1:3) + dt*dr/Rmag
	    if (i .ne. j) dmadelung_matrix(j,i,1:3) = dmadelung_matrix(j,i,1:3) - dt*dr/Rmag

	end do
    end do

    !NO_PARALLEL allocate(t_dmadelung_matrix(N_A,N_A,3))
    !NO_PARALLEL t_dmadelung_matrix = dmadelung_matrix
    !NO_PARALLEL dmadelung_matrix = 0.0_dp
    !NO_PARALLEL call global_sum(t_dmadelung_matrix, dmadelung_matrix)
    !NO_PARALLEL deallocate(t_dmadelung_matrix)

    deallocate(A)

end subroutine add_dmadelung_matrix

!| subroutine add_dmadelung_matrix_dr
!% Calculate derivative of Madelung matrix (interactions of pt charges in 
!%    3-D periodic supercell) dotted with atomic positions.
!% From Nielsen and Martin, Phys. Rev. B {\bf 32}, (1985).
!|
!|I integer N_A: number of atoms
!|I real(dp) a1(3), a2(3), a3(3): periodic supercell vectors
!|I real(dp) A(3,N_A): atomic positions
!|O real(dp) dmadelung_matrix_dr(N_A,N_A): derivative of Madelung matrix
!|    dotted with atomic positions
!|
!| Noam Bernstein 9/12/2001
!| From Nielsen and Martin, PRB _32_, 1985.

subroutine add_dmadelung_matrix_dr(N_A, a1_in, a2_in, a3_in, A_in, component, dmadelung_matrix_dr)
    integer  :: N_A                             !% number of atoms (input) 
    real(dp) :: a1_in(3), a2_in(3), a3_in(3)    !% periodic supercell vectors (input)
    real(dp) :: A_in(3,N_A)                     !% atomic positions (input)
    integer  :: component
    real(dp) :: dmadelung_matrix_dr(N_A,N_A,3)  !% derivative of Madelung matrix dotted with atomic positions (output)

    integer i, j
    integer s_a(3)
    real(dp) b1(3), b2(3), b3(3), h(3,3)
    real(dp) vol
    real(dp) K(3), ri(3), rj(3), dt, tv(3), Kmagsq, Rmag, dr(3)
    real(dp) t_ij(3)

    real(dp) :: a1(3), a2(3), a3(3)
    real(dp), allocatable :: A(:,:)

    integer i_pair
    integer rr, qq

    real(dp) epsilon, fourepsilonsq
    real(dp) recip_space_prefac, K_dep_factor
    real(dp) phase

    allocate(A(3,N_A))
    a1 = a1_in / Bohr
    a2 = a2_in / Bohr
    a3 = a3_in / Bohr
    A = A_in / Bohr

    h(1:3,1) = a1(1:3)
    h(1:3,2) = a2(1:3)
    h(1:3,3) = a3(1:3)
    vol = calc_volume(h)
    call cross_3 (a1, a2, b3)
    call cross_3 (a2, a3, b1)
    call cross_3 (a3, a1, b2)
    b1 = 2.0_dp*PI* b1 / vol
    b2 = 2.0_dp*PI* b2 / vol
    b3 = 2.0_dp*PI* b3 / vol

    !epsilon = 0.7_dp/nn_dist

    epsilon = 1.5_dp / ( 1.0_dp/dsqrt(PI) * 2.0_dp * ((3.0_dp*vol)/(4.0_dp*PI))**(1.0_dp/3.0_dp) )
    !epsilon = 2.25_dp / ( 1.0_dp/dsqrt(PI) * 2.0_dp * ((3.0_dp*vol)/(4.0_dp*PI))**(1.0_dp/3.0_dp) )

    fourepsilonsq = 4.0_dp*epsilon**2

    recip_space_prefac = PI/(vol*epsilon**2)

    do qq=1, n_q

	s_a = q_list(1:3,qq)

	K = s_a(1)*b1 + s_a(2)*b2 + s_a(3)*b3
	Kmagsq = sum(K**2)

	K_dep_factor = recip_space_prefac*dexp(-Kmagsq/fourepsilonsq)/(Kmagsq/fourepsilonsq)

	do i_pair=1, n_pairs
	    i = pair_list(1,i_pair)
	    j = pair_list(2,i_pair)
	    ri = A(1:3,i)
	    rj = A(1:3,j)

	    dr = ri-rj
	    phase = sum(K*(ri-rj))

	    dt = Hartree * K_dep_factor*2.0_dp*dcos(phase)

	    if (component .eq. 0) then
		t_ij(1:3) = dt * (2.0_dp*K(1:3)*K(1:3)/Kmagsq * (Kmagsq/fourepsilonsq + 1.0_dp) - 1.0_dp)
		dmadelung_matrix_dr(i,j,1) = dmadelung_matrix_dr(i,j,1) - sum(t_ij(1:3))
		if (i .ne. j) dmadelung_matrix_dr(j,i,1) = dmadelung_matrix_dr(j,i,1) - sum(t_ij(1:3))
	    else
		t_ij(1:3) = dt * (2.0_dp*K(component)*K(1:3)/Kmagsq * (Kmagsq/fourepsilonsq + 1.0_dp))
		t_ij(component) = t_ij(component) - dt
		dmadelung_matrix_dr(i,j,1:3) = dmadelung_matrix_dr(i,j,1:3) - t_ij(1:3)
		if (i .ne. j) dmadelung_matrix_dr(j,i,1:3) = dmadelung_matrix_dr(j,i,1:3) - t_ij(1:3)
	    end if
	end do

    end do

    tv = 0.0_dp
    do i_pair=1, n_pairs
	i = pair_list(1,i_pair)
	j = pair_list(2,i_pair)
	if (i .eq. j) cycle
	ri = A(1:3,i)
	rj = A(1:3,j)

	dr = ri-tv-rj
	Rmag = sqrt(sum(dr**2))
	if (component .eq. 0) then
	    t_ij(1:3) = Hartree * 0.5_dp*epsilon*2.0_dp*Hprime(epsilon*Rmag)*dr(1:3)*dr(1:3)/Rmag**2
	    dmadelung_matrix_dr(i,j,1) = dmadelung_matrix_dr(i,j,1) - sum(t_ij(1:3))
	    if (i .ne. j) dmadelung_matrix_dr(j,i,1) = dmadelung_matrix_dr(j,i,1) - sum(t_ij(1:3))
	else
	    t_ij(1:3) = Hartree * 0.5_dp*epsilon*2.0_dp*Hprime(epsilon*Rmag)*dr(component)*dr(1:3)/Rmag**2
	    dmadelung_matrix_dr(i,j,1:3) = dmadelung_matrix_dr(i,j,1:3) - t_ij(1:3)
	    if (i .ne. j) dmadelung_matrix_dr(j,i,1:3) = dmadelung_matrix_dr(j,i,1:3) - t_ij(1:3)
	end if
    end do

    do rr=1, n_r
	s_a = r_list(1:3,rr)

	tv = s_a(1)*a1 + s_a(2)*a2 + s_a(3)*a3
	do i_pair=1, n_pairs
	    i = pair_list(1,i_pair)
	    j = pair_list(2,i_pair)
	    ri = A(1:3,i)
	    rj = A(1:3,j)

	    dr = ri-tv-rj
	    Rmag = sqrt(sum(dr**2))
	    if (component .eq. 0) then
		t_ij(1:3) = Hartree * 0.5_dp*epsilon*2.0_dp*Hprime(epsilon*Rmag)*dr(1:3)*dr(1:3)/Rmag**2
		dmadelung_matrix_dr(i,j,1) = dmadelung_matrix_dr(i,j,1) - sum(t_ij(1:3))
		if (i .ne. j) dmadelung_matrix_dr(j,i,1) = dmadelung_matrix_dr(j,i,1) - sum(t_ij(1:3))
	    else
		t_ij(1:3) = Hartree * 0.5_dp*epsilon*2.0_dp*Hprime(epsilon*Rmag)*dr(component)*dr(1:3)/Rmag**2
		dmadelung_matrix_dr(i,j,1:3) = dmadelung_matrix_dr(i,j,1:3) - t_ij(1:3)
		if (i .ne. j) dmadelung_matrix_dr(j,i,1:3) = dmadelung_matrix_dr(j,i,1:3) - t_ij(1:3)
	    end if
	end do
    end do

    if (component .eq. 0) then
	dmadelung_matrix_dr(1:N_A,1:N_A,1) = dmadelung_matrix_dr(1:N_A,1:N_A,1) - &
	    Hartree * 3.0_dp*PI/(vol*epsilon**2)
    else
	dmadelung_matrix_dr(1:N_A,1:N_A,1:3) = dmadelung_matrix_dr(1:N_A,1:N_A,1:3) - &
	    Hartree * PI/(vol*epsilon**2)
    end if

!NO_PARALLEL    if (component .eq. 0) then
!NO_PARALLEL	allocate(t_dmadelung_matrix_dr(N_A,N_A,1))
!NO_PARALLEL	t_dmadelung_matrix_dr(1:N_A,1:N_A,1) = dmadelung_matrix_dr(1:N_A,1:N_A,1)
!NO_PARALLEL	dmadelung_matrix_dr(1:N_A,1:N_A,1) = 0.0_dp
!NO_PARALLEL	call global_sum(t_dmadelung_matrix_dr, dmadelung_matrix_dr)
!NO_PARALLEL	deallocate(t_dmadelung_matrix_dr)
!NO_PARALLEL    else
!NO_PARALLEL	allocate(t_dmadelung_matrix_dr(N_A,N_A,3))
!NO_PARALLEL	t_dmadelung_matrix_dr = dmadelung_matrix_dr
!NO_PARALLEL	dmadelung_matrix_dr = 0.0_dp
!NO_PARALLEL	call global_sum(t_dmadelung_matrix_dr, dmadelung_matrix_dr)
!NO_PARALLEL	deallocate(t_dmadelung_matrix_dr)
!NO_PARALLEL    endif

    deallocate(A)

end subroutine add_dmadelung_matrix_dr

!| subroutine cross_3
!% Calculate cross product of 2 3-vectors.  a(3), b(3): vectors to be crossed, c(3): cross product.
!|
!|I real(dp) a(3), b(3): vectors to be crossed
!|O real(dp) c(3): cross product
!|
!| Noam Bernstein 9/12/2001

subroutine cross_3(a, b, c)
    real(dp) a(3), b(3), c(3)

    c(1) = a(2)*b(3)-a(3)*b(2)
    c(2) = a(3)*b(1)-a(1)*b(3)
    c(3) = a(1)*b(2)-a(2)*b(1)

end subroutine cross_3

!| real(dp) function calc_volume
!| calculate volume of periodic supercell
!| returns volume of supercell
!|
!|I real(dp) a(3,3): matrix of peridic supercell vectors
!|
!| Noam Bernstein 9/12/2001

real(dp) function calc_volume(mat)
    real(dp) mat(3,3)

    calc_volume = abs(det_3_by_3(mat))
end function calc_volume

!| real(dp) function det_3_by_3
!| calculate determinant of 3x3 matrix
!| returns determinant of matrix
!|
!|I real(dp) a(3,3): matrix
!|
!| Noam Bernstein 9/12/2001

real(dp) function det_3_by_3(a)
    real(dp) a(3,3)

    det_3_by_3 = -a(1,3)*a(2,2)*a(3,1) + a(1,2)*a(2,3)*a(3,1) + a(1,3)*a(2,1)*a(3,2) - &
	    a(1,1)*a(2,3)*a(3,2) - a(1,2)*a(2,1)*a(3,3) + a(1,1)*a(2,2)*a(3,3)
end function det_3_by_3

!| real(dp) function Hprime
!| calculate Hprime function for dMadelung_matrix/dr . r
!| returns value of Hprime(x)
!|
!|I real(dp) x: argument
!|
!| Noam Bernstein 9/12/2001

real(dp) function Hprime(x)
    real(dp) x

    Hprime =  (-2.0_dp/sqrt(PI))*exp(-x**2) - erfc(x)/x

end function Hprime

subroutine assign_atom_pairs(N, n_pairs, pair_list)
  integer, intent(in) :: N
  integer, intent(out) :: n_pairs
  integer, allocatable, intent(inout) :: pair_list(:,:)

  integer cur_pair_i, i, j

  n_pairs = N*(N+1)/2

  if (allocated(pair_list)) deallocate(pair_list)
  allocate(pair_list(2,n_pairs))

  cur_pair_i = 1
  do i=1, N
    do j=i, N
      pair_list(1,cur_pair_i) = i
      pair_list(2,cur_pair_i) = j
      cur_pair_i = cur_pair_i + 1
    end do
  end do

end subroutine assign_atom_pairs

end module Ewald_module
