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

! convert from NRL-TB web-site input files to QUIP XML format
module NRL_TB_to_xml_module

use System_module
use TB_Common_module
use TBModel_NRL_TB_defs_module
use fox_wxml

implicit none

contains 

subroutine do_nrl_tb_params_to_xml(fname)
  character(len=*), intent(in) :: fname

  logical is_orthogonal, is_magnetic, has_pair_repulsion, overlap_zero_limit, force_harrison_signs
  character(len=1024) header, label
  integer param_style
  integer n_types
  integer, allocatable :: atomic_num(:), n_orbs(:), n_elecs(:), n_orb_sets(:), orb_set_type(:,:)
  double precision, allocatable :: atomic_mass(:), r_cut(:,:), screen_l(:,:)
  double precision, allocatable :: pair_rep_inner(:,:), pair_rep_outer(:,:)
  double precision :: cutoff
  double precision, allocatable :: lambda_sq(:,:), abcd(:,:,:,:,:)
  double precision, allocatable :: H_coeff(:,:,:,:,:), S_coeff(:,:,:,:,:)

  integer n_elec_s, n_elec_p, n_elec_d
  double precision t_d
  integer i, j, stat
  integer i_mag, n_mag
  character(len=1024) line

  read *, line

  header = trim(line)
  
  if (header(1:1) == 'N') then
    is_orthogonal = .false.
  else if (header(1:1) == 'O') then
    is_orthogonal = .true.
  else
    call system_abort ("confused by first character of header")
  endif

  if (header(2:2) == 'N') then
    is_magnetic = .false.
  else if (header(2:2) == 'M') then
    is_magnetic = .true.
  else
    call system_abort ("confused by first character of header")
  endif

  if (header(3:3) .ne. "0" .or. header(4:4) .ne. "0" .or. header(5:5) .ne. "0") then
    call system_abort( "Unknown parameter format "//trim(header) )
  end if
  if (header(6:6) .eq. "0") then
    has_pair_repulsion = .false.
  else if (header(6:6) .eq. "1") then
    has_pair_repulsion = .true.
  else
    call system_abort ("Unknown parameter format "//trim(header))
  endif
  if (header(7:7) .eq. "0") then
    param_style = 0
  else if (header(7:7) .eq. "1") then
    param_style = 1
  else if (header(7:7) .eq. "2") then
    param_style = 2
  else if (header(7:7) .eq. "3") then
    param_style = 3
  else
    call system_abort ("Unknown parameter format "// trim(header) )
  end if

  select case (param_style)
    case (0)
      overlap_zero_limit = .false.
      force_harrison_signs = .false.
    case (1)
      overlap_zero_limit = .true.
      force_harrison_signs = .false.
    case (2)
      overlap_zero_limit = .false.
      force_harrison_signs = .true.
    case (3)
      overlap_zero_limit = .true.
      force_harrison_signs = .true.
    case default
  end select

  read *, line
  label = trim(line)

!!!!!! read n_types !!!!!!
  read *, n_types

  allocate(atomic_num(n_types))
  allocate(atomic_mass(n_types))
  allocate(n_orbs(n_types))
  allocate(n_elecs(n_types))
  allocate(n_orb_sets(n_types))
  allocate(orb_set_type(n_types,max_n_orb_sets))

  allocate(r_cut(n_types,n_types))
  allocate(screen_l(n_types,n_types))

!!!!!! read non matrix element stuff !!!!!!
  do i=1, n_types
  do j=i, 1, -1
    read *, r_cut(i,j), screen_l(i,j)
    if (i .ne. j) then
      r_cut(j,i) = r_cut(i,j)
      screen_l(j,i) = screen_l(i,j)
    endif
  end do
  end do

  if (has_pair_repulsion) then
    allocate(pair_rep_inner(n_types,n_types))
    allocate(pair_rep_outer(n_types,n_types))
    do i=1, n_types
    do j=i, 1, -1
      read *, pair_rep_inner(i,j), pair_rep_outer(i,j)
      if (i .ne. j) then
	pair_rep_inner(i,j) = pair_rep_inner(j,i)
	pair_rep_outer(i,j) = pair_rep_outer(j,i)
      endif
    end do
    end do
  endif

  cutoff = maxval(r_cut(:,:))
  do i=1, n_types
      read *, n_orbs(i)

      atomic_num(i) = 0
      read (*, *, iostat=stat) atomic_mass(i), atomic_num(i)

      read *, n_elec_s, n_elec_p, n_elec_d
      select case (n_orbs(i))
	case (1)
	  n_orb_sets(i) = 1
	  orb_set_type(i,1) = ORB_S
	  n_elecs(i) = n_elec_s
	case (4)
	  n_orb_sets(i) = 2
	  orb_set_type(i,1) = ORB_S
	  orb_set_type(i,2) = ORB_P
	  n_elecs(i) = n_elec_s + n_elec_p
	case (6)
	  n_orb_sets(i) = 2
	  orb_set_type(i,1) = ORB_S
	  orb_set_type(i,2) = ORB_D
	  n_elecs(i) = n_elec_s + n_elec_d
	case (9)
	  n_orb_sets(i) = 3
	  orb_set_type(i,1) = ORB_S
	  orb_set_type(i,2) = ORB_P
	  orb_set_type(i,3) = ORB_D
	  n_elecs(i) = n_elec_s + n_elec_p + n_elec_d
    end select
  end do

!!  allocate(type_of_atomic_num(maxval(atomic_num(:))))
!!  type_of_atomic_num(:) = 0
!!  do i=1, n_types
!!    type_of_atomic_num(atomic_num(i)) = i
!!  end do

!!!!!!!!! ON-SITE !!!!!!!!!!!!
  if (is_magnetic) then
    n_mag = 2
  else
    n_mag = 1
  endif

  allocate(lambda_sq(n_types,n_mag))
  allocate(abcd(4,3,n_types,n_types,n_mag))
  allocate(H_coeff(4,10,n_types,n_types,n_mag))
  allocate(S_coeff(4,10,n_types,n_types,n_mag))

  do i_mag=1, n_mag

    do i=1, n_types
    j = i
	read *, lambda_sq(i,i_mag)

	read *, abcd(ABCD_A,SPD_S,i,i,i_mag)
	read *, abcd(ABCD_B,SPD_S,i,i,i_mag)
	read *, abcd(ABCD_C,SPD_S,i,i,i_mag)
	read *, abcd(ABCD_D,SPD_S,i,i,i_mag)

	read *, abcd(ABCD_A,SPD_P,i,i,i_mag)
	read *, abcd(ABCD_B,SPD_P,i,i,i_mag)
	read *, abcd(ABCD_C,SPD_P,i,i,i_mag)
	read *, abcd(ABCD_D,SPD_P,i,i,i_mag)

	read *, abcd(ABCD_A,SPD_D,i,i,i_mag)
	read *, abcd(ABCD_B,SPD_D,i,i,i_mag)
	read *, abcd(ABCD_C,SPD_D,i,i,i_mag)
	read *, abcd(ABCD_D,SPD_D,i,i,i_mag)

	read *, t_d
	if (any(orb_set_type(i,1:n_orb_sets(i)) == ORB_D) .and. abcd(ABCD_A,SPD_D,i,i,i_mag) /= t_d) &
	  call system_abort ("NRL-TB No support for E_g T_2g split" // abcd(ABCD_A,SPD_D,i,i,i_mag) // " " // t_d)
	read *, t_d
	if (any(orb_set_type(i,1:n_orb_sets(i)) == ORB_D) .and. abcd(ABCD_B,SPD_D,i,i,i_mag) /= t_d) &
	  call system_abort ("NRL-TB No support for E_g T_2g split" // abcd(ABCD_B,SPD_D,i,i,i_mag) // " " // t_d)
	read *, t_d
	if (any(orb_set_type(i,1:n_orb_sets(i)) == ORB_D) .and. abcd(ABCD_C,SPD_D,i,i,i_mag) /= t_d) &
	  call system_abort ("NRL-TB No support for E_g T_2g split" // abcd(ABCD_C,SPD_D,i,i,i_mag) // " " // t_d)
	read *, t_d
	if (any(orb_set_type(i,1:n_orb_sets(i)) == ORB_D) .and. abcd(ABCD_D,SPD_D,i,i,i_mag) /= t_d) &
	  call system_abort ("NRL-TB No support for E_g T_2g split" // abcd(ABCD_D,SPD_D,i,i,i_mag) // " " // t_d)
    end do

    do i=2, n_types
      do j=1, i-1

	abcd(ABCD_A,SPD_S,i,j,i_mag) = 0.0D0
	read *, abcd(ABCD_B,SPD_S,i,j,i_mag)
	read *, abcd(ABCD_C,SPD_S,i,j,i_mag)
	read *, abcd(ABCD_D,SPD_S,i,j,i_mag)

	abcd(ABCD_A,SPD_P,i,j,i_mag) = 0.0D0
	read *, abcd(ABCD_B,SPD_P,i,j,i_mag)
	read *, abcd(ABCD_C,SPD_P,i,j,i_mag)
	read *, abcd(ABCD_D,SPD_P,i,j,i_mag)

	abcd(ABCD_A,SPD_D,i,j,i_mag) = 0.0D0
	read *, abcd(ABCD_B,SPD_D,i,j,i_mag)
	read *, abcd(ABCD_C,SPD_D,i,j,i_mag)
	read *, abcd(ABCD_D,SPD_D,i,j,i_mag)

	read *, t_d
	if (abcd(ABCD_B,SPD_D,i,j,i_mag) /= t_d) &
	  call system_abort ("NRL-TB No support for E_g T_2g split" // abcd(ABCD_B,SPD_D,i,j,i_mag) // " " // t_d)
	read *, t_d
	if (abcd(ABCD_C,SPD_D,i,j,i_mag) /= t_d) &
	  call system_abort ("NRL-TB No support for E_g T_2g split" // abcd(ABCD_C,SPD_D,i,j,i_mag) // " " // t_d)
	read *, t_d
	if (abcd(ABCD_D,SPD_D,i,j,i_mag) /= t_d) &
	  call system_abort ("NRL-TB No support for E_g T_2g split" // abcd(ABCD_D,SPD_D,i,j,i_mag) // " " // t_d)

	abcd(ABCD_A,SPD_S,j,i,i_mag) = 0.0D0
	read *, abcd(ABCD_B,SPD_S,j,i,i_mag)
	read *, abcd(ABCD_C,SPD_S,j,i,i_mag)
	read *, abcd(ABCD_D,SPD_S,j,i,i_mag)

	abcd(ABCD_A,SPD_P,j,i,i_mag) = 0.0D0
	read *, abcd(ABCD_B,SPD_P,j,i,i_mag)
	read *, abcd(ABCD_C,SPD_P,j,i,i_mag)
	read *, abcd(ABCD_D,SPD_P,j,i,i_mag)

	abcd(ABCD_A,SPD_D,j,i,i_mag) = 0.0D0
	read *, abcd(ABCD_B,SPD_D,j,i,i_mag)
	read *, abcd(ABCD_C,SPD_D,j,i,i_mag)
	read *, abcd(ABCD_D,SPD_D,j,i,i_mag)

	read *, t_d
	if (abcd(ABCD_B,SPD_D,j,i,i_mag) /= t_d) &
	  call system_abort ("NRL-TB No support for E_g T_2g split" // abcd(ABCD_B,SPD_D,j,i,i_mag) // " " // t_d)
	read *, t_d
	if (abcd(ABCD_C,SPD_D,j,i,i_mag) /= t_d) &
	  call system_abort ("NRL-TB No support for E_g T_2g split" // abcd(ABCD_C,SPD_D,j,i,i_mag) // " " // t_d)
	read *, t_d
	if (abcd(ABCD_D,SPD_D,j,i,i_mag) /= t_d) &
	  call system_abort ("NRL-TB No support for E_g T_2g split" // abcd(ABCD_D,SPD_D,j,i,i_mag) // " " // t_d)


      end do
    end do

  !!!!!!!!! HAMILTONIAN !!!!!!!!!!!!

    do i=1, n_types
    j=i
      read *, H_coeff(MC_E,SK_SSS,j,i,i_mag)
      read *, H_coeff(MC_F,SK_SSS,j,i,i_mag)
      read *, H_coeff(MC_FB,SK_SSS,j,i,i_mag)
      read *, H_coeff(MC_G_SQ,SK_SSS,j,i,i_mag)
      read *, H_coeff(MC_E,SK_SPS,j,i,i_mag)
      read *, H_coeff(MC_F,SK_SPS,j,i,i_mag)
      read *, H_coeff(MC_FB,SK_SPS,j,i,i_mag)
      read *, H_coeff(MC_G_SQ,SK_SPS,j,i,i_mag)
      read *, H_coeff(MC_E,SK_PPS,j,i,i_mag)
      read *, H_coeff(MC_F,SK_PPS,j,i,i_mag)
      read *, H_coeff(MC_FB,SK_PPS,j,i,i_mag)
      read *, H_coeff(MC_G_SQ,SK_PPS,j,i,i_mag)
      read *, H_coeff(MC_E,SK_PPP,j,i,i_mag)
      read *, H_coeff(MC_F,SK_PPP,j,i,i_mag)
      read *, H_coeff(MC_FB,SK_PPP,j,i,i_mag)
      read *, H_coeff(MC_G_SQ,SK_PPP,j,i,i_mag)
      read *, H_coeff(MC_E,SK_SDS,j,i,i_mag)
      read *, H_coeff(MC_F,SK_SDS,j,i,i_mag)
      read *, H_coeff(MC_FB,SK_SDS,j,i,i_mag)
      read *, H_coeff(MC_G_SQ,SK_SDS,j,i,i_mag)
      read *, H_coeff(MC_E,SK_PDS,j,i,i_mag)
      read *, H_coeff(MC_F,SK_PDS,j,i,i_mag)
      read *, H_coeff(MC_FB,SK_PDS,j,i,i_mag)
      read *, H_coeff(MC_G_SQ,SK_PDS,j,i,i_mag)
      read *, H_coeff(MC_E,SK_PDP,j,i,i_mag)
      read *, H_coeff(MC_F,SK_PDP,j,i,i_mag)
      read *, H_coeff(MC_FB,SK_PDP,j,i,i_mag)
      read *, H_coeff(MC_G_SQ,SK_PDP,j,i,i_mag)
      read *, H_coeff(MC_E,SK_DDS,j,i,i_mag)
      read *, H_coeff(MC_F,SK_DDS,j,i,i_mag)
      read *, H_coeff(MC_FB,SK_DDS,j,i,i_mag)
      read *, H_coeff(MC_G_SQ,SK_DDS,j,i,i_mag)
      read *, H_coeff(MC_E,SK_DDP,j,i,i_mag)
      read *, H_coeff(MC_F,SK_DDP,j,i,i_mag)
      read *, H_coeff(MC_FB,SK_DDP,j,i,i_mag)
      read *, H_coeff(MC_G_SQ,SK_DDP,j,i,i_mag)
      read *, H_coeff(MC_E,SK_DDD,j,i,i_mag)
      read *, H_coeff(MC_F,SK_DDD,j,i,i_mag)
      read *, H_coeff(MC_FB,SK_DDD,j,i,i_mag)
      read *, H_coeff(MC_G_SQ,SK_DDD,j,i,i_mag)
    end do

    do i=2, n_types
      do j=1, i-1
	read *, H_coeff(MC_E,SK_SSS,j,i,i_mag)
	read *, H_coeff(MC_F,SK_SSS,j,i,i_mag)
	read *, H_coeff(MC_FB,SK_SSS,j,i,i_mag)
	read *, H_coeff(MC_G_SQ,SK_SSS,j,i,i_mag)
	read *, H_coeff(MC_E,SK_SPS,j,i,i_mag)
	read *, H_coeff(MC_F,SK_SPS,j,i,i_mag)
	read *, H_coeff(MC_FB,SK_SPS,j,i,i_mag)
	read *, H_coeff(MC_G_SQ,SK_SPS,j,i,i_mag)
	read *, H_coeff(MC_E,SK_PPS,j,i,i_mag)
	read *, H_coeff(MC_F,SK_PPS,j,i,i_mag)
	read *, H_coeff(MC_FB,SK_PPS,j,i,i_mag)
	read *, H_coeff(MC_G_SQ,SK_PPS,j,i,i_mag)
	read *, H_coeff(MC_E,SK_PPP,j,i,i_mag)
	read *, H_coeff(MC_F,SK_PPP,j,i,i_mag)
	read *, H_coeff(MC_FB,SK_PPP,j,i,i_mag)
	read *, H_coeff(MC_G_SQ,SK_PPP,j,i,i_mag)
	read *, H_coeff(MC_E,SK_SDS,j,i,i_mag)
	read *, H_coeff(MC_F,SK_SDS,j,i,i_mag)
	read *, H_coeff(MC_FB,SK_SDS,j,i,i_mag)
	read *, H_coeff(MC_G_SQ,SK_SDS,j,i,i_mag)
	read *, H_coeff(MC_E,SK_PDS,j,i,i_mag)
	read *, H_coeff(MC_F,SK_PDS,j,i,i_mag)
	read *, H_coeff(MC_FB,SK_PDS,j,i,i_mag)
	read *, H_coeff(MC_G_SQ,SK_PDS,j,i,i_mag)
	read *, H_coeff(MC_E,SK_PDP,j,i,i_mag)
	read *, H_coeff(MC_F,SK_PDP,j,i,i_mag)
	read *, H_coeff(MC_FB,SK_PDP,j,i,i_mag)
	read *, H_coeff(MC_G_SQ,SK_PDP,j,i,i_mag)
	read *, H_coeff(MC_E,SK_DDS,j,i,i_mag)
	read *, H_coeff(MC_F,SK_DDS,j,i,i_mag)
	read *, H_coeff(MC_FB,SK_DDS,j,i,i_mag)
	read *, H_coeff(MC_G_SQ,SK_DDS,j,i,i_mag)
	read *, H_coeff(MC_E,SK_DDP,j,i,i_mag)
	read *, H_coeff(MC_F,SK_DDP,j,i,i_mag)
	read *, H_coeff(MC_FB,SK_DDP,j,i,i_mag)
	read *, H_coeff(MC_G_SQ,SK_DDP,j,i,i_mag)
	read *, H_coeff(MC_E,SK_DDD,j,i,i_mag)
	read *, H_coeff(MC_F,SK_DDD,j,i,i_mag)
	read *, H_coeff(MC_FB,SK_DDD,j,i,i_mag)
	read *, H_coeff(MC_G_SQ,SK_DDD,j,i,i_mag)

	H_coeff(MC_E,SK_SSS,i,j,i_mag) = H_coeff(MC_E,SK_SSS,j,i,i_mag)
	H_coeff(MC_F,SK_SSS,i,j,i_mag) = H_coeff(MC_F,SK_SSS,j,i,i_mag)
	H_coeff(MC_FB,SK_SSS,i,j,i_mag) = H_coeff(MC_FB,SK_SSS,j,i,i_mag)
	H_coeff(MC_G_SQ,SK_SSS,i,j,i_mag) = H_coeff(MC_G_SQ,SK_SSS,j,i,i_mag)

	read *, H_coeff(MC_E,SK_SPS,i,j,i_mag)
	read *, H_coeff(MC_F,SK_SPS,i,j,i_mag)
	read *, H_coeff(MC_FB,SK_SPS,i,j,i_mag)
	read *, H_coeff(MC_G_SQ,SK_SPS,i,j,i_mag)
	H_coeff(MC_E,SK_SPS,i,j,i_mag) = - H_coeff(MC_E,SK_SPS,i,j,i_mag)
	H_coeff(MC_F,SK_SPS,i,j,i_mag) = - H_coeff(MC_F,SK_SPS,i,j,i_mag)
	H_coeff(MC_FB,SK_SPS,i,j,i_mag) = - H_coeff(MC_FB,SK_SPS,i,j,i_mag)

	H_coeff(MC_E,SK_PPS,i,j,i_mag) = H_coeff(MC_E,SK_PPS,j,i,i_mag)
	H_coeff(MC_F,SK_PPS,i,j,i_mag) = H_coeff(MC_F,SK_PPS,j,i,i_mag)
	H_coeff(MC_FB,SK_PPS,i,j,i_mag) = H_coeff(MC_FB,SK_PPS,j,i,i_mag)
	H_coeff(MC_G_SQ,SK_PPS,i,j,i_mag) = H_coeff(MC_G_SQ,SK_PPS,j,i,i_mag)

	H_coeff(MC_E,SK_PPP,i,j,i_mag) = H_coeff(MC_E,SK_PPP,j,i,i_mag)
	H_coeff(MC_F,SK_PPP,i,j,i_mag) = H_coeff(MC_F,SK_PPP,j,i,i_mag)
	H_coeff(MC_FB,SK_PPP,i,j,i_mag) = H_coeff(MC_FB,SK_PPP,j,i,i_mag)
	H_coeff(MC_G_SQ,SK_PPP,i,j,i_mag) = H_coeff(MC_G_SQ,SK_PPP,j,i,i_mag)

	read *, H_coeff(MC_E,SK_SDS,i,j,i_mag)
	read *, H_coeff(MC_F,SK_SDS,i,j,i_mag)
	read *, H_coeff(MC_FB,SK_SDS,i,j,i_mag)
	read *, H_coeff(MC_G_SQ,SK_SDS,i,j,i_mag)

	read *, H_coeff(MC_E,SK_PDS,i,j,i_mag)
	read *, H_coeff(MC_F,SK_PDS,i,j,i_mag)
	read *, H_coeff(MC_FB,SK_PDS,i,j,i_mag)
	read *, H_coeff(MC_G_SQ,SK_PDS,i,j,i_mag)
	H_coeff(MC_E,SK_PDS,i,j,i_mag) = - H_coeff(MC_E,SK_PDS,i,j,i_mag)
	H_coeff(MC_F,SK_PDS,i,j,i_mag) = - H_coeff(MC_F,SK_PDS,i,j,i_mag)
	H_coeff(MC_FB,SK_PDS,i,j,i_mag) = - H_coeff(MC_FB,SK_PDS,i,j,i_mag)

	read *, H_coeff(MC_E,SK_PDP,i,j,i_mag)
	read *, H_coeff(MC_F,SK_PDP,i,j,i_mag)
	read *, H_coeff(MC_FB,SK_PDP,i,j,i_mag)
	read *, H_coeff(MC_G_SQ,SK_PDP,i,j,i_mag)
	H_coeff(MC_E,SK_PDP,i,j,i_mag) = - H_coeff(MC_E,SK_PDP,i,j,i_mag)
	H_coeff(MC_F,SK_PDP,i,j,i_mag) = - H_coeff(MC_F,SK_PDP,i,j,i_mag)
	H_coeff(MC_FB,SK_PDP,i,j,i_mag) = - H_coeff(MC_FB,SK_PDP,i,j,i_mag)

	H_coeff(MC_E,SK_DDS,i,j,i_mag) = H_coeff(MC_E,SK_DDS,j,i,i_mag)
	H_coeff(MC_F,SK_DDS,i,j,i_mag) = H_coeff(MC_F,SK_DDS,j,i,i_mag)
	H_coeff(MC_FB,SK_DDS,i,j,i_mag) = H_coeff(MC_FB,SK_DDS,j,i,i_mag)
	H_coeff(MC_G_SQ,SK_DDS,i,j,i_mag) = H_coeff(MC_G_SQ,SK_DDS,j,i,i_mag)

	H_coeff(MC_E,SK_DDP,i,j,i_mag) = H_coeff(MC_E,SK_DDP,j,i,i_mag)
	H_coeff(MC_F,SK_DDP,i,j,i_mag) = H_coeff(MC_F,SK_DDP,j,i,i_mag)
	H_coeff(MC_FB,SK_DDP,i,j,i_mag) = H_coeff(MC_FB,SK_DDP,j,i,i_mag)
	H_coeff(MC_G_SQ,SK_DDP,i,j,i_mag) = H_coeff(MC_G_SQ,SK_DDP,j,i,i_mag)

	H_coeff(MC_E,SK_DDD,i,j,i_mag) = H_coeff(MC_E,SK_DDD,j,i,i_mag)
	H_coeff(MC_F,SK_DDD,i,j,i_mag) = H_coeff(MC_F,SK_DDD,j,i,i_mag)
	H_coeff(MC_FB,SK_DDD,i,j,i_mag) = H_coeff(MC_FB,SK_DDD,j,i,i_mag)
	H_coeff(MC_G_SQ,SK_DDD,i,j,i_mag) = H_coeff(MC_G_SQ,SK_DDD,j,i,i_mag)

      end do
    end do

  !!!!!!!!! OVERLAP !!!!!!!!!!!!

    do i=1, n_types
    j=i
      read *, S_coeff(MC_E,SK_SSS,j,i,i_mag)
      read *, S_coeff(MC_F,SK_SSS,j,i,i_mag)
      read *, S_coeff(MC_FB,SK_SSS,j,i,i_mag)
      read *, S_coeff(MC_G_SQ,SK_SSS,j,i,i_mag)
      read *, S_coeff(MC_E,SK_SPS,j,i,i_mag)
      read *, S_coeff(MC_F,SK_SPS,j,i,i_mag)
      read *, S_coeff(MC_FB,SK_SPS,j,i,i_mag)
      read *, S_coeff(MC_G_SQ,SK_SPS,j,i,i_mag)
      read *, S_coeff(MC_E,SK_PPS,j,i,i_mag)
      read *, S_coeff(MC_F,SK_PPS,j,i,i_mag)
      read *, S_coeff(MC_FB,SK_PPS,j,i,i_mag)
      read *, S_coeff(MC_G_SQ,SK_PPS,j,i,i_mag)
      read *, S_coeff(MC_E,SK_PPP,j,i,i_mag)
      read *, S_coeff(MC_F,SK_PPP,j,i,i_mag)
      read *, S_coeff(MC_FB,SK_PPP,j,i,i_mag)
      read *, S_coeff(MC_G_SQ,SK_PPP,j,i,i_mag)
      read *, S_coeff(MC_E,SK_SDS,j,i,i_mag)
      read *, S_coeff(MC_F,SK_SDS,j,i,i_mag)
      read *, S_coeff(MC_FB,SK_SDS,j,i,i_mag)
      read *, S_coeff(MC_G_SQ,SK_SDS,j,i,i_mag)
      read *, S_coeff(MC_E,SK_PDS,j,i,i_mag)
      read *, S_coeff(MC_F,SK_PDS,j,i,i_mag)
      read *, S_coeff(MC_FB,SK_PDS,j,i,i_mag)
      read *, S_coeff(MC_G_SQ,SK_PDS,j,i,i_mag)
      read *, S_coeff(MC_E,SK_PDP,j,i,i_mag)
      read *, S_coeff(MC_F,SK_PDP,j,i,i_mag)
      read *, S_coeff(MC_FB,SK_PDP,j,i,i_mag)
      read *, S_coeff(MC_G_SQ,SK_PDP,j,i,i_mag)
      read *, S_coeff(MC_E,SK_DDS,j,i,i_mag)
      read *, S_coeff(MC_F,SK_DDS,j,i,i_mag)
      read *, S_coeff(MC_FB,SK_DDS,j,i,i_mag)
      read *, S_coeff(MC_G_SQ,SK_DDS,j,i,i_mag)
      read *, S_coeff(MC_E,SK_DDP,j,i,i_mag)
      read *, S_coeff(MC_F,SK_DDP,j,i,i_mag)
      read *, S_coeff(MC_FB,SK_DDP,j,i,i_mag)
      read *, S_coeff(MC_G_SQ,SK_DDP,j,i,i_mag)
      read *, S_coeff(MC_E,SK_DDD,j,i,i_mag)
      read *, S_coeff(MC_F,SK_DDD,j,i,i_mag)
      read *, S_coeff(MC_FB,SK_DDD,j,i,i_mag)
      read *, S_coeff(MC_G_SQ,SK_DDD,j,i,i_mag)
    end do

    do i=2, n_types
      do j=1, i-1
	read *, S_coeff(MC_E,SK_SSS,j,i,i_mag)
	read *, S_coeff(MC_F,SK_SSS,j,i,i_mag)
	read *, S_coeff(MC_FB,SK_SSS,j,i,i_mag)
	read *, S_coeff(MC_G_SQ,SK_SSS,j,i,i_mag)
	read *, S_coeff(MC_E,SK_SPS,j,i,i_mag)
	read *, S_coeff(MC_F,SK_SPS,j,i,i_mag)
	read *, S_coeff(MC_FB,SK_SPS,j,i,i_mag)
	read *, S_coeff(MC_G_SQ,SK_SPS,j,i,i_mag)
	read *, S_coeff(MC_E,SK_PPS,j,i,i_mag)
	read *, S_coeff(MC_F,SK_PPS,j,i,i_mag)
	read *, S_coeff(MC_FB,SK_PPS,j,i,i_mag)
	read *, S_coeff(MC_G_SQ,SK_PPS,j,i,i_mag)
	read *, S_coeff(MC_E,SK_PPP,j,i,i_mag)
	read *, S_coeff(MC_F,SK_PPP,j,i,i_mag)
	read *, S_coeff(MC_FB,SK_PPP,j,i,i_mag)
	read *, S_coeff(MC_G_SQ,SK_PPP,j,i,i_mag)
	read *, S_coeff(MC_E,SK_SDS,j,i,i_mag)
	read *, S_coeff(MC_F,SK_SDS,j,i,i_mag)
	read *, S_coeff(MC_FB,SK_SDS,j,i,i_mag)
	read *, S_coeff(MC_G_SQ,SK_SDS,j,i,i_mag)
	read *, S_coeff(MC_E,SK_PDS,j,i,i_mag)
	read *, S_coeff(MC_F,SK_PDS,j,i,i_mag)
	read *, S_coeff(MC_FB,SK_PDS,j,i,i_mag)
	read *, S_coeff(MC_G_SQ,SK_PDS,j,i,i_mag)
	read *, S_coeff(MC_E,SK_PDP,j,i,i_mag)
	read *, S_coeff(MC_F,SK_PDP,j,i,i_mag)
	read *, S_coeff(MC_FB,SK_PDP,j,i,i_mag)
	read *, S_coeff(MC_G_SQ,SK_PDP,j,i,i_mag)
	read *, S_coeff(MC_E,SK_DDS,j,i,i_mag)
	read *, S_coeff(MC_F,SK_DDS,j,i,i_mag)
	read *, S_coeff(MC_FB,SK_DDS,j,i,i_mag)
	read *, S_coeff(MC_G_SQ,SK_DDS,j,i,i_mag)
	read *, S_coeff(MC_E,SK_DDP,j,i,i_mag)
	read *, S_coeff(MC_F,SK_DDP,j,i,i_mag)
	read *, S_coeff(MC_FB,SK_DDP,j,i,i_mag)
	read *, S_coeff(MC_G_SQ,SK_DDP,j,i,i_mag)
	read *, S_coeff(MC_E,SK_DDD,j,i,i_mag)
	read *, S_coeff(MC_F,SK_DDD,j,i,i_mag)
	read *, S_coeff(MC_FB,SK_DDD,j,i,i_mag)
	read *, S_coeff(MC_G_SQ,SK_DDD,j,i,i_mag)

	S_coeff(MC_E,SK_SSS,i,j,i_mag) = S_coeff(MC_E,SK_SSS,j,i,i_mag)
	S_coeff(MC_F,SK_SSS,i,j,i_mag) = S_coeff(MC_F,SK_SSS,j,i,i_mag)
	S_coeff(MC_FB,SK_SSS,i,j,i_mag) = S_coeff(MC_FB,SK_SSS,j,i,i_mag)
	S_coeff(MC_G_SQ,SK_SSS,i,j,i_mag) = S_coeff(MC_G_SQ,SK_SSS,j,i,i_mag)

	read *, S_coeff(MC_E,SK_SPS,i,j,i_mag)
	read *, S_coeff(MC_F,SK_SPS,i,j,i_mag)
	read *, S_coeff(MC_FB,SK_SPS,i,j,i_mag)
	read *, S_coeff(MC_G_SQ,SK_SPS,i,j,i_mag)
	S_coeff(MC_E,SK_SPS,i,j,i_mag) = - S_coeff(MC_E,SK_SPS,i,j,i_mag)
	S_coeff(MC_F,SK_SPS,i,j,i_mag) = - S_coeff(MC_F,SK_SPS,i,j,i_mag)
	S_coeff(MC_FB,SK_SPS,i,j,i_mag) = - S_coeff(MC_FB,SK_SPS,i,j,i_mag)

	S_coeff(MC_E,SK_PPS,i,j,i_mag) = S_coeff(MC_E,SK_PPS,j,i,i_mag)
	S_coeff(MC_F,SK_PPS,i,j,i_mag) = S_coeff(MC_F,SK_PPS,j,i,i_mag)
	S_coeff(MC_FB,SK_PPS,i,j,i_mag) = S_coeff(MC_FB,SK_PPS,j,i,i_mag)
	S_coeff(MC_G_SQ,SK_PPS,i,j,i_mag) = S_coeff(MC_G_SQ,SK_PPS,j,i,i_mag)

	S_coeff(MC_E,SK_PPP,i,j,i_mag) = S_coeff(MC_E,SK_PPP,j,i,i_mag)
	S_coeff(MC_F,SK_PPP,i,j,i_mag) = S_coeff(MC_F,SK_PPP,j,i,i_mag)
	S_coeff(MC_FB,SK_PPP,i,j,i_mag) = S_coeff(MC_FB,SK_PPP,j,i,i_mag)
	S_coeff(MC_G_SQ,SK_PPP,i,j,i_mag) = S_coeff(MC_G_SQ,SK_PPP,j,i,i_mag)

	read *, S_coeff(MC_E,SK_SDS,i,j,i_mag)
	read *, S_coeff(MC_F,SK_SDS,i,j,i_mag)
	read *, S_coeff(MC_FB,SK_SDS,i,j,i_mag)
	read *, S_coeff(MC_G_SQ,SK_SDS,i,j,i_mag)

	read *, S_coeff(MC_E,SK_PDS,i,j,i_mag)
	read *, S_coeff(MC_F,SK_PDS,i,j,i_mag)
	read *, S_coeff(MC_FB,SK_PDS,i,j,i_mag)
	read *, S_coeff(MC_G_SQ,SK_PDS,i,j,i_mag)
	S_coeff(MC_E,SK_PDS,i,j,i_mag) = - S_coeff(MC_E,SK_PDS,i,j,i_mag)
	S_coeff(MC_F,SK_PDS,i,j,i_mag) = - S_coeff(MC_F,SK_PDS,i,j,i_mag)
	S_coeff(MC_FB,SK_PDS,i,j,i_mag) = - S_coeff(MC_FB,SK_PDS,i,j,i_mag)

	read *, S_coeff(MC_E,SK_PDP,i,j,i_mag)
	read *, S_coeff(MC_F,SK_PDP,i,j,i_mag)
	read *, S_coeff(MC_FB,SK_PDP,i,j,i_mag)
	read *, S_coeff(MC_G_SQ,SK_PDP,i,j,i_mag)
	S_coeff(MC_E,SK_PDP,i,j,i_mag) = - S_coeff(MC_E,SK_PDP,i,j,i_mag)
	S_coeff(MC_F,SK_PDP,i,j,i_mag) = - S_coeff(MC_F,SK_PDP,i,j,i_mag)
	S_coeff(MC_FB,SK_PDP,i,j,i_mag) = - S_coeff(MC_FB,SK_PDP,i,j,i_mag)

	S_coeff(MC_E,SK_DDS,i,j,i_mag) = S_coeff(MC_E,SK_DDS,j,i,i_mag)
	S_coeff(MC_F,SK_DDS,i,j,i_mag) = S_coeff(MC_F,SK_DDS,j,i,i_mag)
	S_coeff(MC_FB,SK_DDS,i,j,i_mag) = S_coeff(MC_FB,SK_DDS,j,i,i_mag)
	S_coeff(MC_G_SQ,SK_DDS,i,j,i_mag) = S_coeff(MC_G_SQ,SK_DDS,j,i,i_mag)

	S_coeff(MC_E,SK_DDP,i,j,i_mag) = S_coeff(MC_E,SK_DDP,j,i,i_mag)
	S_coeff(MC_F,SK_DDP,i,j,i_mag) = S_coeff(MC_F,SK_DDP,j,i,i_mag)
	S_coeff(MC_FB,SK_DDP,i,j,i_mag) = S_coeff(MC_FB,SK_DDP,j,i,i_mag)
	S_coeff(MC_G_SQ,SK_DDP,i,j,i_mag) = S_coeff(MC_G_SQ,SK_DDP,j,i,i_mag)

	S_coeff(MC_E,SK_DDD,i,j,i_mag) = S_coeff(MC_E,SK_DDD,j,i,i_mag)
	S_coeff(MC_F,SK_DDD,i,j,i_mag) = S_coeff(MC_F,SK_DDD,j,i,i_mag)
	S_coeff(MC_FB,SK_DDD,i,j,i_mag) = S_coeff(MC_FB,SK_DDD,j,i,i_mag)
	S_coeff(MC_G_SQ,SK_DDD,i,j,i_mag) = S_coeff(MC_G_SQ,SK_DDD,j,i,i_mag)

      end do
    end do
  end do ! i_mag

  lambda_sq = lambda_sq**2
  H_coeff(MC_G_SQ,:,:,:,:) = H_coeff(MC_G_SQ,:,:,:,:)**2
  S_coeff(MC_G_SQ,:,:,:,:) = S_coeff(MC_G_SQ,:,:,:,:)**2

  call print_NRL_TB_xml(trim(fname), header, is_orthogonal, &
    is_magnetic, has_pair_repulsion, overlap_zero_limit, force_harrison_signs, label, &
    n_types, n_mag, atomic_num, n_orbs, n_elecs, n_orb_sets, max_n_orb_sets, orb_set_type, atomic_mass, r_cut, screen_l, &
    pair_rep_inner, pair_rep_outer, cutoff, lambda_sq, abcd, H_coeff, S_coeff)


end subroutine

subroutine print_NRL_TB_xml(fname, header, is_orthogonal, &
    is_magnetic, has_pair_repulsion, overlap_zero_limit, force_harrison_signs, label, &
    n_types, n_mag, atomic_num, n_orbs, n_elecs, n_orb_sets, max_n_orb_sets, orb_set_type, atomic_mass, r_cut, screen_l, &
    pair_rep_inner, pair_rep_outer, cutoff, lambda_sq, abcd, H_coeff, S_coeff)
implicit none
  character(len=*), intent(in) :: fname, header
  logical, intent(in) :: is_orthogonal, is_magnetic, has_pair_repulsion, overlap_zero_limit, force_harrison_signs
  character(len=*), intent(in) :: label
  integer, intent(in) :: n_types, n_mag
  integer, intent(in) :: max_n_orb_sets
  integer, intent(in) :: atomic_num(n_types), n_orbs(n_types), n_elecs(n_types), n_orb_sets(n_types), orb_set_type(n_types,max_n_orb_sets)
  double precision, intent(in) :: atomic_mass(n_types), r_cut(n_types,n_types), screen_l(n_types,n_types)
  double precision, allocatable :: pair_rep_inner(:,:), pair_rep_outer(:,:)
  double precision, intent(in) :: cutoff
  double precision, intent(in) :: lambda_sq(n_types,n_mag), abcd(4,3,n_types,n_types,n_mag)
  double precision, intent(in) :: H_coeff(4,10,n_types,n_types,n_mag), S_coeff(4,10,n_types,n_types,n_mag)

  type(xmlf_t) :: xf

  integer :: i, j, k, i_mag

  call xml_OpenFile(fname,xf)
  ! call xml_AddXMLDeclaration(xf)
  call xml_NewElement(xf,"NRL_TB_params")
  call xml_AddAttribute(xf,"header_str",trim(header))
  call xml_AddAttribute(xf,"label",trim(label))

  call xml_NewElement(xf,"header")
    call xml_AddAttribute(xf,"is_orthogonal","" // is_orthogonal)
    call xml_AddAttribute(xf,"is_magnetic", "" // is_magnetic)
    call xml_AddAttribute(xf,"has_pair_repulsion", "" // has_pair_repulsion)
    call xml_AddAttribute(xf,"overlap_zero_limit", "" // overlap_zero_limit)
    call xml_AddAttribute(xf,"force_harrison_signs", "" // force_harrison_signs)
   call xml_EndElement(xf,"header")

  call xml_NewElement(xf,"n_types")
    call xml_AddAttribute(xf,"v", "" // n_types)
   call xml_EndElement(xf,"n_types")

  call xml_NewElement(xf,"cutoff")
    call xml_AddAttribute(xf,"v", "" // cutoff)
   call xml_EndElement(xf,"cutoff")

  do i=1, n_types
    call xml_NewElement(xf,"per_type_data")
      call xml_AddAttribute(xf,"type", "" // i)
      call xml_AddAttribute(xf,"atomic_num", "" // atomic_num(i))
      call xml_AddAttribute(xf,"atomic_mass", "" // atomic_mass(i))
      call xml_AddAttribute(xf,"n_orbs", "" // n_orbs(i))
      call xml_AddAttribute(xf,"n_elecs", "" // n_elecs(i))
      call xml_AddAttribute(xf,"n_orb_sets", "" // n_orb_sets(i))
      if (is_magnetic) then
	call xml_AddAttribute(xf,"lambda_sq_up", "" // lambda_sq(i,1))
	call xml_AddAttribute(xf,"lambda_sq_down", "" // lambda_sq(i,2))
      else
	call xml_AddAttribute(xf,"lambda_sq", "" // lambda_sq(i,1))
      endif
      call xml_NewElement(xf,"orb_set_type")
	! call xml_AddArray(xf, orb_set_type(i,1:n_orb_sets(i)), '(i0,1x)')
	call xml_AddCharacters(xf, "" // orb_set_type(i,1:n_orb_sets(i)))
      call xml_EndElement(xf,"orb_set_type")
    call xml_EndElement(xf,"per_type_data")
  end do

  do i=1, n_types
  do j=1, n_types
    call xml_NewElement(xf,"per_pair_data")
      call xml_AddAttribute(xf,"type1", "" // i)
      call xml_AddAttribute(xf,"type2", "" // j)
      call xml_AddAttribute(xf,"r_cut", "" // r_cut(i,j))
      call xml_AddAttribute(xf,"screen_l", "" // screen_l(i,j))
      if (has_pair_repulsion) then
	call xml_AddAttribute(xf,"pair_rep_inner", "" // pair_rep_inner(i,j))
	call xml_AddAttribute(xf,"pair_rep_outer", "" // pair_rep_outer(i,j))
      endif
      call xml_NewElement(xf,"abcd")
	do k=SPD_S, SPD_D
	  ! call xml_AddArray(xf,abcd(1:4,k,i,j), '(f20.10,1x)')
	  do i_mag=1, n_mag
	    call xml_AddCharacters(xf, "" // abcd(1:4,k,i,j,i_mag) // " ")
	  end do
	  call xml_AddNewLine(xf)
	end do
      call xml_EndElement(xf,"abcd")
      call xml_NewElement(xf,"H_coeff")
	do k=SK_SSS, SK_DDD
	  ! call xml_AddArray(xf,H_coeff(1:4,k,i,j), '(f20.10,1x)')
	  do i_mag=1, n_mag
	    call xml_AddCharacters(xf, "" // H_coeff(1:4,k,i,j,i_mag) // " ")
	  end do
	  call xml_AddNewLine(xf)
	end do
      call xml_EndElement(xf,"H_coeff")
      call xml_NewElement(xf,"S_coeff")
	do k=SK_SSS, SK_DDD
	  ! call xml_AddArray(xf,S_coeff(1:4,k,i,j), '(f20.10,1x)')
	  do i_mag=1, n_mag
	    call xml_AddCharacters(xf, "" // S_coeff(1:4,k,i,j,i_mag) // " ")
	  end do
	  call xml_AddNewLine(xf)
	end do
      call xml_EndElement(xf,"S_coeff")
    call xml_EndElement(xf,"per_pair_data")
  end do
  end do

  call xml_EndElement(xf,"NRL_TB_params")

  call xml_Close(xf)
end subroutine

end module NRL_TB_to_xml_module

program NRL_TB_to_xml
use NRL_TB_to_xml_module
implicit none

  character(len=1024) fname

  call system_initialise()

  if (command_argument_count() == 0) then
    fname = "tightbind.parms.NRL_TB.xml"
  else if (command_argument_count() == 1) then
    call get_command_argument(1, fname)
  else
    call system_abort("Usage: NRL_TB_to_xml [ output_filename ]")
  endif

  call do_nrl_tb_params_to_xml(fname)

  call system_finalise()
end program


