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

! convert from DFTB web-site input files to QUIP xml format
module DFTB_to_xml_module

use System_module
use TB_Common_module
use fox_wxml

implicit none

contains 

subroutine do_dftb_params_to_xml(fname)
  character(len=*), intent(in) :: fname

  integer n_types
  integer, allocatable :: atomic_num(:), n_orbs(:), n_elecs(:), n_orb_sets(:), orb_set_type(:,:)
  double precision, allocatable :: SK_cutoff(:,:), Vrep_cutoff(:,:)
  double precision, allocatable :: E(:,:), E_spin(:), U(:)
  double precision, allocatable :: H_vals(:,:,:,:), S_vals(:,:,:,:), SK_r(:,:,:)
  double precision, allocatable :: Vrep_efkt(:,:,:), Vrep_r_in(:,:,:,:), Vrep_coeff(:,:,:,:)
  double precision, allocatable :: Vrep_r(:,:,:), Vrep_vals(:,:,:)

  integer max_n_orb_sets

  double precision, allocatable :: SK_dr(:,:), Vrep_maxr(:,:)
  integer, allocatable :: SK_npts(:,:), Vrep_npts(:,:)

  character(len=1024), allocatable :: filename(:,:)

  integer, allocatable :: max_l(:)

  real(dp) :: r

  double precision self_fields(10), SK_fields(20)
  integer n_elec_s, n_elec_p, n_elec_d
  character(len=1024) line

  integer :: n_t_atomic_nums = 200
  integer t_atomic_nums(200)

  integer ti, tj, li, i

  read *, n_types, (t_atomic_nums(i), i=1,min(n_types,n_t_atomic_nums))
  if (n_types > n_t_atomic_nums) then
    call system_abort("too many types in dftb.ini file")
  endif

  allocate(atomic_num(n_types))
  allocate(n_orbs(n_types))
  allocate(n_elecs(n_types))
  allocate(n_orb_sets(n_types))

  atomic_num(1:n_types) = t_atomic_nums(1:n_types)

  allocate(max_l(n_types))

!!!!!! read non matrix element stuff !!!!!!
  read *, (max_l(ti), ti=1,n_types)
  max_n_orb_sets = maxval(max_l)

  allocate(orb_set_type(n_types,max_n_orb_sets))

  do ti=1, n_types
    select case(max_l(ti))
      case (1)
	n_orbs(ti) = 1
	n_orb_sets(ti) = 1
	orb_set_type(ti,1) = ORB_S
      case (2)
	n_orbs(ti) = 4
	n_orb_sets(ti) = 2
	orb_set_type(ti,1) = ORB_S
	orb_set_type(ti,2) = ORB_P
      case (3)
	n_orbs(ti) = 9
	n_orb_sets(ti) = 3
	orb_set_type(ti,1) = ORB_S
	orb_set_type(ti,2) = ORB_P
	orb_set_type(ti,3) = ORB_D
      case default
	write (line, '("Confused by max_l ",2I4)') max_l(i), ti
	call system_abort(line)
    end select
  end do

  allocate(E(max_n_orb_sets,n_types))
  allocate(E_spin(n_types))
  allocate(U(n_types))

  allocate(SK_dr(n_types,n_types))
  allocate(SK_npts(n_types,n_types))
  allocate(Vrep_maxr(n_types,n_types))
  allocate(Vrep_npts(n_types,n_types))

  allocate(filename(n_types,n_types))

  allocate(SK_cutoff(n_types,n_types))
  allocate(Vrep_cutoff(n_types,n_types))

  do tj=1, n_types
  do ti=1, n_types
    read *, filename(ti,tj)
    open (unit=20, file=trim(filename(ti,tj)), status="OLD")
    read (20,*) SK_dr(ti,tj), SK_npts(ti,tj)
    if (ti == tj) then
      read (20,*) self_fields
      do li=1, max_l(ti)
	E(li,ti) = self_fields(4-li)
      end do
      E_spin(ti) = self_fields(4)
      U(ti) = self_fields(5)
      n_elec_s = self_fields(10)
      n_elec_p = self_fields(9)
      n_elec_d = self_fields(8)
      n_elecs(ti) = n_elec_s + n_elec_p + n_elec_d
    endif ! i == j

    SK_cutoff(ti,tj) = SK_dr(ti,tj)*SK_npts(ti,tj)

    do li=1, SK_npts(ti,tj)
      read (20,*) line
    end do

    do while (line /= 'Spline')
      read (20,*) line
    end do
    read (20,*) Vrep_npts(ti,tj), Vrep_maxr(ti,tj)
    Vrep_cutoff(ti,tj) = Vrep_maxr(ti,tj)
  end do
  end do

  allocate(SK_r(maxval(SK_npts)+1,n_types,n_types))
  allocate(H_vals(maxval(SK_npts)+1,10,n_types,n_types))
  allocate(S_vals(maxval(SK_npts)+1,10,n_types,n_types))

  allocate(Vrep_efkt(3,n_types,n_types))
  allocate(Vrep_r_in(maxval(Vrep_npts),2,n_types,n_types))
  allocate(Vrep_coeff(maxval(Vrep_npts),6,n_types,n_types))

  do tj=1, n_types
  do ti=1, n_types
    open (unit=20, file=trim(filename(ti,tj)), status="OLD")
    read (20,*) line ! dr, npts
    if (ti == tj) then
      read (20,*) line ! onsite stuff
    endif
    do li=1, SK_npts(ti,tj)
      read (20, *) H_vals(li, SK_DDS, ti, tj), &
		   H_vals(li, SK_DDP, ti, tj), &
		   H_vals(li, SK_DDD, ti, tj), &
		   H_vals(li, SK_PDS, ti, tj), &
		   H_vals(li, SK_PDP, ti, tj), &
		   H_vals(li, SK_PPS, ti, tj), &
		   H_vals(li, SK_PPP, ti, tj), &
		   H_vals(li, SK_SDS, ti, tj), &
		   H_vals(li, SK_SPS, ti, tj), &
		   H_vals(li, SK_SSS, ti, tj), &
                   S_vals(li, SK_DDS, ti, tj), &
		   S_vals(li, SK_DDP, ti, tj), &
		   S_vals(li, SK_DDD, ti, tj), &
		   S_vals(li, SK_PDS, ti, tj), &
		   S_vals(li, SK_PDP, ti, tj), &
		   S_vals(li, SK_PPS, ti, tj), &
		   S_vals(li, SK_PPP, ti, tj), &
		   S_vals(li, SK_SDS, ti, tj), &
		   S_vals(li, SK_SPS, ti, tj), &
		   S_vals(li, SK_SSS, ti, tj)
      SK_r(li,ti,tj) = (li-1)*SK_dr(ti,tj)
    end do
    H_vals(SK_npts(ti,tj)+1,1:10,ti,tj) = 0.0D0
    S_vals(SK_npts(ti,tj)+1,1:10,ti,tj) = 0.0D0
    SK_r(SK_npts(ti,tj)+1,ti,tj) = SK_npts(ti,tj)*SK_dr(ti,tj)

    do while (line /= 'Spline') 
      read (20,*) line
    end do
    read (20,*) line ! npts, maxr
    read (20,*) Vrep_Efkt(1:3,ti,tj)
    do li=1, Vrep_npts(ti,tj)
      if (li == Vrep_npts(ti,tj)) then
	read(20,*) Vrep_r_in(li,1:2, ti,tj), Vrep_coeff(li,1:6,ti,tj)
      else
	read(20,*) Vrep_r_in(li,1:2, ti,tj), Vrep_coeff(li,1:4,ti,tj)
      endif
    end do
  end do
  end do

  SK_npts = SK_npts + 1

  allocate(Vrep_r(maxval(Vrep_npts+11), n_types, n_types))
  allocate(Vrep_vals(maxval(Vrep_npts+11), n_types, n_types))

  do tj=1, n_types
  do ti=1, n_types
    do li=1, 10
      r = (li-1)/10.0D0*Vrep_r_in(1,1,ti,tj)
      Vrep_r(li,ti,tj) = r
      Vrep_vals(li,ti,tj) = exp(-Vrep_efkt(1,ti,tj)*r+Vrep_efkt(2,ti,tj))+Vrep_efkt(3,ti,tj)
    end do
    do li=1, Vrep_npts(ti,tj)
      Vrep_r(li+10,ti,tj) = Vrep_r_in(li,1,ti,tj)
      Vrep_vals(li+10,ti,tj) = Vrep_coeff(li,1,ti,tj)
    end do
    Vrep_r(11+Vrep_npts(ti,tj),ti,tj) = Vrep_r_in(Vrep_npts(ti,tj),2,ti,tj)
    Vrep_vals(11+Vrep_npts(ti,tj),ti,tj) = 0.0D0
  end do
  end do
  Vrep_npts = Vrep_npts + 11

  deallocate(max_l)

  call print_DFTB_xml(fname, n_types, atomic_num, n_orbs, n_elecs, n_orb_sets, max_n_orb_sets, &
    orb_set_type, E, SK_cutoff, SK_npts, SK_r, H_vals, S_vals, Vrep_cutoff, Vrep_npts, Vrep_r, Vrep_vals)
end subroutine

subroutine print_DFTB_xml(fname, n_types, atomic_num, n_orbs, n_elecs, n_orb_sets, max_n_orb_sets, &
    orb_set_type, E, SK_cutoff, SK_npts, SK_r, H_vals, S_vals, Vrep_cutoff, Vrep_npts, Vrep_r, Vrep_vals)
implicit none
  character(len=*), intent(in) :: fname
  integer, intent(in) :: n_types
  integer, intent(in) :: max_n_orb_sets
  integer, intent(in) :: atomic_num(:), n_orbs(:), n_elecs(:), n_orb_sets(:), &
    orb_set_type(:,:)
  double precision, intent(in) :: E(:,:)
  double precision, intent(in) :: SK_cutoff(:,:)
  integer, intent(in) :: SK_npts(:,:)
  double precision, intent(in) :: SK_r(:,:,:), H_vals(:,:,:,:), S_vals(:,:,:,:)
  double precision, intent(in) :: Vrep_cutoff(:,:)
  integer, intent(in) :: Vrep_npts(:,:)
  double precision, intent(in) :: Vrep_r(:,:,:), Vrep_vals(:,:,:)

  type(xmlf_t) :: xf

  integer :: i, j, k, l
  integer n_non_zero

  call xml_OpenFile(fname,xf)
  call xml_NewElement(xf,"DFTB_params")

  call xml_NewElement(xf,"n_types")
    call xml_AddAttribute(xf,"v", "" // n_types)
   call xml_EndElement(xf,"n_types")

  call xml_NewElement(xf,"cutoff")
    call xml_AddAttribute(xf,"v", "" // max(maxval(SK_cutoff),maxval(Vrep_cutoff)))
   call xml_EndElement(xf,"cutoff")

  call xml_NewElement(xf,"max_n_orb_sets")
    call xml_AddAttribute(xf,"v", "" // max_n_orb_sets)
   call xml_EndElement(xf,"max_n_orb_sets")

  do i=1, n_types
    call xml_NewElement(xf,"per_type_data")
      call xml_AddAttribute(xf,"type", "" // i)
      call xml_AddAttribute(xf,"atomic_num", "" // atomic_num(i))
      call xml_AddAttribute(xf,"n_orbs", "" // n_orbs(i))
      call xml_AddAttribute(xf,"n_elecs", "" // n_elecs(i))
      call xml_AddAttribute(xf,"n_orb_sets", "" // n_orb_sets(i))
      call xml_NewElement(xf,"orb_set_type")
	! call xml_AddArray(xf, orb_set_type(i,1:n_orb_sets(i)), '(i0,1x)')
	call xml_AddCharacters(xf, "" //  orb_set_type(i,1:n_orb_sets(i)))
      call xml_EndElement(xf,"orb_set_type")
      call xml_NewElement(xf,"E")
	! call xml_AddArray(xf, E(1:n_orb_sets(i),i), '(f30.20,1x)')
	call xml_AddCharacters(xf, "" // E(1:n_orb_sets(i),i))
      call xml_EndElement(xf,"E")
    call xml_EndElement(xf,"per_type_data")
  end do

  do i=1, n_types
  do j=1, n_types
print *, "doing per pair data ", i, j
    call xml_NewElement(xf,"per_pair_data")
      call xml_AddAttribute(xf,"type1", "" // i)
      call xml_AddAttribute(xf,"type2", "" // j)
      call xml_AddAttribute(xf,"SK_cutoff", "" // SK_cutoff(i,j))
      call xml_AddAttribute(xf,"Vrep_cutoff", "" // Vrep_cutoff(i,j))
      call xml_AddAttribute(xf,"SK_npts", "" // SK_npts(i,j))
      call xml_AddAttribute(xf,"Vrep_npts", "" // Vrep_npts(i,j))
print *, "done attributes"
      call xml_NewElement(xf,"H_spline")
	do k=1, SK_npts(i,j)
	  call xml_NewElement(xf,"point")
	    mainlog%default_real_precision=3
	    call xml_AddAttribute(xf,"r", "" // SK_r(k,i,j))
	    ! call xml_AddArray(xf,H_vals(k,:,i,j), '(f20.10,1x)')
	    mainlog%default_real_precision=8
	    n_non_zero = 1
	    do l=size(H_vals,2), 1, -1
	      if (H_vals(k,l,i,j) /= 0.0) then
		n_non_zero = l
		exit
	      endif
	    end do
	    call xml_AddCharacters(xf, "" // H_vals(k,1:n_non_zero,i,j))
	    mainlog%default_real_precision=16
	  call xml_EndElement(xf,"point")
	end do
      call xml_EndElement(xf,"H_spline")
print *, "done H_spline"
      call xml_NewElement(xf,"S_spline")
	do k=1, SK_npts(i,j)
	  call xml_NewElement(xf,"point")
	    mainlog%default_real_precision=3
	    call xml_AddAttribute(xf,"r","" // SK_r(k,i,j))
	    mainlog%default_real_precision=8
	    ! call xml_AddArray(xf,S_vals(k,:,i,j), '(f20.10,1x)')
	    n_non_zero = 1
	    do l=size(S_vals,2), 1, -1
	      if (S_vals(k,l,i,j) /= 0.0) then
		n_non_zero = l
		exit
	      endif
	    end do
	    call xml_AddCharacters(xf, "" // S_vals(k,1:n_non_zero,i,j))
	    mainlog%default_real_precision=16
	  call xml_EndElement(xf,"point")
	end do
      call xml_EndElement(xf,"S_spline")
print *, "done S_spline"
      call xml_NewElement(xf,"Vrep_spline")
	do k=1, Vrep_npts(i,j)
	  call xml_NewElement(xf,"point")
	    mainlog%default_real_precision=3
	    call xml_AddAttribute(xf,"r", "" // Vrep_r(k,i,j))
	    mainlog%default_real_precision=8
	    ! call xml_AddArray(xf,Vrep_vals(k:k,i,j), '(f20.10,1x)')
	    call xml_AddCharacters(xf, "" // Vrep_vals(k:k,i,j))
	    mainlog%default_real_precision=16
	  call xml_EndElement(xf,"point")
	end do
      call xml_EndElement(xf,"Vrep_spline")
print *, "done Vrep_spline"
    call xml_EndElement(xf,"per_pair_data")
  end do
  end do

  call xml_EndElement(xf,"DFTB_params")

  call xml_Close(xf)
end subroutine

end module DFTB_to_xml_module

program DFTB_to_xml
use DFTB_to_xml_module
implicit none

  character(len=1024) fname

  call system_initialise()

  if (command_argument_count() == 0) then
    fname = "tightbind.parms.DFTB.xml"
  else if (command_argument_count() == 1) then
    call get_command_argument(1, fname)
  else
    call system_abort("Usage: DFTB_to_xml [ output_filename ]")
  endif

  call do_dftb_params_to_xml(fname)

  call system_finalise()
end program


