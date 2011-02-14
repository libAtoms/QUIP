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
!X IPModel_Brenner module  
!X
!% Module for computing energy, forces and virial with Brenner potential for hydrocarbons. 
!% Ref. D. W. Brenner, Empirical potential for hydrocarbons for use in 
!% simulating the chemical vapor deposition of diamond films.
!% Phys. Rev. B  {\bf 42}, 9458 (1990).
!% The 'IPModel_Brenner' object reads parameters from a 'Brenner_params' XML stanza.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#include "error.inc"

module IPModel_Brenner_module

use libAtoms_module
use mpi_context_module
use QUIP_Common_module

implicit none
private

include 'IPModel_interface.h'

public :: IPModel_Brenner
type IPModel_Brenner
  integer :: n_types = 0         !% Number of atomic types 
  integer, allocatable :: atomic_num(:), type_of_atomic_num(:)  !% Atomic number dimensioned as \texttt{n_types} 

  real(dp) :: cutoff = 0.0_dp

  real(dp), allocatable :: R1(:,:), R2(:,:), De(:,:), S(:,:), beta(:,:)  !% IP parameters
  real(dp), allocatable :: delta(:,:), Re(:,:), a0(:,:,:), c0(:,:,:), d0(:,:,:) !% IP parameters
  real(dp), allocatable :: shift(:,:)                                 

  character(len=FIELD_LENGTH) :: label

end type IPModel_Brenner

logical, private :: parse_in_ip, parse_matched_label
type(IPModel_Brenner), private, pointer :: parse_ip

interface Initialise
  module procedure IPModel_Brenner_Initialise_str
end interface Initialise

interface Finalise
  module procedure IPModel_Brenner_Finalise
end interface Finalise

interface Print
  module procedure IPModel_Brenner_Print
end interface Print

interface Calc
  module procedure IPModel_Brenner_Calc
end interface Calc

contains

subroutine IPModel_Brenner_Initialise_str(this, args_str, param_str)
  type(IPModel_Brenner), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str

  type(Dictionary) :: params

  call Finalise(this)

  call initialise(params)
  this%label=''
  call param_register(params, 'label', '', this%label, help_string="No help yet.  This source file was $LastChangedBy$")
  if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_Brenner_Initialise_str args_str')) then
    call system_abort("IPModel_Brenner_Initialise_str failed to parse label from args_str="//trim(args_str))
  endif
  call finalise(params)

  call IPModel_Brenner_read_params_xml(this, param_str)

end subroutine IPModel_Brenner_Initialise_str

subroutine IPModel_Brenner_Finalise(this)
  type(IPModel_Brenner), intent(inout) :: this

  if (allocated(this%atomic_num)) deallocate(this%atomic_num)
  if (allocated(this%type_of_atomic_num)) deallocate(this%type_of_atomic_num)

  if (allocated(this%R1))    deallocate(this%R1)
  if (allocated(this%R2))    deallocate(this%R2)
  if (allocated(this%De))    deallocate(this%De)
  if (allocated(this%S))     deallocate(this%S)
  if (allocated(this%beta))  deallocate(this%beta)
  if (allocated(this%delta)) deallocate(this%delta)
  if (allocated(this%Re))    deallocate(this%Re)
  if (allocated(this%a0))    deallocate(this%a0)
  if (allocated(this%c0))    deallocate(this%c0)
  if (allocated(this%d0))    deallocate(this%d0)
  if (allocated(this%shift)) deallocate(this%shift)

  this%n_types = 0
  this%label = ''
end subroutine IPModel_Brenner_Finalise

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!X The potential calculator
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!% This routine computes energy, forces and the virial.
!% Derivatives by James Kermode <jrk33@cam.ac.uk>.
subroutine IPModel_Brenner_Calc(this, at, e, local_e, f, virial, local_virial, args_str, mpi, error)
  type(IPModel_Brenner), intent(inout) :: this
  type(Atoms), intent(in) :: at
  real(dp), intent(out), optional :: e, local_e(:) !% \texttt{e} = System total energy, \texttt{local_e} = energy of each atom, vector dimensioned as \texttt{at%N}.  
  real(dp), intent(out), optional :: f(:,:), local_virial(:,:)   !% Forces, dimensioned as \texttt{f(3,at%N)}, local virials, dimensioned as \texttt{local_virial(9,at%N)} 
  real(dp), intent(out), optional :: virial(3,3)   !% Virial
  character(len=*), optional      :: args_str
  type(MPI_Context), intent(in), optional :: mpi
  integer, intent(out), optional :: error

  integer :: i,j,k,m,n,r,s,p, max_neighb
  integer :: ti,tj,tk,tr,ts
  real(dp) :: r_ij, fc_ij, dfc_ij,  &
       Vr1, Va1, Vr, Va, dVr_ij, dVa_ij, G_sum, B(2), Bbar, &
       t1, t2, prefactor
  real(dp), dimension(:), allocatable :: cos_theta, r_rk, fc_rk, dfc_rk, &
       G, dG_dcostheta
  real(dp), dimension(3) :: u_ij, u_rs, grad_s, grad_k, f_s, f_k, f_ij
  real(dp), dimension(:,:), allocatable :: u_rk
  real(dp), pointer :: w_e(:)
  real(dp) :: De_ij, R1_ij, R2_ij, Re_ij, S_ij, beta_ij, delta_ij, shift_ij, w_f
  real(dp) :: a0_ijk, c0_2_ijk, d0_2_ijk, R1_rk, R2_rk, de

  type(Dictionary) :: params
  logical, dimension(:), pointer :: atom_mask_pointer
  logical :: has_atom_mask_name
  character(FIELD_LENGTH) :: atom_mask_name
  real(dp) :: r_scale, E_scale
  logical :: do_rescale_r, do_rescale_E

  INIT_ERROR(error)

  if (present(e)) then
     e = 0.0_dp
  endif

  if (present(local_e)) then
     call check_size('Local_E',local_e,(/at%N/),'IPModel_Brenner_Calc', error)
     local_e = 0.0_dp
  endif

  if (present(f)) then 
     call check_size('Force',f,(/3,at%N/),'IPModel_Brenner_Calc', error)
     f = 0.0_dp
  end if

  if (present(virial)) then
     virial = 0.0_dp
  endif

  if (present(local_virial)) then
     call check_size('Local_virial',local_virial,(/9,at%N/),'IPModel_Brenner_Calc', error)
     local_virial = 0.0_dp
     RAISE_ERROR("IPModel_Brenner_Calc: local_virial calculation requested but not supported yet", error)
  endif

  if (.not. assign_pointer(at, "weight", w_e)) nullify(w_e)

  atom_mask_pointer => null()
  if(present(args_str)) then
     call initialise(params)
     call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, help_string="No help yet.  This source file was $LastChangedBy$")
     call param_register(params, 'r_scale', '1.0',r_scale, has_value_target=do_rescale_r, help_string="Recaling factor for distances. Default 1.0.")
     call param_register(params, 'E_scale', '1.0',E_scale, has_value_target=do_rescale_E, help_string="Recaling factor for energy. Default 1.0.")

     if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='IPModel_Brenner_Calc args_str')) then
        RAISE_ERROR("IPModel_Brenner_Calc failed to parse args_str='"//trim(args_str)//"'", error)
     endif
     call finalise(params)

     if( has_atom_mask_name ) then
        if (.not. assign_pointer(at, trim(atom_mask_name) , atom_mask_pointer)) then
           RAISE_ERROR("IPModel_Brenner_Calc did not find "//trim(atom_mask_name)//" property in the atoms object.",error)
        endif
     else
        atom_mask_pointer => null()
     endif
     if (do_rescale_r .or. do_rescale_E) then
        RAISE_ERROR("IPModel_Brenner_Calc: rescaling of potential with r_scale and E_scale not yet implemented!", error)
     end if

  endif

  ! Find maximum number of neighbours and allocate
  ! temporary storage for bond vectors and angles
  max_neighb = 0
  do i=1,at%N
     n = atoms_n_neighbours(at, i)
     if (n > max_neighb) max_neighb = n
  end do

  allocate(cos_theta(max_neighb))
  allocate(r_rk(max_neighb), u_rk(3,max_neighb))
  allocate(fc_rk(max_neighb), dfc_rk(max_neighb))
  allocate(G(max_neighb), dG_dcostheta(max_neighb))

  do i=1,at%N
    if (present(mpi)) then
       if (mpi%active) then
	 if (mod(i-1, mpi%n_procs) /= mpi%my_proc) cycle
       endif
    endif

     if(associated(atom_mask_pointer)) then
        if(.not. atom_mask_pointer(i)) cycle
     endif

     ti = get_type(this%type_of_atomic_num, at%Z(i))

     ! Loop over neighbour of atom i
     do n=1, atoms_n_neighbours(at, i)
        j = atoms_neighbour(at, i, n,r_ij,cosines=u_ij)
        if (r_ij .feq. 0.0_dp) cycle

        tj = get_type(this%type_of_atomic_num, at%Z(j))

        R1_ij = this%R1(ti,tj)
        R2_ij = this%R2(ti,tj)

        if (r_ij < R2_ij) then
           De_ij = this%De(ti,tj)
           Re_ij = this%Re(ti,tj)
           S_ij  = this%S(ti,tj)
           beta_ij = this%beta(ti,tj)
           delta_ij = this%delta(ti,tj)
           shift_ij = this%shift(ti,tj)
        else
           cycle ! Outside range of cutoff function
        end if

        call brenner_cutoff(r_ij, R1_ij, R2_ij, fc_ij, dfc_ij)

        ! Evaluate attractive and repulsive pair potentials
        Vr1 = De_ij/(S_ij-1.0_dp)* & 
             exp(-sqrt(2.0_dp*S_ij)*beta_ij*(r_ij-Re_ij))
        Va1 = De_ij*S_ij/(S_ij-1.0_dp)* & 
             exp(-sqrt(2.0_dp/S_ij)*beta_ij*(r_ij-Re_ij))

        Vr = Vr1*fc_ij
        Va = Va1*fc_ij

        ! Derivatives of Vr and Va wrt r_ij
        dVr_ij = Vr1*(-sqrt(2.0_dp*S_ij)*beta_ij*fc_ij + dfc_ij)
        dVa_ij = Va1*(-sqrt(2.0_dp/S_ij)*beta_ij*fc_ij + dfc_ij)

        ! Repeat twice, first centred on atom i -> B_ij and derivatives
        ! second time centred on atom j -> B_ji and derivatives
        do p = 1,2

           if (p == 1) then
              r = i; s = j
              tr = ti; ts = tj;
              u_rs = u_ij
           else
              r = j; s = i
              tr = tj; ts = ti;
              u_rs = -u_ij ! reverse bond vector 2nd time
           end if

           ! First we loop over neighbours of atom r (i or j)
           ! to get get total bond order factors B_ij and B_ji.
           ! Store r_rk, u_rk, cos_theta, fc_rk, dfc_rk, G, 
           ! and dG_dcostheta for each neighbour to save having
           ! to recalculate them below when we evaluate forces
           G_sum = 0.0_dp
           do m=1, atoms_n_neighbours(at, r)
              k = atoms_neighbour(at, r, m, r_rk(m), cosines=u_rk(:,m))
              if (k == i .or. k == j) cycle
              if (r_rk(m) .feq. 0.0_dp) cycle

              tk = get_type(this%type_of_atomic_num, at%Z(k))

              R1_rk = this%R1(tr,tk)
              R2_rk = this%R2(tr,tk)

              ! First  time: r_rk = r_ik
              ! Second time: r_rk = r_jk
              if (r_rk(m) > R2_rk) cycle ! Outside range of cutoff function

              a0_ijk = this%a0(ti,tj,tk)
              c0_2_ijk = this%c0(ti,tj,tk)*this%c0(ti,tj,tk)
              d0_2_ijk = this%d0(ti,tj,tk)*this%d0(ti,tj,tk)

              ! first  time: theta_ijk
              ! Second time: theta_jik
              call brenner_cutoff(r_rk(m), R1_rk, R2_rk, fc_rk(m), dfc_rk(m))

              cos_theta(m) = u_rs .dot. u_rk(:,m)

              t1 = (1.0_dp + cos_theta(m)); t2 = t1*t1
              G(m) = a0_ijk*(1.0_dp + c0_2_ijk/d0_2_ijk - &
                   c0_2_ijk/(d0_2_ijk + t2))
              dG_dcostheta(m) = 2.0_dp*a0_ijk*c0_2_ijk*t1/ &
                   ((d0_2_ijk + t2)**2)

              G_sum = G_sum + G(m)*fc_rk(m)
           end do

           ! This is the bond order factor
           B(p) = (1.0_dp + G_sum)**(-delta_ij)

           ! Now we need to loop over neigbours again to work out
           ! derivaties of B_ij or B_ji wrt atom positions. We
           ! get gradient from chain rule and the results that:
           !
           ! grad_j(cos_theta_ijk) = (u_ik - cos_theta*u_ij)/r_ij
           ! grad_k(cos_theta_ijk) = (u_ij - cos_theta*u_ik)/r_ik
           ! grad_i(cos_theta_ijk) = -(grad_j + grad_k)
           !
           ! Use previously cached data for each neighbour

           if (associated(w_e)) then
              w_f = 0.5_dp*(w_e(i)+w_e(j))
           else
              w_f = 1.0_dp
           endif

           if (present(f) .or. present(virial)) then

              prefactor = 0.5_dp*w_f*0.5_dp*Va*delta_ij/ &
                   ((1.0_dp + G_sum)**(1.0_dp + delta_ij))

              do m=1, atoms_n_neighbours(at, r)
                 k = atoms_neighbour(at, r, m)
                 if (k == i .or. k == j) cycle
                 if (r_rk(m) > R2_ij) cycle ! Outside range of cutoff function

                 ! force on atom k depends on r_ik or r_jk, hence second term
                 ! First  time: F_k ~ grad_k(theta_ijk) + grad_k(r_ik)
                 ! Second time: F_k ~ grad_k(theta_jik) + grad_k(r_jk)
                 grad_k = (u_rs - cos_theta(m)*u_rk(:,m))/r_rk(m)
                 f_k = -prefactor*(grad_k*dG_dcostheta(m)*fc_rk(m) + G(m)*dfc_rk(m)*u_rk(:,m))
                 if (present(f)) then
                    f(:,k) = f(:,k) + f_k
                 end if

                 ! First  time: F_j ~ grad_j(theta_ijk)
                 ! Second time: F_i ~ grad_i(theta_jik)
                 grad_s = (u_rk(:,m) - cos_theta(m)*u_rs)/r_ij
                 f_s = -prefactor*(grad_s*dG_dcostheta(m)*fc_rk(m))
                 if (present(f)) then
                    f(:,s) = f(:,s) + f_s
                 end if

                 ! First  time: F_i = -(F_k + F_j)
                 ! Second time: F_j = -(F_k + F_i)
                 if (present(f)) then
                    f(:,r) = f(:,r) - (f_k + f_s)
                 end if

                 if (present(virial)) then
                    virial = virial + ((u_rs*r_ij) .outer. f_s) + &
                         ((u_rk(:,m)*r_rk(m)) .outer. f_k)
                 end if

              end do
           end if
        end do

        ! Finally, average bond order factors B_ij and B_ji
        Bbar = 0.5_dp*(B(1) + B(2)) + shift_ij

        ! Add two-body part of f
        if (present(f) .or. present(virial)) then
           t1 = dVr_ij - Bbar*dVa_ij
           f_ij = 0.5_dp*w_f*t1*u_ij

           if (present(f)) then
              f(:,i) = f(:,i) + f_ij
              f(:,j) = f(:,j) - f_ij
           end if

           ! 2 body contribution to virial
           if (present(virial)) then 
              virial = virial - ((u_ij*r_ij) .outer. f_ij)
           end if

        end if

        if (present(e) .or. present(local_e)) then
           de = 0.5_dp*(Vr - Bbar*Va)
           if (present(local_e)) then
              local_e(i) = local_e(i) + 0.5_dp*de
              local_e(j) = local_e(j) + 0.5_dp*de
           end if
           if (present(e)) then
              e = e + de*w_f
           end if
        end if

     end do
  end do

  if (present(mpi)) then
     if (present(e)) e = sum(mpi, e)
     if (present(local_e)) call sum_in_place(mpi, local_e)
     if (present(f)) call sum_in_place(mpi, f)
     if (present(virial)) call sum_in_place(mpi, virial)
  endif

  deallocate(cos_theta)
  deallocate(r_rk, u_rk)
  deallocate(fc_rk, dfc_rk)
  deallocate(G, dG_dcostheta)

end subroutine IPModel_Brenner_Calc


!% Evaluate cutoff function \texttt{fc} and its derivative \texttt{dfc} wrt r
subroutine brenner_cutoff(r, R1, R2, fc, dfc)
  real(dp), intent(in) :: r, R1, R2
  real(dp), intent(out) :: fc, dfc
  
  if (r < R1) then
     fc = 1.0_dp
     dfc = 0.0_dp !fc constant, so d(fc)/dr = 0
  else
     fc = 0.5_dp*(1.0_dp + cos(PI*(r-R1)/(R2-R1)))
     dfc = -0.5_dp*PI*sin(PI*(r-R1)/(R2-R1))/(R2-R1)
  end if
end subroutine brenner_cutoff




!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!% XML param reader functions.
!% Parameters can be given using an external input file or with an internal string.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

subroutine IPModel_startElement_handler(URI, localname, name, attributes)
  character(len=*), intent(in)   :: URI  
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name 
  type(dictionary_t), intent(in) :: attributes

  integer status
  character(len=1024) :: value

  integer ti, tj, tk, Zi, Zj, Zk

  if (name == 'Brenner_params') then ! new Brenner stanza

    if (parse_in_ip) &
      call system_abort("IPModel_startElement_handler entered Brenner_params with parse_in true. Probably a bug in FoX (4.0.1, e.g.)")

    if (parse_matched_label) return ! we already found an exact match for this label

    call QUIP_FoX_get_value(attributes, 'label', value, status)
    if (status /= 0) value = ''

    if (len(trim(parse_ip%label)) > 0) then ! we were passed in a label
      if (value == parse_ip%label) then ! exact match
        parse_matched_label = .true.
        parse_in_ip = .true.
      else ! no match
        parse_in_ip = .false.
      endif
    else ! no label passed in
      parse_in_ip = .true.
    endif

    if (parse_in_ip) then
      if (parse_ip%n_types /= 0) then
	call finalise(parse_ip)
      endif

      call QUIP_FoX_get_value(attributes, "n_types", value, status)
      if (status /= 0) call system_abort ("IPModel_Brenner_read_params_xml cannot find n_types")
      read (value, *) parse_ip%n_types

      allocate(parse_ip%atomic_num(parse_ip%n_types))
      parse_ip%atomic_num = 0

      allocate(parse_ip%R1(parse_ip%n_types,parse_ip%n_types))
      allocate(parse_ip%R2(parse_ip%n_types,parse_ip%n_types))
      allocate(parse_ip%De(parse_ip%n_types,parse_ip%n_types))
      allocate(parse_ip%S(parse_ip%n_types,parse_ip%n_types))
      allocate(parse_ip%beta(parse_ip%n_types,parse_ip%n_types))
      allocate(parse_ip%delta(parse_ip%n_types,parse_ip%n_types))
      allocate(parse_ip%Re(parse_ip%n_types,parse_ip%n_types))
      allocate(parse_ip%a0(parse_ip%n_types,parse_ip%n_types,parse_ip%n_types))
      allocate(parse_ip%c0(parse_ip%n_types,parse_ip%n_types,parse_ip%n_types))
      allocate(parse_ip%d0(parse_ip%n_types,parse_ip%n_types,parse_ip%n_types))
      allocate(parse_ip%shift(parse_ip%n_types,parse_ip%n_types))
    endif

  elseif (parse_in_ip .and. name == 'per_type_data') then

    call QUIP_FoX_get_value(attributes, "type", value, status)
    if (status /= 0) call system_abort ("IPModel_Brenner_read_params_xml cannot find type")
    read (value, *) ti

    call QUIP_FoX_get_value(attributes, "atomic_num", value, status)
    if (status /= 0) call system_abort ("IPModel_Brenner_read_params_xml cannot find atomic_num")
    read (value, *) parse_ip%atomic_num(ti)


    if (allocated(parse_ip%type_of_atomic_num)) deallocate(parse_ip%type_of_atomic_num)
    allocate(parse_ip%type_of_atomic_num(maxval(parse_ip%atomic_num)))
    parse_ip%type_of_atomic_num = 0
    do ti=1, parse_ip%n_types
      if (parse_ip%atomic_num(ti) > 0) &
	parse_ip%type_of_atomic_num(parse_ip%atomic_num(ti)) = ti
    end do

  elseif (parse_in_ip .and. name == 'per_pair_data') then

    call QUIP_FoX_get_value(attributes, "atnum_i", value, status)
    if (status /= 0) call system_abort ("IPModel_Brenner_read_params_xml cannot find atnum_i")
    read (value, *) Zi
    call QUIP_FoX_get_value(attributes, "atnum_j", value, status)
    if (status /= 0) call system_abort ("IPModel_Brenner_read_params_xml cannot find atnum_j")
    read (value, *) Zj

    ti = get_type(parse_ip%type_of_atomic_num,Zi)
    tj = get_type(parse_ip%type_of_atomic_num,Zj)

    call QUIP_FoX_get_value(attributes, "R1", value, status)
    if (status /= 0) call system_abort ("IPModel_Brenner_read_params_xml cannot find R1")
    read (value, *) parse_ip%R1(ti,tj)
    call QUIP_FoX_get_value(attributes, "R2", value, status)
    if (status /= 0) call system_abort ("IPModel_Brenner_read_params_xml cannot find R2")
    read (value, *) parse_ip%R2(ti,tj)
    call QUIP_FoX_get_value(attributes, "De", value, status)
    if (status /= 0) call system_abort ("IPModel_Brenner_read_params_xml cannot find De")
    read (value, *) parse_ip%De(ti,tj)
    call QUIP_FoX_get_value(attributes, "S", value, status)
    if (status /= 0) call system_abort ("IPModel_Brenner_read_params_xml cannot find S")
    read (value, *) parse_ip%S(ti,tj)
    call QUIP_FoX_get_value(attributes, "beta", value, status)
    if (status /= 0) call system_abort ("IPModel_Brenner_read_params_xml cannot find beta")
    read (value, *) parse_ip%beta(ti,tj)
    call QUIP_FoX_get_value(attributes, "delta", value, status)
    if (status /= 0) call system_abort ("IPModel_Brenner_read_params_xml cannot find delta")
    read (value, *) parse_ip%delta(ti,tj)
    call QUIP_FoX_get_value(attributes, "Re", value, status)
    if (status /= 0) call system_abort ("IPModel_Brenner_read_params_xml cannot find Re")
    read (value, *) parse_ip%Re(ti,tj)
    call QUIP_FoX_get_value(attributes, "shift", value, status)
    if (status /= 0) call system_abort ("IPModel_Brenner_read_params_xml cannot find shift")
    read (value, *) parse_ip%shift(ti,tj)


    if (ti /= tj) then
      parse_ip%R1(tj,ti) = parse_ip%R1(ti,tj)
      parse_ip%R2(tj,ti) = parse_ip%R2(ti,tj)
      parse_ip%De(tj,ti) = parse_ip%De(ti,tj)
      parse_ip%S(tj,ti) = parse_ip%S(ti,tj)
      parse_ip%beta(tj,ti) = parse_ip%beta(ti,tj)
      parse_ip%delta(tj,ti) = parse_ip%delta(ti,tj)
      parse_ip%Re(tj,ti) = parse_ip%Re(ti,tj)
      parse_ip%shift(tj,ti) = parse_ip%shift(ti,tj)
    endif

    parse_ip%cutoff = maxval(parse_ip%R2)

  elseif (parse_in_ip .and. name == 'per_triplet_data') then

    call QUIP_FoX_get_value(attributes, "atnum_c", value, status)
    if (status /= 0) call system_abort ("IPModel_Brenner_read_params_xml cannot find atnum_c")
    read (value, *) Zi
    call QUIP_FoX_get_value(attributes, "atnum_j", value, status)
    if (status /= 0) call system_abort ("IPModel_Brenner_read_params_xml cannot find atnum_j")
    read (value, *) Zj
    call QUIP_FoX_get_value(attributes, "atnum_k", value, status)
    if (status /= 0) call system_abort ("IPModel_Brenner_read_params_xml cannot find atnum_k")
    read (value, *) Zk

    ti = get_type(parse_ip%type_of_atomic_num,Zi)
    tj = get_type(parse_ip%type_of_atomic_num,Zj)
    tk = get_type(parse_ip%type_of_atomic_num,Zk)

    call QUIP_FoX_get_value(attributes, "a0", value, status)
    if (status /= 0) call system_abort ("IPModel_Brenner_read_params_xml cannot find a0")
    read (value, *) parse_ip%a0(ti,tj,tk)
    call QUIP_FoX_get_value(attributes, "c0", value, status)
    if (status /= 0) call system_abort ("IPModel_Brenner_read_params_xml cannot find c0")
    read (value, *) parse_ip%c0(ti,tj,tk)
    call QUIP_FoX_get_value(attributes, "d0", value, status)
    if (status /= 0) call system_abort ("IPModel_Brenner_read_params_xml cannot find d0")
    read (value, *) parse_ip%d0(ti,tj,tk)

    if (tj /= tk) then
      parse_ip%a0(ti,tk,tj) = parse_ip%a0(ti,tj,tk)
      parse_ip%c0(ti,tk,tj) = parse_ip%c0(ti,tj,tk)
      parse_ip%d0(ti,tk,tj) = parse_ip%d0(ti,tj,tk)
    endif
  endif

end subroutine IPModel_startElement_handler

subroutine IPModel_endElement_handler(URI, localname, name)
  character(len=*), intent(in)   :: URI  
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name 

  if (parse_in_ip) then
    if (name == 'Brenner_params') then
      parse_in_ip = .false.
    endif
  endif

end subroutine IPModel_endElement_handler

subroutine IPModel_Brenner_read_params_xml(this, param_str)
  type(IPModel_Brenner), intent(inout), target :: this
  character(len=*), intent(in) :: param_str

  type (xml_t) :: fxml

  if (len(trim(param_str)) <= 0) return

  parse_ip => this
  parse_in_ip = .false.
  parse_matched_label = .false.

  call open_xml_string(fxml, param_str)

  call parse(fxml, &
    startElement_handler = IPModel_startElement_handler, &
    endElement_handler = IPModel_endElement_handler)

  call close_xml_t(fxml)

  if (this%n_types == 0) then
    call system_abort("IPModel_Brenner_read_params_xml parsed file, but n_types = 0")
  endif

end subroutine IPModel_Brenner_read_params_xml

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!% Printing of Brenner potential parameters: number of different types, cutoff radius, atomic numbers, etc.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

subroutine IPModel_Brenner_Print (this, file)
  type(IPModel_Brenner), intent(in) :: this
  type(Inoutput), intent(inout),optional :: file

  integer :: ti, tj, tk

  call Print("IPModel_Brenner : n_types = " // this%n_types // " cutoff = " // this%cutoff, file=file)

  do ti=1, this%n_types
    call Print("IPModel_Brenner : type " // ti // " atomic_num " // this%atomic_num(ti), file=file)
    call verbosity_push_decrement()
    do tj=1, this%n_types
      call Print ("IPModel_Brenner : pair "// ti // " " // tj // " R1 " // this%R1(ti,tj) &
           // " R2 " // this%R2(ti,tj)// " De " // this%De(ti,tj)// " S " // this%S(ti,tj) &
           // " beta " // this%beta(ti,tj)// " delta " // this%delta(ti,tj) &
           // " Re " // this%Re(ti,tj)// " shift " // this%shift(ti,tj), file=file)
      do tk=1,this%n_types
         call Print("IPModel_Brenner : triplet "// ti //" "//tj//" "//tk//" a0 "//this%a0(ti,tj,tk)&
              //" c0 "//this%c0(ti,tj,tj)//" d0 "//this%d0(ti,tj,tk),file=file)
      end do
    end do
    call verbosity_pop()
  end do
    
end subroutine IPModel_Brenner_Print

end module IPModel_Brenner_module
