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
!X IPModel_Tersoff module 
!X  
!% Module for computing energy using the Tersoff potential for C/Si/Ge and alloys
!% (Ref. J. Tersoff,  Phys. Rev. B {\bf 39}, 5566 (1989)).
!% The IPModel_Tersoff object contains all the parameters read from a 'Tersoff_params'
!% XML stanza.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#include "error.inc"

module IPModel_Tersoff_module

use libatoms_module

use mpi_context_module
use QUIP_Common_module

implicit none
private

include 'IPModel_interface.h'

public :: IPModel_Tersoff
type IPModel_Tersoff
  integer :: n_types = 0         !% Number of atomic types 
  integer, allocatable :: atomic_num(:), type_of_atomic_num(:)  !% Atomic number dimensioned as \texttt{n_types} 

  real(dp) :: cutoff = 0.0_dp    !% Cutoff for computing connections

  real(dp), allocatable :: A(:), B(:), lambda(:), mu(:), beta(:)  !% IP parameters dimensioned as \texttt{n_types}
  real(dp), allocatable :: n(:), c(:), d(:), h(:), R(:), S(:)     !% IP parameters dimensioned as \texttt{n_types}
  real(dp), allocatable :: chi(:,:)                               !% IP parameter depending on the pair interaction (mixing hentalpy), dimensioned as \texttt{(n_types,n_types)}

  character(len=FIELD_LENGTH) :: label

end type IPModel_Tersoff

logical, private :: parse_in_ip, parse_matched_label
type(IPModel_Tersoff), private, pointer :: parse_ip

interface Initialise
  module procedure IPModel_Tersoff_Initialise_str
end interface Initialise

interface Finalise
  module procedure IPModel_Tersoff_Finalise
end interface Finalise

interface Print
  module procedure IPModel_Tersoff_Print
end interface Print

interface Calc
  module procedure IPModel_Tersoff_Calc
end interface Calc

contains

subroutine IPModel_Tersoff_Initialise_str(this, args_str, param_str)
  type(IPModel_Tersoff), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str

  type(Dictionary) :: params

  call Finalise(this)

  call initialise(params)
  this%label=''
  call param_register(params, 'label', '', this%label, help_string="No help yet.  This source file was $LastChangedBy$")
  if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_Tersoff_Initialise_str args_str')) then
    call system_abort("IPModel_Tersoff_Initialise_str failed to parse label from args_str="//trim(args_str))
  endif
  call finalise(params)

  call IPModel_Tersoff_read_params_xml(this, param_str)

end subroutine IPModel_Tersoff_Initialise_str

subroutine IPModel_Tersoff_Finalise(this)
  type(IPModel_Tersoff), intent(inout) :: this

  if (allocated(this%atomic_num)) deallocate(this%atomic_num)
  if (allocated(this%type_of_atomic_num)) deallocate(this%type_of_atomic_num)

  if (allocated(this%A)) deallocate(this%A)
  if (allocated(this%B)) deallocate(this%B)
  if (allocated(this%lambda)) deallocate(this%lambda)
  if (allocated(this%mu)) deallocate(this%mu)
  if (allocated(this%beta)) deallocate(this%beta)
  if (allocated(this%n)) deallocate(this%n)
  if (allocated(this%c)) deallocate(this%c)
  if (allocated(this%d)) deallocate(this%d)
  if (allocated(this%h)) deallocate(this%h)
  if (allocated(this%R)) deallocate(this%R)
  if (allocated(this%S)) deallocate(this%S)
  if (allocated(this%chi)) deallocate(this%chi)

  this%n_types = 0
  this%label = ''
end subroutine IPModel_Tersoff_Finalise

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!% The potential calculator. It computes energy, forces and virial.
!% Derivatives are taken from M. Tang, Ph.D. Thesis, MIT 1995.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
subroutine IPModel_Tersoff_Calc(this, at, e, local_e, f, virial, local_virial, args_str, mpi, error)
  type(IPModel_Tersoff), intent(inout) :: this
  type(Atoms), intent(in) :: at
  real(dp), intent(out), optional :: e, local_e(:) !% \texttt{e} = System total energy, \texttt{local_e} = energy of each atom, vector dimensioned as \texttt{at%N}.  
  real(dp), intent(out), optional :: f(:,:), local_virial(:,:)   !% Forces, dimensioned as \texttt{f(3,at%N)}, local virials, dimensioned as \texttt{local_virial(9,at%N)} 
  real(dp), intent(out), optional :: virial(3,3)   !% Virial
  character(len=*), intent(in), optional :: args_str 
  type(MPI_Context), intent(in), optional :: mpi
  integer, intent(out), optional :: error

  real(dp), pointer :: w_e(:)
  integer i, ji, j, ki, k
  integer ti, tj, tk
  real(dp) :: dr_ij(3), dr_ij_mag, dr_ik(3), dr_ik_mag, dr_jk(3), dr_jk_mag
  real(dp) :: dr_ij_diff(3), dr_ik_diff(3)
  real(dp) :: w_f

  real(dp), allocatable :: z_ij(:)
  real(dp) :: AA_ij, BB_ij, lambda_ij, mu_ij, beta_i, n_i, c_i, d_i, h_i, R_ij, S_ij, chi_ij;
  real(dp) :: zeta_ij, b_ij, V_ij_R, V_ij_A, f_C_ij;
  real(dp) :: R_ik, S_ik, cos_theta_ijk;
  real(dp) :: exp_lambda_ij, exp_mu_ij

  real(dp) :: de

  real(dp) :: dr_ij_mag_dr_i(3), dr_ij_mag_dr_j(3)
  real(dp) :: dr_ik_mag_dr_i(3), dr_ik_mag_dr_k(3)
  real(dp) :: dr_jk_mag_dr_j(3), dr_jk_mag_dr_k(3)
  real(dp) :: f_C_ij_d, dV_ij_R_dr_ij_mag 
  real(dp) :: f_C_ik

  real(dp) :: dcos_theta_ijk_dr_ij_mag, dcos_theta_ijk_dr_ik_mag, dcos_theta_ijk_dr_jk_mag
  real(dp) :: dg_dcos_theta_ijk

  real(dp) :: beta_i_db_ij_dz_ij
  real(dp) :: dzeta_ij_dr_ij_mag, dzeta_ij_dr_ik_mag, dzeta_ij_dr_jk_mag
  real(dp) :: db_ij_dr_ij_mag, db_ij_dr_ik_mag, db_ij_dr_jk_mag
  real(dp) :: dV_ij_A_dr_ij_mag, dV_ij_A_dr_ik_mag, dV_ij_A_dr_jk_mag
  real(dp) :: virial_i(3,3)

  type(Dictionary) :: params
  logical, dimension(:), pointer :: atom_mask_pointer
  logical :: has_atom_mask_name
  character(FIELD_LENGTH) :: atom_mask_name

  real(dp) :: r_scale, E_scale
  logical :: do_rescale_r, do_rescale_E

  INIT_ERROR(error)

  if (present(e)) e = 0.0_dp

  if (present(local_e)) then
     call check_size('Local_E',local_e,(/at%N/),'IPModel_Tersoff_Calc', error)
     local_e = 0.0_dp
  endif

  if (present(f)) then 
     call check_size('Force',f,(/3,at%N/),'IPModel_Tersoff_Calc', error)
     f = 0.0_dp
  end if

  if (present(virial)) then
     virial = 0.0_dp
  endif

  if (present(local_virial)) then
     call check_size('Local_virial',local_virial,(/9,at%N/),'IPModel_Tersoff_Calc', error)
     local_virial = 0.0_dp
  endif

  if (.not. assign_pointer(at, "weight", w_e)) nullify(w_e)

  atom_mask_pointer => null()
  if(present(args_str)) then
     call initialise(params)
     call param_register(params, 'atom_mask_name', 'NONE',atom_mask_name,has_value_target=has_atom_mask_name, help_string="No help yet.  This source file was $LastChangedBy$")
     call param_register(params, 'do_rescale_r', 'F',do_rescale_r, help_string="If true, rescale distances by factor r_scale.")
     call param_register(params, 'r_scale', '1.0',r_scale, help_string="Recaling factor for distances. Default 1.0.")
     call param_register(params, 'do_rescale_E', 'F',do_rescale_E, help_string="If true, rescale energy by factor r_scale.")
     call param_register(params, 'E_scale', '1.0',E_scale, help_string="Recaling factor for energy. Default 1.0.")
     if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='IPModel_Tersoff_Calc args_str')) &
     call system_abort("IPModel_Tersoff_Calc failed to parse args_str='"//trim(args_str)//"'")
     call finalise(params)


     if( has_atom_mask_name ) then
        if (.not. assign_pointer(at, trim(atom_mask_name) , atom_mask_pointer)) &
        call system_abort("IPModel_Tersoff_Calc did not find "//trim(atom_mask_name)//" property in the atoms object.")
     else
        atom_mask_pointer => null()
     endif
  endif

  if (do_rescale_r) call print('rescaling distances by factor '//r_scale)
  if (do_rescale_E) call print('rescaling energy by factor '//E_scale)

  do i=1, at%N
    if (present(mpi)) then
       if (mpi%active) then
	 if (mod(i-1, mpi%n_procs) /= mpi%my_proc) cycle
       endif
    endif

     if(associated(atom_mask_pointer)) then
        if(.not. atom_mask_pointer(i)) cycle
     endif

    allocate(z_ij(atoms_n_neighbours(at,i)))
    ti = get_type(this%type_of_atomic_num, at%Z(i))

    ! calculate z_ij = beta_i zeta_ij
    do ji=1, atoms_n_neighbours(at, i)
      j = atoms_neighbour(at, i, ji, dr_ij_mag, cosines = dr_ij)
      if (dr_ij_mag .feq. 0.0_dp) cycle

      if (do_rescale_r) dr_ij_mag = dr_ij_mag*r_scale

      tj = get_type(this%type_of_atomic_num, at%Z(j))

      S_ij = sqrt(this%S(ti)*this%S(tj))
      if (dr_ij_mag < S_ij) then
	beta_i = this%beta(ti)
	R_ij = sqrt(this%R(ti)*this%R(tj))
	c_i = this%c(ti)
	d_i = this%d(ti)
	h_i = this%h(ti)
      else
	cycle
      endif

      zeta_ij = 0.0_dp
      do ki=1, atoms_n_neighbours(at, i)
	if (ki == ji) cycle
	k = atoms_neighbour(at, i, ki, dr_ik_mag, cosines = dr_ik)
	if (dr_ik_mag .feq. 0.0_dp) cycle
        if (do_rescale_r) dr_ik_mag = dr_ik_mag*r_scale

	tk = get_type(this%type_of_atomic_num, at%Z(k))

	R_ik = sqrt(this%R(ti)*this%R(tk))
	S_ik = sqrt(this%S(ti)*this%S(tk))
	cos_theta_ijk = sum(dr_ij*dr_ik)
	zeta_ij = zeta_ij + f_C(dr_ik_mag, R_ik, S_ik)*g(cos_theta_ijk, c_i, d_i, h_i)
      end do
      z_ij(ji) = beta_i*zeta_ij
    end do ! ji (calc z_ij)

    ! calculate energies, forces, etc
    do ji=1, atoms_n_neighbours(at, i)
      j = atoms_neighbour(at, i, ji, dr_ij_mag, dr_ij_diff, dr_ij)
      if (dr_ij_mag .feq. 0.0_dp) cycle

      if (do_rescale_r) then
         dr_ij_mag = dr_ij_mag * r_scale
         dr_ij_diff = dr_ij_diff * r_scale
      end if

      tj = get_type(this%type_of_atomic_num, at%Z(j))

      S_ij = sqrt(this%S(ti)*this%S(tj))
      if (dr_ij_mag < S_ij) then
	n_i = this%n(ti)
	R_ij = sqrt(this%R(ti)*this%R(tj))
	AA_ij = sqrt(this%A(ti)*this%A(tj))
	BB_ij = sqrt(this%B(ti)*this%B(tj))
	lambda_ij = 0.5_dp*(this%lambda(ti)+this%lambda(tj))
	mu_ij = 0.5_dp*(this%mu(ti)+this%mu(tj))
	chi_ij = this%chi(ti,tj)
      else
	cycle
      endif

      b_ij = chi_ij*(1.0_dp + z_ij(ji)**n_i)**(-1.0_dp/(2*n_i))
      f_C_ij = f_C(dr_ij_mag, R_ij, S_ij)
      exp_lambda_ij = exp(-lambda_ij*dr_ij_mag)
      exp_mu_ij = exp(-mu_ij*dr_ij_mag)
      V_ij_R = f_C_ij * (AA_ij * exp_lambda_ij)
      V_ij_A = f_C_ij * (-b_ij*BB_ij * exp_mu_ij)

      if (associated(w_e)) then
	w_f = 0.5_dp*(w_e(i)+w_e(j))
      else
	w_f = 1.0_dp
      endif

      if (present(e) .or. present(local_e)) then
	! factor of 0.5 because Tersoff definition goes over each pair only once
	de = 0.5_dp*(V_ij_R + v_ij_A)
	if (present(local_e)) then
	  local_e(i) = local_e(i) + de
	endif
	if (present(e)) then
	  e = e + de*w_f
	endif
      endif

      if (present(f) .or. present(virial) .or. present(local_virial)) then
	dr_ij_mag_dr_i = -dr_ij
	dr_ij_mag_dr_j = dr_ij

	f_C_ij_d = f_C_d(dr_ij_mag, R_ij, S_ij)
	dV_ij_R_dr_ij_mag = f_C_ij_d*AA_ij*exp_lambda_ij + f_C_ij*AA_ij*(-lambda_ij)*exp_lambda_ij
	if (present(f)) then
	  f(:,i) = f(:,i) - 0.5_dp*w_f * dV_ij_R_dr_ij_mag * dr_ij_mag_dr_i(:)
	  f(:,j) = f(:,j) - 0.5_dp*w_f * dV_ij_R_dr_ij_mag * dr_ij_mag_dr_j(:)
	endif
	if (present(virial) .or. present(local_virial)) virial_i = 0.5_dp*w_f * dV_ij_R_dr_ij_mag * (dr_ij .outer. dr_ij) * dr_ij_mag
        if (present(virial)) virial = virial - virial_i
        if (present(local_virial)) local_virial(:,i) = local_virial(:,i) - reshape(virial_i,(/9/))

	if (z_ij(ji) .fne. 0.0_dp) then
	  beta_i_db_ij_dz_ij = beta_i* ( -0.5_dp*chi_ij * &
	    ( ( 1.0_dp+z_ij(ji)**n_i)**(-(1.0_dp+2.0_dp*n_i)/(2.0_dp*n_i)) ) * &
	    ( z_ij(ji)**(n_i-1.0_dp) ) )
	else
	  beta_i_db_ij_dz_ij = -1.0_dp
	endif

	do ki=1, atoms_n_neighbours(at,i)
	  if (ki == ji) cycle
	  k = atoms_neighbour(at, i, ki, dr_ik_mag, dr_ik_diff, dr_ik)
	  if (dr_ik_mag .feq. 0.0_dp) cycle

          if (do_rescale_r) then
             dr_ik_mag = dr_ik_mag * r_scale
             dr_ik_diff = dr_ik_diff * r_scale
          end if

	  tk = get_type(this%type_of_atomic_num, at%Z(k))

	  dr_ik_mag_dr_i = -dr_ik
	  dr_ik_mag_dr_k = dr_ik

	  cos_theta_ijk = sum(dr_ij*dr_ik)

	  dcos_theta_ijk_dr_ij_mag = 1.0_dp/dr_ik_mag - cos_theta_ijk/dr_ij_mag
	  dcos_theta_ijk_dr_ik_mag = 1.0_dp/dr_ij_mag - cos_theta_ijk/dr_ik_mag

	  dr_jk = dr_ik_diff - dr_ij_diff
	  dr_jk_mag = sqrt(sum(dr_jk*dr_jk))
	  dr_jk = dr_jk/dr_jk_mag

	  dr_jk_mag_dr_j = -dr_jk
	  dr_jk_mag_dr_k = dr_jk

	  dcos_theta_ijk_dr_jk_mag = -dr_jk_mag/(dr_ij_mag*dr_ik_mag)
	  dg_dcos_theta_ijk = dg_dcos_theta(cos_theta_ijk, c_i, d_i, h_i)

	  R_ik = sqrt(this%R(ti)*this%R(tk))
	  S_ik = sqrt(this%S(ti)*this%S(tk))
	  f_C_ik = f_C(dr_ik_mag, R_ik, S_ik)

	  dzeta_ij_dr_ij_mag = f_C_ik * dg_dcos_theta_ijk * dcos_theta_ijk_dr_ij_mag
	  dzeta_ij_dr_jk_mag = f_C_ik * dg_dcos_theta_ijk * dcos_theta_ijk_dr_jk_mag

	  dzeta_ij_dr_ik_mag = f_C_ik * dg_dcos_theta_ijk * dcos_theta_ijk_dr_ik_mag + &
	    f_C_d(dr_ik_mag, R_ik, S_ik) * g(cos_theta_ijk, c_i, d_i, h_i)

	  db_ij_dr_ij_mag = beta_i_db_ij_dz_ij*dzeta_ij_dr_ij_mag
	  db_ij_dr_ik_mag = beta_i_db_ij_dz_ij*dzeta_ij_dr_ik_mag
	  db_ij_dr_jk_mag = beta_i_db_ij_dz_ij*dzeta_ij_dr_jk_mag

	  dV_ij_A_dr_ij_mag = f_C_ij*(-db_ij_dr_ij_mag*BB_ij*exp_mu_ij)
	  dV_ij_A_dr_ik_mag = f_C_ij*(-db_ij_dr_ik_mag*BB_ij*exp_mu_ij)
	  dV_ij_A_dr_jk_mag = f_C_ij*(-db_ij_dr_jk_mag*BB_ij*exp_mu_ij)

	  if (present(f)) then
	    f(:,i) = f(:,i) - w_f*( 0.5_dp * dV_ij_A_dr_ik_mag*dr_ik_mag_dr_i(:) + &
				    0.5_dp * dV_ij_A_dr_ij_mag*dr_ij_mag_dr_i(:) )

	    f(:,j) = f(:,j) - w_f*( 0.5_dp * dV_ij_A_dr_jk_mag*dr_jk_mag_dr_j(:) + &
				    0.5_dp * dV_ij_A_dr_ij_mag*dr_ij_mag_dr_j(:) )

	    f(:,k) = f(:,k) - w_f*( 0.5_dp*dV_ij_A_dr_ik_mag*dr_ik_mag_dr_k(:) + &
				    0.5_dp*dV_ij_A_dr_jk_mag*dr_jk_mag_dr_k(:) )
	  end if

	  if (present(virial) .or. present(local_virial)) virial_i = &
	    w_f*0.5_dp*( dV_ij_A_dr_ij_mag*(dr_ij .outer. dr_ij)*dr_ij_mag + &
            dV_ij_A_dr_ik_mag*(dr_ik .outer. dr_ik)*dr_ik_mag + &
            dV_ij_A_dr_jk_mag*(dr_jk .outer. dr_jk)*dr_jk_mag)

          if (present(virial)) virial = virial - virial_i
          if (present(local_virial)) local_virial(:,i) = local_virial(:,i) - reshape(virial_i,(/9/))

	end do ! ki
	dV_ij_A_dr_ij_mag = f_C_ij_d*(-b_ij*BB_ij*exp_mu_ij) + &
	  f_C_ij*(-b_ij*BB_ij*(-mu_ij)*exp_mu_ij)

	if (present(f)) then
	  f(:,i) = f(:,i) - 0.5_dp*w_f * dV_ij_A_dr_ij_mag * dr_ij_mag_dr_i(:)
	  f(:,j) = f(:,j) - 0.5_dp*w_f * dV_ij_A_dr_ij_mag * dr_ij_mag_dr_j(:)
	endif
	if (present(virial)) virial_i = &
            0.5_dp*w_f * dV_ij_A_dr_ij_mag*(dr_ij .outer. dr_ij)*dr_ij_mag
        if (present(virial)) virial = virial - virial_i
        if (present(local_virial)) local_virial(:,i) = local_virial(:,i) - reshape(virial_i,(/9/))

      endif ! present(f)

    end do ! ji
    deallocate(z_ij)
  end do ! i

  if (present(mpi)) then
     if (present(e)) e = sum(mpi, e)
     if (present(local_e)) call sum_in_place(mpi, local_e)
     if (present(f)) call sum_in_place(mpi, f)
     if (present(virial)) call sum_in_place(mpi, virial)
     if (present(local_virial)) call sum_in_place(mpi, local_virial)
  endif

  if (do_rescale_r) then
     if (present(f)) f = f*r_scale
  end if

  if (do_rescale_E) then
     if (present(e)) e = e*E_scale
     if (present(f)) f = f*E_scale
     if (present(virial)) virial=virial*E_scale
     if (present(local_virial)) local_virial=local_virial*E_scale
  end if

  atom_mask_pointer => null()

end subroutine IPModel_Tersoff_Calc

!% Compute the cutoff function.
function f_C(r, RR, SS)
  real(dp), intent(in) :: r
  real(dp), intent(in) :: RR, SS
  real(dp) :: f_C

  if (r <= RR) then
    f_C = 1.0_dp
  else if (r <= SS) then
    f_C = 0.5+0.5*cos(PI*(r-RR)/(SS-RR))
  else
    f_C = 0.0_dp
  endif
end function f_C

!% Compute the derivative of the cutoff function.
function f_C_d(r, RR, SS)
  real(dp), intent(in) :: r
  real(dp), intent(in) :: RR, SS
  real(dp) :: f_C_d

  if (r <= RR) then
    f_C_d = 0.0_dp
  else if (r <= SS) then
    f_C_d = -0.5*PI/(SS-RR) * sin(PI*(r-RR)/(SS-RR))
  else
    f_C_d = 0.0_dp
  endif
end function f_C_d

!% Compute the angular term g for the bond-order coefficient.
function g(cos_theta, c, d, h)
  real(dp), intent(in) :: cos_theta
  real(dp), intent(in) :: c, d, h
  real(dp) :: g

  g = 1.0 + (c*c)/(d*d) - (c*c)/(d*d+(h-cos_theta)*(h-cos_theta));
end function g

!% Compute the derivative of the angular term for the bond-order coefficient.
function dg_dcos_theta(cos_theta, c, d, h)
  real(dp), intent(in) :: cos_theta
  real(dp), intent(in) :: c, d, h
  real(dp) :: dg_dcos_theta

  real(dp) :: t

  t = d*d + (h-cos_theta)*(h-cos_theta);
  dg_dcos_theta = (-2.0_dp*c*c*(h-cos_theta))/(t*t)

end function dg_dcos_theta

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!X XML param reader functions
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

subroutine IPModel_startElement_handler(URI, localname, name, attributes)
  character(len=*), intent(in)   :: URI  
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name 
  type(dictionary_t), intent(in) :: attributes

  integer status
  character(len=1024) :: value

  integer ti, tj

  if (name == 'Tersoff_params') then ! new Tersoff stanza

    if (parse_in_ip) &
      call system_abort("IPModel_startElement_handler entered Tersoff_params with parse_in true. Probably a bug in FoX (4.0.1, e.g.)")

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
      if (status /= 0) call system_abort ("IPModel_Tersoff_read_params_xml cannot find n_types")
      read (value, *) parse_ip%n_types

      allocate(parse_ip%atomic_num(parse_ip%n_types))
      parse_ip%atomic_num = 0

      allocate(parse_ip%A(parse_ip%n_types))
      allocate(parse_ip%B(parse_ip%n_types))
      allocate(parse_ip%lambda(parse_ip%n_types))
      allocate(parse_ip%mu(parse_ip%n_types))
      allocate(parse_ip%beta(parse_ip%n_types))
      allocate(parse_ip%n(parse_ip%n_types))
      allocate(parse_ip%c(parse_ip%n_types))
      allocate(parse_ip%d(parse_ip%n_types))
      allocate(parse_ip%h(parse_ip%n_types))
      allocate(parse_ip%R(parse_ip%n_types))
      allocate(parse_ip%S(parse_ip%n_types))
      allocate(parse_ip%chi(parse_ip%n_types,parse_ip%n_types))

    endif

  elseif (parse_in_ip .and. name == 'per_type_data') then

    call QUIP_FoX_get_value(attributes, "type", value, status)
    if (status /= 0) call system_abort ("IPModel_Tersoff_read_params_xml cannot find type")
    read (value, *) ti

    call QUIP_FoX_get_value(attributes, "atomic_num", value, status)
    if (status /= 0) call system_abort ("IPModel_Tersoff_read_params_xml cannot find atomic_num")
    read (value, *) parse_ip%atomic_num(ti)

    call QUIP_FoX_get_value(attributes, "A", value, status)
    if (status /= 0) call system_abort ("IPModel_Tersoff_read_params_xml cannot find A")
    read (value, *) parse_ip%A(ti)
    call QUIP_FoX_get_value(attributes, "B", value, status)
    if (status /= 0) call system_abort ("IPModel_Tersoff_read_params_xml cannot find B")
    read (value, *) parse_ip%B(ti)
    call QUIP_FoX_get_value(attributes, "lambda", value, status)
    if (status /= 0) call system_abort ("IPModel_Tersoff_read_params_xml cannot find lambda")
    read (value, *) parse_ip%lambda(ti)
    call QUIP_FoX_get_value(attributes, "mu", value, status)
    if (status /= 0) call system_abort ("IPModel_Tersoff_read_params_xml cannot find mu")
    read (value, *) parse_ip%mu(ti)
    call QUIP_FoX_get_value(attributes, "beta", value, status)
    if (status /= 0) call system_abort ("IPModel_Tersoff_read_params_xml cannot find beta")
    read (value, *) parse_ip%beta(ti)
    call QUIP_FoX_get_value(attributes, "n", value, status)
    if (status /= 0) call system_abort ("IPModel_Tersoff_read_params_xml cannot find n")
    read (value, *) parse_ip%n(ti)
    call QUIP_FoX_get_value(attributes, "c", value, status)
    if (status /= 0) call system_abort ("IPModel_Tersoff_read_params_xml cannot find c")
    read (value, *) parse_ip%c(ti)
    call QUIP_FoX_get_value(attributes, "d", value, status)
    if (status /= 0) call system_abort ("IPModel_Tersoff_read_params_xml cannot find d")
    read (value, *) parse_ip%d(ti)
    call QUIP_FoX_get_value(attributes, "h", value, status)
    if (status /= 0) call system_abort ("IPModel_Tersoff_read_params_xml cannot find h")
    read (value, *) parse_ip%h(ti)
    call QUIP_FoX_get_value(attributes, "R", value, status)
    if (status /= 0) call system_abort ("IPModel_Tersoff_read_params_xml cannot find R")
    read (value, *) parse_ip%R(ti)
    call QUIP_FoX_get_value(attributes, "S", value, status)
    if (status /= 0) call system_abort ("IPModel_Tersoff_read_params_xml cannot find S")
    read (value, *) parse_ip%S(ti)

    parse_ip%cutoff = maxval(parse_ip%S)

    if (allocated(parse_ip%type_of_atomic_num)) deallocate(parse_ip%type_of_atomic_num)
    allocate(parse_ip%type_of_atomic_num(maxval(parse_ip%atomic_num)))
    parse_ip%type_of_atomic_num = 0
    do ti=1, parse_ip%n_types
      if (parse_ip%atomic_num(ti) > 0) &
	parse_ip%type_of_atomic_num(parse_ip%atomic_num(ti)) = ti
    end do

  elseif (parse_in_ip .and. name == 'per_pair_data') then

    call QUIP_FoX_get_value(attributes, "type1", value, status)
    if (status /= 0) call system_abort ("IPModel_Tersoff_read_params_xml cannot find type1")
    read (value, *) ti
    call QUIP_FoX_get_value(attributes, "type2", value, status)
    if (status /= 0) call system_abort ("IPModel_Tersoff_read_params_xml cannot find type2")
    read (value, *) tj

    call QUIP_FoX_get_value(attributes, "chi", value, status)
    if (status /= 0) call system_abort ("IPModel_Tersoff_read_params_xml cannot find chi")
    read (value, *) parse_ip%chi(ti,tj)
    if (ti /= tj) then
      parse_ip%chi(tj,ti) = parse_ip%chi(ti,tj)
    endif
  endif

end subroutine IPModel_startElement_handler

subroutine IPModel_endElement_handler(URI, localname, name)
  character(len=*), intent(in)   :: URI  
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name 

  if (parse_in_ip) then
    if (name == 'Tersoff_params') then
      parse_in_ip = .false.
    endif
  endif

end subroutine IPModel_endElement_handler

subroutine IPModel_Tersoff_read_params_xml(this, param_str)
  type(IPModel_Tersoff), intent(inout), target :: this
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
    call system_abort("IPModel_Tersoff_read_params_xml parsed file, but n_types = 0")
  endif

end subroutine IPModel_Tersoff_read_params_xml

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!% Printing of Tersoff parameters: number of different types, cutoff radius, atomic numbers, ect.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

subroutine IPModel_Tersoff_Print (this, file)
  type(IPModel_Tersoff), intent(in) :: this
  type(Inoutput), intent(inout),optional :: file

  integer :: ti, tj

  call Print("IPModel_Tersoff : n_types = " // this%n_types // " cutoff = " // this%cutoff, file=file)

  do ti=1, this%n_types
    call Print("IPModel_Tersoff : type " // ti // " atomic_num " // this%atomic_num(ti), file=file)
    call verbosity_push_decrement()
    call Print ("IPModel_Tersoff : params A "//this%A(ti) //" B "// this%B(ti)//" lambda "// this%lambda(ti)//" mu "// &
      this%mu(ti)// " beta " // this%beta(ti), file=file)
    call Print ("IPModel_Tersoff : params n "// this%n(ti)//" c "// this%c(ti)//" d "// this%d(ti)//" h "// this%h(ti) &
      // " R "// this%R(ti)//" S " // this%S(ti), file=file)
    do tj=1, this%n_types
      call Print ("IPModel_tersoff : pair "// ti // " " // tj // " chi " // this%chi(ti,tj), file=file)
    end do
    call verbosity_pop()
  end do
    
end subroutine IPModel_Tersoff_Print

end module IPModel_Tersoff_module
