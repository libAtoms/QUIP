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
! H0 X   Public License, version 2, https://www.gnu.org/copyleft/gpl.html
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
! H0 X   https://www.libatoms.org
! H0 X
! H0 X  Additional contributions by
! H0 X    Alessio Comisso, Chiara Gattinoni, and Gianpietro Moras
! H0 X
! H0 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

!X
!X IPModel_ZBL
!X
!% Module for Ziegler-Biersack-Littmark potential for
!% screened nuclear repulsion 
!%
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#include "error.inc"

module IPModel_ZBL_module

use error_module
use system_module, only : dp, inoutput, print, verbosity_push_decrement, verbosity_pop, operator(//)
use units_module
use dictionary_module
use paramreader_module
use linearalgebra_module
use atoms_types_module
use atoms_module
use periodictable_module, only : ElementCovRad

use mpi_context_module
use QUIP_Common_module

implicit none
private

include 'IPModel_interface.h'

public :: IPModel_ZBL
type IPModel_ZBL
  real(dp) :: cutoff = 0.0
  real(dp) :: use_cutoff = 0.0
  real(dp) :: cutoff_width = 0.0
  logical :: shift_cutoff = .false., cutoff_scale_cov_rad = .false.
  real(dp) :: a_pre_exp = 0.46850
  real(dp) :: a_exp = 0.23
  real(dp) :: p_pre_exp_1 = 0.18175
  real(dp) :: p_exp_1 = -3.19980
  real(dp) :: p_pre_exp_2 = 0.50986
  real(dp) :: p_exp_2 = -0.94229
  real(dp) :: p_pre_exp_3 = 0.28022
  real(dp) :: p_exp_3 = -0.40290
  real(dp) :: p_pre_exp_4 = 0.02817
  real(dp) :: p_exp_4 = -0.20162
  real(dp) :: E_scale
  character(len=STRING_LENGTH) :: label

end type IPModel_ZBL

logical, private :: parse_in_ip, parse_matched_label
type(IPModel_ZBL), private, pointer :: parse_ip

interface Initialise
  module procedure IPModel_ZBL_Initialise_str
end interface Initialise

interface Finalise
  module procedure IPModel_ZBL_Finalise
end interface Finalise

interface Print
  module procedure IPModel_ZBL_Print
end interface Print

interface Calc
  module procedure IPModel_ZBL_Calc
end interface Calc

contains

subroutine IPModel_ZBL_Initialise_str(this, args_str, param_str)
  type(IPModel_ZBL), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str

  type(Dictionary) :: params

  call Finalise(this)

  call initialise(params)
  this%label=''
  call param_register(params, 'label', '', this%label, help_string="No help yet.  This source file was $LastChangedBy$")
  if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_ZBL_Initialise_str args_str')) then
    call system_abort("IPModel_ZBL_Initialise_str failed to parse label from args_str="//trim(args_str))
  endif
  call param_register(params, 'E_scale', '1.0', this%E_scale, help_string="E_scale")
  call param_register(params, 'cutoff', '0.0', this%cutoff, help_string="cutoff")
  call param_register(params, 'cutoff_width', '0.0', this%cutoff_width, help_string="smooth cutoff width")
  call param_register(params, 'shift_cutoff', 'F', this%shift_cutoff, help_string="shift value at cutoff to equal 0")
  call param_register(params, 'cutoff_scale_cov_rad', 'F', this%cutoff_scale_cov_rad, help_string="multiply cutoff by sum of covalent radii (but cutoff_width is absolute)")
  call param_register(params, 'a_pre_exp', '0.46850', this%a_pre_exp, help_string="pre-exponential factor for screening parameter")
  call param_register(params, 'a_exp', '0.23', this%a_exp, help_string="exponent of charge of nuclei")
  call param_register(params, 'p_pre_exp_1', '0.18175', this%p_pre_exp_1, help_string="first pre-exponential factor of screening function")
  call param_register(params, 'p_exp_1', '-3.19980', this%p_exp_1, help_string="first exponent of screening function")
  call param_register(params, 'p_pre_exp_2', '0.50986', this%p_pre_exp_2, help_string="second pre-exponential factor of screening function")
  call param_register(params, 'p_exp_2', '-0.94229', this%p_exp_2, help_string="second exponent of screening function")
  call param_register(params, 'p_pre_exp_3', '0.28022', this%p_pre_exp_3, help_string="third pre-exponential factor of screening function")
  call param_register(params, 'p_exp_3', '-0.40290', this%p_exp_3, help_string="third exponent of screening function")
  call param_register(params, 'p_pre_exp_4', '0.02817', this%p_pre_exp_4, help_string="fourth pre-exponential factor of screening function")
  call param_register(params, 'p_exp_4', '-0.20162', this%p_exp_4, help_string="fourth exponent of screening function")

  if (param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_ZBL_Initialise_str args_str')) then
  end if

  ! by default cutoff used is cutoff set
  this%use_cutoff = this%cutoff
  if (this%cutoff_scale_cov_rad) then
       ! if cutoff is actually scale for covalent rad sum, set externally accessible
       ! cutoff to be max possible value for any element
       this%cutoff = 2.0*maxval(ElementCovRad)
  endif

  call finalise(params)

end subroutine IPModel_ZBL_Initialise_str

subroutine IPModel_ZBL_Finalise(this)
  type(IPModel_ZBL), intent(inout) :: this

  this%label = ''
end subroutine IPModel_ZBL_Finalise

subroutine IPModel_ZBL_Calc(this, at, e, local_e, f, virial, local_virial, args_str, mpi, error)
   type(IPModel_ZBL), intent(inout):: this
   type(Atoms), intent(inout)      :: at
   real(dp), intent(out), optional :: e, local_e(:)
   real(dp), intent(out), optional :: f(:,:), local_virial(:,:)   !% Forces, dimensioned as \texttt{f(3,at%N)}, local virials, dimensioned as \texttt{local_virial(9,at%N)} 
   real(dp), intent(out), optional :: virial(3,3)
   character(len=*), optional      :: args_str
   type(MPI_Context), intent(in), optional :: mpi
   integer, intent(out), optional :: error

   integer :: i, ji, j
   logical :: i_is_min_image
   real(dp) :: ke_e2 = 14.3996458521 !% Coulomb's constant x elemntary charge^{2} in eV*A
   real(dp) :: r, rs, dr(3)
   real(dp) :: a, c
   real(dp) :: t_1, t_2, t_3, t_4
   real(dp) :: rs_shifted, t_1_shifted, t_2_shifted, t_3_shifted, t_4_shifted, c_shifted
   real(dp) :: de, de_dr
   real(dp) :: f_cut = 1.0, df_cut = 0.0
   real(dp) :: use_cutoff

   type(Dictionary) :: params
   logical :: has_atom_mask_name
   logical, dimension(:), pointer :: atom_mask_pointer
   character(STRING_LENGTH) :: atom_mask_name

   INIT_ERROR(error)

   if (present(e)) e = 0.0
   if (present(local_e)) then
      call check_size('Local_E',local_e,(/at%N/),'IPModel_ZBL_Calc', error)
      local_e = 0.0
   endif
   if (present(f)) then
      call check_size('Force',f,(/3,at%Nbuffer/),'IPModel_ZBL_Calc', error)
      f = 0.0
   end if
   if (present(virial)) virial = 0.0
   if (present(local_virial)) then
      call check_size('Local_virial',local_virial,(/9,at%Nbuffer/),'IPModel_ZBL_Calc', error)
      local_virial = 0.0
      RAISE_ERROR('IPModel_ZBL_Calc: local_virial calculation requested but not supported yet',error)
   endif
   if (present(local_e)) then
      RAISE_ERROR('IPModel_ZBL_Calc: local_e calculation requested but not supported yet',error)
   end if

   atom_mask_pointer => null()
   if (present(args_str)) then
       call initialise(params)
       call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name,  has_value_target = has_atom_mask_name, help_string="name of atom_mask property")
       if (.not. param_read_line(params, args_str, ignore_unknown=.true., task='IPModel_ZBL_args_str')) then
          RAISE_ERROR("IPModel_ZBL_Calc failed to parse args_str='"//trim(args_str)//"'", error)
       endif
       call finalise(params)
   endif
   if (has_atom_mask_name) then
        if (.not. assign_pointer(at, trim(atom_mask_name) , atom_mask_pointer)) then
           RAISE_ERROR("IPModel_ZBL_Calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
        endif
   endif

   use_cutoff = this%use_cutoff
   do i = 1, at%N
      if (associated(atom_mask_pointer)) then
        if (.not. atom_mask_pointer(i)) cycle
      endif

      i_is_min_image = is_min_image(at,i)

      if (present(mpi)) then
         if (mpi%active) then
            if (mod(i-1, mpi%n_procs) /= mpi%my_proc) cycle
         endif
      endif

      do ji = 1, n_neighbours(at, i)
         j = neighbour(at, i, ji, distance=r, cosines=dr)

         ! if (i < j) cycle

         if (this%cutoff_scale_cov_rad) use_cutoff = this%use_cutoff * (ElementCovRad(at%Z(i))+ElementCovRad(at%Z(j)))

         if (use_cutoff == 0.0 .or. r <= use_cutoff) then
            c = ke_e2*real(at%Z(i))*real(at%Z(j))/r
            a = this%a_pre_exp/(real(at%Z(i))**this%a_exp + real(at%Z(j))**this%a_exp)
            rs = r/a
            t_1 = this%p_pre_exp_1*exp(this%p_exp_1*rs)
            t_2 = this%p_pre_exp_2*exp(this%p_exp_2*rs)
            t_3 = this%p_pre_exp_3*exp(this%p_exp_3*rs)
            t_4 = this%p_pre_exp_4*exp(this%p_exp_4*rs)
            if (use_cutoff > 0.0 .and. this%cutoff_width > 0.0) f_cut = poly_switch(r, use_cutoff, this%cutoff_width)
            de = c*(t_1+t_2+t_3+t_4)
            if (this%shift_cutoff) then
                c_shifted = ke_e2*real(at%Z(i))*real(at%Z(j))/use_cutoff
                rs_shifted = use_cutoff/a
                t_1_shifted = this%p_pre_exp_1*exp(this%p_exp_1*rs_shifted)
                t_2_shifted = this%p_pre_exp_2*exp(this%p_exp_2*rs_shifted)
                t_3_shifted = this%p_pre_exp_3*exp(this%p_exp_3*rs_shifted)
                t_4_shifted = this%p_pre_exp_4*exp(this%p_exp_4*rs_shifted)
                de = de - c_shifted*(t_1_shifted+t_2_shifted+t_3_shifted+t_4_shifted)
            endif
            if (present(e)) then
               e = e + 0.5*de*f_cut
            end if
            if (present(f) .or. present(virial)) then
               if (use_cutoff > 0.0 .and. this%cutoff_width > 0.0) df_cut = dpoly_switch(r, use_cutoff, this%cutoff_width)
               de_dr = -c/r*(t_1+t_2+t_3+t_4) + c/a*(this%p_exp_1*t_1+this%p_exp_2*t_2+this%p_exp_3*t_3+this%p_exp_4*t_4)
               de_dr = de_dr*f_cut + de*df_cut
               if (present(f)) then
                   f(:,i) = f(:,i) + 0.5*de_dr*dr
                   f(:,j) = f(:,j) - 0.5*de_dr*dr
               end if
               if (present(virial)) then
                  virial = virial - 0.5_dp*de_dr*(dr .outer. dr)*r
               endif
            end if
         end if
      end do
   end do

   if(present(e)) e = e*this%E_scale
   if(present(f)) f = f*this%E_scale
   if(present(virial)) virial = virial*this%E_scale
   if(present(local_e)) local_e = local_e*this%E_scale
   
   if (present(mpi)) then
      if (present(e)) e = sum(mpi, e)
      if (present(local_e)) call sum_in_place(mpi, local_e)
      if (present(virial)) call sum_in_place(mpi, virial)
      if (present(f)) call sum_in_place(mpi, f)
   endif

end subroutine IPModel_ZBL_Calc


subroutine IPModel_ZBL_Print(this, file)
  type(IPModel_ZBL), intent(in) :: this
  type(Inoutput), intent(inout),optional :: file

  call Print("IPModel_ZBL : ZBL Potential", file=file)
  call Print("IPModel_ZBL : cutoff = " // this%cutoff, file=file)
  call Print("IPModel_ZBL : cutoff_width = " // this%cutoff_width, file=file)
  call Print("IPModel_ZBL : a_pre_exp = " // this%a_pre_exp, file=file)
  call Print("IPModel_ZBL : a_exp = " // this%a_exp, file=file)
  call Print("IPModel_ZBL : p_pre_exp_1 = " // this%p_pre_exp_1, file=file)
  call Print("IPModel_ZBL : p_exp_1 = " // this%p_exp_1, file=file)
  call Print("IPModel_ZBL : p_pre_exp_2 = " // this%p_pre_exp_2, file=file)
  call Print("IPModel_ZBL : p_exp_2 = " // this%p_exp_2, file=file)
  call Print("IPModel_ZBL : p_pre_exp_3 = " // this%p_pre_exp_3, file=file)
  call Print("IPModel_ZBL : p_exp_3 = " // this%p_exp_3, file=file)
  call Print("IPModel_ZBL : p_pre_exp_4 = " // this%p_pre_exp_4, file=file)
  call Print("IPModel_ZBL : p_exp_4 = " // this%p_exp_4, file=file)

end subroutine IPModel_ZBL_Print

end module IPModel_ZBL_module
