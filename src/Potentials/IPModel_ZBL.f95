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

use mpi_context_module
use QUIP_Common_module

implicit none
private

include 'IPModel_interface.h'

public :: IPModel_ZBL
type IPModel_ZBL
  real(dp) :: cutoff = 0.0
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
  call param_register(params, 'cutoff', '0.0', this%cutoff, help_string="cutoff")
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
   real(dp) :: de, de_dr

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

   do i = 1, at%N
      i_is_min_image = is_min_image(at,i)

      if (present(mpi)) then
         if (mpi%active) then
            if (mod(i-1, mpi%n_procs) /= mpi%my_proc) cycle
         endif
      endif

      do ji = 1, n_neighbours(at, i)
         j = neighbour(at, i, ji, distance=r, cosines=dr)

         if (i < j) cycle

         if (this%cutoff == 0.0 .or. r <= this%cutoff) then
               c = ke_e2*real(at%Z(i))*real(at%Z(j))/r
               a = this%a_pre_exp/(real(at%Z(i))**this%a_exp + real(at%Z(j))**this%a_exp)
               rs = r/a
               t_1 = this%p_pre_exp_1*exp(this%p_exp_1*rs)
               t_2 = this%p_pre_exp_2*exp(this%p_exp_2*rs)
               t_3 = this%p_pre_exp_3*exp(this%p_exp_3*rs)
               t_4 = this%p_pre_exp_4*exp(this%p_exp_4*rs)
            if (present(local_e)) then

            end if
            if (present(e)) then
               de = c*(t_1+t_2+t_3+t_4)
               if (i == j) then
                  e = e + 0.5*de
               else
                  e = e + de
               end if
            end if
            if (present(f)) then
               de_dr = -c/r*(t_1+t_2+t_3+t_4) + c/a*(this%p_exp_1*t_1+this%p_exp_2*t_2+this%p_exp_3*t_3+this%p_exp_4*t_4)
               f(:,i) = f(:,i) + de_dr*dr
               if ( i/=j ) f(:,j) = f(:,j) - de_dr*dr
            end if
         end if
      end do
   end do

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
