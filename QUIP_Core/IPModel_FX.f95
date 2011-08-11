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
!X IPModel_FX 
!X
!% polarisable water model 
!% G. S. Fanourgakis and S. S. Xantheas, Journal of Chemical Physics 128, 074506 (2008)
!% This is a wrapper for a code downloaded from http://www.pnl.gov/science/ttm3f.asp
!% WARNING: it does not deal with periodic boundary conditions
!% UNITS: KCal/mol for energies and Angstroms for distance
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#include "error.inc"

module IPModel_FX_module

use libatoms_module

use mpi_context_module
use QUIP_Common_module

implicit none
private

include 'IPModel_interface.h'

public :: IPModel_FX
type IPModel_FX
   real(dp) :: cutoff = 2.0_dp
   logical :: return_one_body = .false.
   logical :: return_two_body = .false.
   real(dp) :: two_body_weight_roo = 0.0_dp
   real(dp) :: two_body_weight_delta = 0.0_dp
   logical  :: do_two_body_weight = .false.
   type(Spline) :: two_body_weight
end type IPModel_FX


interface Initialise
  module procedure IPModel_FX_Initialise_str
end interface Initialise

interface Finalise
  module procedure IPModel_FX_Finalise
end interface Finalise

interface Print
  module procedure IPModel_FX_Print
end interface Print

interface Calc
  module procedure IPModel_FX_Calc
end interface Calc

contains

subroutine IPModel_FX_Initialise_str(this, args_str, param_str, error)
  type(IPModel_FX), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str
  type(Dictionary)                :: params
  integer, optional, intent(out) :: error

  INIT_ERROR(error)

  call Finalise(this)

  call initialise(params)
  call param_register(params, 'return_two_body', 'F', this%return_two_body, help_string="if set, return the two_body energy and force as the main return data")
  call param_register(params, 'return_one_body', 'F', this%return_one_body, help_string="if set, return the one_body energy and force as the main return data")
  call param_register(params, 'two_body_weight_roo', '0.0', this%two_body_weight_roo, has_value_target=this%do_two_body_weight, help_string="if set, apply weight function to 2-body energy and force based on O-O distance. For a positive two_body_weight_delta, weight is 1 for rOO < two_body_weight_roo-two_body_weight_delta and weight is 0 for rOO > two_body_weight_roo+two_body_weight_delta")
  call param_register(params, 'two_body_weight_delta', '0.25', this%two_body_weight_delta, help_string="width of weighting function for  two_body energy and force based on O-O distance. For weighting to take effect, two_body_weight_roo needs to be explicitly set. For a positive two_body_weight_delta, weight is 1 for rOO < two_body_weight_roo-two_body_weight_delta and weight is 0 for rOO > two_body_weight_roo+two_body_weight_delta")
  
  if(.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_FX_Initialise args_str')) then
     RAISE_ERROR("IPModel_FX_Calc failed to parse args_str='"//trim(args_str)//"'", error)
  endif
  call finalise(params)

  if(this%do_two_body_weight) then
     call initialise(this%two_body_weight, (/this%two_body_weight_roo - this%two_body_weight_delta, this%two_body_weight_roo + this%two_body_weight_delta/), (/1.0_dp, 0.0_dp/), 0.0_dp, 0.0_dp)
  end if

end subroutine IPModel_FX_Initialise_str

subroutine IPModel_FX_Finalise(this)
  type(IPModel_FX), intent(inout) :: this
end subroutine IPModel_FX_Finalise


subroutine IPModel_FX_Calc(this, at, e, local_e, f, virial, local_virial, args_str, mpi, error)
  type(IPModel_FX), intent(inout):: this
  type(Atoms), intent(inout)      :: at
  real(dp), intent(out), optional :: e, local_e(:)
  real(dp), intent(out), optional :: f(:,:), local_virial(:,:)   !% Forces, dimensioned as \texttt{f(3,at%N)}, local virials, dimensioned as \texttt{local_virial(9,at%N)} 
  real(dp), intent(out), optional :: virial(3,3)
  character(len=*), optional      :: args_str
  type(MPI_Context), intent(in), optional :: mpi
  integer, intent(out), optional :: error

  real(dp), dimension(3,at%N) :: RR, dRR
  real(dp) :: energy
  integer, dimension(at%N) :: rindex
  integer :: i, j, k, kk, Oi, Oj

  type(Dictionary)                :: params
  logical :: has_atom_mask_name, do_one_body, do_two_body
  character(FIELD_LENGTH) :: atom_mask_name, one_body_name, two_body_name

  integer, dimension(:,:), allocatable :: water_monomer_index
  real(dp), allocatable :: one_body_energy(:), one_body_force(:,:), two_body_force(:,:)
  real(dp) :: watpos(3,3), wat2pos(3,6), watRR(3,3), watdRR(3,3), wat2RR(3,6), wat2dRR(3,6), wat2_force(3,6)
  real(dp) :: two_body_energy, diff_OiOj(3), watE, wat2E, two_body_energy_ij
  integer :: wat_rindex(3), wat2_rindex(6), watZ(3), wat2Z(6)
   real(dp):: weight, dweight(3), rOiOj

  INIT_ERROR(error)

  if(present(local_e)) then
     RAISE_ERROR('IPModel_FX_Calc: local_e calculation requested but not supported yet.', error)
  end if
  if(present(virial)) then
     RAISE_ERROR('IPModel_FX_Calc: virial calculation requested but not supported yet.', error)
  end if
  if (present(local_virial)) then
     RAISE_ERROR("IPModel_FX_Calc: local_virial calculation requested but not supported yet.", error)
  endif

  if(.not. present(e) .and. .not. present(f)) return ! nothing to do

  if (present(args_str)) then
     call initialise(params)
     call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, help_string="No help yet.  This source file was $LastChangedBy: nb326 $")
     call param_register(params, 'one_body', 'one_body', one_body_name, has_value_target=do_one_body, help_string="compute one-body terms of the cluster expansion and store it using this name")
     call param_register(params, 'two_body', 'two_body', two_body_name, has_value_target=do_two_body, help_string="compute two-body terms of the cluster expansion and store it using this name")

     if(.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_FX_Calc args_str')) then
        RAISE_ERROR("IPModel_FX_Calc failed to parse args_str='"//trim(args_str)//"'",error)
     endif
     call finalise(params)
     if(has_atom_mask_name) then
        RAISE_ERROR('IPModel_FX_Calc: atom_mask_name found, but not supported', error)
     endif
  endif

  
  call nttm3f_readXYZ(at%N/3, at%Z, at%pos, RR, rindex)
  call ttm3f(at%N/3,RR, dRR, energy)

  if(present(e)) e=energy * KCAL_MOL
  if(present(f)) then
     do i=1,at%N
        f(:,i)=-dRR(:,rindex(i)) * KCAL_MOL
     end do
  end if

  ! cluster expansion


  if(do_one_body .or. do_two_body .or. this%return_one_body .or. this%return_two_body) then
     
     allocate(water_monomer_index(3,at%N/3))
     call find_water_monomer(at,water_monomer_index, error)
     call print(at)
     PASS_ERROR(error)

     allocate(one_body_energy(at%N/3))
     allocate(one_body_force(3,at%N))
     ! compute monomer energies and forces
     do i=1,at%N/3
        Oi = water_monomer_index(1,i)
        watpos(:,1) = at%pos(:,Oi)
        watZ(1) = 8 ! Oxygen
        watpos(:,2) = at%pos(:,Oi)+diff_min_image(at, Oi, water_monomer_index(2,i))
        watpos(:,3) = at%pos(:,Oi)+diff_min_image(at, Oi, water_monomer_index(3,i))
        watZ(2:3) = 1 ! Hydrogens

        call nttm3f_readXYZ(1, watZ, watpos, watRR, wat_rindex)
        call ttm3f(1, watRR, watdRR, watE)
        one_body_energy(i) = watE*KCAL_MOL
        do k=1,3
           one_body_force(:,water_monomer_index(k,i)) = -watdRR(:,wat_rindex(k)) * KCAL_MOL
        end do
     end do

     ! if one body terms were asked for, store them
     if(this%return_one_body) then
        if(present(e)) then
           e = sum(one_body_energy)
        end if
        if(present(f)) then
           f = one_body_force
        end if
     end if

     if(do_one_body) then
        call set_value(at%params, trim(one_body_name)//"_energy", sum(one_body_energy))
        call add_property(at, trim(one_body_name)//"_force", one_body_force)
     end if

     ! compute dimer energies and forces
     if(do_two_body .or. this%return_two_body) then
        allocate(two_body_force(3,at%N))
        two_body_energy = 0.0_dp
        two_body_force = 0.0_dp
        do i=1,at%N/3
           do j=i+1,at%N/3
              Oi = water_monomer_index(1,i)
              Oj = water_monomer_index(1,j)
              wat2pos(:,1) = at%pos(:,Oi)
              wat2Z(1) = 8 ! first Oxygen
              wat2pos(:,2) = at%pos(:,Oi) + diff_min_image(at, Oi, water_monomer_index(2,i))
              wat2pos(:,3) = at%pos(:,Oi) + diff_min_image(at, Oi, water_monomer_index(3,i))
              wat2Z(2:3) = 1 ! first Hydrogens
              diff_OiOj = diff_min_image(at, Oi, Oj)
              wat2pos(:,4) = at%pos(:,Oi) + diff_OiOj
              wat2Z(4) = 8 ! second Oxygen
              wat2pos(:,5) = at%pos(:,Oi) + diff_OiOj+diff_min_image(at, Oj, water_monomer_index(2,j))
              wat2pos(:,6) = at%pos(:,Oi) + diff_OiOj+diff_min_image(at, Oj, water_monomer_index(3,j))
              wat2Z(5:6) = 1 ! first Hydrogens

              call nttm3f_readXYZ(2, wat2Z, wat2pos, wat2RR, wat2_rindex)
              call ttm3f(2, wat2RR, wat2dRR, wat2E)

              if(this%do_two_body_weight) then
                 rOiOj = norm(diff_OiOj)
                 weight = spline_value(this%two_body_weight, rOiOj)
                 dweight = spline_deriv(this%two_body_weight, rOiOj)
              else
                 weight = 1.0_dp
                 dweight = 0.0_dp
              end if

              two_body_energy_ij = (wat2E * KCAL_MOL - one_body_energy(i) - one_body_energy(j)) * weight
              two_body_energy = two_body_energy + two_body_energy_ij

              do k=1,6
                 wat2_force(:,k) = -wat2dRR(:,wat2_rindex(k)) * KCAL_MOL
              end do
              do k=1,3
                 kk = water_monomer_index(k, i)
                 two_body_force(:,kk) = two_body_force(:,kk) + (wat2_force(:,k)-one_body_force(:,kk))*weight
                 if(k==1) two_body_force(:,kk) = two_body_force(:,kk) - two_body_energy_ij*dweight*diff_OiOj/rOiOj

                 kk = water_monomer_index(k, j)
                 two_body_force(:,kk) = two_body_force(:,kk) + (wat2_force(:,3+k)-one_body_force(:,kk))*weight
                 if(k==1) two_body_force(:,kk) = two_body_force(:,kk) + two_body_energy_ij*dweight*diff_OiOj/rOiOj
              end do
           end do
        end do

        if(this%return_two_body) then
           if(present(e)) then
              e = two_body_energy
           end if
           if(present(f)) then
              f = two_body_force
           end if
        end if

        ! store two-body terms
        if(do_two_body) then
           call set_value(at%params, trim(two_body_name)//"_energy", two_body_energy)
           call add_property(at, trim(two_body_name)//"_force", two_body_force)
        end if
     end if
  endif

end subroutine IPModel_FX_Calc


subroutine IPModel_FX_Print(this, file)
  type(IPModel_FX), intent(in) :: this
  type(Inoutput), intent(inout),optional :: file
  call Print("IPModel_FX : ", file=file)
  call print("polarisable model for water", file=file)
  call print("by G. S. Fanourgakis and S. S. Xantheas", file=file)
  call print("Journal of Chemical Physics 128, 074506 (2008)", file=file)
  if(this%return_two_body) then
     call print("Returning 2-body term as a result of calc()", file=file)
  end if
  if(this%return_one_body) then
     call print("Returning 1-body term as a result of calc()", file=file)
  end if
  if(this%do_two_body_weight) then
     call print("Two-body term is weighted with parameters x0="//this%two_body_weight_roo//" delta="//this%two_body_weight_delta, file=file)
  end if
end subroutine IPModel_FX_Print


end module IPModel_FX_module
