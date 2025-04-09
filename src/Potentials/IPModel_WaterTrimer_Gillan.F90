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
!X IPModel_WaterTrimer_Gillan
!X
!% 3-body interaction of three water molecules, parametrised by Mike Gillan in 2013
!% 
!% Based on polynomial basis functions of Paesani 2012
!%
!% Energy only
!% 
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#include "error.inc"

module IPModel_WaterTrimer_Gillan_module

use error_module
use system_module, only : dp, inoutput, print, operator(//)
use units_module
use dictionary_module
use paramreader_module
use linearalgebra_module
use spline_module
use quaternions_module
use atoms_types_module
use atoms_module
use topology_module

use mpi_context_module
use QUIP_Common_module

implicit none
private

include 'IPModel_interface.h'

public :: IPModel_WaterTrimer_Gillan
type IPModel_WaterTrimer_Gillan
  real(dp) :: cutoff = 8.0_dp
  character(len=STRING_LENGTH) :: coeff_file
  logical :: OHH_ordercheck = .true.
  type(Spline) :: fcut
  real (dp) :: fcut_delta = 0.5
end type IPModel_WaterTrimer_Gillan

logical, private :: parse_in_ip, parse_matched_label
type(IPModel_WaterTrimer_Gillan), private, pointer :: parse_ip

interface Initialise
  module procedure IPModel_WaterTrimer_Gillan_Initialise_str
end interface Initialise

interface Finalise
  module procedure IPModel_WaterTrimer_Gillan_Finalise
end interface Finalise

interface Print
  module procedure IPModel_WaterTrimer_Gillan_Print
end interface Print

interface Calc
  module procedure IPModel_WaterTrimer_Gillan_Calc
end interface Calc

! module-gobal variables

integer,parameter :: nbas_tot_max=150

integer ibas_code, nbas_tot, nbas_totp
integer natom, nbas_2deg, nbas_3deg
integer iswitch(nbas_tot_max)
real(dp):: drep(nbas_tot_max)
real(dp):: basis(nbas_tot_max)



contains

subroutine IPModel_WaterTrimer_Gillan_Initialise_str(this, args_str, param_str, error)
  type(IPModel_WaterTrimer_Gillan), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str
  type(Dictionary)               :: params
  integer, optional, intent(out) :: error
  integer :: i, nb1
  real(dp) :: rr(3,9)

  INIT_ERROR(error)

  call Finalise(this)

  call initialise(params)
  call param_register(params, 'coeff_file', PARAM_MANDATORY, this%coeff_file, help_string="name of file which contains the polynomial coefficients")

  if(.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_WaterTrimer_Gillan_Initialise args_str')) then
     RAISE_ERROR("IPModel_WaterTrimer_Gillan_Init failed to parse args_str='"//trim(args_str)//"'", error)
  endif
  call finalise(params)

      open(unit=21,file=this%coeff_file)
      read(21,*) ibas_code
!      write(*,121) ibas_code
!  121 format(/'basis code: ',i3)
      read(21,*) nbas_tot
!      write(*,122) nbas_tot
!  122 format(/'total no. of basis functions: ',i5)
      if(nbas_tot .gt. nbas_tot_max) then
!        write(*,901)
!  901   format(/'error: too many basis functions')
        RAISE_ERROR('error: too many basis functions', error)
      endif

!c ... initial call of basis calculator
      call calc_basis(0,ibas_code,nbas_2deg,nbas_3deg,nbas_tot_max,rr,basis)
!      write(*,126) nbas_2deg, nbas_3deg
!  126 format(/'after initial call to basis calculator:'/
!     * 'no. of deg-2 basis functions: ',i5/
!     * 'no. of deg-3 basis functions: ',i5)
      nbas_totp = nbas_2deg + nbas_3deg
!      write(*,127) nbas_totp
!  127 format(/'hence total no. of basis functions: ',i5)
      if(nbas_totp .ne. nbas_tot) then
!        write(*,905)
!  905   format(//'error: conflicting no. of basis functions')
         RAISE_ERROR('error: conflicting no. of basis functions', error)
      endif
!c ... read basis coeffs
!      write(*,131)
!  131 format(/'switch flags, basis coeffs:')
      do i = 1, nbas_tot
        read(21,*) nb1, iswitch(i), drep(i)
!        write(*,132) nb, iswitch(nb), drep(nb)
!  132   format(i5,3x,i3,3x,e15.6)
      enddo
      close(unit=21)

      call initialise(this%fcut, (/this%cutoff - this%fcut_delta, this%cutoff/), (/1.0_dp, 0.0_dp/), 0.0_dp, 0.0_dp)
      

    end subroutine IPModel_WaterTrimer_Gillan_Initialise_str

subroutine IPModel_WaterTrimer_Gillan_Finalise(this)
  type(IPModel_WaterTrimer_Gillan), intent(inout) :: this

  ! Add finalisation code here
  call finalise(this%fcut)


end subroutine IPModel_WaterTrimer_Gillan_Finalise


subroutine IPModel_WaterTrimer_Gillan_Calc(this, at, e, local_e, f, virial, local_virial, args_str, mpi, error)
  type(IPModel_WaterTrimer_Gillan), intent(inout):: this
  type(Atoms), intent(inout)      :: at
  real(dp), intent(out), optional :: e, local_e(:)
  real(dp), intent(out), optional :: f(:,:), local_virial(:,:)   !% Forces, dimensioned as \texttt{f(3,at%N)}, local virials, dimensioned as \texttt{local_virial(9,at%N)} 
  real(dp), intent(out), optional :: virial(3,3)
  character(len=*), optional      :: args_str
  type(MPI_Context), intent(in), optional :: mpi
  integer, intent(out), optional :: error

  integer :: i, w1, w2, w3, w1O, w2O, w3O
  integer, dimension(3,at%N/3) :: water_index

  type(Dictionary)                :: params
  logical :: has_atom_mask_name
  character(STRING_LENGTH) :: atom_mask_name
  real(dp) :: r_scale, E_scale
  logical :: do_rescale_r, do_rescale_E

  real(dp) rr(3,9)
  real(dp) e_trimer

  INIT_ERROR(error)

  if(present(local_e)) then
     RAISE_ERROR('IPModel_WaterTrimer_Gillan_Calc: local_e not implemented', error)
  end if
  if(present(f)) then
     RAISE_ERROR('IPModel_WaterTrimer_Gillan_Calc: forces not implemented', error)
  end if
  if(present(virial)) then
     RAISE_ERROR('IPModel_WaterTrimer_Gillan_Calc: virial not implemented', error)
  end if
  if (present(local_virial)) then
     RAISE_ERROR("IPModel_WaterTrimer_Gillan_Calc: local_virial calculation requested but not supported yet.", error)
  endif

  if(.not. present(e)) return ! nothing to do

  if (present(args_str)) then
     call initialise(params)
     call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, help_string="No help yet.  This source file was $LastChangedBy: nb326 $")
     call param_register(params, 'r_scale', '1.0',r_scale, has_value_target=do_rescale_r, help_string="Recaling factor for distances. Default 1.0.")
     call param_register(params, 'E_scale', '1.0',E_scale, has_value_target=do_rescale_E, help_string="Recaling factor for energy. Default 1.0.")

     if(.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_WaterTrimer_Gillan_Calc args_str')) then
        RAISE_ERROR("IPModel_WaterTrimer_Gillan_Calc failed to parse args_str='"//trim(args_str)//"'",error)
     endif
     call finalise(params)
     if(has_atom_mask_name) then
        RAISE_ERROR('IPModel_WaterTrimer_Gillan_Calc: atom_mask_name found, but not supported', error)
     endif
     if (do_rescale_r .or. do_rescale_E) then
        RAISE_ERROR("IPModel_WaterTrimer_Gillan_Calc: rescaling of potential with r_scale and E_scale not yet implemented!", error)
     end if
  endif

  ! loop through atoms, find oxygens and their closest hydrogens

  call find_water_monomer(at,water_index,this%OHH_ordercheck,error=error)

  e = 0.d0

  ! triple loop through monomers
  do w1 = 1, at%N/3
     do w2 = w1+1, at%N/3
        do w3 = w2+1, at%N/3

           w1O = water_index(1,w1)
           rr(:,1) = at%pos(:,w1O) ! O
           rr(:,2) = rr(:,1) + diff_min_image(at, w1O, water_index(2,w1)) ! H
           rr(:,3) = rr(:,1) + diff_min_image(at, w1O, water_index(3,w1)) ! H

           w2O = water_index(1,w2)
           rr(:,4) = rr(:,1) + diff_min_image(at, w1O, w2O) ! O
           rr(:,5) = rr(:,1) + diff_min_image(at, w1O, water_index(2,w2)) ! H
           rr(:,6) = rr(:,1) + diff_min_image(at, w1O, water_index(3,w2)) ! H

           w3O = water_index(1,w3)
           rr(:,7) = rr(:,1) + diff_min_image(at, w1O, w3O) ! O
           rr(:,8) = rr(:,1) + diff_min_image(at, w1O, water_index(2,w3)) ! H
           rr(:,9) = rr(:,1) + diff_min_image(at, w1O, water_index(3,w3)) ! H

!c ... calculate values of basis functions for current config
           call calc_basis(1,ibas_code,nbas_2deg,nbas_3deg,nbas_tot_max, rr,basis)
           e_trimer = 0.0_dp
           do i = 1, nbas_tot
              if(iswitch(i) .eq. 1) then
                 e_trimer = e_trimer + drep(i)*basis(i)
              endif
           enddo
           e = e+e_trimer* spline_value(this%fcut, distance_min_image(at, w1O, w2O)) * &
                spline_value(this%fcut, distance_min_image(at, w1O,w3O))* &
                spline_value(this%fcut, distance_min_image(at, w2O, w3O))
           

           !        write(*,231) nc, e3bas
           !  231   format(/'config. no. ',i5,' 3-body correction: ',e15.6)
           !        write(13,232) nc, e3bas
           !  232   format(i5,3x,e15.6)
        enddo
     end do
  end do
  !c ... end loop over configs
  

  e = e*HARTREE

  
end subroutine IPModel_WaterTrimer_Gillan_Calc


subroutine IPModel_WaterTrimer_Gillan_Print(this, file)
  type(IPModel_WaterTrimer_Gillan), intent(in) :: this
  type(Inoutput), intent(inout),optional :: file

  call Print("IPModel_WaterTrimer_Gillan : WaterTrimer_Gillan Potential", file=file)
  call Print("IPModel_WaterTrimer_Gillan : cutoff = " // this%cutoff, file=file)
  call Print("IPModel_WaterTrimer_Gillan : cutoff_delta =  " // this%fcut_delta, file=file)
  call Print("IPModel_WaterTrimer_Gillan : WARNING: no multiple periodic images implemented, so minimum-image convention applies to this potential!!!  ", file=file)

end subroutine IPModel_WaterTrimer_Gillan_Print

!c =====================================================================
!c   Subroutine calc_basis: directs calculation to appropriate
!c   version calc_basis_n() according to value of flag ibas_code.
!c   The only versions implemented so far are ibas_code = 0, 1
!c   and any other values will cause the code to stop.
!c ---------------------------------------------------------------------
      subroutine calc_basis(iret,ibas_code,nbas_2deg,nbas_3deg,nbas_tot_max,rr,basis)
      implicit none
      integer iret, ibas_code, nbas_2deg, nbas_3deg, nbas_tot_max
      real*8 rr(3,9)
      real*8 basis(nbas_tot_max)
!c ---------------------------------------------------------------------
      if(ibas_code .eq. 0) then
        call calc_basis_0(iret,nbas_2deg,nbas_3deg,nbas_tot_max,rr,basis)
      elseif(ibas_code .eq. 1) then
        call calc_basis_1(iret,nbas_2deg,nbas_3deg,nbas_tot_max,rr,basis)
      else
!        write(*,901)
!  901   format(//'error: sbrt calc_basis called with flag nbas_2deg ',
!     *   'out of range')
        stop
      endif
!c ---------------------------------------------------------------------
      return
    end subroutine calc_basis
!c =====================================================================

!c =================================================================
!c   Subroutine calc_basis_0: computes deg-2 and deg-3 3-body
!c   basis functions for model where only O-O distances are used
!c   and at most one distance factor is allowed for each
!c   O-O distance.
!c ------------------------------------------------------------------
      subroutine calc_basis_0(iret,nbas_2deg,nbas_3deg,nbas_tot_max, rr,basis)
      implicit none
      integer iret, nbas_2deg, nbas_3deg, nbas_tot_max
      real*8 rr(3,9)
      real*8 basis(nbas_tot_max)
      real*8 kOO, dOO
      real*8 d14, d17, d47, e14, e17, e47
!c ---------------------------------------------------------------------
      nbas_2deg = 1
      nbas_3deg = 1
      if(iret .eq. 0) return
      kOO = 1.3d0
      dOO = 3.d0
!c ...
      d14 = sqrt((rr(1,4) - rr(1,1))**2 + (rr(2,4) - rr(2,1))**2 +           (rr(3,4) - rr(3,1))**2)
      d17 = sqrt((rr(1,7) - rr(1,1))**2 + (rr(2,7) - rr(2,1))**2 +           (rr(3,7) - rr(3,1))**2)
      d47 = sqrt((rr(1,4) - rr(1,7))**2 + (rr(2,4) - rr(2,7))**2 +           (rr(3,4) - rr(3,7))**2)
      e14 = exp(-kOO*(d14 - dOO))
      e17 = exp(-kOO*(d17 - dOO))
      e47 = exp(-kOO*(d47 - dOO))
      basis(1) = e14*e17 + e14*e47 + e17*e47
      basis(2) = e14*e17*e47
!c ---------------------------------------------------------------------
      return
    end subroutine calc_basis_0
!c =====================================================================

!c =================================================================
!c   Subroutine calc_basis_1: computes deg-2 and deg-3 3-body
!c   basis functions for model where O-O, O-H and H-H  distances are
!c   used and at most one distance factor is allowed for each
!c   monomer pair.
!c ------------------------------------------------------------------
      subroutine calc_basis_1(iret,nbas_2deg,nbas_3deg,nbas_tot_max, rr,basis)
      implicit none
      integer iret, nbas_2deg, nbas_3deg, nbas_tot_max
      real*8 rr(3,9)
      real*8 basis(nbas_tot_max)
!c ...
      real*8 kOO, dOO, kOH, dOH, kHH, dHH
      real*8 d14, d15, d16, d17, d18, d19,       d24, d25, d26, d27, d28, d29,       d34, d35, d36, d37, d38, d39,       d47, d48, d49, d57, d58, d59,       d67, d68, d69
      real*8 e14, e15, e16, e17, e18, e19,       e24, e25, e26, e27, e28, e29,       e34, e35, e36, e37, e38, e39,       e47, e48, e49, e57, e58, e59,       e67, e68, e69
!c ---------------------------------------------------------------------
      nbas_2deg = 13
      nbas_3deg = 30
      if(iret .eq. 0) return
!c ---------------------------------------------------------------------
      kOO = 0.6d0
      dOO = 3.d0
      kOH = 1.4d0
      dOH = 3.d0
      kHH = 1.0d0
      dHH = 3.d0
!c ---------------------------------------------------------------------
      d14 = sqrt((rr(1,4) - rr(1,1))**2 +  (rr(2,4) - rr(2,1))**2 +  (rr(3,4) - rr(3,1))**2) 
      d15 = sqrt((rr(1,5) - rr(1,1))**2 +  (rr(2,5) - rr(2,1))**2 +  (rr(3,5) - rr(3,1))**2) 
      d16 = sqrt((rr(1,6) - rr(1,1))**2 +  (rr(2,6) - rr(2,1))**2 +  (rr(3,6) - rr(3,1))**2) 
      d17 = sqrt((rr(1,7) - rr(1,1))**2 +  (rr(2,7) - rr(2,1))**2 +  (rr(3,7) - rr(3,1))**2) 
      d18 = sqrt((rr(1,8) - rr(1,1))**2 +  (rr(2,8) - rr(2,1))**2 +  (rr(3,8) - rr(3,1))**2) 
      d19 = sqrt((rr(1,9) - rr(1,1))**2 +  (rr(2,9) - rr(2,1))**2 +  (rr(3,9) - rr(3,1))**2) 
      d24 = sqrt((rr(1,4) - rr(1,2))**2 +  (rr(2,4) - rr(2,2))**2 +  (rr(3,4) - rr(3,2))**2) 
      d25 = sqrt((rr(1,5) - rr(1,2))**2 +  (rr(2,5) - rr(2,2))**2 +  (rr(3,5) - rr(3,2))**2) 
      d26 = sqrt((rr(1,6) - rr(1,2))**2 +  (rr(2,6) - rr(2,2))**2 +  (rr(3,6) - rr(3,2))**2) 
      d27 = sqrt((rr(1,7) - rr(1,2))**2 +  (rr(2,7) - rr(2,2))**2 +  (rr(3,7) - rr(3,2))**2) 
      d28 = sqrt((rr(1,8) - rr(1,2))**2 +  (rr(2,8) - rr(2,2))**2 +  (rr(3,8) - rr(3,2))**2) 
      d29 = sqrt((rr(1,9) - rr(1,2))**2 +  (rr(2,9) - rr(2,2))**2 +  (rr(3,9) - rr(3,2))**2) 
      d34 = sqrt((rr(1,4) - rr(1,3))**2 +  (rr(2,4) - rr(2,3))**2 +  (rr(3,4) - rr(3,3))**2) 
      d35 = sqrt((rr(1,5) - rr(1,3))**2 +  (rr(2,5) - rr(2,3))**2 +  (rr(3,5) - rr(3,3))**2) 
      d36 = sqrt((rr(1,6) - rr(1,3))**2 +  (rr(2,6) - rr(2,3))**2 +  (rr(3,6) - rr(3,3))**2) 
      d37 = sqrt((rr(1,7) - rr(1,3))**2 +  (rr(2,7) - rr(2,3))**2 +  (rr(3,7) - rr(3,3))**2) 
      d38 = sqrt((rr(1,8) - rr(1,3))**2 +  (rr(2,8) - rr(2,3))**2 +  (rr(3,8) - rr(3,3))**2) 
      d39 = sqrt((rr(1,9) - rr(1,3))**2 +  (rr(2,9) - rr(2,3))**2 +  (rr(3,9) - rr(3,3))**2) 
      d47 = sqrt((rr(1,7) - rr(1,4))**2 +  (rr(2,7) - rr(2,4))**2 +  (rr(3,7) - rr(3,4))**2) 
      d48 = sqrt((rr(1,8) - rr(1,4))**2 +  (rr(2,8) - rr(2,4))**2 +  (rr(3,8) - rr(3,4))**2) 
      d49 = sqrt((rr(1,9) - rr(1,4))**2 +  (rr(2,9) - rr(2,4))**2 +  (rr(3,9) - rr(3,4))**2) 
      d57 = sqrt((rr(1,7) - rr(1,5))**2 +  (rr(2,7) - rr(2,5))**2 +  (rr(3,7) - rr(3,5))**2) 
      d58 = sqrt((rr(1,8) - rr(1,5))**2 +  (rr(2,8) - rr(2,5))**2 +  (rr(3,8) - rr(3,5))**2) 
      d59 = sqrt((rr(1,9) - rr(1,5))**2 +  (rr(2,9) - rr(2,5))**2 +  (rr(3,9) - rr(3,5))**2) 
      d67 = sqrt((rr(1,7) - rr(1,6))**2 +  (rr(2,7) - rr(2,6))**2 +  (rr(3,7) - rr(3,6))**2) 
      d68 = sqrt((rr(1,8) - rr(1,6))**2 +  (rr(2,8) - rr(2,6))**2 +  (rr(3,8) - rr(3,6))**2) 
      d69 = sqrt((rr(1,9) - rr(1,6))**2 +  (rr(2,9) - rr(2,6))**2 +  (rr(3,9) - rr(3,6))**2) 

      e14 = exp(-kOO*(d14 - dOO))
      e15 = exp(-kOH*(d15 - dOH))
      e16 = exp(-kOH*(d16 - dOH))
      e17 = exp(-kOO*(d17 - dOO))
      e18 = exp(-kOH*(d18 - dOH))
      e19 = exp(-kOH*(d19 - dOH))
      e24 = exp(-kOH*(d24 - dOH))
      e25 = exp(-kHH*(d25 - dHH))
      e26 = exp(-kHH*(d26 - dHH))
      e27 = exp(-kOH*(d27 - dOH))
      e28 = exp(-kHH*(d28 - dHH))
      e29 = exp(-kHH*(d29 - dHH))
      e34 = exp(-kOH*(d34 - dOH))
      e35 = exp(-kHH*(d35 - dHH))
      e36 = exp(-kHH*(d36 - dHH))
      e37 = exp(-kOH*(d37 - dOH))
      e38 = exp(-kHH*(d38 - dHH))
      e39 = exp(-kHH*(d39 - dHH))
      e47 = exp(-kOO*(d47 - dOO))
      e48 = exp(-kOH*(d48 - dOH))
      e49 = exp(-kOH*(d49 - dOH))
      e57 = exp(-kOH*(d57 - dOH))
      e58 = exp(-kHH*(d58 - dHH))
      e59 = exp(-kHH*(d59 - dHH))
      e67 = exp(-kOH*(d67 - dOH))
      e68 = exp(-kHH*(d68 - dHH))
      e69 = exp(-kHH*(d69 - dHH))
!c
      basis(1) = e14*e17 + e14*e47 + e17*e47
!c
      basis(2) = e14*e18 + e14*e19 + e24*e47 + e34*e47 + e17*e57 +  e17*e67 + e14*e48 + e14*e49 + e27*e47 + e37*e47 +  e15*e17 + e16*e17
!c
      basis(3) = e14*e27 + e14*e37 + e15*e47 + e16*e47 + e17*e48 +  e17*e49 + e14*e57 + e14*e67 + e18*e47 + e19*e47 +  e17*e24 + e17*e34
!c
      basis(4) = e14*e28 + e14*e29 + e14*e38 + e14*e39 + e25*e47 +  e35*e47 + e26*e47 + e36*e47 + e17*e58 + e17*e68 +  e17*e59 + e17*e69 + e14*e58 + e14*e59 + e14*e68 +  e14*e69 + e28*e47 + e38*e47 + e29*e47 + e39*e47 +  e17*e25 + e17*e26 + e17*e35 + e17*e36
!c
      basis(5) = e15*e18 + e15*e19 + e16*e18 + e16*e19 + e24*e48 +  e34*e48 + e24*e49 + e34*e49 + e27*e57 + e27*e67 +  e37*e57 + e37*e67
!c
      basis(6) = e15*e27 + e16*e27 + e15*e37 + e16*e37 + e15*e48 +  e15*e49 + e16*e48 + e16*e49 + e27*e48 + e37*e48 +  e27*e49 + e37*e49 + e24*e57 + e34*e57 + e24*e67 +  e34*e67 + e18*e57 + e18*e67 + e19*e57 + e19*e67 +  e18*e24 + e19*e24 + e18*e34 + e19*e34
!c
      basis(7) = e15*e28 + e15*e29 + e16*e28 + e16*e29 + e15*e38 +  e15*e39 + e16*e38 + e16*e39 + e25*e48 + e35*e48 +  e25*e49 + e35*e49 + e26*e48 + e36*e48 + e26*e49 +  e36*e49 + e27*e58 + e27*e68 + e37*e58 + e37*e68 +  e27*e59 + e27*e69 + e37*e59 + e37*e69 + e24*e58 +  e24*e59 + e34*e58 + e34*e59 + e24*e68 + e24*e69 +  e34*e68 + e34*e69 + e28*e57 + e38*e57 + e28*e67 +  e38*e67 + e29*e57 + e39*e57 + e29*e67 + e39*e67 +  e18*e25 + e18*e26 + e19*e25 + e19*e26 + e18*e35 +  e18*e36 + e19*e35 + e19*e36
!c
      basis(8) = e15*e57 + e16*e67 + e18*e48 + e19*e49 + e24*e27 +  e34*e37
!c
      basis(9) = e15*e58 + e15*e59 + e16*e68 + e16*e69 + e28*e48 +  e38*e48 + e29*e49 + e39*e49 + e25*e27 + e26*e27 +  e35*e37 + e36*e37 + e24*e28 + e24*e29 + e34*e38 +  e34*e39 + e25*e57 + e35*e57 + e26*e67 + e36*e67 +  e18*e58 + e18*e68 + e19*e59 + e19*e69
!c
      basis(10) = e15*e67 + e16*e57 + e19*e48 + e18*e49 + e27*e34 +  e24*e37
!c
      basis(11) = e15*e68 + e15*e69 + e16*e58 + e16*e59 + e29*e48 +  e39*e48 + e28*e49 + e38*e49 + e27*e35 + e27*e36 +  e25*e37 + e26*e37 + e24*e38 + e24*e39 + e28*e34 +  e29*e34 + e26*e57 + e36*e57 + e25*e67 + e35*e67 +  e18*e59 + e18*e69 + e19*e58 + e19*e68
!c
      basis(12) = e25*e28 + e25*e29 + e26*e28 + e26*e29 + e35*e38 +  e35*e39 + e36*e38 + e36*e39 + e25*e58 + e35*e58 +  e25*e59 + e35*e59 + e26*e68 + e36*e68 + e26*e69 +  e36*e69 + e28*e58 + e28*e68 + e38*e58 + e38*e68 +  e29*e59 + e29*e69 + e39*e59 + e39*e69
!c
      basis(13) = e25*e38 + e25*e39 + e26*e38 + e26*e39 + e28*e35 +  e29*e35 + e28*e36 + e29*e36 + e26*e58 + e36*e58 +  e26*e59 + e36*e59 + e25*e68 + e35*e68 + e25*e69 +  e35*e69 + e28*e59 + e28*e69 + e38*e59 + e38*e69 +  e29*e58 + e29*e68 + e39*e58 + e39*e68
!c
      basis(14) = e14*e17*e47
!c
      basis(15) = e14*e17*e48 + e14*e17*e49 + e14*e27*e47 +  e14*e37*e47 + e15*e17*e47 + e16*e17*e47 +  e14*e18*e47 + e14*e19*e47 + e17*e24*e47 +  e17*e34*e47 + e14*e17*e57 + e14*e17*e67
!c
      basis(16) = e14*e17*e58 + e14*e17*e59 + e14*e17*e68 +  e14*e17*e69 + e14*e28*e47 + e14*e38*e47 +  e14*e29*e47 + e14*e39*e47 + e17*e25*e47 +  e17*e26*e47 + e17*e35*e47 + e17*e36*e47
!c
      basis(17) = e14*e18*e48 + e14*e19*e49 + e24*e27*e47 +  e34*e37*e47 + e15*e17*e57 + e16*e17*e67
!c
      basis(18) = e14*e18*e49 + e14*e19*e48 + e24*e37*e47 +  e27*e34*e47 + e16*e17*e57 + e15*e17*e67
!c
      basis(19) = e14*e18*e57 + e14*e19*e57 + e14*e18*e67 +  e14*e19*e67 + e18*e24*e47 + e18*e34*e47 +  e19*e24*e47 + e19*e34*e47 + e17*e24*e57 +  e17*e24*e67 + e17*e34*e57 + e17*e34*e67 +  e14*e27*e48 + e14*e27*e49 + e14*e37*e48 +  e14*e37*e49 + e15*e27*e47 + e15*e37*e47 +  e16*e27*e47 + e16*e37*e47 + e15*e17*e48 +  e16*e17*e48 + e15*e17*e49 + e16*e17*e49
!c
      basis(20) = e14*e18*e58 + e14*e19*e59 + e14*e18*e68 +  e14*e19*e69 + e24*e28*e47 + e34*e38*e47 +  e24*e29*e47 + e34*e39*e47 + e17*e25*e57 +  e17*e26*e67 + e17*e35*e57 + e17*e36*e67 +  e14*e28*e48 + e14*e29*e49 + e14*e38*e48 +  e14*e39*e49 + e25*e27*e47 + e35*e37*e47 +  e26*e27*e47 + e36*e37*e47 + e15*e17*e58 +  e16*e17*e68 + e15*e17*e59 + e16*e17*e69
!c
      basis(21) = e14*e18*e59 + e14*e19*e58 + e14*e18*e69 +  e14*e19*e68 + e24*e38*e47 + e28*e34*e47 +  e24*e39*e47 + e29*e34*e47 + e17*e26*e57 +  e17*e25*e67 + e17*e36*e57 + e17*e35*e67 +  e14*e29*e48 + e14*e28*e49 + e14*e39*e48 +  e14*e38*e49 + e27*e35*e47 + e25*e37*e47 +  e27*e36*e47 + e26*e37*e47 + e15*e17*e68 +  e16*e17*e58 + e15*e17*e69 + e16*e17*e59
!c
      basis(22) = e14*e27*e57 + e14*e27*e67 + e14*e37*e57 +  e14*e37*e67 + e15*e18*e47 + e15*e19*e47 +  e16*e18*e47 + e16*e19*e47 + e17*e24*e48 +  e17*e34*e48 + e17*e24*e49 + e17*e34*e49
!c
      basis(23) = e14*e27*e58 + e14*e27*e59 + e14*e27*e68 +  e14*e27*e69 + e14*e37*e58 + e14*e37*e59 +  e14*e37*e68 + e14*e37*e69 + e15*e28*e47 +  e15*e38*e47 + e15*e29*e47 + e15*e39*e47 +  e16*e28*e47 + e16*e38*e47 + e16*e29*e47 +  e16*e39*e47 + e17*e25*e48 + e17*e26*e48 +  e17*e35*e48 + e17*e36*e48 + e17*e25*e49 +  e17*e26*e49 + e17*e35*e49 + e17*e36*e49 +  e14*e28*e57 + e14*e29*e57 + e14*e38*e57 +  e14*e39*e57 + e14*e28*e67 + e14*e29*e67 +  e14*e38*e67 + e14*e39*e67 + e18*e25*e47 +  e18*e35*e47 + e18*e26*e47 + e18*e36*e47 +  e19*e25*e47 + e19*e35*e47 + e19*e26*e47 +  e19*e36*e47 + e17*e24*e58 + e17*e24*e68 +  e17*e24*e59 + e17*e24*e69 + e17*e34*e58 +  e17*e34*e68 + e17*e34*e59 + e17*e34*e69
!c
      basis(24) = e14*e28*e58 + e14*e29*e59 + e14*e28*e68 +  e14*e29*e69 + e14*e38*e58 + e14*e39*e59 +  e14*e38*e68 + e14*e39*e69 + e25*e28*e47 +  e35*e38*e47 + e25*e29*e47 + e35*e39*e47 +  e26*e28*e47 + e36*e38*e47 + e26*e29*e47 +  e36*e39*e47 + e17*e25*e58 + e17*e26*e68 +  e17*e35*e58 + e17*e36*e68 + e17*e25*e59 +  e17*e26*e69 + e17*e35*e59 + e17*e36*e69
!c
      basis(25) = e14*e28*e59 + e14*e29*e58 + e14*e28*e69 + &
           e14*e29*e68 + e14*e38*e59 + e14*e39*e58 + &
           e14*e38*e69 + e14*e39*e68 + e25*e38*e47 + &
           e28*e35*e47 + e25*e39*e47 + e29*e35*e47 + &
           e26*e38*e47 + e28*e36*e47 + e26*e39*e47 + &
           e29*e36*e47 + e17*e26*e58 + e17*e25*e68 + &
           e17*e36*e58 + e17*e35*e68 + e17*e26*e59 + &
           e17*e25*e69 + e17*e36*e59 + e17*e35*e69
!c
      basis(26) = e15*e18*e48 + e15*e19*e49 + e16*e18*e48 + &
           e16*e19*e49 + e24*e27*e48 + e34*e37*e48 + &
           e24*e27*e49 + e34*e37*e49 + e15*e27*e57 + &
           e16*e27*e67 + e15*e37*e57 + e16*e37*e67 + &
           e18*e24*e48 + e19*e24*e49 + e18*e34*e48 + &
           e19*e34*e49 + e24*e27*e57 + e34*e37*e57 + &
           e24*e27*e67 + e34*e37*e67 + e15*e18*e57 + &
           e16*e18*e67 + e15*e19*e57 + e16*e19*e67
!c
      basis(27) = e15*e18*e49 + e15*e19*e48 + e16*e18*e49 + &
           e16*e19*e48 + e24*e37*e48 + e27*e34*e48 + &
           e24*e37*e49 + e27*e34*e49 + e16*e27*e57 + &
           e15*e27*e67 + e16*e37*e57 + e15*e37*e67 + &
           e19*e24*e48 + e18*e24*e49 + e19*e34*e48 + &
           e18*e34*e49 + e27*e34*e57 + e24*e37*e57 + &
           e27*e34*e67 + e24*e37*e67 + e15*e18*e67 + &
           e16*e18*e57 + e15*e19*e67 + e16*e19*e57
!c
      basis(28) = e15*e18*e58 + e15*e19*e59 + e16*e18*e68 + &
           e16*e19*e69 + e24*e28*e48 + e34*e38*e48 + &
           e24*e29*e49 + e34*e39*e49 + e25*e27*e57 + &
           e26*e27*e67 + e35*e37*e57 + e36*e37*e67
!c
      basis(29) = e15*e18*e59 + e15*e19*e58 + e16*e18*e69 + &
           e16*e19*e68 + e24*e38*e48 + e28*e34*e48 + &
           e24*e39*e49 + e29*e34*e49 + e26*e27*e57 + &
           e25*e27*e67 + e36*e37*e57 + e35*e37*e67 + &
           e24*e29*e48 + e24*e28*e49 + e34*e39*e48 + &
           e34*e38*e49 + e27*e35*e57 + e25*e37*e57 + &
           e27*e36*e67 + e26*e37*e67 + e15*e18*e68 + &
           e16*e18*e58 + e15*e19*e69 + e16*e19*e59
!c
      basis(30) = e15*e18*e69 + e15*e19*e68 + e16*e18*e59 + &
           e16*e19*e58 + e24*e39*e48 + e29*e34*e48 + &
           e24*e38*e49 + e28*e34*e49 + e27*e36*e57 + &
           e27*e35*e67 + e26*e37*e57 + e25*e37*e67
!c
      basis(31) = e15*e27*e48 + e15*e27*e49 + e16*e27*e48 + &
           e16*e27*e49 + e15*e37*e48 + e15*e37*e49 + &
           e16*e37*e48 + e16*e37*e49 + e18*e24*e57 + &
           e19*e24*e57 + e18*e34*e57 + e19*e34*e57 + &
           e18*e24*e67 + e19*e24*e67 + e18*e34*e67 + &
           e19*e34*e67
!c
      basis(32) = e15*e27*e58 + e15*e27*e59 + e16*e27*e68 + &
           e16*e27*e69 + e15*e37*e58 + e15*e37*e59 + &
           e16*e37*e68 + e16*e37*e69 + e15*e28*e48 + &
           e15*e38*e48 + e15*e29*e49 + e15*e39*e49 + &
           e16*e28*e48 + e16*e38*e48 + e16*e29*e49 + &
           e16*e39*e49 + e25*e27*e48 + e26*e27*e48 + &
           e35*e37*e48 + e36*e37*e48 + e25*e27*e49 + &
           e26*e27*e49 + e35*e37*e49 + e36*e37*e49 + &
           e24*e28*e57 + e24*e29*e57 + e34*e38*e57 + &
           e34*e39*e57 + e24*e28*e67 + e24*e29*e67 + &
           e34*e38*e67 + e34*e39*e67 + e18*e25*e57 + &
           e18*e35*e57 + e18*e26*e67 + e18*e36*e67 + &
           e19*e25*e57 + e19*e35*e57 + e19*e26*e67 + &
           e19*e36*e67 + e18*e24*e58 + e18*e24*e68 + &
           e19*e24*e59 + e19*e24*e69 + e18*e34*e58 + &
           e18*e34*e68 + e19*e34*e59 + e19*e34*e69
!c
      basis(33) = e15*e27*e68 + e15*e27*e69 + e16*e27*e58 + &
           e16*e27*e59 + e15*e37*e68 + e15*e37*e69 + &
           e16*e37*e58 + e16*e37*e59 + e15*e29*e48 + &
           e15*e39*e48 + e15*e28*e49 + e15*e38*e49 + &
           e16*e29*e48 + e16*e39*e48 + e16*e28*e49 + &
           e16*e38*e49 + e27*e35*e48 + e27*e36*e48 + &
           e25*e37*e48 + e26*e37*e48 + e27*e35*e49 + &
           e27*e36*e49 + e25*e37*e49 + e26*e37*e49 + &
           e24*e38*e57 + e24*e39*e57 + e28*e34*e57 + &
           e29*e34*e57 + e24*e38*e67 + e24*e39*e67 + &
           e28*e34*e67 + e29*e34*e67 + e18*e26*e57 + &
           e18*e36*e57 + e18*e25*e67 + e18*e35*e67 + &
           e19*e26*e57 + e19*e36*e57 + e19*e25*e67 + &
           e19*e35*e67 + e18*e24*e59 + e18*e24*e69 + &
           e19*e24*e58 + e19*e24*e68 + e18*e34*e59 + &
           e18*e34*e69 + e19*e34*e58 + e19*e34*e68
!c
      basis(34) = e15*e28*e57 + e15*e29*e57 + e16*e28*e67 + &
           e16*e29*e67 + e15*e38*e57 + e15*e39*e57 + &
           e16*e38*e67 + e16*e39*e67 + e18*e25*e48 + &
           e18*e35*e48 + e19*e25*e49 + e19*e35*e49 + &
           e18*e26*e48 + e18*e36*e48 + e19*e26*e49 + &
           e19*e36*e49 + e24*e27*e58 + e24*e27*e68 + &
           e34*e37*e58 + e34*e37*e68 + e24*e27*e59 + &
           e24*e27*e69 + e34*e37*e59 + e34*e37*e69
!c
      basis(35) = e15*e28*e58 + e15*e29*e59 + e16*e28*e68 + &
           e16*e29*e69 + e15*e38*e58 + e15*e39*e59 + &
           e16*e38*e68 + e16*e39*e69 + e25*e28*e48 + &
           e35*e38*e48 + e25*e29*e49 + e35*e39*e49 + &
           e26*e28*e48 + e36*e38*e48 + e26*e29*e49 + &
           e36*e39*e49 + e25*e27*e58 + e26*e27*e68 + &
           e35*e37*e58 + e36*e37*e68 + e25*e27*e59 + &
           e26*e27*e69 + e35*e37*e59 + e36*e37*e69 + &
           e24*e28*e58 + e24*e29*e59 + e34*e38*e58 + &
           e34*e39*e59 + e24*e28*e68 + e24*e29*e69 + &
           e34*e38*e68 + e34*e39*e69 + e25*e28*e57 + &
           e35*e38*e57 + e26*e28*e67 + e36*e38*e67 + &
           e25*e29*e57 + e35*e39*e57 + e26*e29*e67 + &
           e36*e39*e67 + e18*e25*e58 + e18*e26*e68 + &
           e19*e25*e59 + e19*e26*e69 + e18*e35*e58 + &
           e18*e36*e68 + e19*e35*e59 + e19*e36*e69
!c
      basis(36) = e15*e28*e59 + e15*e29*e58 + e16*e28*e69 + &
           e16*e29*e68 + e15*e38*e59 + e15*e39*e58 + &
           e16*e38*e69 + e16*e39*e68 + e25*e38*e48 + &
           e28*e35*e48 + e25*e39*e49 + e29*e35*e49 + &
           e26*e38*e48 + e28*e36*e48 + e26*e39*e49 + &
           e29*e36*e49 + e26*e27*e58 + e25*e27*e68 + &
           e36*e37*e58 + e35*e37*e68 + e26*e27*e59 + &
           e25*e27*e69 + e36*e37*e59 + e35*e37*e69 + &
           e24*e29*e58 + e24*e28*e59 + e34*e39*e58 + &
           e34*e38*e59 + e24*e29*e68 + e24*e28*e69 + &
           e34*e39*e68 + e34*e38*e69 + e28*e35*e57 + &
           e25*e38*e57 + e28*e36*e67 + e26*e38*e67 + &
           e29*e35*e57 + e25*e39*e57 + e29*e36*e67 + &
           e26*e39*e67 + e18*e25*e68 + e18*e26*e58 + &
           e19*e25*e69 + e19*e26*e59 + e18*e35*e68 + &
           e18*e36*e58 + e19*e35*e69 + e19*e36*e59
!c
      basis(37) = e15*e28*e67 + e15*e29*e67 + e16*e28*e57 + &
           e16*e29*e57 + e15*e38*e67 + e15*e39*e67 + &
           e16*e38*e57 + e16*e39*e57 + e19*e25*e48 + &
           e19*e35*e48 + e18*e25*e49 + e18*e35*e49 + &
           e19*e26*e48 + e19*e36*e48 + e18*e26*e49 + &
           e18*e36*e49 + e27*e34*e58 + e27*e34*e68 + &
           e24*e37*e58 + e24*e37*e68 + e27*e34*e59 + &
           e27*e34*e69 + e24*e37*e59 + e24*e37*e69
!c
      basis(38) = e15*e28*e68 + e15*e29*e69 + e16*e28*e58 + &
           e16*e29*e59 + e15*e38*e68 + e15*e39*e69 + &
           e16*e38*e58 + e16*e39*e59 + e25*e29*e48 + &
           e35*e39*e48 + e25*e28*e49 + e35*e38*e49 + &
           e26*e29*e48 + e36*e39*e48 + e26*e28*e49 + &
           e36*e38*e49 + e27*e35*e58 + e27*e36*e68 + &
           e25*e37*e58 + e26*e37*e68 + e27*e35*e59 + &
           e27*e36*e69 + e25*e37*e59 + e26*e37*e69 + &
           e24*e38*e58 + e24*e39*e59 + e28*e34*e58 + &
           e29*e34*e59 + e24*e38*e68 + e24*e39*e69 + &
           e28*e34*e68 + e29*e34*e69 + e26*e28*e57 + &
           e36*e38*e57 + e25*e28*e67 + e35*e38*e67 + &
           e26*e29*e57 + e36*e39*e57 + e25*e29*e67 + &
           e35*e39*e67 + e18*e25*e59 + e18*e26*e69 + &
           e19*e25*e58 + e19*e26*e68 + e18*e35*e59 + &
           e18*e36*e69 + e19*e35*e58 + e19*e36*e68
!c
      basis(39) = e15*e28*e69 + e15*e29*e68 + e16*e28*e59 + &
           e16*e29*e58 + e15*e38*e69 + e15*e39*e68 + &
           e16*e38*e59 + e16*e39*e58 + e25*e39*e48 + &
           e29*e35*e48 + e25*e38*e49 + e28*e35*e49 + &
           e26*e39*e48 + e29*e36*e48 + e26*e38*e49 + &
           e28*e36*e49 + e27*e36*e58 + e27*e35*e68 + &
           e26*e37*e58 + e25*e37*e68 + e27*e36*e59 + &
           e27*e35*e69 + e26*e37*e59 + e25*e37*e69 + &
           e24*e39*e58 + e24*e38*e59 + e29*e34*e58 + &
           e28*e34*e59 + e24*e39*e68 + e24*e38*e69 + &
           e29*e34*e68 + e28*e34*e69 + e28*e36*e57 + &
           e26*e38*e57 + e28*e35*e67 + e25*e38*e67 + &
           e29*e36*e57 + e26*e39*e57 + e29*e35*e67 + &
           e25*e39*e67 + e18*e25*e69 + e18*e26*e59 + &
           e19*e25*e68 + e19*e26*e58 + e18*e35*e69 + &
           e18*e36*e59 + e19*e35*e68 + e19*e36*e58
!c
      basis(40) = e25*e28*e58 + e25*e29*e59 + e26*e28*e68 + &
           e26*e29*e69 + e35*e38*e58 + e35*e39*e59 + &
           e36*e38*e68 + e36*e39*e69
!c
      basis(41) = e25*e28*e59 + e25*e29*e58 + e26*e28*e69 + &
           e26*e29*e68 + e35*e38*e59 + e35*e39*e58 + &
           e36*e38*e69 + e36*e39*e68 + e25*e38*e58 + &
           e28*e35*e58 + e25*e39*e59 + e29*e35*e59 + &
           e26*e38*e68 + e28*e36*e68 + e26*e39*e69 + &
           e29*e36*e69 + e26*e28*e58 + e25*e28*e68 + &
           e36*e38*e58 + e35*e38*e68 + e26*e29*e59 + &
           e25*e29*e69 + e36*e39*e59 + e35*e39*e69
!c
      basis(42) = e25*e28*e69 + e25*e29*e68 + e26*e28*e59 + &
           e26*e29*e58 + e35*e38*e69 + e35*e39*e68 + &
           e36*e38*e59 + e36*e39*e58 + e25*e39*e58 + &
           e29*e35*e58 + e25*e38*e59 + e28*e35*e59 + &
           e26*e39*e68 + e29*e36*e68 + e26*e38*e69 + &
           e28*e36*e69 + e28*e36*e58 + e28*e35*e68 + &
           e26*e38*e58 + e25*e38*e68 + e29*e36*e59 + &
           e29*e35*e69 + e26*e39*e59 + e25*e39*e69
!c
      basis(43) = e25*e38*e69 + e25*e39*e68 + e26*e38*e59 + &
           e26*e39*e58 + e28*e35*e69 + e29*e35*e68 + &
           e28*e36*e59 + e29*e36*e58
!c ---------------------------------------------------------------------
      return
    end subroutine calc_basis_1
    
  end module IPModel_WaterTrimer_Gillan_module
  


