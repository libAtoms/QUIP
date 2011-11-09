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
!X IPModel_WaterDimer_Gillan
!X
!% Interaction of two water dimers, parametrised by Mike Gillan in 2011
!% 
!% only expected to be accurate at long range (> 5 A). Parameters are
!% fitted to MP2 AVTZ energies
!%
!% Energy only
!% 
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#include "error.inc"

module IPModel_WaterDimer_Gillan_module

use libatoms_module

use mpi_context_module
use QUIP_Common_module

implicit none
private

include 'IPModel_interface.h'

public :: IPModel_WaterDimer_Gillan
type IPModel_WaterDimer_Gillan
  real(dp) :: cutoff = 0.0_dp
  character(len=20) :: fname_d, fname_q
  real(dp) :: two_body_weight_roo = 0.0_dp
  real(dp) :: two_body_weight_delta = 0.0_dp
  logical  :: do_two_body_weight = .false.
  type(Spline) :: two_body_weight
end type IPModel_WaterDimer_Gillan

logical, private :: parse_in_ip, parse_matched_label
type(IPModel_WaterDimer_Gillan), private, pointer :: parse_ip

interface Initialise
  module procedure IPModel_WaterDimer_Gillan_Initialise_str
end interface Initialise

interface Finalise
  module procedure IPModel_WaterDimer_Gillan_Finalise
end interface Finalise

interface Print
  module procedure IPModel_WaterDimer_Gillan_Print
end interface Print

interface Calc
  module procedure IPModel_WaterDimer_Gillan_Calc
end interface Calc

contains

subroutine IPModel_WaterDimer_Gillan_Initialise_str(this, args_str, param_str, error)
  type(IPModel_WaterDimer_Gillan), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str
  type(Dictionary)                :: params
  integer, optional, intent(out) :: error
  integer:: init3d
  real(dp) :: theta, r1, r2, f1, f2, f3

  INIT_ERROR(error)

  call Finalise(this)


  call initialise(params)
  call param_register(params, 'two_body_weight_roo', '0.0', this%two_body_weight_roo, has_value_target=this%do_two_body_weight, help_string="if set, apply weight function to 2-body energy and force based on O-O distance. For a positive two_body_weight_delta, weight is 1 for rOO < two_body_weight_roo-two_body_weight_delta and weight is 0 for rOO > two_body_weight_roo+two_body_weight_delta")
  call param_register(params, 'two_body_weight_delta', '0.25', this%two_body_weight_delta, help_string="width of weighting function for  two_body energy and force based on O-O distance. For weighting to take effect, two_body_weight_roo needs to be explicitly set. For a positive two_body_weight_delta, weight is 1 for rOO < two_body_weight_roo-two_body_weight_delta and weight is 0 for rOO > two_body_weight_roo+two_body_weight_delta")
  
  if(.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_WaterDimer_Gillan_Initialise args_str')) then
     RAISE_ERROR("IPModel_WaterDimer_Gillan_Init failed to parse args_str='"//trim(args_str)//"'", error)
  endif
  call finalise(params)


  if(this%do_two_body_weight) then
     call initialise(this%two_body_weight, (/this%two_body_weight_roo - this%two_body_weight_delta, this%two_body_weight_roo + this%two_body_weight_delta/), (/1.0_dp, 0.0_dp/), 0.0_dp, 0.0_dp)
  end if

  ! read interpolation tables
  init3d = 0

  this%fname_d = "dip_grid_1-1331"
  this%fname_q = "quad_grid_1-1331"

  call lin3d_2(init3d,theta,r1,r2,f1,f2, this%fname_d)
  call lin3d_3(init3d,theta,r1,r2,f1,f2,f3,this%fname_q)


end subroutine IPModel_WaterDimer_Gillan_Initialise_str


subroutine IPModel_WaterDimer_Gillan_Finalise(this)
  type(IPModel_WaterDimer_Gillan), intent(inout) :: this

  ! Add finalisation code here

end subroutine IPModel_WaterDimer_Gillan_Finalise


subroutine IPModel_WaterDimer_Gillan_Calc(this, at, e, local_e, f, virial, local_virial, args_str, mpi, error)
   type(IPModel_WaterDimer_Gillan), intent(inout):: this
   type(Atoms), intent(inout)      :: at
   real(dp), intent(out), optional :: e, local_e(:)
   real(dp), intent(out), optional :: f(:,:), local_virial(:,:)   !% Forces, dimensioned as \texttt{f(3,at%N)}, local virials, dimensioned as \texttt{local_virial(9,at%N)} 
   real(dp), intent(out), optional :: virial(3,3)
   character(len=*), optional      :: args_str
   type(MPI_Context), intent(in), optional :: mpi
   integer, intent(out), optional :: error

   real(dp) :: dip_frac, pol_o, pol_h, c6oo, c6oh, c6hh
   real(dp) :: atpos(2,3,3), axis(3)
   real(dp) :: e_dqtot, e_ind, e_disp, elong_tot
   integer :: nmolec, init3d
   type(Quaternion) :: quat

   real(dp):: weight, dweight(3), rOiOj

   INIT_ERROR(error)


   ! Parameters

   dip_frac = 0.4_dp
   pol_o = 7.0_dp
   pol_h = 2.0_dp
   c6oo = 25.0_dp
   c6oh = 4.5_dp
   c6hh = 2.0_dp
   nmolec = 2

   ! Rotate atoms so that rOO is along the z axis and O1 is at the origin

   !get rotation axis
   axis = (at%pos(:,4)-at%pos(:,1)) .cross. (/0.0_dp,0.0_dp,1.0_dp/)
   atpos(1,1,:) = at%pos(:,1)-at%pos(:,1)
   atpos(1,2,:) = at%pos(:,2)-at%pos(:,1)
   atpos(1,3,:) = at%pos(:,3)-at%pos(:,1)
   atpos(2,1,:) = at%pos(:,4)-at%pos(:,1)
   atpos(2,2,:) = at%pos(:,5)-at%pos(:,1)
   atpos(2,3,:) = at%pos(:,6)-at%pos(:,1)

   if(axis .fne. (/0.0_dp,0.0_dp,0.0_dp/)) then
      quat = orientation(at%pos(:,4)-at%pos(:,1), axis, (/0.0_dp,0.0_dp,1.0_dp/), axis)

      call rotate(atpos(1,1,:), quat)
      call rotate(atpos(1,2,:), quat)
      call rotate(atpos(1,3,:), quat)
      call rotate(atpos(2,1,:), quat)
      call rotate(atpos(2,2,:), quat)
      call rotate(atpos(2,3,:), quat)
   end if

   ! make that call
   call h2o_dimer_far(dip_frac,pol_o,pol_h,c6oo,c6oh,c6hh,atpos,e_dqtot,e_ind,e_disp,elong_tot)

   if(this%do_two_body_weight) then
      rOiOj = norm(diff_min_image(at, 1, 4))
      weight = spline_value(this%two_body_weight, rOiOj)
      dweight = spline_deriv(this%two_body_weight, rOiOj)
   else
      weight = 1.0_dp
      dweight = 0.0_dp
   end if

   if (present(e)) e = elong_tot*HARTREE*weight
   if (present(local_e)) then
      call check_size('Local_E',local_e,(/at%N/),'IPModel_WaterDimer_Gillan_Calc', error)
      local_e = 0.0_dp
   endif
   if (present(f)) then
      RAISE_ERROR('Forces not implemented', error)
   end if
   if (present(virial)) virial = 0.0_dp
   if (present(local_virial)) then
      call check_size('Local_virial',local_virial,(/9,at%Nbuffer/),'IPModel_WaterDimer_Gillan_Calc', error)
      local_virial = 0.0_dp
   endif


end subroutine IPModel_WaterDimer_Gillan_Calc


subroutine IPModel_WaterDimer_Gillan_Print(this, file)
  type(IPModel_WaterDimer_Gillan), intent(in) :: this
  type(Inoutput), intent(inout),optional :: file

  call Print("IPModel_WaterDimer_Gillan : WaterDimer_Gillan Potential", file=file)
  call Print("IPModel_WaterDimer_Gillan : cutoff = " // this%cutoff, file=file)
  if(this%do_two_body_weight) then
     call print("Two-body term is weighted with parameters x0="//this%two_body_weight_roo//" delta="//this%two_body_weight_delta, file=file)
  end if

end subroutine IPModel_WaterDimer_Gillan_Print

! =====================================================================
!   sbrt computes long-range 2-body energy of H2O dimer, using simple
!   distributed multipole model plus simple models for induction
!   and dispersion.
! ---------------------------------------------------------------------
      subroutine h2o_dimer_far(dip_frac,pol_o,pol_h,c6oo,c6oh,c6hh,atpos,e_dqtot,e_ind,e_disp,elong_tot)
      implicit none
! ... variables in argument list
      real(dp):: dip_frac, pol_o, pol_h, c6oo, c6oh, c6hh
      real(dp):: atpos(2,3,3)
      real(dp):: e_dqtot, e_ind, e_disp, elong_tot
! ... local variables
      real*8 pi, bohr
      parameter (pi = 3.14159265359, bohr = 0.52917720859)
      integer init3d
      integer mo, na, ns, i, j, k, m
      real*8 f1, f2, f3
      real*8 r1, r2, theta, theta_deg, r1sq, r2sq, scprod
      real*8 sinfun, cosfun, tanfun, cotfun
      real*8 mux, muz, qxx, qzz, qxz, rmu_scal
      real*8 xmagsq, zmagsq, xnorm, znorm
      real*8 bvec_0(2,3), bvec_bohr_0(2,3), bvec_r(2,2,3)
      real*8 rlocvec(3,3)
      real*8 univec(2,3)
      real*8 dipvec_0(3)
      real*8 quadten_0(9), quadten_d(9)
      real*8 ddvec_0(3,3)
      real*8 ddvec_r(2,3,3)
      real*8 quadten_r(2,9)
      real*8 dv1o(3), dv1h1(3), dv1h2(3), dv2o(3), dv2h1(3),dv2h2(3)
      real*8 qt1(9) ,qt2(9), sepoo(3)
      real*8 rroo(3), rroh1(3), rroh2(3), rrh1o(3), rrh1h1(3), rrh1h2(3), rrh2o(3), rrh2h1(3), rrh2h2(3)
      real*8 dd_int_oo, dd_int_oh1, dd_int_oh2, dd_int_h1o, dd_int_h1h1, dd_int_h1h2, dd_int_h2o, dd_int_h2h1,  dd_int_h2h2, dd_int
      real*8 dq_int_oo, dq_int_h1o, dq_int_h2o, dq_int
      real*8 qd_int_oo, qd_int_oh1, qd_int_oh2, qd_int, qq_int
      real*8 e_ind_1o, e_ind_1h1, e_ind_1h2,  e_ind_2o, e_ind_2h1, e_ind_2h2
      real*8 d6oo, d6oh1, d6oh2, d6h1o, d6h1h1, d6h1h2, d6h2o, d6h2h1, d6h2h2
      real*8 evec1(3), evec2(3), evec3(3), evec4(3), evect(3)
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!   begin loop over monomers
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      do mo = 1, 2
! ... make the two bond vectors bvec_r drawn from O to H1 and O to H2 
! ... for present monomer
!        write(*,111)
!  111   format(/'bond lengths r1, r2 and bond angle theta ',
!     *   'for monomers:')
!        write(*,112) mo
!  112   format('monomer no:',i5)
        r1sq = 0.d0
        r2sq = 0.d0
        scprod = 0.d0
        do i = 1, 3
          bvec_r(mo,1,i) = atpos(mo,2,i) - atpos(mo,1,i)
          r1sq = r1sq + bvec_r(mo,1,i)*bvec_r(mo,1,i)
          bvec_r(mo,2,i) = atpos(mo,3,i) - atpos(mo,1,i)
          r2sq = r2sq + bvec_r(mo,2,i)*bvec_r(mo,2,i)
          scprod = scprod + bvec_r(mo,1,i)*bvec_r(mo,2,i)
        enddo
        r1 = sqrt(r1sq)
        r2 = sqrt(r2sq)
        theta = acos(scprod/(r1*r2))
        theta_deg = (180.d0/pi)*theta
! ... unit vectors along local x, y and z axes
        do i = 1, 3
          univec(1,i) = bvec_r(mo,1,i)/r1
          univec(2,i) = bvec_r(mo,2,i)/r2
        enddo
        xmagsq = 0.d0
        zmagsq = 0.d0
        do i = 1, 3
          rlocvec(1,i) = univec(1,i) - univec(2,i)
          xmagsq = xmagsq + rlocvec(1,i)*rlocvec(1,i)
          rlocvec(3,i) = univec(1,i) + univec(2,i)
          zmagsq = zmagsq + rlocvec(3,i)*rlocvec(3,i)
        enddo
        xnorm = sqrt(xmagsq)
        znorm = sqrt(zmagsq)
        do i = 1, 3
          rlocvec(1,i) = rlocvec(1,i)/xnorm
          rlocvec(3,i) = rlocvec(3,i)/znorm
        enddo
        rlocvec(2,1) = rlocvec(3,2)*rlocvec(1,3) - rlocvec(3,3)*rlocvec(1,2)
        rlocvec(2,2) = rlocvec(3,3)*rlocvec(1,1) - rlocvec(3,1)*rlocvec(1,3)
        rlocvec(2,3) = rlocvec(3,1)*rlocvec(1,2) - rlocvec(3,2)*rlocvec(1,1)
!        write(*,121)
!  121   format(/'local x, y and z unit vectors:')
!        do ns = 1, 3
!          write(*,122) (rlocvec(ns,i), i = 1, 3)
!  122     format(3x,3f12.6)
!        enddo
! ... bond vectors in local frame
        sinfun = sin(0.5d0*theta)
        cosfun = cos(0.5d0*theta)
        tanfun = sinfun/cosfun
        cotfun = cosfun/sinfun
        bvec_0(1,1) = r1*sinfun
        bvec_0(1,2) = 0.d0
        bvec_0(1,3) = r1*cosfun
        bvec_0(2,1) = -r2*sinfun
        bvec_0(2,2) = 0.d0
        bvec_0(2,3) = r2*cosfun
        do ns = 1, 2
          do i = 1, 3
            bvec_bohr_0(ns,i) = bvec_0(ns,i)/bohr
          enddo
        enddo
! ... for present monomer, compute components of
! ... dipole vector and quadrupole tensor in local frame.
! ... Local frame is defined thus:
! ... z-axis is along bisector of H-O-H angle, with bond vectors drawn 
! ... from O to H having positive z-components
! ... x-axis is in molecular plane and is perpendicular to z-axis,
! ... with bond vector drawn from O to H1 having positive x-component
! ... bond vector from O to H2 having negative x-component.
! ... y-axis is perpendicular to x- and z-axes, in direction such that
! ... x-, y- and z-axes form a normal right-handed set.
! ... dipole vector
        init3d = 1
        call lin3d_2(init3d,theta_deg,r1,r2,mux,muz,repeat(' ',20))
        dipvec_0(1) = mux
        dipvec_0(2) = 0.d0
        dipvec_0(3) = muz
!        write(*,131)
!  131   format('x-, y- and z-components of dipole vector:')
!        write(*,132) (dipvec_0(i), i = 1, 3)
!  132   format(3x,3f15.6)
! ... compute components of quadrupole tensor
        init3d = 1
        call lin3d_3(init3d,theta_deg,r1,r2,qxx,qzz,qxz,repeat(' ',20))
        quadten_0(1) = qxx
        quadten_0(2) = 0.d0
        quadten_0(3) = qxz
        quadten_0(4) = 0.d0
        quadten_0(5) = -qxx - qzz
        quadten_0(6) = 0.d0
        quadten_0(7) = qxz
        quadten_0(8) = 0.d0
        quadten_0(9) = qzz
!        write(*,141)
!  141   format('xx, xy, xz, yx, yy, yz, zx, zy, zz components of ',
!     *     'quadrupole tensor:')
!        do i = 1, 3
!          write(*,142) (quadten_0(3*(i-1)+j), j = 1, 3)
!  142     format(3x,3f15.6)
!        enddo
! ... distribute dipole vector between O and H atoms
        ddvec_0(1,1) = (1.d0 - dip_frac)*dipvec_0(1)
        ddvec_0(1,2) = 0.d0
        ddvec_0(1,3) = (1.d0 - dip_frac)*dipvec_0(3)
        ddvec_0(2,1) = 0.5d0*dip_frac*(dipvec_0(1) +  tanfun*dipvec_0(3))
        ddvec_0(2,2) = 0.d0
        ddvec_0(2,3) = 0.5d0*dip_frac*(cotfun*dipvec_0(1) +  dipvec_0(3))
        ddvec_0(3,1) = 0.5d0*dip_frac*(dipvec_0(1) -  tanfun*dipvec_0(3))
        ddvec_0(3,2) = 0.d0
        ddvec_0(3,3) = 0.5d0*dip_frac*(-cotfun*dipvec_0(1) +  dipvec_0(3))
!        write(*,151)
!  151   format('components of dipole vectors on O, H1 and H2: ')
!        do i = 1, 3
!          write(*,152) (ddvec_0(i,j), j = 1, 3)
!  152     format(3x,3f15.6)
!        enddo
! ... adjust quadrupole tensor to compensate for distributed dipole
        do ns = 1, 2
          rmu_scal = 0.d0
          do i = 1, 3
            rmu_scal = rmu_scal +   bvec_bohr_0(ns,i)*ddvec_0(ns+1,i)
          enddo
          do i = 1, 3
            k = 3*(i-1) + i
            quadten_0(k) = quadten_0(k) + rmu_scal
          enddo
          do i = 1, 3
            do j = 1, 3
              k = 3*(i-1) + j
            quadten_0(k) = quadten_0(k) -   1.5d0*(bvec_bohr_0(ns,i)*ddvec_0(ns+1,j) + ddvec_0(ns+1,i)*bvec_bohr_0(ns,j))
            enddo
          enddo
        enddo
!        write(*,161)
!  161   format('quad tensor adj to compensate for distrib dip: ')
!        do i = 1, 3
!          write(*,162) (quadten_0(3*(i-1)+j), j = 1, 3)
!  162     format(3x,3f15.6)
!        enddo
! ... distributed dipole vectors in lab frame
        do ns = 1, 3
          do i = 1, 3
            ddvec_r(mo,ns,i) = 0.d0
            do k = 1, 3
              ddvec_r(mo,ns,i) = ddvec_r(mo,ns,i) + ddvec_0(ns,k)*rlocvec(k,i)
            enddo
          enddo
        enddo
!        write(*,171)
!  171   format(/'distributed dipole vectors in lab frame:')
!        do ns = 1, 3
!          write(*,172) (ddvec_r(mo,ns,i), i = 1, 3)
!  172     format(3x,3f15.6)
!        enddo
! ... quadrupole tensor in lab frame
        do i = 1, 3
          do j = 1, 3
            m = 3*(i-1) + j
            quadten_d(m) = 0.d0
            do k = 1, 3
              ns = 3*(i-1) + k
              quadten_d(m) = quadten_d(m) +  quadten_0(ns)*rlocvec(k,j)
            enddo
          enddo
        enddo
        do i = 1, 3
          do j = 1, 3
            m = 3*(i-1) + j
            quadten_r(mo,m) = 0.d0
            do k = 1, 3
              ns = 3*(k-1) + j
              quadten_r(mo,m) = quadten_r(mo,m) +   rlocvec(k,i)*quadten_d(ns)
            enddo
          enddo
        enddo
!        write(*,181)
!  181   format(/'quadrupole tensor in lab frame:')
!        do i = 1, 3
!          write(*,182) (quadten_r(mo,3*(i-1)+j), j = 1, 3)
!  182     format(3x,3f15.6)
!        enddo
      enddo
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!   end loop over monomers
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!
! ... calculate 1st-order electrostatic interactions .........
!      write(*,211)
!  211 format(/75('-'))
!      write(*,212)
!  212 format(/'d-d, d-q, q-d and q-q interactions: (a.u.)')
! ... temp arrays for distributed dipoles and quadrupoles
      do i = 1, 3
        dv1o(i) = ddvec_r(1,1,i)
        dv1h1(i) = ddvec_r(1,2,i)
        dv1h2(i) = ddvec_r(1,3,i)
        dv2o(i) = ddvec_r(2,1,i)
        dv2h1(i) = ddvec_r(2,2,i)
        dv2h2(i) = ddvec_r(2,3,i)
      enddo
      do k = 1, 9
        qt1(k) = quadten_r(1,k)
        qt2(k) = quadten_r(2,k)
      enddo
! ... vector posn of O in 2nd monomer relative to O in 1st monomer
      do i = 1, 3
        sepoo(i) = atpos(2,1,i) - atpos(1,1,i)
      enddo
! ... make interatomic separation vectors
      do i = 1, 3
        rroo(i) = sepoo(i)/bohr
        rroh1(i) = (bvec_r(2,1,i) + sepoo(i))/bohr
        rroh2(i) = (bvec_r(2,2,i) + sepoo(i))/bohr
        rrh1o(i) = (-bvec_r(1,1,i) + sepoo(i))/bohr
        rrh1h1(i) = (bvec_r(2,1,i) - bvec_r(1,1,i) + sepoo(i))/bohr
        rrh1h2(i) = (bvec_r(2,2,i) - bvec_r(1,1,i) + sepoo(i))/bohr
        rrh2o(i) = (-bvec_r(1,2,i) + sepoo(i))/bohr
        rrh2h1(i) = (bvec_r(2,1,i) - bvec_r(1,2,i) + sepoo(i))/bohr
        rrh2h2(i) = (bvec_r(2,2,i) - bvec_r(1,2,i) + sepoo(i))/bohr
      enddo
! ... dipole-dipole interactions
      call dip_dip(dv1o,dv2o,rroo,dd_int_oo)
      call dip_dip(dv1o,dv2h1,rroh1,dd_int_oh1)
      call dip_dip(dv1o,dv2h2,rroh2,dd_int_oh2)
      call dip_dip(dv1h1,dv2o,rrh1o,dd_int_h1o)
      call dip_dip(dv1h1,dv2h1,rrh1h1,dd_int_h1h1)
      call dip_dip(dv1h1,dv2h2,rrh1h2,dd_int_h1h2)
      call dip_dip(dv1h2,dv2o,rrh2o,dd_int_h2o)
      call dip_dip(dv1h2,dv2h1,rrh2h1,dd_int_h2h1)
      call dip_dip(dv1h2,dv2h2,rrh2h2,dd_int_h2h2)
      dd_int = dd_int_oo + dd_int_oh1 + dd_int_oh2 + dd_int_h1o + dd_int_h1h1 + dd_int_h1h2 + dd_int_h2o + dd_int_h2h1 +   dd_int_h2h2
! ... dipole-quadrupole interactions
      call dip_quad(dv1o,qt2,rroo,dq_int_oo)
      call dip_quad(dv1h1,qt2,rrh1o,dq_int_h1o)
      call dip_quad(dv1h2,qt2,rrh2o,dq_int_h2o)
      dq_int = dq_int_oo + dq_int_h1o + dq_int_h2o
! ... quadrupole-dipole interactions
      call quad_dip(qt1,dv2o,rroo,qd_int_oo)
      call quad_dip(qt1,dv2h1,rroh1,qd_int_oh1)
      call quad_dip(qt1,dv2h2,rroh2,qd_int_oh2)
      qd_int = qd_int_oo + qd_int_oh1 + qd_int_oh2
! ... quadrupole-quadrupole interaction
      call quad_quad(qt1,qt2,rroo,qq_int)
!      write(*,221) dd_int, dq_int, qd_int, qq_int
!  221 format(3x,4e15.6)
      e_dqtot = dd_int + dq_int + qd_int + qq_int
! --- induction energy ------------------------------------------------
! ... pol of mono2/O from dip and quad on mono1
      call dip_field(rroo,dv1o,evec1)
      call dip_field(rrh1o,dv1h1,evec2)
      call dip_field(rrh2o,dv1h2,evec3)
      call quad_field(rroo,qt1,evec4)
      do i = 1, 3
        evect(i) = evec1(i) + evec2(i) + evec3(i) + evec4(i)
      enddo
      e_ind_2o =  -0.5d0*pol_o*(evect(1)**2 + evect(2)**2 + evect(3)**2)
! ... pol of mono2/H1 from dip and quad on mono1
      call dip_field(rroh1,dv1o,evec1)
      call dip_field(rrh1h1,dv1h1,evec2)
      call dip_field(rrh2h1,dv1h2,evec3)
      call quad_field(rroh1,qt1,evec4)
      do i = 1, 3
        evect(i) = evec1(i) + evec2(i) + evec3(i) + evec4(i)
      enddo
      e_ind_2h1 = -0.5d0*pol_h*(evect(1)**2 + evect(2)**2 + evect(3)**2)
! ... pol of mono2/H2 from dip and quad on mono1
      call dip_field(rroh2,dv1o,evec1)
      call dip_field(rrh1h2,dv1h1,evec2)
      call dip_field(rrh2h2,dv1h2,evec3)
      call quad_field(rroh2,qt1,evec4)
      do i = 1, 3
        evect(i) = evec1(i) + evec2(i) + evec3(i) + evec4(i)
      enddo
      e_ind_2h2 =  -0.5d0*pol_h*(evect(1)**2 + evect(2)**2 + evect(3)**2)
! ... pol of mono1/O from dip and quad on mono2
      call dip_field(rroo,dv2o,evec1)
      call dip_field(rroh1,dv2h1,evec2)
      call dip_field(rroh2,dv2h2,evec3)
      call quad_field(rroo,qt2,evec4)
      do i = 1, 3
        evect(i) = evec1(i) + evec2(i) + evec3(i) - evec4(i)
      enddo
      e_ind_1o =  -0.5d0*pol_o*(evect(1)**2 + evect(2)**2 + evect(3)**2)
! ... pol of mono1/H1 from dip and quad on mono2
      call dip_field(rrh1o,dv2o,evec1)
      call dip_field(rrh1h1,dv2h1,evec2)
      call dip_field(rrh1h2,dv2h2,evec3)
      call quad_field(rrh1o,qt2,evec4)
      do i = 1, 3
        evect(i) = evec1(i) + evec2(i) + evec3(i) - evec4(i)
      enddo
      e_ind_1h1 = -0.5d0*pol_h*(evect(1)**2 + evect(2)**2 + evect(3)**2)
! ... pol of mono1/H2 from dip and quad on mono2
      call dip_field(rrh2o,dv2o,evec1)
      call dip_field(rrh2h1,dv2h1,evec2)
      call dip_field(rrh2h2,dv2h2,evec3)
      call quad_field(rrh2o,qt2,evec4)
      do i = 1, 3
        evect(i) = evec1(i) + evec2(i) + evec3(i) - evec4(i)
      enddo
      e_ind_1h2 = -0.5d0*pol_h*(evect(1)**2 + evect(2)**2 + evect(3)**2)
      e_ind = e_ind_1o + e_ind_1h1 + e_ind_1h2 + e_ind_2o + e_ind_2h1 + e_ind_2h2
! --- dispersion energy ----------------------------------------------
      d6oo = (rroo(1)**2 + rroo(2)**2 + rroo(3)**2)**3
      d6oh1 = (rroh1(1)**2 + rroh1(2)**2 + rroh1(3)**2)**3
      d6oh2 = (rroh2(1)**2 + rroh2(2)**2 + rroh2(3)**2)**3
      d6h1o = (rrh1o(1)**2 + rrh1o(2)**2 + rrh1o(3)**2)**3
      d6h1h1 = (rrh1h1(1)**2 + rrh1h1(2)**2 + rrh1h1(3)**2)**3
      d6h1h2 = (rrh1h2(1)**2 + rrh1h2(2)**2 + rrh1h2(3)**2)**3
      d6h2o = (rrh2o(1)**2 + rrh2o(2)**2 + rrh2o(3)**2)**3
      d6h2h1 = (rrh2h1(1)**2 + rrh2h1(2)**2 + rrh2h1(3)**2)**3
      d6h2h2 = (rrh2h2(1)**2 + rrh2h2(2)**2 + rrh2h2(3)**2)**3
      e_disp = -(c6oo/d6oo + c6oh/d6oh1 + c6oh/d6oh2 + c6oh/d6h1o + c6hh/d6h1h1 + c6hh/d6h1h2 + c6oh/d6h2o + c6hh/d6h2h1 + c6hh/d6h2h2)
! ... output components of and total 2-body energy
      elong_tot = e_dqtot + e_ind + e_disp
!      write(*,231)
!  231 format('1st-order e-stat, induct, disp, total energy (a.u.):')
!      write(*,232) e_dqtot, e_ind, e_disp, elong_tot
!  232 format(3x,3e15.6,3x,e15.6)
! ---------------------------------------------------------------------
      return
    end subroutine h2o_dimer_far
! =====================================================================

! =====================================================================
!   sbrt lin3d_2: performs 3D linear interpolation to evaluate two
!   functions f1, f2 of 3 variables, the values of the functions
!   being given on a regular 3D grid. Sbrt is called first with
!   flag init3d = 0, whereupon it reads in the grid data from
!   a file, which contains also the numbers of grid points and
!   the upper and lower limits of x, y and z. On subsequent call
!   with init3d != 0, sbrt uses linear interpolation to evaluate
!   f1 and f2 at requested point (x,y,z).
! ---------------------------------------------------------------------
      subroutine lin3d_2(init3d,x,y,z,f1,f2,fname)
      implicit none
      character*20 fname
      integer ngrid_mx, ngrid3_mx
      parameter (ngrid_mx = 20, ngrid3_mx = ngrid_mx*ngrid_mx*ngrid_mx)
      integer init3d
      real*8 x, y, z, f1, f2
      integer igotdata
      integer ix, iy, iz, nx, ny, nz
      integer nxsup, nysup, nzsup
      integer itab_print
      real*8 xinf, xsup, yinf, ysup, zinf, zsup
      real*8 deltax, deltay, deltaz
      real*8 f1tab(ngrid_mx,ngrid_mx,ngrid_mx)
      real*8 f2tab(ngrid_mx,ngrid_mx,ngrid_mx)
      real*8 xi, eta, zeta, px, py, pz, qx, qy, qz
      save
      data igotdata /0/
      data itab_print /0/
      data nxsup, nysup, nzsup / 3 * 0 /
      data xinf, xsup, yinf, ysup, zinf, zsup / 6 * 0.d0 /
      data deltax, deltay, deltaz / 3 * 0.d0 /
      data f1tab / ngrid3_mx * 0.d0 /
      data f2tab / ngrid3_mx * 0.d0 /
! ---------------------------------------------------------------------
      if(init3d .eq. 0) then
        igotdata = 1
        open(unit=21,file=fname)
        read(21,*) nxsup, xinf, xsup
        read(21,*) nysup, yinf, ysup
        read(21,*) nzsup, zinf, zsup
        deltax = (xsup - xinf)/float(nxsup - 1)
        deltay = (ysup - yinf)/float(nysup - 1)
        deltaz = (zsup - zinf)/float(nzsup - 1)
! .....................................................................
        do nx = 1, nxsup
          do ny = 1, nysup
            do nz = 1, nzsup
              read(21,*) ix, iy, iz,  f1tab(ix,iy,iz), f2tab(ix,iy,iz)
            enddo
          enddo
        enddo
        close(unit=21)
      else
        xi = (x - xinf)/deltax
        if((x - xinf) .lt. 0.d0) then
          ix = 0
        elseif(((x - xinf) .ge. 0.d0) .and. ((x - xsup) .le. 0.d0)) then
          ix = int(xi)
        else 
          ix = nxsup - 2
        endif
        px = xi - float(ix)
        qx = 1.d0 - px
        ix = ix + 1
        eta = (y - yinf)/deltay
        if((y - yinf) .lt. 0.d0) then
          iy = 0
        elseif(((y - yinf) .ge. 0.d0) .and. ((y - ysup) .le. 0.d0)) then
          iy = int(eta)
        else 
          iy = nysup - 2
        endif
        py = eta - float(iy)
        qy = 1.d0 - py
        iy = iy + 1
        zeta = (z - zinf)/deltaz
        if((z - zinf) .lt. 0.d0) then
          iz = 0
        elseif(((z - zinf) .ge. 0.d0) .and. ((z - zsup) .le. 0.d0)) then
          iz = int(zeta)
        else 
          iz = nzsup - 2
        endif
        pz = zeta - float(iz)
        qz = 1.d0 - pz
        iz = iz + 1
        f1 =      qx*qy*qz*f1tab(ix,iy,iz) + qx*qy*pz*f1tab(ix,iy,iz+1)    + qx*py*qz*f1tab(ix,iy+1,iz) + qx*py*pz*f1tab(ix,iy+1,iz+1)   + px*qy*qz*f1tab(ix+1,iy,iz) + px*qy*pz*f1tab(ix+1,iy,iz+1)   + px*py*qz*f1tab(ix+1,iy+1,iz) + px*py*pz*f1tab(ix+1,iy+1,iz+1)
        f2 =      qx*qy*qz*f2tab(ix,iy,iz) + qx*qy*pz*f2tab(ix,iy,iz+1)    + qx*py*qz*f2tab(ix,iy+1,iz) + qx*py*pz*f2tab(ix,iy+1,iz+1)   + px*qy*qz*f2tab(ix+1,iy,iz) + px*qy*pz*f2tab(ix+1,iy,iz+1)   + px*py*qz*f2tab(ix+1,iy+1,iz) + px*py*pz*f2tab(ix+1,iy+1,iz+1)
      endif
! ---------------------------------------------------------------------
      return
    end subroutine lin3d_2
! =====================================================================


! =====================================================================
!   sbrt lin3d_3: performs 3D linear interpolation to evaluate three
!   functions f1, f2, f3 of 3 variables, the values of the functions
!   being given on a regular 3D grid. Sbrt is called first with
!   flag init3d = 0, whereupon it reads in the grid data from
!   a file, which contains also the numbers of grid points and
!   the upper and lower limits of x, y and z. On subsequent call
!   with init3d != 0, sbrt uses linear interpolation to evaluate
!   f1, f2 and f3 at requested point (x,y,z).
! ---------------------------------------------------------------------
      subroutine lin3d_3(init3d,x,y,z,f1,f2,f3,fname)
      implicit none
      character*20 fname
      integer ngrid_mx, ngrid3_mx
      parameter (ngrid_mx = 20, ngrid3_mx = ngrid_mx*ngrid_mx*ngrid_mx)
      integer init3d
      real*8 x, y, z, f1, f2, f3
      integer igotdata
      integer itab_print
      integer ix, iy, iz, nx, ny, nz
      integer nxsup, nysup, nzsup
      real*8 xinf, xsup, yinf, ysup, zinf, zsup
      real*8 deltax, deltay, deltaz
      real*8 f1tab(ngrid_mx,ngrid_mx,ngrid_mx)
      real*8 f2tab(ngrid_mx,ngrid_mx,ngrid_mx)
      real*8 f3tab(ngrid_mx,ngrid_mx,ngrid_mx)
      real*8 xi, eta, zeta, px, py, pz, qx, qy, qz
      save
      data igotdata /0/
      data itab_print /0/
      data nxsup, nysup, nzsup / 3 * 0 /
      data xinf, xsup, yinf, ysup, zinf, zsup / 6 * 0.d0 /
      data deltax, deltay, deltaz / 3 * 0.d0 /
      data f1tab / ngrid3_mx * 0.d0 /
      data f2tab / ngrid3_mx * 0.d0 /
      data f3tab / ngrid3_mx * 0.d0 /
! ---------------------------------------------------------------------
      if(init3d .eq. 0) then
        igotdata = 1
        open(unit=21,file=fname)
        read(21,*) nxsup, xinf, xsup
        read(21,*) nysup, yinf, ysup
        read(21,*) nzsup, zinf, zsup
        deltax = (xsup - xinf)/float(nxsup - 1)
        deltay = (ysup - yinf)/float(nysup - 1)
        deltaz = (zsup - zinf)/float(nzsup - 1)
! .....................................................................
        do nx = 1, nxsup
          do ny = 1, nysup
            do nz = 1, nzsup
              read(21,*) ix, iy, iz,  f1tab(ix,iy,iz), f2tab(ix,iy,iz), f3tab(ix,iy,iz)
            enddo
          enddo
        enddo
        close(unit=21)
      else
        xi = (x - xinf)/deltax
        if((x - xinf) .lt. 0.d0) then
          ix = 0
        elseif(((x - xinf) .ge. 0.d0) .and. ((x - xsup) .le. 0.d0)) then
          ix = int(xi)
        else 
          ix = nxsup - 2
        endif
        px = xi - float(ix)
        qx = 1.d0 - px
        ix = ix + 1
        if((ix .lt. 1) .or. (ix .gt. nxsup)) then
          write(*,931) ix
  931     format(//'error: index ix out of range: ',i10)
          stop
        endif
        eta = (y - yinf)/deltay
        if((y - yinf) .lt. 0.d0) then
          iy = 0
        elseif(((y - yinf) .ge. 0.d0) .and. ((y - ysup) .le. 0.d0)) then
          iy = int(eta)
        else 
          iy = nysup - 2
        endif
        py = eta - float(iy)
        qy = 1.d0 - py
        iy = iy + 1
        if((iy .lt. 1) .or. (iy .gt. nysup)) then
          write(*,932) iy
  932     format(//'error: index iy out of range: ',i10)
          stop
        endif
        zeta = (z - zinf)/deltaz
        if((z - zinf) .lt. 0.d0) then
          iz = 0
        elseif(((z - zinf) .ge. 0.d0) .and. ((z - zsup) .le. 0.d0)) then
          iz = int(zeta)
        else 
          iz = nzsup - 2
        endif
        pz = zeta - float(iz)
        qz = 1.d0 - pz
        iz = iz + 1
        if((iz .lt. 1) .or. (iz .gt. nzsup)) then
          write(*,933) iz
  933     format(//'error: index iz out of range: ',i10)
          stop
        endif
        f1 =      qx*qy*qz*f1tab(ix,iy,iz) + qx*qy*pz*f1tab(ix,iy,iz+1)    + qx*py*qz*f1tab(ix,iy+1,iz) + qx*py*pz*f1tab(ix,iy+1,iz+1)   + px*qy*qz*f1tab(ix+1,iy,iz) + px*qy*pz*f1tab(ix+1,iy,iz+1)   + px*py*qz*f1tab(ix+1,iy+1,iz) + px*py*pz*f1tab(ix+1,iy+1,iz+1)
        f2 =      qx*qy*qz*f2tab(ix,iy,iz) + qx*qy*pz*f2tab(ix,iy,iz+1)    + qx*py*qz*f2tab(ix,iy+1,iz) + qx*py*pz*f2tab(ix,iy+1,iz+1)   + px*qy*qz*f2tab(ix+1,iy,iz) + px*qy*pz*f2tab(ix+1,iy,iz+1)   + px*py*qz*f2tab(ix+1,iy+1,iz) + px*py*pz*f2tab(ix+1,iy+1,iz+1)
        f3 =      qx*qy*qz*f3tab(ix,iy,iz) + qx*qy*pz*f3tab(ix,iy,iz+1)    + qx*py*qz*f3tab(ix,iy+1,iz) + qx*py*pz*f3tab(ix,iy+1,iz+1)   + px*qy*qz*f3tab(ix+1,iy,iz) + px*qy*pz*f3tab(ix+1,iy,iz+1)   + px*py*qz*f3tab(ix+1,iy+1,iz) + px*py*pz*f3tab(ix+1,iy+1,iz+1)
      endif
! ---------------------------------------------------------------------
      return
    end subroutine lin3d_3
! =====================================================================


! ============================================================
!   sbrt dip_dip: evaluates dipole-dipole interaction energy
! ------------------------------------------------------------
      subroutine dip_dip(dv1,dv2,rr,vdd)
      implicit none
      real*8 vdd
      real*8 dv1(*), dv2(*), rr(*)
      real*8 d1_d2, d1_r, d2_r, r
! ... scalar product of dipoles
      d1_d2 = dv1(1)*dv2(1) + dv1(2)*dv2(2) + dv1(3)*dv2(3)
! ... scalar products of dipoles and separation vector
      d1_r = dv1(1)*rr(1) + dv1(2)*rr(2) + dv1(3)*rr(3)
      d2_r = dv2(1)*rr(1) + dv2(2)*rr(2) + dv2(3)*rr(3)
! ... distance
      r = sqrt(rr(1)*rr(1) + rr(2)*rr(2) + rr(3)*rr(3))
! ... dipole-dipole interaction energy
      vdd = -(3.d0*d1_r*d2_r/(r*r) - d1_d2)/(r*r*r)
      return
    end subroutine dip_dip
! ============================================================

! ============================================================
!   sbrt dip_quad: evaluates dipole-quadrupole interaction energy
! ------------------------------------------------------------
      subroutine dip_quad(dv,qt,rr,vdq)
      implicit none
      real*8 vdq
      real*8 dv(*), qt(*), rr(*)
      integer i, j, k
      real*8 d_r, r_q_r, r_q_d, r
! ... scalar product of dipole vector and separation vector
      d_r = 0.d0
      do i = 1, 3
        d_r = d_r + dv(i)*rr(i)
      enddo
! ... scalar product separation with quad tensor with separation
      r_q_r = 0.d0
      do i = 1, 3
        do j = 1, 3
          k = 3*(i-1) + j
          r_q_r = r_q_r + rr(i)*qt(k)*rr(j)
        enddo
      enddo
! ... scalar product separation with quad tensor with dipole
      r_q_d = 0.d0
      do i = 1, 3
        do j = 1, 3
          k = 3*(i-1) + j
          r_q_d = r_q_d + rr(i)*qt(k)*dv(j)
        enddo
      enddo
! ... distance
      r = sqrt(rr(1)*rr(1) + rr(2)*rr(2) + rr(3)*rr(3))
! ... dipole-quadrupole interaction energy
      vdq = 5.d0*d_r*r_q_r/(r**7) - 2.d0*r_q_d/(r**5)
! ...
      return
    end subroutine dip_quad
! =============================================================

! ============================================================
!   sbrt quad_dip: evaluates quadrupole-dipole interaction energy
! ------------------------------------------------------------
      subroutine quad_dip(qt,dv,rr,vqd)
      implicit none
      real*8 vqd
      real*8 dv(*), qt(*), rr(*)
      integer i, j, k
      real*8 d_r, r_q_r, r_q_d, r
! ... scalar product of dipole vector and separation vector
      d_r = 0.d0
      do i = 1, 3
        d_r = d_r + dv(i)*rr(i)
      enddo
! ... scalar product separation with quad tensor with separation
      r_q_r = 0.d0
      do i = 1, 3
        do j = 1, 3
          k = 3*(i-1) + j
          r_q_r = r_q_r + rr(i)*qt(k)*rr(j)
        enddo
      enddo
! ... scalar product separation with quad tensor with dipole
      r_q_d = 0.d0
      do i = 1, 3
        do j = 1, 3
          k = 3*(i-1) + j
          r_q_d = r_q_d + rr(i)*qt(k)*dv(j)
        enddo
      enddo
! ... distance
      r = sqrt(rr(1)*rr(1) + rr(2)*rr(2) + rr(3)*rr(3))
! ... quadrupole-dipole interaction energy
      vqd = -5.d0*d_r*r_q_r/(r**7) + 2.d0*r_q_d/(r**5)
! ...
      return
    end subroutine quad_dip
! =============================================================

! =============================================================
!   sbrt quad_quad: evaluates quadrupole-quadrupole interaction energy
! -------------------------------------------------------------
      subroutine quad_quad(qt1,qt2,rr,vqq)
      implicit none
      real*8 vqq
      real*8 qt1(*), qt2(*), rr(*)
      integer i, j, k
      real*8 tr_q_q, r_q1_r, r_q2_r, r_q_q_r, s1, s2, r
! ... trace of qt1*qt2
      tr_q_q = 0.d0
      do i = 1, 3
        do j = 1, 3
          k = 3*(i-1) + j
          tr_q_q = tr_q_q + qt1(k)*qt2(k)
        enddo
      enddo
! ... scalar products r*qt1*r and r*qt2*r
      r_q1_r = 0.d0
      r_q2_r = 0.d0
      do i = 1, 3
        do j = 1, 3
          k = 3*(i-1) + j
          r_q1_r = r_q1_r + rr(i)*qt1(k)*rr(j)
          r_q2_r = r_q2_r + rr(i)*qt2(k)*rr(j)
        enddo
      enddo
! ... scalar product r*qt1*qt2*r
      r_q_q_r = 0.d0
      do i = 1, 3
        s1 = 0.d0
        s2 = 0.d0
        do j = 1, 3
          k = 3*(i-1) + j
          s1 = s1 + qt1(k)*rr(j)
          s2 = s2 + qt2(k)*rr(j)
        enddo
        r_q_q_r = r_q_q_r + s1*s2
      enddo
! ... distance
      r = sqrt(rr(1)*rr(1) + rr(2)*rr(2) + rr(3)*rr(3))
! ... quad-quad interaction energy
      vqq = 2.d0*tr_q_q/(3.d0*(r**5)) - 20.d0*r_q_q_r/(3.d0*(r**7)) + 35.d0*r_q1_r*r_q2_r/(3.d0*(r**9))
! ...
      return
    end subroutine quad_quad
! ===============================================================

! ===============================================================
!   sbrt dip_field: computes electric field vector evec2 at vector
!   position pos2 due to electric dipole vector mu1 at origin
! ---------------------------------------------------------------
      subroutine dip_field(pos2,mu1,evec2)
      implicit none
      real*8 pos2(*), mu1(*), evec2(*)
      integer n
      real*8 drr, mudotr, drr3, drr5
! ---------------------------------------------------------------
      drr = 0.d0
      mudotr = 0.d0
      do n = 1, 3
        drr = drr + pos2(n)*pos2(n)
        mudotr = mudotr + mu1(n)*pos2(n)
      enddo
      drr = sqrt(drr)
      drr3 = drr*drr*drr
      drr5 = drr*drr*drr3
      do n = 1, 3
        evec2(n) = 3.d0*mudotr*pos2(n)/drr5 - mu1(n)/drr3
      enddo
! ---------------------------------------------------------------
      return
    end subroutine dip_field
! ===============================================================

! ===============================================================
!   sbrt quad_field: computes electric field evec2 at vector 
!   position pos2 due to electric quadrupole tensor q1 at origin
! ---------------------------------------------------------------
      subroutine quad_field(pos2,q1,evec2)
      implicit none
      real*8 pos2(*), q1(*), evec2(*)
      integer k, m, n
      real*8 drr, drr5, drr7, rdqdr
      real*8 qdotr(3)
! ---------------------------------------------------------------
      drr = 0.d0
      rdqdr = 0.d0
      do m = 1, 3
        drr = drr + pos2(m)*pos2(m)
        qdotr(m) = 0.d0
        do n = 1, 3
          k = 3*(m-1) + n
          rdqdr = rdqdr + pos2(m)*q1(k)*pos2(n)
          qdotr(m) = qdotr(m) + q1(k)*pos2(n)
        enddo
      enddo
      drr = sqrt(drr)
      drr5 = drr*drr*drr*drr*drr
      drr7 = drr*drr*drr5
      do n = 1, 3
        evec2(n) = 5.d0*rdqdr*pos2(n)/drr7 - 2.d0*qdotr(n)/drr5
      enddo
! ----------------------------------------------------------------
      return
    end subroutine quad_field
! ================================================================



end module IPModel_WaterDimer_Gillan_module
