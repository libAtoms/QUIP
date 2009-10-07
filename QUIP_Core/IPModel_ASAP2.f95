!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X     QUIP: quantum mechanical and interatomic potential simulation package
!X     
!X     Portions written by Noam Bernstein, while working at the
!X     Naval Research Laboratory, Washington DC. 
!X
!X     Portions written by Gabor Csanyi, Copyright 2006-2007.   
!X
!X     When using this software,  please cite the following reference:
!X
!X     reference
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X IPModel_ASAP2
!X
!% Reimplementation of ASAP potential:
!% P. Tangney and S. Scandolo,
!% An ab initio parametrized interatomic force field for silica
!% J. Chem. Phys, 117, 8898 (2002). 
!%
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
module IPModel_ASAP2_module

use libatoms_module

use mpi_context_module
use QUIP_Common_module

implicit none
private

include 'IPModel_interface.h'

logical, private :: asap_initialised = .false.

public :: IPModel_ASAP2
type IPModel_ASAP2
  integer :: n_types = 0
  real(dp) :: betapol, tolpol, yukalpha, yuksmoothlength
  integer :: maxipol, pred_order
  integer, allocatable :: atomic_num(:), type_of_atomic_num(:)
  real(dp), allocatable, dimension(:) :: pol, z
  real(dp), allocatable, dimension(:,:) :: D_ms, gamma_ms, R_ms, B_pol, C_pol

  real(dp) :: cutoff_coulomb, cutoff_ms

  character(len=FIELD_LENGTH) :: label
  type(mpi_context) :: mpi
  logical :: initialised, tdip_sr

end type IPModel_ASAP2

logical :: parse_in_ip, parse_matched_label
type(IPModel_ASAP2), pointer :: parse_ip

interface Initialise
  module procedure IPModel_ASAP2_Initialise_str
end interface Initialise

interface Finalise
  module procedure IPModel_ASAP2_Finalise
end interface Finalise

interface Print
  module procedure IPModel_ASAP2_Print
end interface Print

interface Calc
  module procedure IPModel_ASAP2_Calc
end interface Calc

contains


subroutine IPModel_ASAP2_Initialise_str(this, args_str, param_str, mpi)
  type(IPModel_ASAP2), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str
  type(mpi_context), intent(in), optional :: mpi

  type(Dictionary) :: params

  this%initialised = .false.

  call Finalise(this)

  call initialise(params)
  this%label=''
  call param_register(params, 'label', '', this%label)
  if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_ASAP2_Initialise_str args_str')) then
    call system_abort("IPModel_ASAP2_Initialise_str failed to parse label from args_str="//trim(args_str))
  endif
  call finalise(params)

  call IPModel_ASAP2_read_params_xml(this, param_str)
  this%initialised = .true.
  
  if (present(mpi)) this%mpi = mpi

end subroutine IPModel_ASAP2_Initialise_str

subroutine IPModel_ASAP2_Finalise(this)
  type(IPModel_ASAP2), intent(inout) :: this

  this%initialised = .false.

  if (allocated(this%atomic_num)) deallocate(this%atomic_num)
  if (allocated(this%type_of_atomic_num)) deallocate(this%type_of_atomic_num)
  if (allocated(this%pol)) deallocate(this%pol)
  if (allocated(this%z)) deallocate(this%z)
  if (allocated(this%D_ms)) deallocate(this%D_ms)
  if (allocated(this%gamma_ms)) deallocate(this%gamma_ms)
  if (allocated(this%R_ms)) deallocate(this%R_ms)
  if (allocated(this%B_pol)) deallocate(this%B_pol)
  if (allocated(this%C_pol)) deallocate(this%C_pol)
  this%n_types = 0
  this%label = ''
end subroutine IPModel_ASAP2_Finalise

!% Smooth cutoff function from D.J. Cole \emph{et al.}, J. Chem. Phys. {\bf 127}, 204704 (2007).
subroutine smooth_cutoff(x,R,D,fc,dfc_dx)

  real(dp) x,R,D,fc,dfc_dx
  real(dp), parameter :: pi = dacos(-1.0d0)

  if (x .lt. (R-D)) then
     fc = 1.0d0
     dfc_dx = 0.0d0
     return
  else if (x .gt. (R+D)) then
     fc = 0.0d0
     dfc_dx = 0.0d0
     return
  else
     fc = 1.0d0 - (x-R+D)/(2.0d0*D) + 1.0d0/(2.0d0*pi)*dsin((pi/D)*(x-R+D))
     dfc_dx = 1.0d0/(2.0d0*D)* (dcos((pi/D)*(x-R+D)) - 1.0d0)
     return
  end if
end subroutine smooth_cutoff


!% Charge-charge interactions, screened by Yukawa function
subroutine asap_rs_charges(this, at, e, local_e, f, virial, efield)
   type(IPModel_ASAP2), intent(inout):: this
   type(Atoms), intent(inout)      :: at
   real(dp), intent(out), optional :: e, local_e(:)
   real(dp), intent(out), optional :: f(:,:)
   real(dp), intent(out), optional :: virial(3,3)
   real(dp), intent(out), optional :: efield(:,:)

   integer i, j, m, ti, tj
   real(dp) :: r_ij, u_ij(3), zv2, gamjir, gamjir3, gamjir2, fc, dfc_dr
   real(dp) :: de, dforce, expfactor
   logical :: i_is_min_image

   call system_timer('asap_rs_charges')

   do i=1, at%n
      if (allocated(at%connect%is_min_image)) then
         i_is_min_image = at%connect%is_min_image(i)
      else
         i_is_min_image = is_min_image(at, i)
      end if
      ti = get_type(this%type_of_atomic_num, at%Z(i))
      do m = 1, atoms_n_neighbours(at, i)
         
         j = atoms_neighbour(at, i, m, distance=r_ij, cosines=u_ij, max_dist=(this%cutoff_coulomb*BOHR))
         if (j <= 0) cycle
         if (r_ij .feq. 0.0_dp) cycle
         if (i < j .and. i_is_min_image) cycle

         r_ij = r_ij/BOHR
         tj = get_type(this%type_of_atomic_num, at%Z(j))
         zv2 = this%z(ti)*this%z(tj)

         gamjir = zv2/r_ij
         gamjir3 = gamjir/(r_ij**2.0_dp)
         gamjir2 = zv2/(r_ij**2.0_dp)
         expfactor = exp(-this%yukalpha*r_ij)

         call smooth_cutoff(r_ij, this%cutoff_coulomb-this%yuksmoothlength, &
              this%yuksmoothlength, fc, dfc_dr)

         de = gamjir

         if (present(e) .or. present(local_e)) then
            if (present(e)) then
               if (i_is_min_image) then
                  e = e + de*expfactor*fc
               else
                  e = e + 0.5_dp*de*expfactor*fc
               end if
            end if
            if (present(local_e)) then
               local_e(i) = local_e(i) + 0.5_dp*de*expfactor*fc
               if (i_is_min_image) local_e(j) = local_e(j) + 0.5_dp*de*expfactor*fc
            end if
         end if

         if (present(f) .or. present(virial) .or. present(efield)) then
            dforce = gamjir3*expfactor*fc*r_ij + de*(this%yukalpha*fc - dfc_dr)*expfactor

            if (present(f)) then
               f(:,i) = f(:,i) - dforce*u_ij
               if (i_is_min_image) f(:,j) = f(:,j) + dforce*u_ij
            end if

            if (present(virial)) then
               if (i_is_min_image) then
                  virial = virial + dforce*(u_ij .outer. u_ij)*r_ij
               else
                  virial = virial + 0.5_dp*dforce*(u_ij .outer. u_ij)*r_ij
               end if
            end if

            if (present(efield)) then
               efield(:,i) = efield(:,i) - gamjir3*expfactor*fc*u_ij*r_ij/this%z(ti)
               if (i_is_min_image) efield(:,j) = efield(:,j) + gamjir3*expfactor*fc*u_ij*r_ij/this%z(tj)
            end if
         end if
      end do
   end do

   call system_timer('asap_rs_charges')

end subroutine asap_rs_charges


!% Charge-dipole and dipole-dipole interactions, screened by Yukawa function
subroutine asap_rs_dipoles(this, at, dip, e, local_e, f, virial, efield)
   type(IPModel_ASAP2), intent(inout):: this
   type(Atoms), intent(inout)      :: at
   real(dp), intent(in)            :: dip(:,:)
   real(dp), intent(out), optional :: e, local_e(:)
   real(dp), intent(out), optional :: f(:,:)
   real(dp), intent(out), optional :: virial(3,3)
   real(dp), intent(out), optional :: efield(:,:)

   integer i, j, m, ti, tj, k
   real(dp) :: r_ij, u_ij(3), gamjir3, gamjir2, fc, dfc_dr
   real(dp) :: expfactor, dipi(3), dipj(3), qj, qi, pp, pri, prj
   real(dp) :: de_ind, de_dd, de_qd, dfqdip(3), dfdipdip(3), factor1, dist3, dist5
   real(dp) :: const1, const2, factork, de_sr, df_sr(3), gij, dgijdrij, bij, cij
   logical :: i_is_min_image, tpoli, tpolj, qipj, qjpi, pipj

   call system_timer('asap_rs_dipoles')

   do i=1, at%n
      if (allocated(at%connect%is_min_image)) then
         i_is_min_image = at%connect%is_min_image(i)
      else
         i_is_min_image = is_min_image(at, i)
      end if
      ti = get_type(this%type_of_atomic_num, at%Z(i))

      qi = this%z(ti)
      dipi = dip(:,i)
      tpoli = abs(this%pol(ti)) > 0.0_dp

      ! Induced contribution to energy
      if ((present(e) .or. present(local_e)) .and. tpoli) then
         de_ind = 0.5_dp*(dipi .dot. dipi)/this%pol(ti)
         if (present(e))       e = e + de_ind
         if (present(local_e)) local_e(i) = local_e(i) + de_ind
      end if

      do m = 1, atoms_n_neighbours(at, i)
         
         j = atoms_neighbour(at, i, m, distance=r_ij, diff=u_ij, max_dist=(this%cutoff_coulomb*BOHR))
         if (j <= 0) cycle
         if (r_ij .feq. 0.0_dp) cycle
         if (i < j .and. i_is_min_image) cycle

         r_ij = r_ij/BOHR
         u_ij = u_ij/BOHR
         tj = get_type(this%type_of_atomic_num, at%Z(j))

         qj = this%z(tj)
         dipj = dip(:,j)
         tpolj = abs(this%pol(tj)) > 0.0_dp

         qipj = (abs(qi) > 0.0_dp) .and. tpolj
         qjpi = (abs(qj) > 0.0_dp) .and. tpoli
         pipj = tpoli .and. tpolj

         if (.not. (pipj .or. qipj .or. qjpi)) cycle
        
         gamjir2 = 1.0_dp/r_ij**2.0_dp
         gamjir3 = gamjir2/r_ij
         factor1 = 3.0_dp*gamjir2*gamjir3

         expfactor = exp(-this%yukalpha*r_ij)

         call smooth_cutoff(r_ij, this%cutoff_coulomb-this%yuksmoothlength, &
              this%yuksmoothlength, fc, dfc_dr)

         pp = 0.0_dp
         if (pipj)  pp  = dipi .dot. dipj
         pri = 0.0_dp
         if (tpoli) pri = dipi .dot. u_ij

         prj = 0.0_dp
         if (tpolj) prj = dipj .dot. u_ij

         if (present(efield)) then
            efield(:,i) = efield(:,i) + &
                 (3.0_dp*prj*u_ij*gamjir2 - dipj)*gamjir2*dsqrt(gamjir2)*expfactor*fc
            if (i_is_min_image) &
                 efield(:,j) = efield(:,j) + &
                 (3.0_dp*pri*u_ij*gamjir2 - dipi)*gamjir2*dsqrt(gamjir2)*expfactor*fc
         end if

         if (present(e) .or. present(local_e) .or. present(virial) .or. present(f)) then

            de_dd = (pp - 3.0_dp*pri*prj*gamjir2)*gamjir3
            de_qd = -(qi*prj - qj*pri)*gamjir3

            de_sr = 0.0_dp
            df_sr = 0.0_dp
            if (this%tdip_sr) then

               bij = this%b_pol(ti, tj)
               cij = this%c_pol(ti, tj)
         
               dist3 = 1.0_dp/(r_ij**3.0_dp)
               dist5 = 1.0_dp/(r_ij**5.0_dp)

               gij = 0.0_dp
               factork = cij*exp(-bij*r_ij)
               do k=1,4
                  gij = gij + factork 
                  factork = factork*bij*r_ij/float(k)
               enddo
               gij = gij + factork
               dgijdrij = -bij*factork

               const1 = gij*dist3
               const2 = (qj*pri-qi*prj) &
                    * (r_ij*dgijdrij-3.0_dp*gij)*dist5

               de_sr = -(qi*prj - qj*pri)*gij*dist3

               df_sr = ((qj*dipi-qi*dipj)*const1 + u_ij*const2) &
                    *expfactor*fc  &
                    -de_sr*(this%yukalpha*fc - dfc_dr)*expfactor/r_ij*u_ij
            end if

            if (present(e)) then
               if (i_is_min_image) then
                  e = e + (de_dd + de_qd + de_sr)*expfactor*fc
               else
                  e = e + 0.5_dp*(de_dd + de_qd + de_sr)*expfactor*fc
               end if
            end if

            if (present(local_e)) then
               local_e(i) = local_e(i) + 0.5_dp*(de_dd + de_qd + de_sr)*expfactor*fc
               if (i_is_min_image) local_e(j) = local_e(j) + 0.5_dp*(de_dd + de_qd + de_sr)*expfactor*fc
            end if

            if (present(f) .or. present(virial)) then
               dfqdip = .0_dp
               if (qipj .or. qjpi) then
                  dfqdip = ((qj*dipi(:) - qi*dipj(:) &
                       - 3.0_dp*(qj*pri-qi*prj)*u_ij*gamjir2)*gamjir3)*expfactor*fc &
                       - de_qd*(this%yukalpha*fc - dfc_dr)*expfactor/r_ij*u_ij
               end if

               dfdipdip = 0.0_dp
               if (pipj) then
                  dfdipdip = (-(pp*u_ij + dipj(:)*pri+dipi(:)*prj &
                       - 5.0_dp*prj*pri*u_ij(:)*gamjir2)*factor1)*expfactor*fc &
                       - de_dd*(this%yukalpha*fc - dfc_dr)*expfactor/r_ij*u_ij(:)
               end if

               if (present(f)) then
                  f(:,i) = f(:,i) + dfqdip + dfdipdip + df_sr
                  if (i_is_min_image) f(:,j) = f(:,j) - (dfqdip + dfdipdip + df_sr)
               end if

               if (present(virial)) then
                  if (i_is_min_image) then
                     virial = virial - ((dfqdip+dfdipdip+df_sr) .outer. u_ij)
                  else
                     virial = virial - 0.5_dp*((dfqdip+dfdipdip+df_sr) .outer. u_ij)
                  end if
               end if

            end if

         end if

      end do
   end do

   call system_timer('asap_rs_dipoles')

end subroutine asap_rs_dipoles


subroutine asap_short_range_dipole_moments(this, at, dip_sr)
  type(IPModel_ASAP2), intent(inout) :: this
  type(Atoms), intent(inout) :: at
  real(dp), dimension(:,:), intent(out) :: dip_sr

  integer :: i, ti, m, j, tj, k
  real(dp) :: r_ij, u_ij(3), qj, bij, cij, dist3, dist5, gij, factork, expfactor, fc, dfc_dr
  integer, parameter :: nk = 4

  call system_timer('asap_short_range_dipole_moments')

  dip_sr = 0.0_dp
  do i=1, at%n
     ti = get_type(this%type_of_atomic_num, at%Z(i))
     if (.not. (abs(this%pol(ti)) > 0.0_dp)) cycle

     do m = 1, atoms_n_neighbours(at, i)
         
         j = atoms_neighbour(at, i, m, distance=r_ij, diff=u_ij, max_dist=(this%cutoff_ms*BOHR))
         if (j <= 0) cycle
         if (r_ij .feq. 0.0_dp) cycle

         r_ij = r_ij/BOHR
         u_ij = u_ij/BOHR

         tj = get_type(this%type_of_atomic_num, at%Z(j))
         qj = this%z(tj)
         bij = this%b_pol(ti, tj)
         cij = this%c_pol(ti, tj)

         dist3 = r_ij**3.0_dp
         dist5 = r_ij**5.0_dp
         
         gij = 0.0_dp
         factork = cij*exp(-bij*r_ij)
         do k=1,nk
            gij = gij + factork 
            factork = factork*bij*r_ij/real(k,dp)
         enddo
         gij = gij + factork
         
         expfactor = exp(-this%yukalpha*r_ij)
         call smooth_cutoff(r_ij, this%cutoff_coulomb-this%yuksmoothlength, this%yuksmoothlength, fc, dfc_dr)
         dip_sr(:,i) = dip_sr(:,i) - this%pol(ti)*qj*u_ij*gij/dist3*expfactor*fc
  
      end do
  end do

  call system_timer('asap_short_range_dipole_moments')

end subroutine asap_short_range_dipole_moments

!% Morse-stretch potential, defined by
!% \begin{displaymath}
!% U_ij = D_ij \left[ e^{\gamma_{ij}\left( 1 - r_{ij}/r^0_{ij} \right)} 
!%                  - 2 e^{\left( 1 - r_{ij}/r^0_{ij} \right)} \right]
!% \end{displaymath}
subroutine asap_morse_stretch(this, at, e, local_e, f, virial)
   type(IPModel_ASAP2), intent(inout):: this
   type(Atoms), intent(inout)      :: at
   real(dp), intent(out), optional :: e, local_e(:)
   real(dp), intent(out), optional :: f(:,:)
   real(dp), intent(out), optional :: virial(3,3)

   integer i, j, m, ti, tj
   real(dp) :: r_ij, u_ij(3), dms, gammams, rms
   real(dp) :: exponentms, factorms, phi, de
   real(dp) :: dforce
   real(dp) :: elimitij(this%n_types, this%n_types)
   logical :: i_is_min_image

   call system_timer('asap_morse_stretch')

   ! Evaluate potential at cutoff. Will be subtracted from total energy.
   elimitij = 0.0_dp
   do ti=1,this%n_types
      do tj=1,this%n_types
         phi = exp(this%gamma_ms(ti,tj)*(1.0_dp - this%cutoff_ms/this%r_ms(ti,tj)))
         elimitij(ti,tj) = this%d_ms(ti,tj)*(phi - 2.0_dp*sqrt(phi))
      end do
   end do
 
   do i=1, at%n
      if (allocated(at%connect%is_min_image)) then
         i_is_min_image = at%connect%is_min_image(i)
      else
         i_is_min_image = is_min_image(at, i)
      end if
      ti = get_type(this%type_of_atomic_num, at%Z(i))
      do m = 1, atoms_n_neighbours(at, i)
         
         j = atoms_neighbour(at, i, m, distance=r_ij, cosines=u_ij, max_dist=(this%cutoff_ms*BOHR))
         if (j <= 0) cycle
         if (r_ij .feq. 0.0_dp) cycle
         if (i < j .and. i_is_min_image) cycle

         r_ij = r_ij/BOHR
         tj = get_type(this%type_of_atomic_num, at%Z(j))

         dms = this%d_ms(ti,tj)
         gammams = this%gamma_ms(ti,tj)
         rms = this%r_ms(ti, tj)

         exponentms = gammams*(1.0_dp-r_ij/rms)
         factorms = gammams/rms
         phi = exp(exponentms)

         if (present(e) .or. present(local_e)) then
            de = dms*(phi-2.0_dp*sqrt(phi)) - elimitij(ti, tj)

            if (present(e)) then
               if (i_is_min_image) then
                  e = e + de
               else
                  e = e + 0.5_dp*de
               end if
            end if
            if (present(local_e)) then
               local_e(i) = local_e(i) + 0.5_dp*de
               if (i_is_min_image) local_e(j) = local_e(j) + 0.5_dp*de
            end if
         end if

         if (present(f) .or. present(virial)) then
            dforce  = -dms*(factorms*phi - factorms*dsqrt(phi))

            if (present(f)) then
               f(:,i) = f(:,i) + dforce*u_ij
               if (i_is_min_image) f(:,j) = f(:,j) - dforce*u_ij
            end if
            if (present(virial)) then
               if (i_is_min_image) then
                  virial = virial - dforce*(u_ij .outer. u_ij)*r_ij
               else
                  virial = virial - 0.5_dp*dforce*(u_ij .outer. u_ij)*r_ij
               end if
            end if
         end if
      end do
   end do

   call system_timer('asap_morse_stretch')

end subroutine asap_morse_stretch

subroutine IPModel_ASAP2_prepare_atoms(this, at)
   type(IPModel_ASAP2), intent(inout):: this
   type(Atoms), intent(inout)      :: at

   ! Add properties to Atoms object - it's not safe to do this from within calc()
   if (.not. has_property(at, 'efield')) call add_property(at, 'efield', 0.0_dp, n_cols=3)
   if (.not. has_property(at, 'dipoles')) call add_property(at, 'dipoles', 0.0_dp, n_cols=3)
   if (.not. has_property(at, 'efield_old1')) call add_property(at, 'efield_old1', 0.0_dp, n_cols=3)
   if (.not. has_property(at, 'efield_old2')) call add_property(at, 'efield_old2', 0.0_dp, n_cols=3)
   if (.not. has_property(at, 'efield_old3')) call add_property(at, 'efield_old3', 0.0_dp, n_cols=3)

end subroutine IPModel_ASAP2_prepare_atoms


subroutine IPModel_ASAP2_Calc(this, at, e, local_e, f, virial, args_str)
   type(IPModel_ASAP2), intent(inout):: this
   type(Atoms), intent(inout)      :: at
   real(dp), intent(out), optional :: e, local_e(:)
   real(dp), intent(out), optional :: f(:,:)
   real(dp), intent(out), optional :: virial(3,3)
   character(len=*), optional, intent(in) :: args_str

   type(Dictionary) :: params
   logical :: save_efield, save_dipoles, restart
   real(dp), allocatable, target :: theefield(:,:), thedipoles(:,:) 
   real(dp), allocatable :: efield_charge(:,:), efield_dipole(:,:), dip_sr(:,:)
   real(dp), pointer :: efield(:,:), dipoles(:,:)
   real(dp), pointer, dimension(:,:) :: efield_old1, efield_old2, efield_old3
   real(dp) :: diff, diff_old
   integer :: n_efield_old
   integer :: i, npol, ti, vv

   real, parameter :: difftol = 500.0_dp

   call system_timer('asap_calc')
   vv = current_verbosity()

   call initialise(params)
   call param_register(params, 'save_efield', 'T', save_efield)
   call param_register(params, 'save_dipoles', 'T', save_dipoles)
   call param_register(params, 'restart', 'F', restart)
   if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_ASAP2_Calc args_str')) then
      call system_abort("IPModel_ASAP2_Calc failed to parse args_str="//trim(args_str))
   endif
   call finalise(params)

   if (present(e)) e = 0.0_dp
   if (present(f)) f = 0.0_dp
   if (present(local_e)) local_e = 0.0_dp
   if (present(virial)) virial = 0.0_dp

   allocate(efield_charge(3,at%n), efield_dipole(3,at%n))

   ! Assign pointers
   if (save_efield) then
      if (.not. assign_pointer(at, 'efield', efield)) &
           call system_abort('IPModel_ASAP2_calc failed to assign pointer to "efield" property')
   else
      allocate(theefield(3,at%n))
      efield => theefield
   end if

   if (save_dipoles) then
      if (.not. assign_pointer(at, 'dipoles', dipoles)) &
           call system_abort('IPModel_ASAP2_calc failed to assign pointer to "dipoles" property')
   else
      allocate(thedipoles(3,at%n))
      dipoles => thedipoles
   end if

   if (.not. assign_pointer(at, 'efield_old1', efield_old1)) &
        call system_abort('IPModel_ASAP2_calc failed to assign pointer to "efield_old1" property')
   if (.not. assign_pointer(at, 'efield_old2', efield_old2)) &
        call system_abort('IPModel_ASAP2_calc failed to assign pointer to "efield_old2" property')
   if (.not. assign_pointer(at, 'efield_old3', efield_old3)) &
        call system_abort('IPModel_ASAP2_calc failed to assign pointer to "efield_old3" property')

   if (.not. get_value(at%params, 'n_efield_old', n_efield_old)) n_efield_old = 0
      

   if (restart) then
      efield_old1(:,:) = 0.0_dp
      efield_old2(:,:) = 0.0_dp
      efield_old3(:,:) = 0.0_dp
   end if

   efield = 0.0_dp
   dipoles = 0.0_dp

   ! Extrapolate from the old electric efield
   if (this%pred_order == 0) then
      if (n_efield_old >= 1) efield = efield_old1(:,:)

   else if (this%pred_order == 1) then

      if (n_efield_old >= 2) &
           efield = 2.0_dp*efield_old1(:,:) - efield_old2(:,:)

      efield_old2(:,:) = efield_old1(:,:)
   else if (this%pred_order == 2) then
        
      if (n_efield_old >= 3) &
           efield = 3.0_dp*efield_old1(:,:) - 3.0_dp*efield_old2(:,:) + efield_old3(:,:)

      efield_old3(:,:) = efield_old2(:,:)
      efield_old2(:,:) = efield_old1(:,:)
   end if

   call set_value(at%params, 'n_efield_old', min(n_efield_old+1,3))

   efield_dipole = 0.0_dp
   efield_charge = 0.0_dp

   if (maxval(abs(this%z)) > 0.0_dp) then
      call asap_rs_charges(this, at, e, local_e, f, virial, efield_charge)
   end if

   allocate(dip_sr(3,at%n))
   dip_sr = 0.0_dp
   if (this%tdip_sr) call asap_short_range_dipole_moments(this, at, dip_sr)

   if (maxval(abs(this%pol)) > 0.0_dp) then
      ! Self-consistent determination of dipole moments
      diff_old = 1.0_dp
      npol = 1
      call system_timer('asap_self_consistent_dipoles')
      do 
         ! Mix current and previous total efields
         if (npol == 1) then
            efield = efield_dipole + efield_charge
         else
            efield = this%betapol*efield_dipole + &
                 (1.0_dp - this%betapol)*efield_old1(:,:) + efield_charge
         end if

         ! Calculate dipole moment in response to total efield
         do i=1,at%n
            ti = get_type(this%type_of_atomic_num, at%Z(i))
            dipoles(:,i) = efield(:,i)*this%pol(ti) + dip_sr(:,i)
         end do

         ! Calculate new efield and measure of convergence
         efield_old1(:,:) = efield_dipole
         efield_dipole = 0.0_dp
         call asap_rs_dipoles(this, at, dipoles, efield=efield_dipole)

         diff = 0.0_dp
         do i=1,at%n
            ti = get_type(this%type_of_atomic_num, at%Z(i))
            diff = diff + ((efield_dipole(:,i) - efield_old1(:,i)) .dot. (efield_dipole(:,i) - efield_old1(:,i)))* &
                 this%pol(ti)*this%pol(ti)
         end do
         diff = sqrt(diff/at%n)
         
         if (vv >= NORMAL) then
            write (line,'("Polarisation iteration : ",i5,3e16.8)') npol, diff_old, diff
            call print(line, NORMAL)
         end if

         if (diff > difftol) &
              call system_abort('IPModel_ASAP2_calc: Polarisation diverges - diff='//diff)
      
         if (abs(diff - diff_old) < this%tolpol) exit

         diff_old = diff
         npol = npol + 1
         if (npol >= this%maxipol) &
              call system_abort('IPModel_ASAP2_calc: Polarisation not converged in '//this%maxipol//' steps - diff='//diff)
      end do
      call system_timer('asap_self_consistent_dipoles')

      ! Compute final energy, local energies, forces, virial and electric efield
      efield = efield_charge
      call asap_rs_dipoles(this, at, dipoles, e, local_e, f, virial, efield)
   end if

   ! Finally, add the short-range contribution
   call asap_morse_stretch(this, at, e, local_e, f, virial)

      ! Unit conversion
   if (present(e)) e = e*HARTREE
   if (present(local_e)) local_e = local_e*HARTREE
   if (present(f)) f = f*(HARTREE/BOHR)
   if (present(virial)) virial = virial*HARTREE

   if (allocated(theefield)) deallocate(theefield)
   if (allocated(thedipoles)) deallocate(thedipoles)
   if (allocated(efield_charge)) deallocate(efield_charge)
   if (allocated(efield_dipole)) deallocate(efield_dipole)
   if (allocated(dip_sr)) deallocate(dip_sr)

   call system_timer('asap_calc')

end subroutine IPModel_ASAP2_Calc


subroutine IPModel_ASAP2_Print(this, file)
  type(IPModel_ASAP2), intent(in) :: this
  type(Inoutput), intent(inout),optional :: file

  integer :: ti, tj

  call Print("IPModel_ASAP2 : ASAP2 Potential", file=file)
  call Print("IPModel_ASAP2 : n_types = " // this%n_types, file=file)
  call Print("IPModel_ASAP2 : betapol = "//this%betapol//" maxipol = "//this%maxipol//" tolpol = "//this%tolpol//" pred_order = "//this%pred_order, file=file)
  call Print("IPModel_ASAP2 : yukalpha = "//this%yukalpha//" yuksmoothlength = "//this%yuksmoothlength, file=file)

  do ti=1, this%n_types
    call Print ("IPModel_ASAP2 : type " // ti // " atomic_num " // this%atomic_num(ti), file=file)
    call Print ("IPModel_ASAP2 : pol = "//this%pol(ti), file=file)
    call Print ("IPModel_ASAP2 : z   = "//this%z(ti), file=file)
   call verbosity_push_decrement()
    do tj =1,this%n_types
       call Print ("IPModel_ASAP2 : pair interaction ti tj " // ti // " " // tj // " Zi Zj " // this%atomic_num(ti) //&
            " " // this%atomic_num(tj), file=file)
       call Print ("IPModel_ASAP2 : pair " // this%D_ms(ti,tj) // " " // this%gamma_ms(ti,tj) // " " &
            // this%R_ms(ti,tj) // " " // this%B_pol(ti,tj) // " " // this%C_pol(ti, tj), file=file)
    end do
   call verbosity_pop()
  end do

end subroutine IPModel_ASAP2_Print

subroutine IPModel_ASAP2_read_params_xml(this, param_str)
  type(IPModel_ASAP2), intent(inout), target :: this
  character(len=*), intent(in) :: param_str

  type(xml_t) :: fxml

  if (len(trim(param_str)) <= 0) return

  parse_in_ip = .false. 
  parse_matched_label = .false.
  parse_ip => this

  call open_xml_string(fxml, param_str)
  call parse(fxml,  &
    startElement_handler = IPModel_startElement_handler, &
    endElement_handler = IPModel_endElement_handler)
  call close_xml_t(fxml)

  if (this%n_types == 0) then
    call system_abort("IPModel_ASAP2_read_params_xml parsed file, but n_types = 0")
  endif

end subroutine IPModel_ASAP2_read_params_xml

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!% XML param reader functions
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
subroutine IPModel_startElement_handler(URI, localname, name, attributes)
  character(len=*), intent(in)   :: URI
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name
  type(dictionary_t), intent(in) :: attributes

  integer :: status
  character(len=FIELD_LENGTH) :: value

  integer ti, tj, Zi, Zj

  if (name == 'ASAP_params') then ! new ASAP stanza

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

      call QUIP_FoX_get_value(attributes, 'n_types', value, status)
      if (status == 0) then
        read (value, *), parse_ip%n_types
      else
        call system_abort("Can't find n_types in ASAP_params")
      endif

      allocate(parse_ip%atomic_num(parse_ip%n_types))
      parse_ip%atomic_num = 0

      allocate(parse_ip%pol(parse_ip%n_types))
      parse_ip%pol = 0.0_dp
      allocate(parse_ip%z(parse_ip%n_types))

      allocate(parse_ip%D_ms(parse_ip%n_types,parse_ip%n_types))
      parse_ip%D_ms = 0.0_dp
      allocate(parse_ip%gamma_ms(parse_ip%n_types,parse_ip%n_types))
      parse_ip%gamma_ms = 0.0_dp
      allocate(parse_ip%R_ms(parse_ip%n_types,parse_ip%n_types))
      parse_ip%R_ms = 0.0_dp
      allocate(parse_ip%B_pol(parse_ip%n_types,parse_ip%n_types))
      parse_ip%B_pol = 0.0_dp
      allocate(parse_ip%C_pol(parse_ip%n_types,parse_ip%n_types))
      parse_ip%C_pol = 0.0_dp

      call QUIP_FoX_get_value(attributes, "cutoff_coulomb", value, status)
      if (status /= 0) call system_abort ("IPModel_ASAP2_read_params_xml cannot find cutoff_coulomb")
      read (value, *) parse_ip%cutoff_coulomb

      call QUIP_FoX_get_value(attributes, "cutoff_ms", value, status)
      if (status /= 0) call system_abort ("IPModel_ASAP2_read_params_xml cannot find cutoff_ms")
      read (value, *) parse_ip%cutoff_ms

      call QUIP_FoX_get_value(attributes, "betapol", value, status)
      if (status == 0) read (value, *) parse_ip%betapol

      call QUIP_FoX_get_value(attributes, "maxipol", value, status)
      if (status == 0) read (value, *) parse_ip%maxipol

      call QUIP_FoX_get_value(attributes, "tolpol", value, status)
      if (status == 0) read (value, *) parse_ip%tolpol

      call QUIP_FoX_get_value(attributes, "pred_order", value, status)
      if (status == 0) read (value, *) parse_ip%pred_order

      call QUIP_FoX_get_value(attributes, "yukalpha", value, status)
      if (status == 0) read (value, *) parse_ip%yukalpha

      call QUIP_FoX_get_value(attributes, "yuksmoothlength", value, status)
      if (status == 0) read (value, *) parse_ip%yuksmoothlength

      parse_ip%tdip_sr = .true.
      call QUIP_FoX_get_value(attributes, "tdip_sr", value, status)
      if (status == 0) read (value, *), parse_ip%tdip_sr

    endif

  elseif (parse_in_ip .and. name == 'per_type_data') then

    call QUIP_FoX_get_value(attributes, "type", value, status)
    if (status /= 0) call system_abort ("IPModel_ASAP2_read_params_xml cannot find type")
    read (value, *) ti

    call QUIP_FoX_get_value(attributes, "atomic_num", value, status)
    if (status /= 0) call system_abort ("IPModel_ASAP2_read_params_xml cannot find atomic_num")
    read (value, *) parse_ip%atomic_num(ti)

    call QUIP_FoX_get_value(attributes, "pol", value, status)
    if (status /= 0) call system_abort ("IPModel_ASAP2_read_params_xml cannot find pol")
    read (value, *) parse_ip%pol(ti)

    call QUIP_FoX_get_value(attributes, "z", value, status)
    if (status /= 0) call system_abort ("IPModel_ASAP2_read_params_xml cannot find z")
    read (value, *) parse_ip%z(ti)

    if (allocated(parse_ip%type_of_atomic_num)) deallocate(parse_ip%type_of_atomic_num)
    allocate(parse_ip%type_of_atomic_num(maxval(parse_ip%atomic_num)))
    parse_ip%type_of_atomic_num = 0
    do ti=1, parse_ip%n_types
      if (parse_ip%atomic_num(ti) > 0) &
        parse_ip%type_of_atomic_num(parse_ip%atomic_num(ti)) = ti
    end do

  elseif (parse_in_ip .and. name == 'per_pair_data') then

    call QUIP_FoX_get_value(attributes, "atnum_i", value, status)
    if (status /= 0) call system_abort ("IPModel_SW_read_params_xml cannot find atnum_i")
    read (value, *) Zi
    call QUIP_FoX_get_value(attributes, "atnum_j", value, status)
    if (status /= 0) call system_abort ("IPModel_SW_read_params_xml cannot find atnum_j")
    read (value, *) Zj

    ti = get_type(parse_ip%type_of_atomic_num,Zi)
    tj = get_type(parse_ip%type_of_atomic_num,Zj)

    call QUIP_FoX_get_value(attributes, "D_ms", value, status)
    if (status /= 0) call system_abort ("IPModel_ASAP2_read_params_xml cannot find D_ms")
    read (value, *) parse_ip%D_ms(ti,tj)
    call QUIP_FoX_get_value(attributes, "gamma_ms", value, status)
    if (status /= 0) call system_abort ("IPModel_ASAP2_read_params_xml cannot find gamma_ms")
    read (value, *) parse_ip%gamma_ms(ti,tj)
    call QUIP_FoX_get_value(attributes, "R_ms", value, status)
    if (status /= 0) call system_abort ("IPModel_ASAP2_read_params_xml cannot find R_ms")
    read (value, *) parse_ip%R_ms(ti,tj)
    call QUIP_FoX_get_value(attributes, "B_pol", value, status)
    if (status /= 0) call system_abort ("IPModel_ASAP2_read_params_xml cannot find B_pol")
    read (value, *) parse_ip%B_pol(ti,tj)
    call QUIP_FoX_get_value(attributes, "C_pol", value, status)
    if (status /= 0) call system_abort ("IPModel_ASAP2_read_params_xml cannot find C_pol")
    read (value, *) parse_ip%C_pol(ti,tj)

    if (ti /= tj) then
      parse_ip%D_ms(tj,ti) = parse_ip%D_ms(ti,tj)
      parse_ip%gamma_ms(tj,ti) = parse_ip%gamma_ms(ti,tj)
      parse_ip%R_ms(tj,ti) = parse_ip%R_ms(ti,tj)
      parse_ip%B_pol(tj,ti) = parse_ip%B_pol(ti,tj)
      parse_ip%C_pol(tj,ti) = parse_ip%C_pol(ti,tj)
    endif

  endif

end subroutine IPModel_startElement_handler

subroutine IPModel_endElement_handler(URI, localname, name)
  character(len=*), intent(in)   :: URI
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name

  if (parse_in_ip) then
    if (name == 'ASAP2_params') then
      parse_in_ip = .false.
    end if
  endif

end subroutine IPModel_endElement_handler

end module IPModel_ASAP2_module
