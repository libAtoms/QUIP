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

#include "error.inc"

module Yukawa_module

use libatoms_module
use functions_module
use QUIP_Common_module

implicit none
private

public :: yukawa_charges, yukawa_dipoles

real(dp), parameter :: sqrt_2 = 1.41421356237309504880168872421_dp
real(dp), parameter :: sqrt_pi = 1.77245385090551602729816748334_dp

contains

!% Charge-charge interactions, screened by Yukawa function
subroutine yukawa_charges(at, charge, cutoff_coulomb, yukalpha, yuksmoothlength, &
     e, local_e, f, virial, efield, mpi, atom_mask_name, source_mask_name, pseudise, grid_size, error)
   type(Atoms), intent(inout)      :: at
   real(dp), dimension(:), intent(in) :: charge
   real(dp), intent(in) :: cutoff_coulomb, yukalpha, yuksmoothlength
   real(dp), intent(out), optional :: e, local_e(:)
   real(dp), intent(out), optional :: f(:,:)
   real(dp), intent(out), optional :: virial(3,3)
   real(dp), intent(out), optional :: efield(:,:)
   type(MPI_Context), intent(in), optional :: mpi
   character(len=*), optional, intent(in) :: atom_mask_name, source_mask_name
   logical, optional, intent(in) :: pseudise
   real(dp), optional, intent(in) :: grid_size
   integer, intent(out), optional :: error

   real(dp) :: erf_val, erf_deriv
   integer i, j, m
   real(dp) :: r_ij, u_ij(3), zv2, gamjir, gamjir3, gamjir2, fc, dfc_dr
   real(dp) :: de, dforce, expfactor, defield(3)
   logical :: i_is_min_image, j_is_min_image, do_pseudise
   real(dp) :: private_virial(3,3), private_e, pseudise_sigma
   real(dp), allocatable :: private_f(:,:), private_efield(:,:), private_local_e(:)
   logical, pointer, dimension(:) :: atom_mask, source_mask

   INIT_ERROR(error)
   call system_timer('yukawa_charges')

   do_pseudise = optional_default(.false., pseudise)

   atom_mask => null()
   if (present(atom_mask_name)) then
      if (trim(atom_mask_name) /= '') then
         call assign_property_pointer(at, atom_mask_name, atom_mask, error)
         PASS_ERROR(error)
      end if
   end if

   source_mask => null()
   if (present(source_mask_name)) then
      if (trim(source_mask_name) /= '') then
         call assign_property_pointer(at, source_mask_name, source_mask, error)
         PASS_ERROR(error)
      end if
   end if

   !$omp parallel default(none) shared(mpi, charge, at, e, local_e, f, virial, efield, atom_mask, source_mask, cutoff_coulomb, yukalpha, yuksmoothlength, grid_size, do_pseudise) private(i, j, m, r_ij, u_ij, zv2, gamjir, gamjir3, gamjir2, fc, dfc_dr, de, dforce, expfactor, i_is_min_image, j_is_min_image, private_virial, private_e, private_f, private_local_e, private_efield, erf_val, erf_deriv, pseudise_sigma)

   if (present(e)) private_e = 0.0_dp
   if (present(local_e)) then
      allocate(private_local_e(at%N))
      private_local_e = 0.0_dp
   endif
   if (present(f)) then
      allocate(private_f(3,at%N))
      private_f = 0.0_dp
   endif
   if (present(efield)) then
      allocate(private_efield(3,at%N))
      private_efield = 0.0_dp
   end if
   if (present(virial)) private_virial = 0.0_dp

   !$omp do schedule(runtime)
   do i=1, at%n
      if (present(mpi)) then
	 if (mpi%active) then
	    if (mod(i-1, mpi%n_procs) /= mpi%my_proc) cycle
	 endif
      endif

      if (associated(atom_mask)) then
         if (.not. atom_mask(i)) cycle
         i_is_min_image = .false.
      else
         if (allocated(at%connect%is_min_image)) then
            i_is_min_image = at%connect%is_min_image(i)
         else
            i_is_min_image = is_min_image(at, i)
         end if
      end if

      do m = 1, atoms_n_neighbours(at, i)
         
         j = atoms_neighbour(at, i, m, distance=r_ij, cosines=u_ij, max_dist=(cutoff_coulomb*BOHR))
         if (j <= 0) cycle
         
         if (associated(source_mask)) then
            if (.not. source_mask(j)) cycle
         end if

         if (r_ij .feq. 0.0_dp) then
            if (present(grid_size)) then 
               r_ij = grid_size
            else
               cycle
            end if
         end if

         if (allocated(at%connect%is_min_image)) then
            j_is_min_image = at%connect%is_min_image(j)
         else
            j_is_min_image = is_min_image(at, j)
         end if

         if (i < j .and. i_is_min_image .and. j_is_min_image) cycle

         r_ij = r_ij/BOHR
         zv2 = charge(i)*charge(j)

         gamjir = zv2/r_ij
         gamjir3 = gamjir/(r_ij**2.0_dp)
         gamjir2 = zv2/(r_ij**2.0_dp)
         expfactor = exp(-yukalpha*r_ij)

         call smooth_cutoff(r_ij, cutoff_coulomb-yuksmoothlength, yuksmoothlength, fc, dfc_dr)

         de = gamjir

         if (do_pseudise) then
            pseudise_sigma = (ElementCovRad(at%z(i)) + ElementCovRad(at%z(j)))/(2.0_dp*BOHR)
            erf_val = 1.0_dp
            erf_deriv = 0.0_dp
            ! pseudise if sigma > 0, correction for r >= 9s is < 1e-16
            if (pseudise_sigma > 0.0_dp .and. r_ij < 9.0_dp*pseudise_sigma) then
               erf_val = derf(r_ij/(sqrt_2 * pseudise_sigma))
               erf_deriv = sqrt_2/(sqrt_pi*pseudise_sigma)*exp(-r_ij*r_ij/(2.0_dp*pseudise_sigma*pseudise_sigma))
            end if
         end if

         if (present(e) .or. present(local_e)) then
            if (present(e)) then
               if (i_is_min_image .and. j_is_min_image) then
                  if (do_pseudise) then ! .and. pseudise_sigma > 0.0_dp .and. r_ij < 9.0_dp*pseudise_sigma) then
                     private_e = private_e + de*expfactor*fc*erf_val
                  else
                     private_e = private_e + de*expfactor*fc
                  end if
               else
                  if (do_pseudise) then ! .and. pseudise_sigma > 0.0_dp .and. r_ij < 9.0_dp*pseudise_sigma) then
                     private_e = private_e + 0.5_dp*de*expfactor*fc*erf_val
                  else
                     private_e = private_e + 0.5_dp*de*expfactor*fc
                  end if
               end if
            end if
            if (present(local_e)) then
               if (do_pseudise) then ! .and. pseudise_sigma > 0.0_dp .and. r_ij < 9.0_dp*pseudise_sigma) then
                  private_local_e(i) = private_local_e(i) + 0.5_dp*de*expfactor*fc*erf_val
                  if (i_is_min_image .and. j_is_min_image) private_local_e(j) = private_local_e(j) + 0.5_dp*de*expfactor*fc*erf_val
               else
                  private_local_e(i) = private_local_e(i) + 0.5_dp*de*expfactor*fc
                  if (i_is_min_image .and. j_is_min_image) private_local_e(j) = private_local_e(j) + 0.5_dp*de*expfactor*fc
               end if
            end if
         end if

         if (present(f) .or. present(virial) .or. present(efield)) then
            dforce = gamjir3*expfactor*fc*r_ij + de*(yukalpha*fc - dfc_dr)*expfactor

            if (do_pseudise) then ! .and. pseudise_sigma > 0.0_dp .and. r_ij < 9.0_dp*pseudise_sigma) then
               dforce = dforce*erf_val - de*expfactor*fc*erf_deriv
            end if

            if (present(f)) then
               private_f(:,i) = private_f(:,i) - dforce*u_ij
               if (i_is_min_image .and. j_is_min_image) private_f(:,j) = private_f(:,j) + dforce*u_ij
            end if

            if (present(virial)) then
               if (i_is_min_image .and. j_is_min_image) then
                  private_virial = private_virial + dforce*(u_ij .outer. u_ij)*r_ij
               else
                  private_virial = private_virial + 0.5_dp*dforce*(u_ij .outer. u_ij)*r_ij
               end if
            end if

            if (present(efield)) then
               defield = gamjir3*expfactor*fc*u_ij*r_ij
               if (do_pseudise) defield = defield*erf_val
               private_efield(:,i) = private_efield(:,i) - defield/charge(i)
               if (i_is_min_image .and. j_is_min_image) private_efield(:,j) = private_efield(:,j) + defield/charge(j)
            end if
         end if

      end do
   end do

   if (present(mpi)) then
      if (mpi%active) then
	 if (present(e)) private_e = sum(mpi, private_e) 
	 if (present(local_e)) call sum_in_place(mpi, private_local_e)
	 if (present(f)) call sum_in_place(mpi, private_f)
	 if (present(virial)) call sum_in_place(mpi, private_virial)
	 if (present(efield)) call sum_in_place(mpi, private_efield)
      end if
   end if

   !$omp critical
   if (present(e)) e = e + private_e
   if (present(local_e)) local_e = local_e + private_local_e
   if (present(f)) f = f + private_f
   if (present(virial)) virial = virial + private_virial
   if (present(efield)) efield = efield + private_efield
   !$omp end critical 

   if (allocated(private_f)) deallocate(private_f)
   if (allocated(private_local_e)) deallocate(private_local_e)
   if (allocated(private_efield)) deallocate(private_efield)

   !$omp end parallel

   call system_timer('yukawa_charges')

 end subroutine yukawa_charges


!% Charge-dipole and dipole-dipole interactions, screened by Yukawa function
subroutine yukawa_dipoles(at, charge, dip, cutoff_coulomb, yukalpha, yuksmoothlength, pol, b_pol, c_pol, &
     type_of_atomic_num, tdip_sr, e, local_e, f, virial, efield, mpi, atom_mask_name, source_mask_name, pseudise, grid_size, error)
#ifdef _OPENMP
   use omp_lib
#endif

   type(Atoms), intent(inout)      :: at
   real(dp), intent(in)            :: charge(:), dip(:,:)
   real(dp), intent(in) :: cutoff_coulomb, yukalpha, yuksmoothlength, pol(:), b_pol(:,:), c_pol(:,:)
   integer, intent(in) :: type_of_atomic_num(:)
   logical, intent(in) :: tdip_sr
   real(dp), intent(out), optional :: e, local_e(:)
   real(dp), intent(out), optional :: f(:,:)
   real(dp), intent(out), optional :: virial(3,3)
   real(dp), intent(out), optional :: efield(:,:)
   type(MPI_Context), intent(in), optional :: mpi
   character(len=*), optional, intent(in) :: atom_mask_name, source_mask_name
   logical, optional, intent(in) :: pseudise
   real(dp), optional, intent(in) :: grid_size
   integer, intent(out), optional :: error

   integer i, j, m, ti, tj, k
   real(dp) :: r_ij, u_ij(3), gamjir3, gamjir2, fc, dfc_dr, de, pseudise_sigma, dforce(3), erf_val, erf_deriv
   real(dp) :: expfactor, dipi(3), dipj(3), qj, qi, pp, pri, prj, defield_i(3), defield_j(3)
   real(dp) :: de_ind, de_dd, de_qd, dfqdip(3), dfdipdip(3), factor1, dist3, dist5
   real(dp) :: const1, const2, factork, de_sr, df_sr(3), gij, dgijdrij, bij, cij
   logical :: i_is_min_image, j_is_min_image, tpoli, tpolj, qipj, qjpi, pipj, do_pseudise

   real(dp) :: private_virial(3,3), private_e
   real(dp), allocatable :: private_f(:,:), private_local_e(:), private_efield(:,:)
   logical, pointer, dimension(:) :: atom_mask, source_mask

   INIT_ERROR(error)
   call system_timer('yukawa_dipoles')

   do_pseudise = optional_default(.false., pseudise)

   atom_mask => null()
   if (present(atom_mask_name)) then
      if (trim(atom_mask_name) /= '') then
         call assign_property_pointer(at, atom_mask_name, atom_mask, error)
         PASS_ERROR(error)
      end if
   end if

   source_mask => null()
   if (present(source_mask_name)) then
      if (trim(source_mask_name) /= '') then
         call assign_property_pointer(at, source_mask_name, source_mask, error)
         PASS_ERROR(error)
      end if
   end if

   !$omp parallel default(none) shared(mpi, at, charge, dip, e, local_e, f, virial, efield, type_of_atomic_num, cutoff_coulomb, yukalpha, yuksmoothlength, pol, b_pol, c_pol, tdip_sr, do_pseudise, atom_mask, source_mask) private(i, j, m, ti, tj, k, r_ij, u_ij, gamjir3, gamjir2, fc, dfc_dr, expfactor, dipi, dipj, qj, qi, pp, pri, prj, de_ind, de_dd, de_qd, dfqdip, dfdipdip, factor1, dist3, dist5, const1, const2, factork, de_sr, df_sr, gij, dgijdrij, bij, cij, i_is_min_image, j_is_min_image, tpoli, tpolj, qipj, qjpi, pipj, private_e, private_local_e, private_virial, private_f, private_efield, pseudise_sigma)

   if (present(e)) private_e = 0.0_dp
   if (present(local_e)) then
      allocate(private_local_e(at%N))
      private_local_e = 0.0_dp
   endif
   if (present(f)) then
      allocate(private_f(3,at%N))
      private_f = 0.0_dp
   endif
   if (present(efield)) then
      allocate(private_efield(3,at%N))
      private_efield = 0.0_dp
   end if
   if (present(virial)) private_virial = 0.0_dp

   !$omp do schedule(runtime)
   do i=1, at%n
      if (present(mpi)) then
	 if (mpi%active) then
	    if (mod(i-1, mpi%n_procs) /= mpi%my_proc) cycle
	 endif
      endif

      if (associated(atom_mask)) then
         if (.not. atom_mask(i)) cycle
         i_is_min_image = .false.
      else
         if (allocated(at%connect%is_min_image)) then
            i_is_min_image = at%connect%is_min_image(i)
         else
            i_is_min_image = is_min_image(at, i)
         end if
      end if
      ti = get_type(type_of_atomic_num, at%Z(i))

      qi = charge(i)
      dipi = dip(:,i)
      tpoli = abs(pol(ti)) > 0.0_dp

      ! Induced contribution to energy
      if ((present(e) .or. present(local_e)) .and. tpoli) then
         de_ind = 0.5_dp*(dipi .dot. dipi)/pol(ti)
         if (present(e))       private_e = private_e + de_ind
         if (present(local_e)) private_local_e(i) = private_local_e(i) + de_ind
      end if

      do m = 1, atoms_n_neighbours(at, i)
         
         j = atoms_neighbour(at, i, m, distance=r_ij, diff=u_ij, max_dist=(cutoff_coulomb*BOHR))
         if (j <= 0) cycle

         if (associated(source_mask)) then
            if (.not. source_mask(j)) cycle
         end if

         if (r_ij .feq. 0.0_dp) then
            if (present(grid_size)) then 
               r_ij = grid_size
            else
               cycle
            end if
         end if

         if (allocated(at%connect%is_min_image)) then
            j_is_min_image = at%connect%is_min_image(j)
         else
            j_is_min_image = is_min_image(at, j)
         end if

         if (i < j .and. i_is_min_image .and. j_is_min_image) cycle

         r_ij = r_ij/BOHR
         u_ij = u_ij/BOHR
         tj = get_type(type_of_atomic_num, at%Z(j))

         qj = charge(j)
         dipj = dip(:,j)
         tpolj = abs(pol(tj)) > 0.0_dp

         qipj = (abs(qi) > 0.0_dp) .and. tpolj
         qjpi = (abs(qj) > 0.0_dp) .and. tpoli
         pipj = tpoli .and. tpolj

         if (.not. (pipj .or. qipj .or. qjpi)) cycle
        
         gamjir2 = 1.0_dp/r_ij**2.0_dp
         gamjir3 = gamjir2/r_ij
         factor1 = 3.0_dp*gamjir2*gamjir3

         expfactor = exp(-yukalpha*r_ij)

         call smooth_cutoff(r_ij, cutoff_coulomb-yuksmoothlength, yuksmoothlength, fc, dfc_dr)

         pp = 0.0_dp
         if (pipj)  pp  = dipi .dot. dipj
         pri = 0.0_dp
         if (tpoli) pri = dipi .dot. u_ij

         prj = 0.0_dp
         if (tpolj) prj = dipj .dot. u_ij
         
         if (do_pseudise) then
            pseudise_sigma = (ElementCovRad(at%z(i)) + ElementCovRad(at%z(j)))/(2.0_dp*BOHR)
            erf_val = 1.0_dp
            erf_deriv = 0.0_dp
            ! pseudise if sigma > 0, correction for r >= 9s is < 1e-16
            if (pseudise_sigma > 0.0_dp .and. r_ij < 9.0_dp*pseudise_sigma) then
               erf_val = derf(r_ij/(sqrt_2 * pseudise_sigma))
               erf_deriv = sqrt_2/(sqrt_pi*pseudise_sigma)*exp(-r_ij*r_ij/(2.0_dp*pseudise_sigma*pseudise_sigma))
            end if
         end if

         if (present(efield)) then
            defield_i = (3.0_dp*prj*u_ij*gamjir2 - dipj)*gamjir2*dsqrt(gamjir2)*expfactor*fc
            defield_j = (3.0_dp*pri*u_ij*gamjir2 - dipi)*gamjir2*dsqrt(gamjir2)*expfactor*fc
            if (do_pseudise) then
               defield_i = defield_i*erf_val
               defield_j = defield_j*erf_val
            end if
            private_efield(:,i) = private_efield(:,i) + defield_i
            if (i_is_min_image .and. j_is_min_image) private_efield(:,j) = private_efield(:,j) + defield_j
         end if

         if (present(e) .or. present(local_e) .or. present(virial) .or. present(f)) then

            de_dd = (pp - 3.0_dp*pri*prj*gamjir2)*gamjir3
            de_qd = -(qi*prj - qj*pri)*gamjir3

            de_sr = 0.0_dp
            df_sr = 0.0_dp
            if (tdip_sr) then

               bij = b_pol(ti, tj)
               cij = c_pol(ti, tj)
         
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
                    -de_sr*(yukalpha*fc - dfc_dr)*expfactor/r_ij*u_ij
            end if

            de = de_dd + de_qd + de_sr

            if (present(e)) then
               if (do_pseudise) then
                  if (i_is_min_image .and. j_is_min_image) then
                     private_e = private_e + de*expfactor*fc*erf_val
                  else
                     private_e = private_e + 0.5_dp*de*expfactor*fc*erf_val
                  end if
               else
                  if (i_is_min_image .and. j_is_min_image) then
                     private_e = private_e + de*expfactor*fc
                  else
                     private_e = private_e + 0.5_dp*de*expfactor*fc
                  end if
               end if
            end if

            if (present(local_e)) then
               private_local_e(i) = private_local_e(i) + 0.5_dp*de*expfactor*fc
               if (i_is_min_image .and. j_is_min_image) private_local_e(j) = private_local_e(j) + 0.5_dp*de*expfactor*fc
            end if

            if (present(f) .or. present(virial)) then
               dfqdip = .0_dp
               if (qipj .or. qjpi) then
                  dfqdip = ((qj*dipi(:) - qi*dipj(:) &
                       - 3.0_dp*(qj*pri-qi*prj)*u_ij*gamjir2)*gamjir3)*expfactor*fc &
                       - de_qd*(yukalpha*fc - dfc_dr)*expfactor/r_ij*u_ij
               end if

               dfdipdip = 0.0_dp
               if (pipj) then
                  dfdipdip = (-(pp*u_ij + dipj(:)*pri+dipi(:)*prj &
                       - 5.0_dp*prj*pri*u_ij(:)*gamjir2)*factor1)*expfactor*fc &
                       - de_dd*(yukalpha*fc - dfc_dr)*expfactor/r_ij*u_ij(:)
               end if

               dforce = dfqdip + dfdipdip + df_sr

               if (do_pseudise) then
                  dforce = dforce*erf_val + de*expfactor*fc*erf_deriv*u_ij/r_ij
               end if
               
               if (present(f)) then
                  private_f(:,i) = private_f(:,i) + dforce
                  if (i_is_min_image .and. j_is_min_image) private_f(:,j) = private_f(:,j) - dforce
               end if

               if (present(virial)) then
                  if (i_is_min_image .and. j_is_min_image) then
                     private_virial = private_virial - (dforce .outer. u_ij)
                  else
                     private_virial = private_virial - 0.5_dp*(dforce .outer. u_ij)
                  end if
               end if
            end if
         end if

      end do
   end do
   !$omp end do

   if (present(mpi)) then
      if (mpi%active) then
	 if (present(e)) private_e = sum(mpi, private_e) 
	 if (present(local_e)) call sum_in_place(mpi, private_local_e)
	 if (present(f)) call sum_in_place(mpi, private_f)
	 if (present(virial)) call sum_in_place(mpi, private_virial)
	 if (present(efield)) call sum_in_place(mpi, private_efield)
      end if
   end if

   !$omp critical
   if (present(e)) e = e + private_e
   if (present(f)) f = f + private_f
   if (present(local_e)) local_e = local_e + private_local_e
   if (present(virial)) virial = virial + private_virial
   if (present(efield)) efield = efield + private_efield
   !$omp end critical 
   
   if (allocated(private_f)) deallocate(private_f)
   if (allocated(private_local_e)) deallocate(private_local_e)
   if (allocated(private_efield)) deallocate(private_efield)

   !$omp end parallel   

   call system_timer('yukawa_dipoles')

end subroutine yukawa_dipoles

end module Yukawa_module
