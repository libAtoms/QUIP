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

module IPEwald_module

use error_module
use system_module, only : dp, optional_default, PRINT_ANAL, operator(//)
use units_module
use Atoms_module
use linearalgebra_module
use functions_module

implicit none

real(dp), parameter :: reciprocal_time_by_real_time = 1.0_dp / 3.0_dp

private
public :: Ewald_calc, Ewald_corr_calc, Direct_Coulomb_Calc, DSF_Coulomb_calc

contains

  ! Ewald routine
  ! input: atoms object, has to have charge property
  ! input, optional: ewald_error (controls speed and ewald_error)
  ! input, optional: use_ewald_cutoff (forces original cutoff to be used)
  ! output: energy, force, virial

  ! procedure to determine optimal Ewald parameters:
  ! Optimization of the Ewald sum for large systems, Mol. Simul. 13 (1994), no. 1, 1-9.

  subroutine Ewald_calc(at_in, charge, e, f, virial, ewald_error, use_ewald_cutoff, smooth_coulomb_cutoff, error)

    type(Atoms), intent(in), target    :: at_in
    real(dp), dimension(:), intent(in) :: charge

    real(dp), intent(out), optional                    :: e
    real(dp), dimension(:,:), intent(out), optional    :: f
    real(dp), dimension(3,3), intent(out), optional    :: virial
    real(dp), intent(in), optional                     :: ewald_error
    logical, intent(in), optional                      :: use_ewald_cutoff
    real(dp), intent(in), optional                     :: smooth_coulomb_cutoff
    integer, intent(out), optional                     :: error

    integer  :: i, j, k, n, n1, n2, n3, not_needed !for reciprocal force
    integer, dimension(3) :: nmax !how many reciprocal vectors are to be taken

    logical :: my_use_ewald_cutoff

    real(dp) :: r_ij, erfc_ar, arg, my_ewald_error, alpha, kmax, kmax2, prefac, infac, two_alpha_over_sqrt_pi, v, &
    & ewald_precision, ewald_cutoff, my_cutoff, my_smooth_coulomb_cutoff, smooth_arg, smooth_f, dsmooth_f

    real(dp), dimension(3) :: force, u_ij, a, b, c, h
    real(dp), dimension(3,3) :: identity3x3, k3x3
    real(dp), dimension(:,:,:,:), allocatable :: coskr, sinkr 
    real(dp), dimension(:,:,:,:), allocatable :: k_vec   ! reciprocal vectors
    real(dp), dimension(:,:,:,:), allocatable :: force_factor
    real(dp), dimension(:,:,:), allocatable   :: energy_factor
    real(dp), dimension(:,:,:), allocatable   :: mod2_k  !square of length of reciprocal vectors

    type(Atoms), target :: my_at
    type(Atoms), pointer :: at

    INIT_ERROR(error)

    call check_size('charge',charge,at_in%N,'IPEwald',error)

    identity3x3 = 0.0_dp
    call add_identity(identity3x3)

    ! Set up Ewald calculation
    my_ewald_error = optional_default(1e-06_dp,ewald_error) * 4.0_dp * PI * EPSILON_0 ! convert eV to internal units
    my_use_ewald_cutoff = optional_default(.true.,use_ewald_cutoff) ! can choose between optimal Ewald 
    my_smooth_coulomb_cutoff = optional_default(0.0_dp, smooth_coulomb_cutoff) ! default is not to use smooth Coulomb

    a = at_in%lattice(:,1); b = at_in%lattice(:,2); c = at_in%lattice(:,3)
    v = cell_volume(at_in)

    h(1) = v / norm(b .cross. c)
    h(2) = v / norm(c .cross. a)
    h(3) = v / norm(a .cross. b)

    ewald_precision = -log(my_ewald_error)
    ewald_cutoff = sqrt(ewald_precision/PI) * reciprocal_time_by_real_time**(1.0_dp/6.0_dp) * &
    & minval(sqrt( sum(at_in%lattice(:,:)**2,dim=1) )) / at_in%N**(1.0_dp/6.0_dp)

    call print('Ewald cutoff = '//ewald_cutoff,PRINT_ANAL)

    if( my_use_ewald_cutoff .and. (ewald_cutoff > at_in%cutoff) ) then
        my_at = at_in
        call set_cutoff(my_at,ewald_cutoff)
        call calc_connect(my_at)
        at => my_at
    else
        at => at_in
    endif

    if( my_use_ewald_cutoff ) then
       my_cutoff = ewald_cutoff
    else
       my_cutoff = at_in%cutoff
    endif

    if( my_cutoff < my_smooth_coulomb_cutoff ) then
       RAISE_ERROR('Cutoff='//my_cutoff//' is smaller than the smooth region specified by smooth_coulomb_cutoff='//my_smooth_coulomb_cutoff, error)
    endif
         
    alpha = sqrt(ewald_precision) / my_cutoff 
    call print('Ewald alpha = '//alpha,PRINT_ANAL)

    kmax = 2.0_dp * ewald_precision / my_cutoff
    kmax2 = kmax**2
    call print('Ewald kmax = '//kmax,PRINT_ANAL)

    nmax = nint( kmax * h / 2.0_dp / PI )
    call print('Ewald nmax = '//nmax,PRINT_ANAL)

    two_alpha_over_sqrt_pi = 2.0_dp * alpha / sqrt(PI)

    prefac = 4.0_dp * PI / v
    infac  = - 1.0_dp / (4.0_dp * alpha**2.0_dp) 

    allocate( k_vec(3,-nmax(3):nmax(3),-nmax(2):nmax(2),0:nmax(1)), &
            & mod2_k(-nmax(3):nmax(3),-nmax(2):nmax(2),0:nmax(1) ),  &
            & force_factor(-nmax(3):nmax(3),-nmax(2):nmax(2),0:nmax(1),3), & 
            & energy_factor(-nmax(3):nmax(3),-nmax(2):nmax(2),0:nmax(1) ) )

    allocate( coskr(-nmax(3):nmax(3),-nmax(2):nmax(2),0:nmax(1),at%N), &
            & sinkr(-nmax(3):nmax(3),-nmax(2):nmax(2),0:nmax(1),at%N) )

    k_vec = 0.0_dp
    mod2_k = 0.0_dp
    force_factor = 0.0_dp
    energy_factor = 0.0_dp

    coskr = 0.0_dp
    sinkr = 0.0_dp

    not_needed = ( 2*nmax(3) + 1 ) * nmax(2) + nmax(3) + 1 ! lot of symmetries for k and -k so count only

    n = 0
    do n1 = 0, nmax(1)
       do n2 = -nmax(2), nmax(2)
           do n3 = -nmax(3), nmax(3)

              n = n + 1

              if( n>not_needed ) then
                  k_vec(:,n3,n2,n1) = ( at%g(1,:)*n1 + at%g(2,:)*n2 + at%g(3,:)*n3 ) * 2.0_dp * PI
                  mod2_k(n3,n2,n1)  = normsq( k_vec(:,n3,n2,n1) )

                  force_factor(n3,n2,n1,:) = 1.0_dp/mod2_k(n3,n2,n1) * &
                  & exp(mod2_k(n3,n2,n1)*infac) * k_vec(:,n3,n2,n1)
                  energy_factor(n3,n2,n1)  = 1.0_dp/mod2_k(n3,n2,n1) * exp(infac*mod2_k(n3,n2,n1))
              endif

           enddo
       enddo
    enddo

    do i = 1, at%N
       n = 0
       do n1 = 0, nmax(1)
          do n2 = -nmax(2), nmax(2)
             do n3 = -nmax(3), nmax(3)
             
                n = n + 1

                if( (n>not_needed) .and. ( mod2_k(n3,n2,n1)<kmax2 ) ) then
                   arg = dot_product(at%pos(:,i), k_vec(:,n3,n2,n1))
                   coskr(n3,n2,n1,i) = cos(arg)*charge(i)
                   sinkr(n3,n2,n1,i) = sin(arg)*charge(i)

                endif
             enddo
          enddo
       enddo
    enddo

    if(present(e)) e = 0.0_dp
    if(present(f)) f  = 0.0_dp
    if(present(virial)) virial  = 0.0_dp

    do i=1,at%N
       !Loop over neighbours
       do n = 1, n_neighbours(at,i)
          j = neighbour(at,i,n,distance=r_ij,cosines=u_ij) ! nth neighbour of atom i
          if( r_ij > my_cutoff )  cycle

          if( r_ij < my_smooth_coulomb_cutoff ) then
             smooth_arg = r_ij * PI / my_smooth_coulomb_cutoff / 2.0_dp
             smooth_f = ( 1.0_dp - sin(smooth_arg) ) / r_ij
             dsmooth_f = cos(smooth_arg) * PI / my_smooth_coulomb_cutoff / 2.0_dp
          else
             smooth_f = 0.0_dp
             dsmooth_f = 0.0_dp
          endif
           
          erfc_ar = erfc(r_ij*alpha)/r_ij

          if( present(e) ) e = e + 0.5_dp * charge(i)*charge(j)* ( erfc_ar - smooth_f )

          if( present(f) .or. present(virial) ) then
              force(:) = charge(i)*charge(j) * &
              & ( two_alpha_over_sqrt_pi * exp(-(r_ij*alpha)**2) + erfc_ar - smooth_f - dsmooth_f) / r_ij * u_ij(:)

              if(present(f)) then
                 f(:,i) = f(:,i) - force(:) 
              endif

              if (present(virial)) virial = virial + 0.5_dp * (force .outer. u_ij) * r_ij
          endif
 
      enddo
    enddo
             
    ! reciprocal energy
    if(present(e)) e = e + sum((sum(coskr,dim=4)**2 + sum(sinkr,dim=4)**2)*energy_factor) * prefac &
    & - sum(charge**2) * alpha / sqrt(PI) - PI / ( 2.0_dp * alpha**2 * v ) * sum(charge)**2

    ! reciprocal force
    if( present(f) ) then
        do i = 1, at%N
           do j = 1, at%N

              if( i<=j ) cycle

              force = (/( sum(force_factor(:,:,:,k) * &
              & ( sinkr(:,:,:,j)*coskr(:,:,:,i) - sinkr(:,:,:,i)*coskr(:,:,:,j) ) ), k=1,3 )/) * prefac * 2.0_dp

              !force acting on atom j by atom i
              f(:,i) = f(:,i) - force(:)
              f(:,j) = f(:,j) + force(:)

           enddo
        enddo
    endif

    ! reciprocal contribution to virial
    if(present(virial)) then
       n = 0
       do n1 = 0, nmax(1)
          do n2 = -nmax(2), nmax(2)
             do n3 = -nmax(3), nmax(3)
             
                n = n + 1

                if( (n>not_needed) .and. ( mod2_k(n3,n2,n1)<kmax2 ) ) then
   
                   k3x3 = k_vec(:,n3,n2,n1) .outer. k_vec(:,n3,n2,n1)
                   virial = virial &
                   & + ( identity3x3 - 2*(-infac + 1.0_dp/mod2_k(n3,n2,n1))*k3x3 ) * &
                   & (sum(coskr(n3,n2,n1,:))**2 + sum(sinkr(n3,n2,n1,:))**2) * energy_factor(n3,n2,n1) * &
                   & prefac

                endif
             enddo
          enddo
       enddo
    endif

    if(present(virial)) virial = virial - identity3x3 * sum(charge)**2 * PI / v / alpha**2 / 2

   ! if(present(e)) e = e / ( 4.0_dp * PI * EPSILON_0 ) ! convert from internal units to eV
   ! if(present(f)) f = f / ( 4.0_dp * PI * EPSILON_0 ) ! convert from internal units to eV/A
   ! if(present(virial)) virial = virial / ( 4.0_dp * PI * EPSILON_0 )

    if(present(e)) e = e * HARTREE*BOHR ! convert from internal units to eV
    if(present(f)) f = f * HARTREE*BOHR ! convert from internal units to eV/A
    if(present(virial)) virial = virial * HARTREE*BOHR


    deallocate( coskr, sinkr )
    deallocate( k_vec, mod2_k, force_factor, energy_factor )
    if (associated(at,my_at)) call finalise(my_at)

  endsubroutine Ewald_calc

  subroutine Ewald_corr_calc(at_in,charge, e,f,virial,cutoff,error)

    type(Atoms), intent(in), target    :: at_in
    real(dp), dimension(:), intent(in) :: charge

    real(dp), intent(out), optional                    :: e
    real(dp), dimension(:,:), intent(out), optional    :: f
    real(dp), dimension(3,3), intent(out), optional    :: virial
    real(dp), intent(in), optional                     :: cutoff
    integer, intent(out), optional                     :: error

    integer  :: i, j, n

    real(dp) :: my_cutoff, r_ij, de
    real(dp), dimension(3) :: force, u_ij

    real(dp), dimension(:), pointer :: my_charge

    type(Atoms), target :: my_at
    type(Atoms), pointer :: at => null()

    INIT_ERROR(error)
    call check_size('charge',charge,(/at_in%N/),'Ewald_corr_calc',error)

    my_cutoff = optional_default(at_in%cutoff,cutoff)

    if( present(cutoff) .and. (my_cutoff > at_in%cutoff) ) then
        my_at = at_in
        call set_cutoff(my_at,cutoff)
        call calc_connect(my_at)
        at => my_at
    else
        at => at_in
    endif
         
    if( present(e) ) e = 0.0_dp
    if( present(f) ) f = 0.0_dp
    if( present(virial) ) virial = 0.0_dp

    do i = 1, at%N
       !Loop over neighbours
       do n = 1, n_neighbours(at,i)
          j = neighbour(at,i,n,distance=r_ij,cosines=u_ij) ! nth neighbour of atom i
          if( r_ij > my_cutoff )  cycle
           
          de = 0.5_dp * ( cos(r_ij*PI/my_cutoff) + 1.0_dp ) / r_ij

          if( present(e) ) e = e + 0.5_dp * de * charge(i)*charge(j)

          if( present(f) .or. present(virial) ) then
              force = charge(i)*charge(j) * &
              & ( -de - 0.5*PI*sin(r_ij*PI/my_cutoff)/my_cutoff ) / r_ij * u_ij

              if(present(f)) then
                 f(:,i) = f(:,i) + force
              endif

              if (present(virial)) virial = virial - 0.5_dp * (force .outer. u_ij) * r_ij
          endif
 
      enddo
    enddo
             
    if(present(e)) e = e * HARTREE*BOHR ! convert from internal units to eV
    if(present(f)) f = f * HARTREE*BOHR ! convert from internal units to eV/A
    if(present(virial)) virial = virial * HARTREE*BOHR

    if(associated(at,my_at)) call finalise(my_at)
    at => null()

  endsubroutine Ewald_corr_calc

  subroutine Direct_Coulomb_calc(at_in,charge, e,f,virial,local_e,cutoff,error)

    type(Atoms), intent(in), target    :: at_in
    real(dp), dimension(:), intent(in) :: charge

    real(dp), intent(out), optional                    :: e
    real(dp), dimension(:,:), intent(out), optional    :: f
    real(dp), dimension(3,3), intent(out), optional    :: virial
    real(dp), dimension(:), intent(out), optional      :: local_e
    real(dp), intent(in), optional                     :: cutoff
    integer, intent(out), optional                     :: error

    integer  :: i, j, n

    real(dp) :: my_cutoff, r_ij, de
    real(dp), dimension(3) :: force, u_ij

    real(dp), dimension(:), pointer :: my_charge

    type(Atoms), target :: my_at
    type(Atoms), pointer :: at => null()

    INIT_ERROR(error)
    call check_size('charge',charge,(/at_in%N/),'Ewald_corr_calc',error)

    my_cutoff = optional_default(at_in%cutoff,cutoff)

    if( present(cutoff) .and. (my_cutoff > at_in%cutoff) ) then
        my_at = at_in
        call set_cutoff(my_at,cutoff)
        call calc_connect(my_at)
        at => my_at
    else
        at => at_in
    endif
         
    if( present(e) ) e = 0.0_dp
    if( present(f) ) f = 0.0_dp
    if( present(virial) ) virial = 0.0_dp
    if( present(local_e) ) local_e = 0.0_dp

    do i = 1, at%N
       !Loop over neighbours
       do n = 1, n_neighbours(at,i)
          j = neighbour(at,i,n,distance=r_ij,cosines=u_ij) ! nth neighbour of atom i
          if( r_ij > my_cutoff )  cycle
           
          de = 0.5_dp * charge(i)*charge(j) / r_ij

          if( present(e) ) e = e + de
          if( present(local_e) ) local_e(i) = local_e(i) + de

          if( present(f) .or. present(virial) ) then
              force = - de / r_ij * u_ij

              if(present(f)) then
                 f(:,i) = f(:,i) + force
                 f(:,j) = f(:,j) - force
              endif

              if (present(virial)) virial = virial - (force .outer. u_ij) * r_ij
          endif
 
      enddo
    enddo
             
    if(present(e)) e = e * HARTREE*BOHR ! convert from internal units to eV
    if(present(f)) f = f * HARTREE*BOHR ! convert from internal units to eV/A
    if(present(virial)) virial = virial * HARTREE*BOHR

    if(associated(at,my_at)) call finalise(my_at)
    at => null()

  endsubroutine Direct_Coulomb_calc

  subroutine DSF_Coulomb_calc(at_in,charge, alpha, e, f, virial, local_e, e_potential, e_field, cutoff, error)

    type(Atoms), intent(in), target    :: at_in
    real(dp), dimension(:), intent(in) :: charge
    real(dp), intent(in)              :: alpha

    real(dp), intent(out), optional                    :: e
    real(dp), dimension(:,:), intent(out), optional    :: f
    real(dp), dimension(3,3), intent(out), optional    :: virial
    real(dp), dimension(:), intent(out), optional      :: local_e
    real(dp), dimension(:), intent(out), optional      :: e_potential
    real(dp), dimension(:,:), intent(out), optional    :: e_field
    real(dp), intent(in), optional                     :: cutoff
    integer, intent(out), optional                     :: error

    integer  :: i, j, n

    real(dp) :: my_cutoff, r_ij, phi_i, e_i, v_ij, dv_ij, v_cutoff, dv_cutoff, two_alpha_over_square_root_pi

    real(dp), dimension(3) :: u_ij, dphi_i, dphi_ij
    real(dp), dimension(3,3) :: dphi_ij_outer_r_ij

    type(Atoms), target :: my_at
    type(Atoms), pointer :: at => null()

    INIT_ERROR(error)
    call check_size('charge',charge,(/at_in%N/),'DSF_Coulomb_calc',error)

    my_cutoff = optional_default(at_in%cutoff,cutoff)

    two_alpha_over_square_root_pi = 2.0_dp * alpha / sqrt(pi)

    v_cutoff = erfc(alpha * my_cutoff) / my_cutoff
    dv_cutoff = ( v_cutoff + two_alpha_over_square_root_pi*exp(-(alpha*my_cutoff)**2) ) / my_cutoff

    if( present(cutoff) .and. (my_cutoff > at_in%cutoff) ) then
        my_at = at_in
        call set_cutoff(my_at,cutoff)
        call calc_connect(my_at)
        at => my_at
    else
        at => at_in
    endif
         
    if( present(e) ) e = 0.0_dp
    if( present(f) ) f = 0.0_dp
    if( present(virial) ) virial = 0.0_dp
    if( present(local_e) ) local_e = 0.0_dp
    if( present(e_potential) ) e_potential = 0.0_dp
    if( present(e_field) ) e_field = 0.0_dp

    do i = 1, at%N
       !Loop over neighbours

       phi_i = 0.0_dp
       dphi_i = 0.0_dp
       dphi_ij_outer_r_ij = 0.0_dp


       do n = 1, n_neighbours(at,i)
          j = neighbour(at, i, n, distance=r_ij, cosines=u_ij, max_dist=my_cutoff) ! nth neighbour of atom i
          if (j <= 0) cycle
          if (r_ij .feq. 0.0_dp) cycle
          
          v_ij = erfc(alpha*r_ij) / r_ij 
          phi_i = phi_i + charge(j) * ( v_ij - v_cutoff + dv_cutoff * (r_ij - my_cutoff) )

          if( present(f) .or. present(virial) .or. present(e_field) ) then
             dv_ij = ( v_ij + two_alpha_over_square_root_pi * exp(-(alpha*r_ij)**2) ) / r_ij
             dphi_ij = charge(j) * ( dv_ij - dv_cutoff ) * u_ij
             dphi_i = dphi_i + dphi_ij

             if(present(virial) ) dphi_ij_outer_r_ij = dphi_ij_outer_r_ij + ( dphi_ij .outer. u_ij ) * r_ij
          endif
       enddo
       
       if( present(e) .or. present(local_e) ) then
          e_i = 0.5_dp * phi_i * charge(i)
          if( present(e) ) e = e + e_i
          if( present(local_e) ) local_e(i) = e_i
       endif

       if(present(e_potential)) e_potential(i) = phi_i

       if(present(f)) f(:,i) = f(:,i) - dphi_i * charge(i)

       if (present(virial)) virial = virial + 0.5_dp * charge(i) * dphi_ij_outer_r_ij

       if(present(e_field)) e_field(:,i) = dphi_i
    enddo
             
    if(present(e)) e = e * HARTREE*BOHR ! convert from internal units to eV
    if(present(local_e)) local_e = local_e * HARTREE*BOHR ! convert from internal units to eV
    if(present(f)) f = f * HARTREE*BOHR ! convert from internal units to eV/A
    if(present(virial)) virial = virial * HARTREE*BOHR

    if(associated(at,my_at)) call finalise(my_at)
    at => null()

  endsubroutine DSF_Coulomb_calc

endmodule IPEwald_module
