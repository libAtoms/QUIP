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
!X  IPModel_SW_VP  module 
!X 
!%  Danny potential J ChemPhys 127, 204704(2007) for SiO_2. 
!%  
!%  Only the case of fourfold coordinated silcon and twofold coordinated
!%  is considered. (-->Fixed charges).    
!%
!%  The IPModel_SW_VP object contains all the parameters read 
!%  from an 'SW_VP_params' XML stanza.
!X
!X  Contribution from Anke Butenuth
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#include "error.inc"

module IPModel_SW_VP_module

use libAtoms_module

use mpi_context_module
use QUIP_Common_module
use Functions_module
implicit none

private

include 'IPModel_interface.h'

public :: IPModel_SW_VP
type IPModel_SW_VP
  integer :: n_types = 0         !% Number of atomic types 
  integer, allocatable :: atomic_num(:), type_of_atomic_num(:)  !% Atomic number dimensioned as \texttt{n_types} 

  real(dp) :: cutoff = 0.0_dp

  
  real(dp), allocatable :: a(:,:), AA(:,:), BB(:,:), p(:,:), q(:,:), sigma(:,:), eps2(:,:), C_0(:,:), C_1(:,:), b(:,:), D_SiO(:,:), C_OOprime(:,:), D_OOprime(:,:) !% IP parameters
  real(dp), allocatable :: lambda(:,:,:), d1(:,:,:), d2(:,:,:), gamma1(:,:,:), gamma2(:,:,:), eps3(:,:,:), costheta0(:,:,:) !% IP parameters


  character(len=STRING_LENGTH) :: label

end type IPModel_SW_VP

logical, private :: parse_in_ip, parse_matched_label
type(IPModel_SW_VP), private, pointer :: parse_ip

interface Initialise
  module procedure IPModel_SW_VP_Initialise_str
end interface Initialise

interface Finalise
  module procedure IPModel_SW_VP_Finalise
end interface Finalise

interface Print
  module procedure IPModel_SW_VP_Print
end interface Print

interface Calc
  module procedure IPModel_SW_VP_Calc
end interface Calc

contains

subroutine IPModel_SW_VP_Initialise_str(this, args_str, param_str)
  type(IPModel_SW_VP), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str

  type(Dictionary) :: params

  call Finalise(this)

  call initialise(params)
  this%label=''
  call param_register(params, 'label', '', this%label, help_string="No help yet.  This source file was $LastChangedBy: nb326 $")
  if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_SW_VP_Initialise_str args_str')) then
    call system_abort("IPModel_SW_VP_Initialise_str failed to parse label from args_str="//trim(args_str))
  endif
  call finalise(params)

  call IPModel_SW_VP_read_params_xml(this, param_str)

end subroutine IPModel_SW_VP_Initialise_str

subroutine IPModel_SW_VP_Finalise(this)
  type(IPModel_SW_VP), intent(inout) :: this

  if (allocated(this%atomic_num)) deallocate(this%atomic_num)
  if (allocated(this%type_of_atomic_num)) deallocate(this%type_of_atomic_num)

  if (allocated(this%a)) deallocate(this%a)
  if (allocated(this%AA)) deallocate(this%AA)
  if (allocated(this%BB)) deallocate(this%BB)
  if (allocated(this%p)) deallocate(this%p)
  if (allocated(this%q)) deallocate(this%q)
  if (allocated(this%sigma)) deallocate(this%sigma)
  if (allocated(this%eps2)) deallocate(this%eps2)
  if (allocated(this%lambda)) deallocate(this%lambda)
  if (allocated(this%gamma1)) deallocate(this%gamma1)
  if (allocated(this%gamma2)) deallocate(this%gamma2)
  if (allocated(this%eps3)) deallocate(this%eps3)
  if (allocated(this%costheta0)) deallocate(this%costheta0)
  if (allocated(this%d1)) deallocate(this%d1)
  if (allocated(this%d2)) deallocate(this%d2)
  if (allocated(this%C_0)) deallocate(this%C_0)
  if (allocated(this%C_1)) deallocate(this%C_1)
  if (allocated(this%b)) deallocate(this%b)
  if (allocated(this%D_SiO)) deallocate(this%D_SiO)
  if (allocated(this%C_OOprime)) deallocate(this%C_OOprime)
  if (allocated(this%D_OOprime)) deallocate(this%D_OOprime)

  this%n_types = 0
  this%label = ''
end subroutine IPModel_SW_VP_Finalise

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!% The potential calculator. It computes energy, forces and virial.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

subroutine IPModel_SW_VP_Calc(this, at, e, local_e, f, virial, local_virial, args_str, mpi, error)
  type(IPModel_SW_VP), intent(inout) :: this
  type(Atoms), intent(in) :: at
  real(dp), intent(out), optional :: e, local_e(:) !% \texttt{e} = System total energy, \texttt{local_e} = energy of each atom, vector dimensioned as \texttt{at%N}.  
  real(dp), intent(out), optional :: f(:,:), local_virial(:,:)   !% Forces, dimensioned as \texttt{f(3,at%N)}, local virials, dimensioned as \texttt{local_virial(9,at%N)} 
  real(dp), intent(out), optional :: virial(3,3)   !% Virial
  character(len=*), optional      :: args_str
  type(MPI_Context), intent(in), optional :: mpi
  integer, intent(out), optional :: error

  real(dp), pointer :: w_e(:)
  integer i, ji, j, ki, k
  real(dp) :: drij(3), drij_mag, drik(3), drik_mag, drij_dot_drik 
  real(dp) :: w_f

  integer ti, tj, tk

  real(dp) :: drij_dri(3), drij_drj(3), drik_dri(3), drik_drk(3) 
  real(dp) :: dcos_ijk_dri(3), dcos_ijk_drj(3), dcos_ijk_drk(3) 
  real(dp) :: virial_i(3,3) 

  real(dp) :: de, de_dr, de_drij, de_drik, de_dcos_ijk 
  real(dp) :: cur_cutoff
 

  integer :: n_neigh_i

  type(Dictionary)                :: params
  logical :: has_atom_mask_name, atom_mask_exclude_all
  character(STRING_LENGTH) :: atom_mask_name
  logical, dimension(:), pointer :: atom_mask_pointer

#ifdef _OPENMP
  real(dp) :: private_virial(3,3), private_e
  real(dp), allocatable :: private_f(:,:), private_local_e(:), private_local_virial(:,:)
#endif

  INIT_ERROR(error)

  call print("IPModel_SW_VP_Calc starting ", PRINT_ANAL)
  if (present(e)) e = 0.0_dp

  if (present(local_e)) then
     call check_size('Local_E',local_e,(/at%N/),'IPModel_SW_VP_Calc', error)
     local_e = 0.0_dp
  endif

  if (present(f)) then 
     call check_size('Force',f,(/3,at%N/),'IPModel_SW_VP_Calc', error)
     f = 0.0_dp
  end if

  if (present(virial)) virial = 0.0_dp

  if (present(local_virial)) then
     call check_size('Local_virial',local_virial,(/9,at%N/),'IPModel_GAP_Calc', error)
     local_virial = 0.0_dp
  endif

  atom_mask_pointer => null()
  if (present(args_str)) then
     call initialise(params)
     call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, help_string="Name of logical property used to mask atoms")
     call param_register(params, 'atom_mask_exclude_all', 'F', atom_mask_exclude_all, help_string="If true, exclude contributions made by atoms with atom_mask==.false. on their neighbour with atom_mask==.true., i.e. behave as if excluded atoms were not there at all")

     if(.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_SW_VP_Calc args_str')) then
        RAISE_ERROR("IPModel_SW_VP_Calc failed to parse args_str='"//trim(args_str)//"'",error)
     endif
     call finalise(params)
     if(has_atom_mask_name) then
        if (.not. assign_pointer(at, trim(atom_mask_name) , atom_mask_pointer)) &
        call system_abort("IPModel_SW_VP_Calc did not find "//trim(atom_mask_name)//" property in the atoms object.")
     else
        atom_mask_pointer => null()
     endif
  endif

  if (.not.assign_pointer(at,"weight", w_e)) nullify(w_e)

#ifdef _OPENMP
      !$omp parallel default(none) private(i, ji, j, ki, k, drij, drij_mag, drik, drik_mag, drij_dot_drik, w_f, ti, tj, tk, drij_dri, drij_drj, drik_dri, drik_drk, dcos_ijk_dri, dcos_ijk_drj, dcos_ijk_drk, de, de_dr, de_drij, de_drik, de_dcos_ijk, cur_cutoff, private_virial, private_e, private_f, private_local_e, private_local_virial, n_neigh_i, virial_i) shared(this,at,e,local_e,f,virial,local_virial,mpi,atom_mask_pointer,atom_mask_exclude_all,w_e)

  if (present(e)) private_e = 0.0_dp
  if (present(local_e)) then
    allocate(private_local_e(at%N))
    private_local_e = 0.0_dp
  endif
  if (present(f)) then
    allocate(private_f(3,at%N))
    private_f = 0.0_dp
  endif
  if (present(virial)) private_virial = 0.0_dp
  if (present(local_virial)) then
    allocate(private_local_virial(9,at%N))
    private_local_virial = 0.0_dp
  endif

!$omp do
#endif
  do i=1, at%N 
   if (present(mpi)) then
       if (mpi%active) then
         if (mod(i-1, mpi%n_procs) /= mpi%my_proc) cycle
       endif
    endif

    if(associated(atom_mask_pointer)) then
       if(.not. atom_mask_pointer(i)) cycle
    endif

    ti = get_type(this%type_of_atomic_num, at%Z(i)) 

    if (current_verbosity() >= PRINT_ANAL) call print ("IPModel_SW_VP_Calc i " // i // " " // atoms_n_neighbours(at,i), PRINT_ANAL)
    cur_cutoff = maxval(this%a(ti,:)*this%sigma(ti,:))   
    n_neigh_i = atoms_n_neighbours(at, i)


    do ji=1, n_neigh_i 
      j = atoms_neighbour(at, i, ji, drij_mag, cosines = drij, max_dist = cur_cutoff) 
      
      if(associated(atom_mask_pointer) .and. atom_mask_exclude_all) then
         if(.not. atom_mask_pointer(j)) cycle
      endif

      if (j <= 0) cycle  
      if (drij_mag .feq. 0.0_dp) cycle   

      tj = get_type(this%type_of_atomic_num, at%Z(j))
      
      if (drij_mag/this%sigma(ti,tj)  > this%a(ti,tj)) cycle 
      

      if (current_verbosity() >= PRINT_ANAL) call print ("IPModel_SW_VP_Calc i j " // i // " " // j, PRINT_ANAL)

      if (associated(w_e)) then
	w_f = 0.5_dp*(w_e(i)+w_e(j))  ! 
      else
	w_f = 1.0_dp
      endif

      if (present(e) .or. present(local_e)) then
	! factor of 0.5 because SW_VP definition goes over each pair only once
	
        ! select functional form of pair potentials depending on species : oxygen -oxygen
        if ((ti .eq. 2) .and. (tj  .eq. 2)) then  ! 
          de = 0.5_dp*this%eps2(ti,tj)*f2OO(this, drij_mag, ti, tj)

        elseif ((ti .eq. 2) .and. (tj .eq. 3)) then  ! Oxygen Silicon
          de = 0.5_dp*this%eps2(ti,tj)*f2SiO(this, drij_mag, ti, tj) 
       
 
        elseif ((ti .eq. 3) .and. (tj .eq.2 )) then ! Calls the Silicon-Oxygen function again  
          de = 0.5_dp*this%eps2(ti,tj)*f2SiO(this, drij_mag, ti, tj)

        elseif ((ti .eq. 3) .and. (tj .eq. 3)) then ! Silicon -Silicon
          de = 0.5_dp*this%eps2(ti,tj)*f2SiSi(this, drij_mag, ti, tj)  
       
        else 
          de= 0.0_dp
        endif
       

	
        if (present(local_e)) then
#ifdef _OPENMP
	  private_local_e(i) = private_local_e(i) + de
#else
	  local_e(i) = local_e(i) + de 
#endif
	endif
	if (present(e)) then
#ifdef _OPENMP
	  private_e = private_e + de*w_f
#else
	  e = e + de*w_f 
#endif
	endif
      endif

      if (present(f) .or. present(virial) .or. present(local_virial)) then

        
        !select pair potentials depending on species : oxygen-oxygen
        if ((ti .eq. 2) .and. (tj .eq. 2)) then  
        de_dr = 0.5_dp*this%eps2(ti,tj)*df2OO_dr(this, drij_mag, ti, tj) 
        
        elseif ((ti .eq. 3) .and. (tj .eq. 2)) then  ! oxygen-silicon
        de_dr = 0.5_dp*this%eps2(ti,tj)* df2SiO_dr(this, drij_mag, ti, tj) 

        elseif ((ti.eq. 2) .and. (tj .eq. 3 )) then ! Calls the silicon-oxygen function again  
        de_dr = 0.5_dp*this%eps2(ti,tj)* df2SiO_dr(this, drij_mag, ti, tj)

        elseif ((ti .eq. 3 ) .and. (tj .eq. 3)) then ! silicon-silicon
        de_dr = 0.5_dp*this%eps2(ti,tj)*df2SiSi_dr(this, drij_mag, ti, tj) 
        else 
        de_dr = 0.0_dp  
        endif
        
       

 
	if (present(f)) then
#ifdef _OPENMP
	  private_f(:,i) = private_f(:,i) + de_dr*w_f*drij
	  private_f(:,j) = private_f(:,j) - de_dr*w_f*drij
#else
	  f(:,i) = f(:,i) + de_dr*w_f*drij 
	  f(:,j) = f(:,j) - de_dr*w_f*drij 
#endif
	endif

        if(present(virial) .or. present(local_virial)) virial_i = de_dr*w_f*(drij .outer. drij)*drij_mag 



	if (present(virial)) then
#ifdef _OPENMP
	  private_virial = private_virial - virial_i
#else
	  virial = virial - virial_i
#endif
	endif
        if (present(local_virial)) then
#ifdef _OPENMP
           private_local_virial(:,i) = private_local_virial(:,i) - reshape(virial_i, (/9/)) 
#else
           local_virial(:,i) = local_virial(:,i) - reshape(virial_i, (/9/))
#endif
        endif

      endif

      if (associated(w_e)) then
	w_f = w_e(i)
      else
	w_f = 1.0_dp
      endif
      enddo !  closes j loop
      enddo !  closes i loop 
      
      ! New loop for three body terms, because (see parameters/ip.params.SW_VP.xml) cutoff for three body terms is in SiSiO case larger then the  SiSi cutoff

      !$omp do
      do i=1, at%N 

        if (present(mpi)) then
           if (mpi%active) then
             if (mod(i-1, mpi%n_procs) /= mpi%my_proc) cycle
           endif
        endif

        if(associated(atom_mask_pointer)) then
           if(.not. atom_mask_pointer(i)) cycle
        endif

    ti = get_type(this%type_of_atomic_num, at%Z(i)) ! 


   ! if (current_verbosity() >= PRINT_ANAL) call print ("IPModel_SW_VP_Calc i " // i // " " // atoms_n_neighbours(at,i), PRINT_ANAL)
    cur_cutoff = maxval(this%a(ti,:)*this%sigma(ti,:)) 
    n_neigh_i = atoms_n_neighbours(at, i)


      do ji=1, n_neigh_i 
      j = atoms_neighbour(at, i, ji, drij_mag, cosines = drij, max_dist = cur_cutoff) 

      if(associated(atom_mask_pointer) .and. atom_mask_exclude_all) then
         if(.not. atom_mask_pointer(j)) cycle
      endif

      if (j <= 0) cycle  
      if (drij_mag .feq. 0.0_dp) cycle 

      tj = get_type(this%type_of_atomic_num, at%Z(j)) 

    
        
      do ki=1, n_neigh_i 
         
	if (ki <= ji) cycle  
	
        k = atoms_neighbour(at, i, ki, drik_mag, cosines = drik, max_dist=cur_cutoff) 

        if(associated(atom_mask_pointer) .and. atom_mask_exclude_all) then
           if(.not. atom_mask_pointer(k)) cycle
        endif

	if (k <= 0) cycle
	if (drik_mag .feq. 0.0_dp) cycle 

	tk = get_type(this%type_of_atomic_num, at%Z(k))
        if (drik_mag > this%d2(ti, tj, tk)) cycle 
        if (drij_mag > this%d1(ti, tj, tk)) cycle

	drij_dot_drik = sum(drij*drik) 
	if (present(e) .or. present(local_e)) then
	  de = this%eps3(ti,tj,tk)*f3(this, drij_mag, drik_mag, drij_dot_drik, ti, tj, tk) 
        
      
	  if (present(local_e)) then
#ifdef _OPENMP
	    private_local_e(i) = private_local_e(i) + de
#else
	    local_e(i) = local_e(i) + de
#endif
	  endif
	  if (present(e)) then
#ifdef _OPENMP
	    private_e = private_e + de*w_f 
#else
	    e = e + de*w_f 
#endif
	  endif
	endif
	if (present(f) .or. present(virial)) then
	  call df3_dr(this, drij_mag, drik_mag, drij_dot_drik, ti, tj, tk, de_drij, de_drik, de_dcos_ijk)  
	  drij_dri = drij 
	  drij_drj = -drij !grad.|vec(r_i) -vec(r_j)|= -( vec(r_i) -vec(r_j))/|vec(r_i) -vec(r_j)| = -drij 
	  drik_dri = drik
	  drik_drk = -drik

! dcos_ijk = drij_dot_drik
! (ri_x - rj_x)*(ri_x - rk_x) + (ri_y - rj_y)*(ri_y-rk_y) / (rij rik)
!
! d/dri
! 
! ((rij rik)(rij_x + rik_x) - (rij . rik)(rij rik_x / rik + rik rij_x/rij)) / (rij rik)^2
! 
! rij rik ( rij_x + rik_x) -  (rij.rik)(rij rik_x / rik + rik rij_x/rij)
! ---------------------------------------------------------------------
!                       rij^2 rik^2
!                       
! rij_x + rik_x     cos_ijk (rij rik_x / rik + rik rij_x /rij)
! -------------  - ------------------------------------------
! rij rik                    rij rik
! 
! rij_x + rik_x
! ------------- - cos_ijk (rik_x / rik^2 + rij_x / rij^2)
!    rij rik
!    
!    
! rhatij_x/rik + rhatik_x/rij - cos_ijk(rhatik_x/rik + rhatij_x/rij)
! 
! 
! d/drj
! 
! ((rij rik)(-rik_x) - (rij.rik)(rik (-rij_x)/rij)) / (rij rik)^2
! 
!  -rik_x     (rij.rik) ( rik rij_x/rij)
! -------  + ---------------------------
! rij rik           rij^2 rik^2
! 
! -rhatik_x/rij + cos_ijk rhatij_x/rij
!
! d/drij
!
! (ri_x - rj_x)*(ri_x - rk_x) + (ri_y - rj_y)*(ri_y-rk_y) / (rij rik)
!
! (rij rik) rik_x - (rij . rik) (rij_x/rij rik)
! ---------------------------------------------
!                    (rij rik)^2
!
! rik_x      (rij.rik) rij_x 
! ------   - ---------------
! rij rik    rij^3 rik
!
! rhatik_x      (rhatij . rhatik) rhatij_x
! --------   -  ------------------------
! rij           rij
!
! rhatik_x      cos_ijk * rhatij_x
! --------   -  ------------------
!   rij                 rij

	  dcos_ijk_dri = drij/drik_mag + drik/drij_mag - drij_dot_drik * (drik/drik_mag + drij/drij_mag) 
	  dcos_ijk_drj = -drik/drij_mag + drij_dot_drik * drij/drij_mag 
	  dcos_ijk_drk = -drij/drik_mag + drij_dot_drik * drik/drik_mag 




	  if (present(f)) then
#ifdef _OPENMP
	    private_f(:,i) = private_f(:,i) + w_f*this%eps3(ti,tj,tk)*(de_drij*drij_dri(:) + de_drik*drik_dri(:) + &
						       de_dcos_ijk * dcos_ijk_dri(:))
	    private_f(:,j) = private_f(:,j) + w_f*this%eps3(ti,tj,tk)*(de_drij*drij_drj(:) + de_dcos_ijk*dcos_ijk_drj(:))
	    private_f(:,k) = private_f(:,k) + w_f*this%eps3(ti,tj,tk)*(de_drik*drik_drk(:) + de_dcos_ijk*dcos_ijk_drk(:)) 
#else
	    f(:,i) = f(:,i) + w_f*this%eps3(ti,tj,tk)*(de_drij*drij_dri(:) + de_drik*drik_dri(:) + &
						       de_dcos_ijk * dcos_ijk_dri(:)) 
	    f(:,j) = f(:,j) + w_f*this%eps3(ti,tj,tk)*(de_drij*drij_drj(:) + de_dcos_ijk*dcos_ijk_drj(:))  
	    f(:,k) = f(:,k) + w_f*this%eps3(ti,tj,tk)*(de_drik*drik_drk(:) + de_dcos_ijk*dcos_ijk_drk(:)) 
#endif
	  end if

          if( present(virial) .or. present(local_virial) ) virial_i = &
              w_f*this%eps3(ti,tj,tk)*( &
	      de_drij*(drij .outer. drij)*drij_mag + de_drik*(drik .outer. drik)*drik_mag + &
	      de_dcos_ijk * ((drik .outer. drij) - drij_dot_drik * (drij .outer. drij)) + &
	      de_dcos_ijk * ((drij .outer. drik) - drij_dot_drik * (drik .outer. drik)) ) 

	  if (present(virial)) then
#ifdef _OPENMP
	    private_virial = private_virial - virial_i
#else
	    virial = virial - virial_i
#endif
	  end if
	  if (present(local_virial)) then
#ifdef _OPENMP
	    private_local_virial(:,i) = private_local_virial(:,i) - reshape(virial_i,(/9/))
#else
	    local_virial(:,i) = local_virial(:,i) - reshape(virial_i,(/9/))
#endif
	  end if
	endif

      end do 

    end do 
  end do 

#ifdef _OPENMP
!$omp critical
  if (present(e)) e = e + private_e
  if (present(f)) f = f + private_f
  if (present(local_e)) local_e = local_e + private_local_e
  if (present(virial)) virial = virial + private_virial
  if (present(local_virial)) local_virial = local_virial + private_local_virial
!$omp end critical 

  if(allocated(private_f)) deallocate(private_f)
  if(allocated(private_local_e)) deallocate(private_local_e)
  if(allocated(private_local_virial)) deallocate(private_local_virial)

!$omp end parallel
#endif

  if (present(mpi)) then
     if (present(e)) e = sum(mpi, e) 
     if (present(local_e)) call sum_in_place(mpi, local_e)
     if (present(f)) call sum_in_place(mpi, f)
     if (present(virial)) call sum_in_place(mpi, virial)
     if (present(local_virial)) call sum_in_place(mpi, local_virial)
  endif

end subroutine IPModel_SW_VP_Calc




function f3(this, rij_i, rik_i, cos_ijk, ti, tj, tk) 
  type(IPModel_SW_VP), intent(in) :: this
  real(dp), intent(in) :: rij_i, rik_i, cos_ijk
  integer, intent(in) :: ti, tj, tk
  real(dp) :: f3

  real(dp) :: rij, rik
  !We give all cosines as parameters from xml file

  rij = rij_i 
  rik = rik_i 

  if (rij_i >= this%d1(ti,tj,tk) .or. rik_i >= this%d2(ti,tj,tk)) then 
     f3 = 0.0_dp
    return
  endif
  
  f3 = this%lambda(ti,tj,tk)*exp( this%gamma1(ti,tj,tk)/(rij-this%d1(ti,tj,tk)) + this%gamma2(ti,tj,tk)/(rik-this%d2(ti,tj,tk)) )* &
     (cos_ijk - this%costheta0(ti,tj,tk))**2

 

end function f3

subroutine df3_dr(this, rij_i, rik_i, cos_ijk, ti, tj, tk, de_drij, de_drik, de_dcos_ijk) 
  type(IPModel_SW_VP), intent(in) :: this
  real(dp), intent(in) :: rij_i, rik_i, cos_ijk
  integer, intent(in) :: ti, tj, tk
  real(dp), intent(out) :: de_drij, de_drik, de_dcos_ijk

  real(dp) :: rij, rik 
  real(dp) :: expf, cosf


  rij = rij_i
  rik = rik_i

  if (rij_i >= this%d1(ti,tj,tk) .or. rik_i >= this%d2(ti,tj,tk)) then ! Anke_new
    de_drij = 0.0_dp
    de_drik = 0.0_dp
    de_dcos_ijk = 0.0_dp
    return
  endif 

  
  expf = this%lambda(ti,tj,tk)*exp(this%gamma1(ti,tj,tk)/(rij - this%d1(ti,tj,tk)) + this%gamma2(ti,tj,tk)/(rik - this%d2(ti,tj,tk))) 
  cosf = (cos_ijk -this%costheta0(ti,tj,tk)) 

  de_drij = -expf * this%gamma1(ti,tj,tk)/(rij-this%d1(ti,tj,tk))**2 * cosf**2    
  de_drik = -expf * this%gamma2(ti,tj,tk)/(rik-this%d2(ti,tj,tk))**2 * cosf**2    
  de_dcos_ijk = expf * 2.0_dp * cosf  





end subroutine df3_dr




!************Begin: Anke modification**********




function f2SiSi(this, ri, ti, tj) 
  type(IPModel_SW_VP), intent(in) :: this
  real(dp), intent(in) :: ri
  integer, intent(in) :: ti, tj
  
  real(dp) :: f2SiSi


!  if (ri >= this%a(ti,tj)) then
!    f2SiSi = 0.0_dp
!    return
!  endif

!  f2SiSi = this%AA(ti,tj)*(1.0_dp+3.2_dp*qi*qj)*faqj*faqi*(this%BB(ti,tj) * r**(-this%p(ti,tj)) - r**(-this%q(ti,tj))) * exp(1.0_dp/(r-this%a(ti,tj)))
! for now, fixed coordination. If Si is fourfold coordinated, in Si0_2 is faq(Si) = 0.
! 
! Potential has to be modified if one wants Si-Si interactions as well
 
  f2SiSi = 0.0_dp

end function f2SiSi

function df2SiSi_dr(this, ri, ti, tj) 
  type(IPModel_SW_VP), intent(in) :: this
  real(dp), intent(in) :: ri
  integer, intent(in) :: ti, tj
  real(dp) :: df2SiSi_dr

  real(dp) :: expf, dexpf_dr, powf, dpowf_dr


 ! if (r >= this%a(ti,tj)) then
 !   df2SiSi_dr = 0.0_dp

 !   return
 ! endif


 ! expf = exp(this%sigma(ti,tj)/(r-this%a(ti,tj)))
 ! dexpf_dr = -expf*this%sigma(ti,tj)/(r-this%a(ti,tj))**2

 ! powf = this%BB(ti,tj)*(r**(-this%p(ti,tj))) - r**(-this%q(ti,tj))
 ! dpowf_dr = -this%p(ti,tj)*this%BB(ti,tj)*(r**(-this%p(ti,tj)-1)) + this%q(ti,tj)*r**(-this%q(ti,tj)-1)

   df2SiSi_dr = 0.0_dp

end function df2SiSi_dr 


function f2SiO(this, r, ti, tj)
  type(IPModel_SW_VP), intent(in) :: this
  real(dp), intent(in) :: r
  integer, intent(in) :: ti, tj
  real(dp) :: f2SiO
   
  real(dp) :: expf
  real(dp) :: f2SiO_at_cutoff, df2SiO_dr_at_cutoff

! Cutoff sigma*a for longrange non-coulombic terms
  
  if (r >= (this%sigma(ti,tj)*this%a(ti,tj))) then
    f2SiO = 0.0_dp
    return
  endif
 


  expf = exp(-r/this%a(ti,tj))
  f2SiO_at_cutoff = (this%C_0(ti,tj)-this%C_1(ti,tj)*1.6_dp)*(this%sigma(ti,tj)*this%a(ti,tj))**(-9.0_dp)-this%D_SiO(ti,tj)*(this%sigma(ti,tj)*this%a(ti,tj))**(-4.0_dp)*exp(-this%sigma(ti,tj))
  df2SiO_dr_at_cutoff= -9.0_dp*(this%C_0(ti, tj)-this%C_1(ti,tj)*1.6_dp)*(this%sigma(ti,tj)*this%a(ti,tj)) **(-10.0_dp) &
                         + 4.0_dp*this%D_SiO(ti,tj)*(this%sigma(ti,tj)*this%a(ti,tj)) **(-5.0_dp)* exp(-this%sigma(ti,tj)) &
                         + this%D_SiO(ti,tj)*(this%sigma(ti,tj)*this%a(ti,tj) )**(-4.0_dp)*exp(-this%sigma(ti,tj))/this%a(ti,tj)
  


  f2SiO = (this%C_0(ti,tj)-this%C_1(ti,tj)*1.6_dp)*r**(-9.0_dp) - this%D_SiO(ti,tj)*r**(-4.0_dp)*expf  - f2SiO_at_cutoff - (r - this%sigma(ti,tj)*this%a(ti,tj))*df2SiO_dr_at_cutoff



end function f2SiO



function df2SiO_dr(this, r, ti, tj)
  type(IPModel_SW_VP), intent(in) :: this
  real(dp), intent(in) :: r
  integer, intent(in) :: ti, tj
  real(dp) :: df2SiO_dr
  real(dp) :: expf
  real(dp) :: df2SiO_dr_at_cutoff
    
  if (r >= (this%sigma(ti,tj)*this%a(ti,tj))) then
    df2SiO_dr = 0.0_dp
    return
  endif



  expf = exp(-r/this%a(ti,tj))
  df2SiO_dr_at_cutoff= -9.0_dp*(this%C_0(ti, tj)-this%C_1(ti,tj)*1.6_dp)*(this%sigma(ti,tj)*this%a(ti,tj)) **(-10.0_dp) &
                         + 4.0_dp*this%D_SiO(ti,tj)*(this%sigma(ti,tj)*this%a(ti,tj)) **(-5.0_dp)* exp(-this%sigma(ti,tj)) &
                         + this%D_SiO(ti,tj)*(this%sigma(ti,tj)*this%a(ti,tj) )**(-4.0_dp)*exp(-this%sigma(ti,tj))/this%a(ti,tj)


  df2SiO_dr = -9.0_dp*(this%C_0(ti, tj)-this%C_1(ti,tj)*1.6_dp)*r**(-10.0_dp) + 4.0_dp*this%D_SiO(ti,tj)*r**(-5.0_dp)*expf + this%D_SiO(ti,tj)*r**(-4.0_dp)*expf/this%a(ti,tj) - df2SiO_dr_at_cutoff



end function df2SiO_dr


function f2OO(this, r, ti, tj)
  type(IPModel_SW_VP), intent(in) :: this
  real(dp), intent(in) :: r
  integer, intent(in) :: ti, tj
  real(dp) :: f2OO
  real(dp) :: expf
  real(dp) :: f2OO_at_cutoff, df2OO_dr_at_cutoff

  if (r >= (this%sigma(ti,tj)*this%a(ti,tj))) then
    f2OO = 0.0_dp
    return
  endif


  expf = exp(-r/this%a(ti,tj))
  f2OO_at_cutoff= (this%C_OOprime(ti, tj))*(this%sigma(ti,tj)*this%a(ti,tj))**(-7.0_dp) - this%D_OOprime(ti, tj)*(this%sigma(ti,tj)*this%a(ti,tj))**(-4.0_dp)*exp(-this%sigma(ti,tj))
 ! Cutoff is at r = sigma*a, thus can be altered via parameter list 
  df2OO_dr_at_cutoff=  -7.0_dp*(this%C_OOprime(ti,tj))*(this%sigma(ti,tj)*this%a(ti,tj) )**(-8.0_dp) + 4.0_dp*this%D_OOprime(ti,tj)*(this%sigma(ti,tj)*this%a(ti,tj) )**(-5.0_dp)* exp(-this%sigma(ti,tj)) &
                       +  this%D_OOprime(ti,tj)*(this%sigma(ti,tj)*this%a(ti,tj) )**(-4.0_dp)*exp(-this%sigma(ti,tj))/this%a(ti,tj)

  f2OO = (this%C_OOprime(ti, tj))*r**(-7.0_dp) - this%D_OOprime(ti, tj)*r**(-4.0_dp)*expf  - f2OO_at_cutoff &
        - (r - this%sigma(ti,tj)*this%a(ti,tj))*df2OO_dr_at_cutoff
 

end function f2OO


function df2OO_dr(this, r, ti, tj)
  type(IPModel_SW_VP), intent(in) :: this
  real(dp), intent(in) :: r
  integer, intent(in) :: ti, tj
  real(dp) :: df2OO_dr, df2OO_dr_at_cutoff 

  real(dp) :: expf
 
  if (r >= (this%sigma(ti,tj)*this%a(ti,tj))) then
   df2OO_dr = 0.0_dp
    return
  endif


  expf = exp(-r/this%a(ti,tj))

  df2OO_dr_at_cutoff = -7.0_dp*(this%C_OOprime(ti,tj))*(this%sigma(ti,tj)*this%a(ti,tj) )**(-8.0_dp) + 4.0_dp*this%D_OOprime(ti,tj)*(this%sigma(ti,tj)*this%a(ti,tj) )**(-5.0_dp)* exp(-this%sigma(ti,tj)) &
                       +  this%D_OOprime(ti,tj)*(this%sigma(ti,tj)*this%a(ti,tj) )**(-4.0_dp)*exp(-this%sigma(ti,tj))/this%a(ti,tj)

  df2OO_dr = -7.0_dp*(this%C_OOprime(ti,tj))*r**(-8.0_dp) + 4.0_dp*this%D_OOprime(ti,tj)*r**(-5.0_dp)*expf +  this%D_OOprime(ti,tj)*r**(-4.0_dp)*expf/this%a(ti,tj)- df2OO_dr_at_cutoff 



end function df2OO_dr




!*************End: Anke modification**********




!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!% XML param reader functions.
!% An example for the input file.
!%  
!% <SW_VP_potential>
!%  <Potential label="SW_VP_Potential" init_args="Sum init_args_pot1={IP SW_VP} init_args_pot2={IP Coulomb}"/>
!%
!%   <Coulomb_params n_types="3" cutoff="6.0" method="ewald" ewald_error="1.0e-3" label="default_ewald">
!%   <!-- <per_type_data type="1" atomic_num="1" charge="1.0"/>
!%      <per_type_data type="2" atomic_num="8" charge="-0.8"/>
!%      <per_type_data type="3" atomic_num="14" charge="1.6"/> charge can be directly read from extended xyz file -->
!%   </Coulomb_params>
!%
!%
!%  <SW_VP_params n_types="3" label="PRB_31_plus_H">
!%  <comment> Danny potential J ChemPhys 127, 204704(2007).The potential has to be called as a pot_init_args='Sum init_args_pot1={IP SW_VP} init_args_pot2={IP Coulomb}'</comment>
!%  <per_type_data type="1" atomic_num="1" />
!%  <per_type_data type="2" atomic_num="8" />
!%  <per_type_data type="3" atomic_num= "14"/>
!%  <per_pair_data atnum_i="1" atnum_j="1" AA="0.0" BB="0.0"
!%        p="0" q="0" a="1.0" sigma="1.0" eps="0.0" C_0="0.0" C_1="0.0" D_SiO="0" C_OOprime="0" D_OOprime="0" />
!%  <per_pair_data atnum_i="1" atnum_j="8" AA="0.0" BB="0.0"
!%        p="0" q="0" a="1.0" sigma="1.0" eps="0.0" C_0="0.0" C_1="0.0" D_SiO="0" C_OOprime="0" D_OOprime="0" />
!%  <per_pair_data atnum_i="1" atnum_j="14" AA="8.581214" BB="0.0327827"
!%        p="4" q="0" a="1.25" sigma="2.537884" eps="0.000" C_0="0.0" C_1="0.0" D_SiO="0" C_OOprime="0" D_OOprime="0"/>
!%  <per_pair_data atnum_i="8" atnum_j="8" AA="0.0" BB="0.0"
!%        p="0" q="0" a="4.43" sigma="1.24" eps="14.39" C_0="0.0" C_1="0.0" D_SiO="0" C_OOprime="51.692" D_OOprime="1.536"/>
!%  <per_pair_data atnum_i="8" atnum_j="14" AA="0.0" BB="0.0"
!%        p="0" q="0" a="4.43" sigma="1.24" eps="14.39" C_0="14.871" C_1="2.178" D_SiO="3.072" C_OOprime="0.0" D_OOprime="0.0"/>
!%  <per_pair_data atnum_i="14" atnum_j="14" AA="7.049556277" BB="0.6022245584"
!%        p="4" q="0" a="1.80" sigma="2.0951" eps="0" C_0="0.0" C_1="0.0" D_SiO="0" C_OOprime="0.0" D_OOprime="0.0"/>
!%
!%  <!-- triplet terms: atnum_c is the center atom, neighbours j and k -->
!%  <per_triplet_data atnum_c="1"  atnum_j="1"  atnum_k="1"
!%        lambda="21.0" gamma1="1.20" gamma2="1.20" d1="1.80" d2="1.80" costheta0="-0.333" eps="0.0" />
!%  <per_triplet_data atnum_c="1"  atnum_j="1"  atnum_k="8"
!%        lambda="21.0" gamma1="1.20" gamma2="1.20" d1="1.80" d2="1.80" costheta0="-0.333" eps="0.0" />
!%  <per_triplet_data atnum_c="1"  atnum_j="1"  atnum_k="14"
!%        lambda="21.0" gamma1="1.20" gamma2="1.20" d1="1.80" d2="1.80" costheta0="-0.333" eps="0.0" />
!%  <per_triplet_data atnum_c="1"  atnum_j="8"  atnum_k="8"
!%        lambda="21.0" gamma1="1.20" gamma2="1.20" d1="1.80" d2="1.80" costheta0="-0.333" eps="0.0" />
!%  <per_triplet_data atnum_c="1"  atnum_j="8"  atnum_k="14"
!%        lambda="21.0" gamma1="1.20" gamma2="1.20" d1="1.80" d2="1.80" costheta0="-0.333" eps="0.0" />
!%  <per_triplet_data atnum_c="1"  atnum_j="14" atnum_k="14"
!%        lambda="21.0" gamma1="1.20" gamma2="1.20" d1="1.80" d2="1.80" costheta0="-0.333" eps="0.0" />
!%  <per_triplet_data atnum_c="8" atnum_j="1" atnum_k="1"
!%        lambda="3.46" gamma1="0.48" gamma2="0.48" d1="1.24" d2="1.24"  costheta0="-0.5373"  eps="0.0" />
!%  <per_triplet_data atnum_c="8" atnum_j="1" atnum_k="8"
!%        lambda="3.46" gamma1="0.48" gamma2="0.48" d1="1.24" d2="1.24"  costheta0="-0.5373"  eps="0.0" />
!%  <per_triplet_data atnum_c="8" atnum_j="1" atnum_k="14"
!%        lambda="0.521" gamma1="1" gamma2="1" d1="2.6" d2="2.6"  costheta0="-0.5373"  eps="14.39" />
!%  <per_triplet_data atnum_c="8" atnum_j="8" atnum_k="8"
!%        lambda="3.46" gamma1="0.48" gamma2="0.48" d1="1.24" d2="1.24"  costheta0="-0.5373"  eps="0.0" />
!%  <per_triplet_data atnum_c="8" atnum_j="8" atnum_k="14"
!%        lambda="3.46" gamma1="0.48" gamma2="0.48" d1="1.24" d2="1.24"  costheta0="-0.5373"  eps="0.0" />
!%  <per_triplet_data atnum_c="8" atnum_j="14" atnum_k="14"
!%        lambda="1.4" gamma1="1" gamma2="1" d1="2.6" d2="2.6"  costheta0="-0.777"  eps="14.39" />
!%  <per_triplet_data atnum_c="14" atnum_j="1"  atnum_k="1"
!%        lambda="21.0" gamma1="1.20" gamma2="1.20" d1="1.80" d2="1.80" costheta0="-0.333" eps="0.0" />
!%  <per_triplet_data atnum_c="14" atnum_j="1"  atnum_k="8"
!%        lambda="21.0" gamma1="1.20" gamma2="1.20" d1="1.80" d2="1.80" costheta0="-0.333" eps="0.0" />
!%  <per_triplet_data atnum_c="14" atnum_j="1"  atnum_k="14"
!%        lambda="21.0" gamma1="1.20" gamma2="1.20" d1="1.80" d2="1.80" costheta0="-0.333" eps="0.0" />
!%  <per_triplet_data atnum_c="14" atnum_j="8"  atnum_k="8"
!%        lambda="0.350" gamma1="1" gamma2="1" d1="2.6" d2="2.6" costheta0="-0.333" eps="14.39" />
!%  <!--In original paper, the c=14,j=14,k=8 term is listed, hence d2 and d1 and gamma1 and gamma2 are interchanged  -->
!%  <per_triplet_data atnum_c="14" atnum_j="8" atnum_k="14"
!%        lambda="3.164" gamma1="4.06" gamma2="0.52" d1="3.981" d2="2.933"  costheta0="-0.333" eps="14.39" /> 
!%  <per_triplet_data atnum_c="14" atnum_j="14" atnum_k="14"
!%        lambda="3.164" gamma1="2.51" gamma2="2.51" d1="3.771" d2="3.771" costheta0="-0.333"  eps="14.39" />
!%  </SW_VP_params>
!%
!%  <constraints N="1">
!% <!-- All OH-bond lengths have to be constrained to avoid coulomb clash between OH for surfaces-->
!%  </constraints>
!% </SW_VP_potential>
!%
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

subroutine IPModel_startElement_handler(URI, localname, name, attributes)
  character(len=*), intent(in)   :: URI  
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name 
  type(dictionary_t), intent(in) :: attributes

  integer status
  character(len=1024) :: value
  integer ti, tj, tk, Zi, Zj, Zk

  if (name == 'SW_VP_params') then ! new SW_VP stanza
    call print("startElement_handler SW_VP_params", PRINT_NERD)
    if (parse_in_ip) &
      call system_abort("IPModel_startElement_handler entered SW_VP_params with parse_in true. Probably a bug in FoX (4.0.1, e.g.)")

    if (parse_matched_label) then 
       call print("SW_VP_params startElement_handler bailing because we already matched our label", PRINT_NERD)
       return ! we already found an exact match for this label
    end if

    call QUIP_FoX_get_value(attributes, 'label', value, status)
    if (status /= 0) value = ''

    call print("SW_VP_params startElement_handler found xml label '"//trim(value)//"'", PRINT_NERD)

    if (len(trim(parse_ip%label)) > 0) then ! we were passed in a label
      call print("SW_VP_params startElement_handler was passed in label '"//trim(parse_ip%label)//"'", PRINT_NERD)
      if (value == parse_ip%label) then ! exact match
        parse_matched_label = .true.
        parse_in_ip = .true.
      else ! no match
	call print("SW_VP_params startElement_handler got label didn't match", PRINT_NERD)
        parse_in_ip = .false.
      endif
    else ! no label passed in
      call print("SW_VP_params startElement_handler was not passed in a label", PRINT_NERD)
      parse_in_ip = .true.
    endif

    if (parse_in_ip) then
      if (parse_ip%n_types /= 0) then
	call print("SW_VP_params startElement_handler finalising old data, restarting to parse new section", PRINT_NERD)
        call finalise(parse_ip)
      endif

      call QUIP_FoX_get_value(attributes, "n_types", value, status)
      if (status /= 0) call system_abort ("IPModel_SW_VP_read_params_xml cannot find n_types")
      read (value, *) parse_ip%n_types

      allocate(parse_ip%atomic_num(parse_ip%n_types))
      parse_ip%atomic_num = 0
      allocate(parse_ip%a(parse_ip%n_types,parse_ip%n_types))
      parse_ip%a = 0.0_dp
      allocate(parse_ip%AA(parse_ip%n_types,parse_ip%n_types))
      allocate(parse_ip%BB(parse_ip%n_types,parse_ip%n_types))
      allocate(parse_ip%p(parse_ip%n_types,parse_ip%n_types))
      allocate(parse_ip%q(parse_ip%n_types,parse_ip%n_types))
      allocate(parse_ip%sigma(parse_ip%n_types,parse_ip%n_types))
      parse_ip%sigma = 0.0_dp
      allocate(parse_ip%eps2(parse_ip%n_types,parse_ip%n_types))
      allocate(parse_ip%C_0(parse_ip%n_types,parse_ip%n_types))
      allocate(parse_ip%C_1(parse_ip%n_types,parse_ip%n_types))
      allocate(parse_ip%D_SiO(parse_ip%n_types,parse_ip%n_types))
      allocate(parse_ip%C_OOprime(parse_ip%n_types,parse_ip%n_types))
      allocate(parse_ip%D_OOprime(parse_ip%n_types,parse_ip%n_types))
      
      allocate(parse_ip%lambda(parse_ip%n_types,parse_ip%n_types,parse_ip%n_types))
      allocate(parse_ip%gamma1(parse_ip%n_types,parse_ip%n_types,parse_ip%n_types))
      allocate(parse_ip%gamma2(parse_ip%n_types,parse_ip%n_types,parse_ip%n_types))
      allocate(parse_ip%eps3(parse_ip%n_types,parse_ip%n_types,parse_ip%n_types))
      allocate(parse_ip%d1(parse_ip%n_types,parse_ip%n_types,parse_ip%n_types))
      allocate(parse_ip%d2(parse_ip%n_types,parse_ip%n_types,parse_ip%n_types))
      allocate(parse_ip%costheta0(parse_ip%n_types,parse_ip%n_types,parse_ip%n_types))
    endif

  elseif (parse_in_ip .and. name == 'per_type_data') then

    call QUIP_FoX_get_value(attributes, "type", value, status)
    if (status /= 0) call system_abort ("IPModel_SW_VP_read_params_xml cannot find type")
    read (value, *) ti

    call QUIP_FoX_get_value(attributes, "atomic_num", value, status)
    if (status /= 0) call system_abort ("IPModel_SW_VP_read_params_xml cannot find atomic_num")
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
    if (status /= 0) call system_abort ("IPModel_SW_VP_read_params_xml cannot find atnum_i")
    read (value, *) Zi
    call QUIP_FoX_get_value(attributes, "atnum_j", value, status)
    if (status /= 0) call system_abort ("IPModel_SW_VP_read_params_xml cannot find atnum_j")
    read (value, *) Zj

    ti = get_type(parse_ip%type_of_atomic_num,Zi)
    tj = get_type(parse_ip%type_of_atomic_num,Zj)

    call QUIP_FoX_get_value(attributes, "AA", value, status)
    if (status /= 0) call system_abort ("IPModel_SW_VP_read_params_xml cannot find AA")
    read (value, *) parse_ip%AA(ti,tj)
    call QUIP_FoX_get_value(attributes, "BB", value, status)
    if (status /= 0) call system_abort ("IPModel_SW_VP_read_params_xml cannot find BB")
    read (value, *) parse_ip%BB(ti,tj)
    call QUIP_FoX_get_value(attributes, "p", value, status)
    if (status /= 0) call system_abort ("IPModel_SW_VP_read_params_xml cannot find p")
    read (value, *) parse_ip%p(ti,tj)
    call QUIP_FoX_get_value(attributes, "q", value, status)
    if (status /= 0) call system_abort ("IPModel_SW_VP_read_params_xml cannot find q")
    read (value, *) parse_ip%q(ti,tj)
    call QUIP_FoX_get_value(attributes, "a", value, status)
    if (status /= 0) call system_abort ("IPModel_SW_VP_read_params_xml cannot find a")
    read (value, *) parse_ip%a(ti,tj)
    call QUIP_FoX_get_value(attributes, "sigma", value, status)
    if (status /= 0) call system_abort ("IPModel_SW_VP_read_params_xml cannot find sigma")
    read (value, *) parse_ip%sigma(ti,tj)
    call QUIP_FoX_get_value(attributes, "eps", value, status)
    if (status /= 0) call system_abort ("IPModel_SW_VP_read_params_xml cannot find eps")
    read (value, *) parse_ip%eps2(ti,tj)
    call QUIP_FoX_get_value(attributes, "C_0", value, status)
    if (status /= 0) call system_abort ("IPModel_SW_VP_read_params_xml cannot find C_0")
    read (value, *) parse_ip%C_0(ti,tj)
    call QUIP_FoX_get_value(attributes, "C_1", value, status)
    if (status /= 0) call system_abort ("IPModel_SW_VP_read_params_xml cannot find C_1")
    read (value, *) parse_ip%C_1(ti,tj)
    call QUIP_FoX_get_value(attributes, "D_SiO", value, status)
    if (status /= 0) call system_abort ("IPModel_SW_VP_read_params_xml cannot find D_SiO")
    read (value, *) parse_ip%D_SiO(ti,tj)
    call QUIP_FoX_get_value(attributes, "C_OOprime", value, status)
    if (status /= 0) call system_abort ("IPModel_SW_VP_read_params_xml cannot find C_OOprime")
    read (value, *) parse_ip%C_OOprime(ti,tj)
    call QUIP_FoX_get_value(attributes, "D_OOprime", value, status)
    if (status /= 0) call system_abort ("IPModel_SW_VP_read_params_xml cannot find D_OOprime")
    read (value, *) parse_ip%D_OOprime(ti,tj)
    
    if (ti /= tj) then
      parse_ip%AA(tj,ti) = parse_ip%AA(ti,tj)
      parse_ip%BB(tj,ti) = parse_ip%BB(ti,tj)
      parse_ip%p(tj,ti) = parse_ip%p(ti,tj)
      parse_ip%q(tj,ti) = parse_ip%q(ti,tj)
      parse_ip%a(tj,ti) = parse_ip%a(ti,tj)
      parse_ip%sigma(tj,ti) = parse_ip%sigma(ti,tj)
      parse_ip%eps2(tj,ti) = parse_ip%eps2(ti,tj)
      parse_ip%C_0(tj,ti) = parse_ip%C_0(ti,tj)
      parse_ip%C_1(tj,ti) = parse_ip%C_1(ti,tj)
      parse_ip%D_SiO(tj,ti) = parse_ip%D_SiO(ti,tj)
      parse_ip%C_OOprime(tj,ti) = parse_ip%C_OOprime(ti,tj)
      parse_ip%D_OOprime(tj,ti) = parse_ip%D_OOprime(ti,tj)
    endif

    parse_ip%cutoff = maxval(parse_ip%a*parse_ip%sigma)

  elseif (parse_in_ip .and. name == 'per_triplet_data') then

    call QUIP_FoX_get_value(attributes, "atnum_c", value, status)
    if (status /= 0) call system_abort ("IPModel_SW_VP_read_params_xml cannot find atnum_c")
    read (value, *) Zi
    call QUIP_FoX_get_value(attributes, "atnum_j", value, status)
    if (status /= 0) call system_abort ("IPModel_SW_VP_read_params_xml cannot find atnum_j")
    read (value, *) Zj            
    call QUIP_FoX_get_value(attributes, "atnum_k", value, status)
    if (status /= 0) call system_abort ("IPModel_SW_VP_read_params_xml cannot find atnum_k")
    read (value, *) Zk

    ti = get_type(parse_ip%type_of_atomic_num,Zi)
    tj = get_type(parse_ip%type_of_atomic_num,Zj)
    tk = get_type(parse_ip%type_of_atomic_num,Zk)

    call QUIP_FoX_get_value(attributes, "lambda", value, status)
    if (status /= 0) call system_abort ("IPModel_SW_VP_read_params_xml cannot find lambda")
    read (value, *) parse_ip%lambda(ti,tj,tk)
    call QUIP_FoX_get_value(attributes, "gamma1", value, status)
    if (status /= 0) call system_abort ("IPModel_SW_VP_read_params_xml cannot find gamma1")
    read (value, *) parse_ip%gamma1(ti,tj,tk)
    call QUIP_FoX_get_value(attributes, "gamma2", value, status)
    if (status /= 0) call system_abort ("IPModel_SW_VP_read_params_xml cannot find gamma2")
    read (value, *) parse_ip%gamma2(ti,tj,tk)
    call QUIP_FoX_get_value(attributes, "d1", value, status)
    if (status /= 0) call system_abort ("IPModel_SW_VP_read_params_xml cannot find d1")
    read (value, *) parse_ip%d1(ti,tj,tk)  
    call QUIP_FoX_get_value(attributes, "d2", value, status)
    if (status /= 0) call system_abort ("IPModel_SW_VP_read_params_xml cannot find d2")
    read (value, *) parse_ip%d2(ti,tj,tk) 
    call QUIP_FoX_get_value(attributes, "costheta0", value, status)
    if (status /= 0) call system_abort ("IPModel_SW_VP_read_params_xml cannot find costheta0")
    read (value, *) parse_ip%costheta0(ti,tj,tk)  
    call QUIP_FoX_get_value(attributes, "eps", value, status)
    if (status /= 0) call system_abort ("IPModel_SW_VP_read_params_xml cannot find eps")
    read (value, *) parse_ip%eps3(ti,tj,tk)


    if (tj /= tk) then
      parse_ip%lambda(ti,tk,tj) = parse_ip%lambda(ti,tj,tk)
      parse_ip%gamma1(ti,tk,tj) = parse_ip%gamma2(ti,tj,tk)     ! ti is central atom for tj /=tk, the legs change and hence gamma1 and gamma2 interchange
      parse_ip%gamma2(ti,tk,tj) = parse_ip%gamma1(ti,tj,tk)     !
      parse_ip%d1(ti,tk,tj) = parse_ip%d2(ti,tj,tk)             !  ti is central atom for tj /=tk, the legs change and hence d1 and d2 interchange
      parse_ip%d2(ti,tk,tj) = parse_ip%d1(ti,tj,tk)             ! 
      parse_ip%eps3(ti,tk,tj) = parse_ip%eps3(ti,tj,tk)
      parse_ip%costheta0(ti,tk,tj) = parse_ip%costheta0(ti,tj,tk)
    endif
  endif

end subroutine IPModel_startElement_handler

subroutine IPModel_endElement_handler(URI, localname, name)
  character(len=*), intent(in)   :: URI  
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name 

  if (parse_in_ip) then
    if (name == 'SW_VP_params') then
      call print("endElement_handler SW_VP_params", PRINT_NERD)
      parse_in_ip = .false.
    endif
  endif

end subroutine IPModel_endElement_handler

subroutine IPModel_SW_VP_read_params_xml(this, param_str)
  type(IPModel_SW_VP), intent(inout), target :: this
  character(len=*) :: param_str

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
    call system_abort("IPModel_SW_VP_read_params_xml parsed file, but n_types = 0")
  endif

    end subroutine IPModel_SW_VP_read_params_xml


!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!X printing
!% Printing of SW_VP parameters: number of different types, cutoff radius, atomic numbers, ect.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


subroutine IPModel_SW_VP_Print (this, file)
  type(IPModel_SW_VP), intent(in) :: this
  type(Inoutput), intent(inout),optional :: file

  integer ti, tj, tk

  call Print("IPModel_SW_VP : Combined Stillinger-Weber/Vashishta potential", file=file)
  call Print("IPModel_SW_VP : n_types = " // this%n_types // " cutoff = " // this%cutoff, file=file)

  do ti=1, this%n_types
    call Print ("IPModel_SW_VP : type "// ti // " atomic_num " // this%atomic_num(ti), file=file)
    call verbosity_push_decrement()
    do tj=1, this%n_types
      call Print ("IPModel_SW_VP : pair interaction ti tj " // ti // " " // tj // " Zi Zj " // this%atomic_num(ti) // &
	" " // this%atomic_num(tj), file=file)
      call Print ("IPModel_SW_VP : pair " // this%AA(ti,tj) // " " // this%BB(ti,tj) // " " // this%p(ti,tj) // " " // &
	this%q(ti,tj) // " " // this%a(ti,tj), file=file)
      call Print ("IPModel_SW_VP :      " // this%sigma(ti,tj) // " " // this%eps2(ti,tj), file=file)
      do tk=1, this%n_types
	call Print ("IPModel_SW_VP : triplet interaction ti tj " // ti // " " // tj // " " // tk // &
	  " Zi Zj Zk " // this%atomic_num(ti) // " " // this%atomic_num(tj) // " " // this%atomic_num(tk), file=file)
	call Print ("IPModel_SW_VP : triplet " // this%lambda(ti,tj,tk) // " " // this%gamma1(ti,tj,tk) // " " // &
	  this%eps3(ti,tj,tk), file=file)
      end do
    end do
    call verbosity_pop()
  end do

end subroutine IPModel_SW_VP_Print

end module IPModel_SW_VP_module
