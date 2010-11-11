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
!X IPModel_GAP module  
!X
!% Module for Gaussian Approximation Potential.
!%
!% The IPModel_GAP object contains all the parameters read from a
!% 'GAP_params' XML stanza.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#include "error.inc"

module IPModel_GAP_module

use libatoms_module

use mpi_context_module
use QUIP_Common_module

#ifdef HAVE_GP_PREDICT
use descriptors_module
use gp_predict_module
#endif

implicit none

private 

include 'IPModel_interface.h'

! this stuff is here for now, but it should live somewhere else eventually
! lower down in the GP

public :: IPModel_GAP

type IPModel_GAP

  real(dp) :: cutoff = 0.0_dp                                  !% Cutoff for computing connection.

  ! bispectrum parameters
  integer :: j_max = 0
  real(dp) :: z0 = 0.0_dp
  integer :: n_species = 0                                       !% Number of atomic types.
  integer, dimension(:), allocatable :: Z
  real(dp), dimension(116) :: z_eff = 0.0_dp
  real(dp), dimension(116) :: w_Z = 1.0_dp
  real(dp) :: e0 = 0.0_dp
  real(dp) :: f0 = 0.0_dp

  ! qw parameters
  integer :: qw_l_max = 0
  integer :: qw_f_n = 0
  logical :: qw_do_q = .false.
  logical :: qw_do_w = .false.
  real(dp), allocatable :: qw_cutoff(:)
  integer, allocatable :: qw_cutoff_f(:)
  real(dp), allocatable :: qw_cutoff_r1(:)
  real(dp), dimension(:,:), allocatable :: pca_matrix

  logical :: do_pca = .false.

  character(len=256) :: coordinates             !% Coordinate system used in GAP database

  character(len=FIELD_LENGTH) :: label

#ifdef HAVE_GP_PREDICT
  type(gp) :: my_gp
#endif
  logical :: initialised = .false.
  type(extendable_str) :: command_line

end type IPModel_GAP

logical, private :: parse_in_ip, parse_matched_label, parse_in_ip_done
integer, private :: parse_n_row, parse_cur_row

type(IPModel_GAP), private, pointer :: parse_ip
type(extendable_str), save :: parse_cur_data

interface Initialise
  module procedure IPModel_GAP_Initialise_str
end interface Initialise

interface Finalise
  module procedure IPModel_GAP_Finalise
end interface Finalise

interface Print
  module procedure IPModel_GAP_Print
end interface Print

interface Calc
  module procedure IPModel_GAP_Calc
end interface Calc

contains

subroutine IPModel_GAP_Initialise_str(this, args_str, param_str)
  type(IPModel_GAP), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str
  type(Dictionary) :: params

  call Finalise(this)

  ! now initialise the potential
#ifndef HAVE_GP_PREDICT
  call system_abort('IPModel_GAP_Initialise_str: compiled without HAVE_GP_PREDICT')
#else
  
  call initialise(params)
  this%label=''

  call param_register(params, 'label', '', this%label)
  if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_SW_Initialise_str args_str')) &
  call system_abort("IPModel_GAP_Initialise_str failed to parse label from args_str="//trim(args_str))
  call finalise(params)

  call IPModel_GAP_read_params_xml(this, param_str)
  !call gp_read_binary(this%my_gp, this%datafile)
  call gp_read_xml(this%my_gp, param_str,label=this%label)


  if (trim(this%coordinates) == 'qw') then
     this%cutoff = maxval(this%qw_cutoff)
  endif
#endif  

end subroutine IPModel_GAP_Initialise_str

subroutine IPModel_GAP_Finalise(this)
  type(IPModel_GAP), intent(inout) :: this
#ifdef HAVE_GP_PREDICT
  if (allocated(this%qw_cutoff)) deallocate(this%qw_cutoff)
  if (allocated(this%qw_cutoff_f)) deallocate(this%qw_cutoff_f)
  if (allocated(this%qw_cutoff_r1)) deallocate(this%qw_cutoff_r1)

  if (allocated(this%Z)) deallocate(this%Z)

  if (this%my_gp%initialised) call finalise(this%my_gp)


  this%cutoff = 0.0_dp
  this%j_max = 0
  this%z0 = 0.0_dp
  this%n_species = 0
  this%z_eff = 0.0_dp
  this%w_Z = 1.0_dp
  this%qw_l_max = 0
  this%qw_f_n = 0
  this%qw_do_q = .false.
  this%qw_do_w = .false.

  this%coordinates = ''

  this%label = ''
  this%initialised = .false.
#endif

  call finalise(this%command_line)

end subroutine IPModel_GAP_Finalise

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!% The potential calculator: this routine computes energy, forces and the virial.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

subroutine IPModel_GAP_Calc(this, at, e, local_e, f, virial, local_virial, args_str, mpi, error)
  type(IPModel_GAP), intent(in) :: this
  type(Atoms), intent(in) :: at
  real(dp), intent(out), optional :: e, local_e(:) !% \texttt{e} = System total energy, \texttt{local_e} = energy of each atom, vector dimensioned as \texttt{at%N}.  
  real(dp), intent(out), optional :: f(:,:), local_virial(:,:)   !% Forces, dimensioned as \texttt{f(3,at%N)}, local virials, dimensioned as \texttt{local_virial(9,at%N)} 
  real(dp), intent(out), optional :: virial(3,3)   !% Virial
  character(len=*), intent(in), optional :: args_str 
  type(MPI_Context), intent(in), optional :: mpi
  integer, intent(out), optional :: error

#ifdef HAVE_GP_PREDICT
  real(dp), pointer :: w_e(:)
  real(dp) :: e_i, f_gp, f_gp_k, water_monomer_energy, water_dimer_energy
  real(dp), dimension(:), allocatable   :: local_e_in, w
  real(dp), dimension(:,:,:), allocatable   :: virial_in
  real(dp), dimension(:,:), allocatable   :: vec
  real(dp), dimension(:,:,:), allocatable   :: jack
  integer, dimension(:,:), allocatable :: water_monomer_index
  integer :: d, i, j, k, n, nei_max, jn, iAo, iBo, n_water_pair

  real(dp), dimension(:,:), allocatable :: covariance

  integer, dimension(3) :: shift
  type(Dictionary) :: params
  logical, dimension(:), pointer :: atom_mask_pointer
  logical :: has_atom_mask_name
  character(FIELD_LENGTH) :: atom_mask_name


  type(fourier_so4), save :: f_hat
  type(grad_fourier_so4), save :: df_hat
  type(bispectrum_so4), save :: bis
  type(grad_bispectrum_so4), save :: dbis

  type(fourier_so3), save :: f3_hat
  type(grad_fourier_so3), save :: df3_hat
  type(qw_so3), save :: qw
  type(grad_qw_so3), save :: dqw

  !$omp threadprivate(f_hat,df_hat,bis,dbis)  
  !$omp threadprivate(f3_hat,df3_hat,qw,dqw)  

  INIT_ERROR(error)

  if (present(e)) then
     e = 0.0_dp
  endif

  if (present(local_e)) then
     call check_size('Local_E',local_e,(/at%N/),'IPModel_GAP_Calc', error)
     local_e = 0.0_dp
  endif

  if (present(f)) then 
     call check_size('Force',f,(/3,at%N/),'IPModel_GAP_Calc', error)
     f = 0.0_dp
  end if

  if (present(virial)) then
     virial = 0.0_dp
  endif

  if (present(local_virial)) then
     call check_size('Local_virial',local_virial,(/9,at%N/),'IPModel_GAP_Calc', error)
     local_virial = 0.0_dp
  endif

  if(present(e) .or. present(local_e) ) then
     allocate(local_e_in(at%N))
     local_e_in = 0.0_dp
  endif

  if (present(virial) .or. present(local_virial)) then
     allocate(virial_in(3,3,at%N))
     virial_in = 0.0_dp
  endif

  if (.not. assign_pointer(at, "weight", w_e)) nullify(w_e)

  atom_mask_pointer => null()
  if(present(args_str)) then
     call initialise(params)
     call param_register(params, 'atom_mask_name', 'NONE',atom_mask_name,has_value_target=has_atom_mask_name)
     if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='IPModel_GAP_Calc args_str')) &
     call system_abort("IPModel_GAP_Calc failed to parse args_str='"//trim(args_str)//"'")
     call finalise(params)


     if( has_atom_mask_name ) then
        if (.not. assign_pointer(at, trim(atom_mask_name) , atom_mask_pointer)) &
        call system_abort("IPModel_GAP_Calc did not find "//trim(atom_mask_name)//" propery in the atoms object.")
     else
        atom_mask_pointer => null()
     endif
  endif

  allocate(w(at%N))
  select case(trim(this%coordinates))
  case('water_monomer','water_dimer')
     w = 1.0_dp
  case default
     do i = 1, at%N
        w(i) = this%w_Z(at%Z(i))
     enddo
  endselect


  select case(trim(this%coordinates))
  case('water_monomer')
     d = 3
  case('water_dimer')
     d = WATER_DIMER_D
  case('bispectrum')
     d = j_max2d(this%j_max)
     call cg_initialise(this%j_max,2)
  case('qw')
     d = (this%qw_l_max / 2) * this%qw_f_n
     if (this%qw_do_q .and. this%qw_do_w) d = d * 2
  case default
     call system_abort('IPModel_GAP_Calc: coordinates = '//trim(this%coordinates)//' unknown')
  endselect


  nei_max = 0
  do i = 1, at%N
     if( nei_max < (atoms_n_neighbours(at,i)+1) ) nei_max = atoms_n_neighbours(at,i)+1
  enddo

  select case(trim(this%coordinates))
  case('water_monomer')
     if(mod(at%N,3) /= 0) call system_abort('IPModel_GAP_Calc: number of atoms '//at%N//' cannot be divided by 3 - cannot be pure water.')

     allocate(vec(d,at%N/3),water_monomer_index(3,at%N/3))
     call find_water_monomer(at,water_monomer_index)
     do i = 1, at%N/3
        vec(:,i) = water_monomer(at,water_monomer_index(:,i))
     enddo

  case('water_dimer')
     if(mod(at%N,3) /= 0) call system_abort('IPModel_GAP_Calc: number of atoms '//at%N//' cannot be divided by 3 - cannot be pure water.')

     n_water_pair = 0

     do i = 1, at%N
        if(at%Z(i) == 8) then
           do n = 1, atoms_n_neighbours(at,i)
              j = atoms_neighbour(at,i,n)
              if(at%Z(j) == 8) n_water_pair = n_water_pair + 1
           enddo
        endif
     enddo

     n_water_pair = n_water_pair / 2 ! Water dimers were double counted

     allocate(vec(d,n_water_pair),water_monomer_index(3,at%N/3))
     call find_water_monomer(at,water_monomer_index)
     k = 0
     do i = 1, at%N/3
        iAo = water_monomer_index(1,i)
        do n = 1, atoms_n_neighbours(at,iAo)
           iBo = atoms_neighbour(at,iAo,n)
           if(at%Z(iBo) == 8) then
              j = find_in_array(water_monomer_index(1,:),iBo)
              if( i < j ) then
                 k = k + 1
                 vec(:,k) = water_dimer(at,water_monomer_index(:,i),water_monomer_index(:,j),this%cutoff)
              endif
           endif
        enddo
     enddo
  case default
     allocate(vec(d,at%N))
     vec = 0.0_dp
     if(present(f) .or. present(virial) .or. present(local_virial)) then
        allocate(jack(d,3*nei_max,at%N))
        jack = 0.0_dp
     endif
  endselect

!$omp parallel 

  if (trim(this%coordinates) == 'bispectrum') then
     call initialise(f_hat,this%j_max,this%z0,this%cutoff)
     if(present(f).or.present(virial) .or. present(local_virial)) call initialise(df_hat,this%j_max,this%z0,this%cutoff)
  elseif (trim(this%coordinates) == 'qw') then
     call initialise(f3_hat, this%qw_l_max, this%qw_cutoff, this%qw_cutoff_f, this%qw_cutoff_r1)
     call initialise(qw, this%qw_l_max, this%qw_f_n, do_q = this%qw_do_q, do_w = this%qw_do_w)
     if (present(f) .or. present(virial) .or. present(local_virial)) then
        call initialise(df3_hat, this%qw_l_max, this%qw_cutoff, this%qw_cutoff_f, this%qw_cutoff_r1)
        call initialise(dqw, this%qw_l_max, this%qw_f_n, do_q = this%qw_do_q, do_w = this%qw_do_w)
     endif
  endif

!$omp do private(n)
  do i = 1, at%N
     if (present(mpi)) then
	if (mpi%active) then
	  if (mod(i-1, mpi%n_procs) /= mpi%my_proc) cycle
	endif
      endif

     if(associated(atom_mask_pointer)) then
        if(.not. atom_mask_pointer(i)) cycle
     endif


     if (trim(this%coordinates) == 'bispectrum') then
        call fourier_transform(f_hat,at,i,w)
        call calc_bispectrum(bis,f_hat)
        call bispectrum2vec(bis,vec(:,i))

        if(this%do_pca) vec(:,i) = matmul(vec(:,i),this%pca_matrix)

        if(present(f).or.present(virial) .or. present(local_virial)) then
           do n = 0, atoms_n_neighbours(at,i)
              call fourier_transform(df_hat,at,i,n,w)
              call calc_bispectrum(dbis,f_hat,df_hat)
              call bispectrum2vec(dbis,jack(:,3*n+1:3*(n+1),i))

              if(this%do_pca) jack(:,3*n+1:3*(n+1),i) = matmul(transpose(this%pca_matrix), jack(:,3*n+1:3*(n+1),i))

           enddo
        endif
     elseif (trim(this%coordinates) == 'qw') then
        call fourier_transform(f3_hat, at, i)
        call calc_qw(qw, f3_hat)
        call qw2vec(qw, vec(:,i))
        if (present(f) .or. present(virial) .or. present(local_virial)) then
           do n = 0, atoms_n_neighbours(at, i)
              call fourier_transform(df3_hat, at, i, n)
              call calc_qw(dqw, f3_hat, df3_hat)
              call qw2vec(dqw, jack(:,3*n+1:3*(n+1),i))
           enddo
        endif
     endif

  enddo
!$omp end do 


  if (trim(this%coordinates) == 'bispectrum') then
     call finalise(f_hat)
     call finalise(df_hat)
     call finalise(bis)
     call finalise(dbis)
  elseif (trim(this%coordinates) == 'qw') then
     call finalise(f3_hat)
     call finalise(df3_hat)
     call finalise(qw)
     call finalise(dqw)
  endif

!$omp end parallel


  select case(trim(this%coordinates))
  case('water_monomer','water_dimer')

  case default
     allocate(covariance(this%my_gp%n,at%N))
     call gp_precompute_covariance(this%my_gp,vec,covariance,at%Z,mpi)
  endselect

    
  select case(trim(this%coordinates))
  case('water_monomer','water_dimer')
     do i = 1, size(vec,2)
        if(present(e)) then
           call gp_predict(gp_data=this%my_gp, mean=water_monomer_energy,x_star=vec(:,i))
           local_e_in(1) = local_e_in(1) + water_monomer_energy + this%e0
        endif
     enddo

  case default
!$omp parallel do private(k,f_gp,f_gp_k,n,j,jn,shift)
     do i = 1, at%N
        if (present(mpi)) then
           if (mpi%active) then
              if (mod(i-1, mpi%n_procs) /= mpi%my_proc) cycle
           endif
        endif

        if(associated(atom_mask_pointer)) then
           if(.not. atom_mask_pointer(i)) cycle
        endif

        if(present(e) .or. present(local_e)) then

           call gp_predict(gp_data=this%my_gp, mean=local_e_in(i),x_star=vec(:,i),Z=at%Z(i),c_in=covariance(:,i))

           local_e_in(i) = local_e_in(i) + this%e0
        endif

        if(present(f).or.present(virial) .or. present(local_virial)) then
           do k = 1, 3
              f_gp = 0.0_dp
       

              call gp_predict(gp_data=this%my_gp, mean=f_gp_k,x_star=vec(:,i),x_prime_star=jack(:,k,i),Z=at%Z(i),c_in=covariance(:,i))

       
              if( present(f) ) f(k,i) = f(k,i) - f_gp_k
              if( present(virial) .or. present(local_virial) ) virial_in(:,k,i) = virial_in(:,k,i) - f_gp_k*at%pos(:,i)

              do n = 1, atoms_n_neighbours(at,i)
                 j = atoms_neighbour(at,i,n,jn=jn,shift=shift)
       

                 call gp_predict(gp_data=this%my_gp,mean=f_gp_k,x_star=vec(:,i),x_prime_star=jack(:,n*3+k,i),Z=at%Z(i),c_in=covariance(:,i))

!$omp critical
                 if( present(f) ) f(k,j) = f(k,j) - f_gp_k
                 if( present(virial) .or. present(local_virial) ) virial_in(:,k,j) = virial_in(:,k,j) - f_gp_k*( at%pos(:,j) + matmul(at%lattice,shift) )
!$omp end critical
              enddo
       
           enddo
        endif
     enddo
  endselect

  if (present(mpi)) then
     if(present(f)) call sum_in_place(mpi,f)
     if(present(virial) .or. present(local_virial)) call sum_in_place(mpi,virial_in)
     if(present(e) .or. present(local_e) ) call sum_in_place(mpi,local_e_in)
  endif

  if(present(e)) e = sum(local_e_in)
  if(present(local_e)) local_e = local_e_in
  if(present(virial)) virial = sum(virial_in,dim=3)
  if(present(local_virial)) then
     do i = 1, at%N
        local_virial(:,i) = reshape(virial_in(:,:,i),(/9/))
     enddo
  endif

  if(allocated(local_e_in)) deallocate(local_e_in)
  if(allocated(virial_in)) deallocate(virial_in)
  if(allocated(w)) deallocate(w)
  if(allocated(covariance)) deallocate(covariance)
  if(allocated(vec)) deallocate(vec)
  if(allocated(jack)) deallocate(jack)
  if(allocated(water_monomer_index)) deallocate(water_monomer_index)
  atom_mask_pointer => null()

#endif

end subroutine IPModel_GAP_Calc

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!% XML param reader functions.
!% An example for XML stanza is given below, please notice that
!% they are simply dummy parameters for testing purposes, with no physical meaning.
!%
!%> <GAP_params datafile="file" label="default">
!%> </GAP_params>
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

subroutine IPModel_startElement_handler(URI, localname, name, attributes)
  character(len=*), intent(in)   :: URI  
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name 
  type(dictionary_t), intent(in) :: attributes

  integer :: status
  character(len=1024) :: value

  integer :: ti, ri

  if(name == 'GAP_params') then ! new GAP stanza
     
     if(parse_in_ip) &
        call system_abort("IPModel_startElement_handler entered GAP_params with parse_in true. Probably a bug in FoX (4.0.1, e.g.)")
     
     if(parse_matched_label) return ! we already found an exact match for this label
     
     call QUIP_FoX_get_value(attributes, 'label', value, status)
     if(status /= 0) value = ''
     
     if(len(trim(parse_ip%label)) > 0) then ! we were passed in a label
        if(value == parse_ip%label) then ! exact match
 	   parse_matched_label = .true.
 	   parse_in_ip = .true.
        else ! no match
 	   parse_in_ip = .false.
        endif
     else ! no label passed in
        parse_in_ip = .true.
     endif

    if(parse_in_ip) then
       if(parse_ip%initialised) call finalise(parse_ip)
    endif

  elseif(parse_in_ip .and. name == 'GAP_data') then

     call QUIP_FoX_get_value(attributes, 'n_species', value, status)
     if(status == 0) then
        read (value, *) parse_ip%n_species
     else
        call system_abort('IPModel_GAP_read_params_xml cannot find n_species')
     endif

     call QUIP_FoX_get_value(attributes, 'e0', value, status)
     if(status == 0) then
        read (value, *) parse_ip%e0
     else
        call system_abort('IPModel_GAP_read_params_xml cannot find e0')
     endif

     call QUIP_FoX_get_value(attributes, 'do_pca', value, status)
     if(status == 0) then
        read (value, *) parse_ip%do_pca
     else
        parse_ip%do_pca = .false.
     endif

     call QUIP_FoX_get_value(attributes, 'f0', value, status)
     if(status == 0) then
        read (value, *) parse_ip%f0
     else
        call system_abort('IPModel_GAP_read_params_xml cannot find f0')
     endif

     call QUIP_FoX_get_value(attributes, 'coordinates', value, status)
     if(status == 0) then
        read (value, *) parse_ip%coordinates
     else
        call system_abort('IPModel_GAP_read_params_xml cannot find coordinates')
     endif

     allocate( parse_ip%Z(parse_ip%n_species) )

  elseif(parse_in_ip .and. name == 'water_monomer_params') then

     call QUIP_FoX_get_value(attributes, 'cutoff', value, status)
     if(status == 0) then
        read (value, *) parse_ip%cutoff
     else
        call system_abort('IPModel_GAP_read_params_xml cannot find cutoff')
     endif

  elseif(parse_in_ip .and. name == 'water_dimer_params') then

     call QUIP_FoX_get_value(attributes, 'cutoff', value, status)
     if(status == 0) then
        read (value, *) parse_ip%cutoff
     else
        call system_abort('IPModel_GAP_read_params_xml cannot find cutoff')
     endif

  elseif(parse_in_ip .and. name == 'bispectrum_so4_params') then

     call QUIP_FoX_get_value(attributes, 'cutoff', value, status)
     if(status == 0) then
        read (value, *) parse_ip%cutoff
     else
        call system_abort('IPModel_GAP_read_params_xml cannot find cutoff')
     endif

     call QUIP_FoX_get_value(attributes, 'j_max', value, status)
     if(status == 0) then
        read (value, *) parse_ip%j_max
     else
        call system_abort('IPModel_GAP_read_params_xml cannot find j_max')
     endif

     call QUIP_FoX_get_value(attributes, 'z0', value, status)
     if(status == 0) then
        read (value, *) parse_ip%z0
     else
        call system_abort('IPModel_GAP_read_params_xml cannot find z0')
     endif

  elseif(parse_in_ip .and. name == 'qw_so3_params') then

     call QUIP_FoX_get_value(attributes, 'l_max', value, status)
     if(status == 0) then
        read (value, *) parse_ip%qw_l_max
     else
        call system_abort('IPModel_GAP_read_params_xml cannot find l_max')
     endif

     call QUIP_FoX_get_value(attributes, 'n_radial', value, status)
     if(status == 0) then
        read (value, *) parse_ip%qw_f_n
     else
        call system_abort('IPModel_GAP_read_params_xml cannot find n_radial')
     endif

     call QUIP_FoX_get_value(attributes, 'do_q', value, status)
     if(status == 0) then
        read (value, *) parse_ip%qw_do_q
     else
        call system_abort('IPModel_GAP_read_params_xml cannot find do_q')
     endif

     call QUIP_FoX_get_value(attributes, 'do_w', value, status)
     if(status == 0) then
        read (value, *) parse_ip%qw_do_w
     else
        call system_abort('IPModel_GAP_read_params_xml cannot find do_w')
     endif

     allocate(parse_ip%qw_cutoff(parse_ip%qw_f_n), parse_ip%qw_cutoff_f(parse_ip%qw_f_n), parse_ip%qw_cutoff_r1(parse_ip%qw_f_n))

  elseif(parse_in_ip .and. name == 'radial_function') then

     call QUIP_FoX_get_value(attributes, 'i', value, status)
     if(status == 0) then
        read (value, *) ri
     else
        call system_abort('IPModel_GAP_read_params_xml cannot find i')
     endif

     call QUIP_FoX_get_value(attributes, 'cutoff', value, status)
     if(status == 0) then
        read (value, *) parse_ip%qw_cutoff(ri)
     else
        call system_abort('IPModel_GAP_read_params_xml cannot find cutoff')
     endif

     call QUIP_FoX_get_value(attributes, 'cutoff_type', value, status)
     if(status == 0) then
        read (value, *) parse_ip%qw_cutoff_f(ri)
     else
        call system_abort('IPModel_GAP_read_params_xml cannot find cutoff_type')
     endif

     call QUIP_FoX_get_value(attributes, 'cutoff_r1', value, status)
     if(status == 0) then
        read (value, *) parse_ip%qw_cutoff_r1(ri)
     else
        call system_abort('IPModel_GAP_read_params_xml cannot find cutoff_r1')
     endif

  elseif(parse_in_ip .and. name == 'PCA_matrix') then

     call QUIP_FoX_get_value(attributes, 'n', value, status)
     if(status == 0) then
        read (value, *) parse_n_row
     else
        call system_abort('IPModel_GAP_read_params_xml cannot find n')
     endif
     allocate(parse_ip%pca_matrix(parse_n_row,parse_n_row))

  elseif(parse_in_ip .and. name == 'row') then

     call QUIP_FoX_get_value(attributes, 'i', value, status)
     if(status == 0) then
        read (value, *) parse_cur_row
     else
        call system_abort('IPModel_GAP_read_params_xml cannot find i')
     endif

     call zero(parse_cur_data)

  elseif(parse_in_ip .and. name == 'per_type_data') then

     call QUIP_FoX_get_value(attributes, 'i', value, status)
     if(status == 0) then
        read (value, *) ti
     else
        call system_abort('IPModel_GAP_read_params_xml cannot find i')
     endif

     call QUIP_FoX_get_value(attributes, 'atomic_num', value, status)
     if(status == 0) then
        read (value, *) parse_ip%Z(ti)
     else
        call system_abort('IPModel_GAP_read_params_xml cannot find atomic_num')
     endif

     call QUIP_FoX_get_value(attributes, 'weight', value, status)
     if(status == 0) then
        read (value, *) parse_ip%w_Z(parse_ip%Z(ti))
     else
        call system_abort('IPModel_GAP_read_params_xml cannot find weight')
     endif

  elseif(parse_in_ip .and. name == 'command_line') then
      call zero(parse_cur_data)

  endif

end subroutine IPModel_startElement_handler

subroutine IPModel_endElement_handler(URI, localname, name)
  character(len=*), intent(in)   :: URI  
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name 

  character(len=100*parse_n_row) :: val

  if (parse_in_ip) then
    if(name == 'GAP_params') then
       parse_in_ip = .false.
       parse_in_ip_done = .true.
    elseif(name == 'GAP_data') then

    elseif(name == 'bispectrum_so4_params') then

    elseif(name == 'water_monomer_params') then

    elseif(name == 'water_dimer_params') then

    elseif(name == 'qw_so3_params') then

    elseif(name == 'radial_function') then

    elseif(name == 'per_type_data') then

    elseif(name == 'row') then

       val = string(parse_cur_data)
       read(val,*) parse_ip%pca_matrix(:,parse_cur_row)

    elseif(name == 'command_line') then
       parse_ip%command_line = parse_cur_data
    end if
  endif

end subroutine IPModel_endElement_handler

subroutine IPModel_characters_handler(in)
   character(len=*), intent(in) :: in

   if(parse_in_ip) then
     call concat(parse_cur_data, in, keep_lf=.false.)
   endif

end subroutine IPModel_characters_handler

subroutine IPModel_GAP_read_params_xml(this, param_str)
  type(IPModel_GAP), intent(inout), target :: this
  character(len=*), intent(in) :: param_str

  type(xml_t) :: fxml

  if (len(trim(param_str)) <= 0) then
     call system_abort('IPModel_GAP_read_params_xml: invalid param_str length '//len(trim(param_str)) )
  else
     parse_in_ip = .false.
     parse_in_ip_done = .false.
     parse_matched_label = .false.
     parse_ip => this
     call initialise(parse_cur_data)

     call open_xml_string(fxml, param_str)
     call parse(fxml,  &
       startElement_handler = IPModel_startElement_handler, &
       endElement_handler = IPModel_endElement_handler, &
       characters_handler = IPModel_characters_handler)
     call close_xml_t(fxml)

     if(.not. parse_in_ip_done) &
     call  system_abort('IPModel_GAP_read_params_xml: could not initialise GAP potential. No GAP_params present?')
     this%initialised = .true.
  endif

end subroutine IPModel_GAP_read_params_xml

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!% Printing of GAP parameters: number of different types, cutoff radius, atomic numbers, etc.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

subroutine IPModel_GAP_Print (this, file)
  type(IPModel_GAP), intent(in) :: this
  type(Inoutput), intent(inout),optional :: file
  integer :: i

#ifdef HAVE_GP_PREDICT
  call Print("IPModel_GAP : Gaussian Approximation Potential", file=file)
  call Print("IPModel_GAP : cutoff = "//this%cutoff, file=file)
  call Print("IPModel_GAP : j_max = "//this%j_max, file=file)
  call Print("IPModel_GAP : z0 = "//this%z0, file=file)
  call Print("IPModel_GAP : n_species = "//this%n_species, file=file)

  do i = 1, this%n_species
     call Print("IPModel_GAP : Z = "//this%Z(i), file=file)
     call Print("IPModel_GAP : z_eff = "//this%z_eff(this%Z(i)), file=file)
     call Print("IPModel_GAP : delta = "//this%my_gp%delta(i), file=file)
     call Print("IPModel_GAP : theta = "//this%my_gp%theta(:,i), file=file)
  enddo

  call Print("IPModel_GAP : command_line = "//string(this%command_line))
#endif

end subroutine IPModel_GAP_Print

end module IPModel_GAP_module
