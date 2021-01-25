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
!X IPModel_vdW module  
!X
!% Module for Gaussian Approximation Potential.
!%
!% The IPModel_vdW object contains all the parameters read from a
!% 'vdW_params' XML stanza.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#include "error.inc"

module IPModel_vdW_module

use error_module
use system_module, only : dp, inoutput, string_to_int, reallocate, system_timer , print, PRINT_NERD
use units_module
use dictionary_module
use periodictable_module, only: total_elements
use extendable_str_module
use paramreader_module
use linearalgebra_module
use atoms_types_module
use atoms_module

use mpi_context_module
use QUIP_Common_module

#ifdef HAVE_GAP
use descriptors_module
use gp_predict_module
#endif

implicit none

private 

include 'IPModel_interface.h'

#ifdef GAP_VERSION
   integer, parameter :: gap_version = GAP_VERSION
#else
   integer, parameter :: gap_version = 0
#endif

! this stuff is here for now, but it should live somewhere else eventually
! lower down in the GP

public :: IPModel_vdW

type IPModel_vdW

  real(dp) :: cutoff = 0.0_dp                                  !% Cutoff for computing connection.

  real(dp) :: E_scale = 0.0_dp                                 !% scale factor for the potential 
  real(dp) :: sr = 0.0_dp
  real(dp) :: d = 0.0_dp
  real(dp) :: vdW_cutoff = 0.0_dp                              
  real(dp) :: vdW_cutoff_buffer = 0.0_dp                             

  character(len=STRING_LENGTH) :: label

#ifdef HAVE_GAP
  type(gpSparse) :: my_gp
  type(descriptor), dimension(:), allocatable :: my_descriptor
#endif
  real(dp), dimension(total_elements) :: e0 = 0.0_dp
  integer(dp), dimension(total_elements) :: map = -1
  integer :: n_types = 0
  integer, dimension(:), allocatable :: Z
  real(dp), dimension(:), allocatable :: r0
  real(dp), dimension(:), allocatable :: alpha0
  real(dp), dimension(:), allocatable :: c6
  real(dp), dimension(:,:), allocatable :: c6_pair


  logical :: initialised = .false.
  type(extendable_str) :: command_line
  integer :: xml_version

end type IPModel_vdW

logical, private :: parse_in_ip, parse_in_vdW_data, parse_matched_label, parse_in_ip_done

type(IPModel_vdW), private, pointer :: parse_ip
type(extendable_str), save :: parse_cur_data

interface Initialise
  module procedure IPModel_vdW_Initialise_str
end interface Initialise

interface Finalise
  module procedure IPModel_vdW_Finalise
end interface Finalise

interface Print
  module procedure IPModel_vdW_Print
end interface Print

interface Calc
  module procedure IPModel_vdW_Calc
end interface Calc

contains

subroutine IPModel_vdW_Initialise_str(this, args_str, param_str)
  type(IPModel_vdW), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str
  type(Dictionary) :: params

  integer :: n, m, i_coordinate

  call Finalise(this)

  ! now initialise the potential
#ifndef HAVE_GAP
  call system_abort('IPModel_vdW_Initialise_str: must be compiled with HAVE_GAP')
#else

  call initialise(params)
  this%label=''

  call param_register(params, 'label', '', this%label, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'E_scale', '1.0', this%E_scale, help_string="rescaling factor for the potential")

  if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_SW_Initialise_str args_str')) &
  call system_abort("IPModel_vdW_Initialise_str failed to parse label from args_str="//trim(args_str))
  call finalise(params)

  call IPModel_vdW_read_params_xml(this, param_str)
  call gp_readXML(this%my_gp, param_str,label=trim(this%label))

  do n = 1, this%n_types
     this%map(this%Z(n)) = n
  enddo

  allocate(this%c6_pair(this%n_types,this%n_types))
  do n = 1, this%n_types
     do m = 1, this%n_types
        this%c6_pair(m,n) = 2.0_dp * ( this%alpha0(n)*this%alpha0(m)*this%c6(n)*this%c6(m) ) / &
           ( this%alpha0(n)**2 * this%c6(m) + this%alpha0(m)**2 * this%c6(n) )
     enddo
  enddo

  allocate(this%my_descriptor(this%my_gp%n_coordinate))

  this%cutoff = 0.0_dp
  do i_coordinate = 1, this%my_gp%n_coordinate
     call concat(this%my_gp%coordinate(i_coordinate)%descriptor_str," xml_version="//this%xml_version)
     call initialise(this%my_descriptor(i_coordinate),string(this%my_gp%coordinate(i_coordinate)%descriptor_str))
     this%cutoff = max(this%cutoff,cutoff(this%my_descriptor(i_coordinate)))
  enddo

#endif  

end subroutine IPModel_vdW_Initialise_str

subroutine IPModel_vdW_Finalise(this)
  type(IPModel_vdW), intent(inout) :: this
#ifdef HAVE_GAP

  if (this%my_gp%initialised) call finalise(this%my_gp)
  if(allocated(this%Z)) deallocate(this%Z)
  if(allocated(this%r0)) deallocate(this%r0)
  if(allocated(this%alpha0)) deallocate(this%alpha0)
  if(allocated(this%c6)) deallocate(this%c6)
  if(allocated(this%c6_pair)) deallocate(this%c6_pair)
  this%n_types = 0
  this%map = -1

  this%cutoff = 0.0_dp
  this%sr = 0.0_dp
  this%d = 0.0_dp

  this%label = ''
  this%initialised = .false.
#endif

  call finalise(this%command_line)

end subroutine IPModel_vdW_Finalise

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!% The potential calculator: this routine computes energy, forces and the virial.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

subroutine IPModel_vdW_Calc(this, at, e, local_e, f, virial, local_virial, args_str, mpi, error)
  type(IPModel_vdW), intent(inout) :: this
  type(Atoms), intent(inout) :: at
  real(dp), intent(out), optional :: e, local_e(:) !% \texttt{e} = System total energy, \texttt{local_e} = energy of each atom, vector dimensioned as \texttt{at%N}.  
  real(dp), intent(out), optional :: f(:,:), local_virial(:,:)   !% Forces, dimensioned as \texttt{f(3,at%N)}, local virials, dimensioned as \texttt{local_virial(9,at%N)} 
  real(dp), intent(out), optional :: virial(3,3)   !% Virial
  character(len=*), intent(in), optional :: args_str 
  type(MPI_Context), intent(in), optional :: mpi
  integer, intent(out), optional :: error

#ifdef HAVE_GAP
  type v_data
     real(dp) :: v
     integer :: n_lo, n_up
     integer, dimension(:), allocatable :: j
     real(dp), dimension(:,:), allocatable :: deriv
     real(dp), dimension(:,:), allocatable :: diff
  endtype v_data

  character(STRING_LENGTH) :: atom_mask_name
  
  logical :: has_atom_mask_name, do_rescale_r, do_rescale_E, do_select_descriptor, do_derivative
  logical, dimension(:), pointer :: atom_mask_pointer
  logical, dimension(:), allocatable :: mpi_local_mask

  integer :: only_descriptor, i_coordinate, d, i, j, k, i_desc, n, m

  real(dp) :: r_scale, E_scale, v_i, v_i_cutoff, r_ij, &
     r6_ij, dr6_ij_dr_ij, scaled_buffer, f_cut_ij, df_cut_ij_dr_ij, &
     c6_eff, r0_ij, d_sR_r0_ij, exp_factor, f_damp, e_i, &
     df_damp_dr_ij, df_damp_dr0_ij, df_damp_dvi, df_damp_dvj, dc6_eff_dvi, dc6_eff_dvj, &
     de_i_dr_ij
  real(dp), dimension(3) :: v_deriv, u_ij, f_ij, de_i_drk
  real(dp), dimension(:), allocatable     :: grad_predict, r0_eff, dr0_eff_dv
  real(dp), dimension(:), allocatable     :: local_e_in
  real(dp), dimension(:,:), allocatable   :: f_in
  real(dp), dimension(:,:,:), allocatable :: virial_in

  type(Dictionary) :: params
  type(descriptor_data), target :: my_descriptor_data
  type(extendable_str) :: my_args_str
  type(v_data), dimension(:), allocatable :: local_v

  !real(dp) :: e_i, e_i_cutoff
  !integer :: d, i, j, n, m, i_coordinate, i_pos0


  !real(dp), dimension(3) :: pos, f_gp
  !real(dp), dimension(3,3) :: virial_i
  !type(Dictionary) :: params
  !logical :: has_atom_mask_name
  !character(STRING_LENGTH) :: atom_mask_name, calc_local_gap_variance, calc_energy_per_coordinate
  !real(dp) :: r_scale, E_scale

  !real(dp) :: gap_variance_i_cutoff
  !real(dp), dimension(:), allocatable :: gap_variance, local_gap_variance_in
  !real(dp), dimension(:), pointer :: local_gap_variance_pointer
  !real(dp), dimension(:,:), allocatable :: gap_variance_gradient_in
  !real(dp), dimension(:,:), pointer :: gap_variance_gradient_pointer
  !real(dp) :: gap_variance_regularisation
  !logical :: do_rescale_r, do_rescale_E, do_gap_variance, print_gap_variance, do_local_gap_variance, do_energy_per_coordinate
  !integer :: only_descriptor
  !logical :: do_select_descriptor
  !logical :: mpi_parallel_descriptor

  !type(descriptor_data) :: my_descriptor_data
  !type(extendable_str) :: my_args_str

  INIT_ERROR(error)

  if (present(e)) then
     e = 0.0_dp
  endif

  if (present(local_e)) then
     call check_size('Local_E',local_e,(/at%N/),'IPModel_vdW_Calc', error)
     local_e = 0.0_dp
  endif

  if (present(f)) then 
     call check_size('Force',f,(/3,at%N/),'IPModel_vdW_Calc', error)
     f = 0.0_dp
  end if

  if (present(virial)) then
     virial = 0.0_dp
  endif

  if (present(local_virial)) then
     call check_size('Local_virial',local_virial,(/9,at%N/),'IPModel_vdW_Calc', error)
     local_virial = 0.0_dp
  endif

  do_derivative = present(f) .or. present(virial) .or. present(local_virial)

  atom_mask_pointer => null()
  has_atom_mask_name = .false.
  atom_mask_name = ""
  only_descriptor = 0

  if( any(this%map(at%Z) == -1) ) then
     RAISE_ERROR("IPModel_vdW_Calc: atoms object contains atom types for which vdW parameters were not specified",error)
  endif

   call initialise(params)
   
   call param_register(params, 'atom_mask_name', 'NONE',atom_mask_name,has_value_target=has_atom_mask_name, &
   help_string="Name of a logical property in the atoms object. For atoms where this property is true, energies, forces, virials etc. are " // &
    "calculated")
   call param_register(params, 'r_scale', '1.0',r_scale, has_value_target=do_rescale_r, help_string="Rescaling factor for distances. Default 1.0.")
   call param_register(params, 'E_scale', '1.0',E_scale, has_value_target=do_rescale_E, help_string="Rescaling factor for energy. Default 1.0.")
   call param_register(params, 'only_descriptor', '0', only_descriptor, has_value_target=do_select_descriptor, help_string="Only select a single coordinate")

   if(present(args_str)) then
     if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='IPModel_vdW_Calc args_str')) &
       call system_abort("IPModel_vdW_Calc failed to parse args_str='"//trim(args_str)//"'")
     call finalise(params)

     if( has_atom_mask_name ) then
        if (.not. assign_pointer(at, trim(atom_mask_name) , atom_mask_pointer)) &
            call system_abort("IPModel_vdW_Calc did not find "//trim(atom_mask_name)//" property in the atoms object.")
     endif
     if (do_rescale_r .or. do_rescale_E) then
        RAISE_ERROR("IPModel_vdW_Calc: rescaling of potential at the calc() stage with r_scale and E_scale not yet implemented!", error)
     end if

     my_args_str = trim(args_str)
  else
     ! call parser to set defaults
     if (.not. param_read_line(params,"",ignore_unknown=.true.,task='IPModel_vdW_Calc args_str')) &
       call system_abort("IPModel_vdW_Calc failed to parse args_str='"//trim(args_str)//"'")
     call finalise(params)
     call initialise(my_args_str)
  endif

  if( present(mpi) ) then
     if(mpi%active) then
        if(has_atom_mask_name) then
           RAISE_ERROR("IPModel_vdW: atom_mask_name "//trim(atom_mask_name)//" present while running MPI version. &
              The use of atom_mask_name is intended for serial-compiled code called from an external parallel code, such as LAMMPS",error)
        endif

        if( has_property(at,"mpi_local_mask") ) then
           RAISE_ERROR("IPModel_vdW: mpi_local_mask property already present", error)
        endif

        allocate(mpi_local_mask(at%N))
        call add_property_from_pointer(at,'mpi_local_mask',mpi_local_mask,error=error)
        call concat(my_args_str," atom_mask_name=mpi_local_mask")
     endif
  endif

  call concat(my_args_str," xml_version="//this%xml_version)

  allocate(local_v(at%N))
  do i = 1, at%N
     local_v(i)%v = this%e0(at%Z(i))
  enddo

  do i_coordinate = 1, this%my_gp%n_coordinate

     if (do_select_descriptor .and. (this%my_gp%n_coordinate > 1)) then
        if (i_coordinate /= only_descriptor) then
           call print("vdW label="//trim(this%label)//" skipping coordinate "//i_coordinate, PRINT_NERD)
           cycle
        end if
     end if

     if(mpi%active) call descriptor_MPI_setup(this%my_descriptor(i_coordinate),at,mpi,mpi_local_mask,error)

     d = descriptor_dimensions(this%my_descriptor(i_coordinate))

     if( do_derivative ) then
        if (allocated(grad_predict)) deallocate(grad_predict)
        allocate(grad_predict(d))
     end if     
     call calc(this%my_descriptor(i_coordinate),at,my_descriptor_data, &
        do_descriptor=.true.,do_grad_descriptor=present(f) .or. present(virial) .or. present(local_virial), args_str=trim(string(my_args_str)), error=error)
     PASS_ERROR(error)
     do i_desc = 1, size(my_descriptor_data%x)
        if( size(my_descriptor_data%x(i_desc)%ci) /= 1 ) then
           RAISE_ERROR("IPModel_vdW_Calc: descriptor is not local and atomic",error)
        endif
        i = my_descriptor_data%x(i_desc)%ci(1)

        if ( do_derivative ) then
           local_v(i)%n_lo = lbound(my_descriptor_data%x(i_desc)%ii,1)
           local_v(i)%n_up = ubound(my_descriptor_data%x(i_desc)%ii,1)
           allocate(local_v(i)%deriv(3,local_v(i)%n_lo:local_v(i)%n_up), &
              local_v(i)%j(local_v(i)%n_lo:local_v(i)%n_up), &
              local_v(i)%diff(3,local_v(i)%n_lo:local_v(i)%n_up))
           do n = local_v(i)%n_lo, local_v(i)%n_up
              local_v(i)%diff(:,n) = my_descriptor_data%x(i_desc)%pos(:,n) - &
                 my_descriptor_data%x(i_desc)%pos(:,local_v(i)%n_lo)
           enddo
           local_v(i)%j = my_descriptor_data%x(i_desc)%ii
           local_v(i)%deriv = 0.0_dp
        endif
     enddo

!$omp parallel default(none) private(i_desc,n,i,grad_predict,v_i,v_i_cutoff,v_deriv) &
!$omp shared(this,d,i_coordinate,my_descriptor_data,f,virial,local_virial,local_v,do_derivative)
!$omp do schedule(dynamic)
     do i_desc = 1, size(my_descriptor_data%x)
        if(present(f) .or. present(virial) .or. present(local_virial)) then
           call reallocate(grad_predict,d,zero=.true.)
           v_i = gp_predict(this%my_gp%coordinate(i_coordinate), &
              xStar=my_descriptor_data%x(i_desc)%data(:), &
              gradPredict = grad_predict )
        else
           v_i =  gp_predict(this%my_gp%coordinate(i_coordinate), xStar=my_descriptor_data%x(i_desc)%data(:))
        endif

        v_i_cutoff = v_i * my_descriptor_data%x(i_desc)%covariance_cutoff 
        i = my_descriptor_data%x(i_desc)%ci(1)
        local_v(i)%v = local_v(i)%v + v_i_cutoff

        if( do_derivative ) then
           do n = lbound(my_descriptor_data%x(i_desc)%ii,1), ubound(my_descriptor_data%x(i_desc)%ii,1)
              if( .not. my_descriptor_data%x(i_desc)%has_grad_data(n) ) cycle
              v_deriv = matmul( grad_predict,my_descriptor_data%x(i_desc)%grad_data(:,:,n)) * my_descriptor_data%x(i_desc)%covariance_cutoff + &
                 v_i * my_descriptor_data%x(i_desc)%grad_covariance_cutoff(:,n)
              local_v(i)%deriv(:,n) = local_v(i)%deriv(:,n) + v_deriv
           enddo
        endif
     enddo
!$omp end do
     if(allocated(grad_predict)) deallocate(grad_predict)
!$omp end parallel
     call finalise(my_descriptor_data)
  enddo

  ! Has to be allocated as it's in the reduction clause.
  allocate(local_e_in(at%N))
  local_e_in = 0.0_dp

  allocate(f_in(3,at%N))
  f_in = 0.0_dp

  allocate(virial_in(3,3,at%N))
  virial_in = 0.0_dp

  allocate(r0_eff(at%N))
  if( do_derivative ) then
     allocate(dr0_eff_dv(at%N))
  endif

  do i = 1, at%N
     r0_eff(i) = local_v(i)%v**(1.0_dp/3.0_dp) * this%r0(this%map(at%Z(i)))
     if( do_derivative ) then
        dr0_eff_dv(i) = r0_eff(i) / local_v(i)%v / 3.0_dp
     endif
  enddo


!$omp parallel default(none) private(i,j,k,n,m,r_ij,u_ij,r6_ij,dr6_ij_dr_ij,scaled_buffer,f_cut_ij,df_cut_ij_dr_ij, &
!$omp c6_eff,r0_ij,d_sR_r0_ij,exp_factor,f_damp,e_i,df_damp_dr_ij,df_damp_dr0_ij,df_damp_dvi,df_damp_dvj, &
!$omp dc6_eff_dvi,dc6_eff_dvj,de_i_dr_ij,f_ij,de_i_drk) &
!$omp shared(this,at,e,f,virial,do_derivative,local_v,r0_eff,dr0_eff_dv) &
!$omp reduction(+:local_e_in,f_in,virial_in)
!$omp do schedule(dynamic)
  do i = 1, at%N
     do n = 1, n_neighbours(at,i)
        j = neighbour(at, i, n, distance=r_ij, cosines=u_ij)

        if( r_ij > this%vdW_cutoff ) cycle

        r6_ij = 1.0_dp / r_ij**6
        if( do_derivative ) dr6_ij_dr_ij = -6.0_dp * r6_ij / r_ij

        if( r_ij > this%vdW_cutoff - this%vdW_cutoff_buffer ) then
           scaled_buffer = ( r_ij - this%vdW_cutoff + this%vdW_cutoff_buffer ) / this%vdW_cutoff_buffer
           f_cut_ij = 1.0_dp - 3.0_dp*scaled_buffer**2 + 2.0_dp*scaled_buffer**3

           if( do_derivative ) then
              df_cut_ij_dr_ij = -6.0_dp*(scaled_buffer+scaled_buffer**2) / this%vdW_cutoff_buffer
              dr6_ij_dr_ij = dr6_ij_dr_ij*f_cut_ij + r6_ij*df_cut_ij_dr_ij
           endif
           r6_ij = r6_ij*f_cut_ij
        endif

        c6_eff = local_v(i)%v * local_v(j)%v * this%c6_pair(this%map(at%Z(j)),this%map(at%Z(i)))
        r0_ij = r0_eff(i) + r0_eff(j)

        d_sR_r0_ij = this%d / this%sR / r0_ij
        exp_factor = exp( -d_sR_r0_ij*r_ij + this%d )
        f_damp = 1.0_dp / ( 1.0_dp + exp_factor )

        e_i = -0.5*f_damp*c6_eff*r6_ij
        local_e_in(i) = local_e_in(i) + e_i

        if(do_derivative) then
           df_damp_dr_ij = f_damp**2 * exp_factor * d_sR_r0_ij
           df_damp_dr0_ij = -df_damp_dr_ij * r_ij / r0_ij

           df_damp_dvi = df_damp_dr0_ij * dr0_eff_dv(i)
           df_damp_dvj = df_damp_dr0_ij * dr0_eff_dv(j)

           dc6_eff_dvi = local_v(j)%v * this%c6_pair(this%map(at%Z(j)),this%map(at%Z(i)))
           dc6_eff_dvj = local_v(i)%v * this%c6_pair(this%map(at%Z(j)),this%map(at%Z(i)))

           de_i_dr_ij = 0.5_dp*(f_damp*c6_eff*dr6_ij_dr_ij + df_damp_dr_ij*c6_eff*r6_ij)
           f_ij = de_i_dr_ij*u_ij
           if(present(f)) then
              f_in(:,i) = f_in(:,i) - f_ij
              f_in(:,j) = f_in(:,j) + f_ij
           endif
           if(present(virial)) then
              virial_in(:,:,j) = virial_in(:,:,j) + ( f_ij .outer. u_ij ) * r_ij
           endif

           do m = local_v(i)%n_lo, local_v(i)%n_up
              k = local_v(i)%j(m)
              de_i_drk = (df_damp_dvi*c6_eff + f_damp*dc6_eff_dvi)*r6_ij* &
                 local_v(i)%deriv(:,m)
              if(present(f)) f_in(:,k) = f_in(:,k) + de_i_drk
              if(present(virial)) then
                 virial_in(:,:,k) = virial_in(:,:,k) + (de_i_drk .outer. local_v(i)%diff(:,m))
              endif
           enddo
        endif

     enddo
  enddo
!$omp end do
!$omp end parallel 

  if(present(f)) f = f_in
  if(present(e)) e = sum(local_e_in)
  if(present(local_e)) local_e = local_e_in
  if(present(virial)) virial = sum(virial_in,dim=3)

  if(present(local_virial)) then
     do i = 1, at%N
        local_virial(:,i) = reshape(virial_in(:,:,i),(/9/))
     enddo
  endif

  if(allocated(local_e_in)) deallocate(local_e_in)
  if(allocated(f_in)) deallocate(f_in)
  if(allocated(virial_in)) deallocate(virial_in)

  if (present(mpi)) then
     if( mpi%active ) then
        if(present(f)) call sum_in_place(mpi,f)
        if(present(virial)) call sum_in_place(mpi,virial)
        if(present(local_virial)) call sum_in_place(mpi,local_virial)
        if(present(e)) e = sum(mpi,e)
        if(present(local_e) ) call sum_in_place(mpi,local_e)

        call remove_property(at,'mpi_local_mask', error=error)
        deallocate(mpi_local_mask)
     endif
  endif

!  if(present(e)) then
!     if( associated(atom_mask_pointer) ) then
!        e = e + sum(this%e0(at%Z),mask=atom_mask_pointer)
!     else
!        e = e + sum(this%e0(at%Z))
!     endif
!  endif
!
!  if(present(local_e)) then
!     if( associated(atom_mask_pointer) ) then
!        where (atom_mask_pointer) local_e = local_e + this%e0(at%Z)
!     else
!        local_e = local_e + this%e0(at%Z)
!     endif
!  endif

  if(present(f)) f = this%E_scale * f
  if(present(e)) e = this%E_scale * e
  if(present(local_e)) local_e = this%E_scale * local_e
  if(present(virial)) virial = this%E_scale * virial
  if(present(local_virial)) local_virial = this%E_scale * local_virial
  
  atom_mask_pointer => null()
  call finalise(my_args_str)

  if(allocated(local_v)) then
     do i = 1, at%N
        local_v(i)%v = 0.0_dp
        if(allocated(local_v(i)%deriv)) deallocate(local_v(i)%deriv)
        if(allocated(local_v(i)%j)) deallocate(local_v(i)%j)
        if(allocated(local_v(i)%diff)) deallocate(local_v(i)%diff)
     enddo
     deallocate(local_v)
  endif
  if( allocated(r0_eff) ) deallocate(r0_eff)
  if( allocated(dr0_eff_dv) ) deallocate(dr0_eff_dv)

!        if(present(e) .or. present(local_e)) then
!
!           e_i_cutoff = e_i * my_descriptor_data%x(i)%covariance_cutoff / size(my_descriptor_data%x(i)%ci)
!           call print("vdWDEBUG ci="//my_descriptor_data%x(i)%ci//" e_i="//e_i//" e_i_cutoff="//e_i_cutoff, PRINT_NERD)
!
!           do n = 1, size(my_descriptor_data%x(i)%ci)
!              local_e_in( my_descriptor_data%x(i)%ci(n) ) = local_e_in( my_descriptor_data%x(i)%ci(n) ) + e_i_cutoff
!           enddo
!        endif
!
!        if(present(f) .or. present(virial) .or. present(local_virial)) then
!           i_pos0 = lbound(my_descriptor_data%x(i)%ii,1)
!
!           do n = lbound(my_descriptor_data%x(i)%ii,1), ubound(my_descriptor_data%x(i)%ii,1)
!              if( .not. my_descriptor_data%x(i)%has_grad_data(n) ) cycle
!              j = my_descriptor_data%x(i)%ii(n)
!              pos = my_descriptor_data%x(i)%pos(:,n)
!              f_gp = matmul( gradPredict,my_descriptor_data%x(i)%grad_data(:,:,n)) * my_descriptor_data%x(i)%covariance_cutoff + &
!              e_i * my_descriptor_data%x(i)%grad_covariance_cutoff(:,n)
!              if( present(f) ) then
!                 f_in(:,j) = f_in(:,j) - f_gp
!              endif
!              if( do_local_gap_variance ) then
!                 gap_variance_gradient_in(:,j) = gap_variance_gradient_in(:,j) + & 
!                    matmul( grad_variance_estimate, my_descriptor_data%x(i)%grad_data(:,:,n)) * my_descriptor_data%x(i)%covariance_cutoff**2 + &
!                    2.0_dp * gap_variance(i) * my_descriptor_data%x(i)%covariance_cutoff * my_descriptor_data%x(i)%grad_covariance_cutoff(:,n)
!              endif
!              if( present(virial) .or. present(local_virial) ) then
!                 virial_i = ((pos-my_descriptor_data%x(i)%pos(:,i_pos0)) .outer. f_gp)
!                 virial_in(:,:,j) = virial_in(:,:,j) - virial_i
!
!                 !virial_i = (pos .outer. f_gp) / size(my_descriptor_data%x(i)%ci)
!                 !do m = 1, size(my_descriptor_data%x(i)%ci)
!                 !   virial_in(:,:,my_descriptor_data%x(i)%ci(m)) = virial_in(:,:,my_descriptor_data%x(i)%ci(m)) - virial_i
!                 !enddo
!              endif
!           enddo
!        endif
!     enddo loop_over_descriptor_instances
!     call system_timer('IPModel_vdW_Calc_gp_predict')

#endif

end subroutine IPModel_vdW_Calc

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!% XML param reader functions.
!% An example for XML stanza is given below, please notice that
!% they are simply dummy parameters for testing purposes, with no physical meaning.
!%
!%> <vdW_params datafile="file" label="default">
!%> </vdW_params>
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

subroutine IPModel_startElement_handler(URI, localname, name, attributes)
  character(len=*), intent(in)   :: URI  
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name 
  type(dictionary_t), intent(in) :: attributes

  integer :: status
  character(len=1024) :: value

  integer :: ri, Z

  if(name == 'vdW_params') then ! new vdW stanza
     
     if(parse_in_ip) &
        call system_abort("IPModel_startElement_handler entered vdW_params with parse_in true. Probably a bug in FoX (4.0.1, e.g.)")
     
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
        parse_ip%label = trim(value) ! if we found a label, AND didn't have one originally, pass it back to the object.
     endif

     if(parse_in_ip) then
        if(parse_ip%initialised) call finalise(parse_ip)
     endif


     call QUIP_FoX_get_value(attributes, 'gap_version', value, status)
     if( (status == 0) ) then
        parse_ip%xml_version = string_to_int(value)
        if( parse_ip%xml_version > gap_version ) &
        call system_abort( &
           'Database was created with a later version of the code.' // &
           'Version of code used to generate the database is '//trim(value)//'.'// &
           'Version of current code is '//gap_version//'. Please update your code.')
     else
        parse_ip%xml_version = 0
     endif

     call QUIP_FoX_get_value(attributes, 'sR', value, status)
     if(status == 0) then
        read (value, *) parse_ip%sR
     else
        call system_abort('IPModel_vdW_read_params_xml cannot find sR')
     endif

     call QUIP_FoX_get_value(attributes, 'd', value, status)
     if(status == 0) then
        read (value, *) parse_ip%d
     else
        call system_abort('IPModel_vdW_read_params_xml cannot find d')
     endif

     call QUIP_FoX_get_value(attributes, 'cutoff', value, status)
     if(status == 0) then
        read (value, *) parse_ip%vdW_cutoff
     else
        call system_abort('IPModel_vdW_read_params_xml cannot find cutoff')
     endif

     call QUIP_FoX_get_value(attributes, 'cutoff_buffer', value, status)
     if(status == 0) then
        read (value, *) parse_ip%vdW_cutoff_buffer
     else
        call system_abort('IPModel_vdW_read_params_xml cannot find cutoff_buffer')
     endif

     ! For good measure
     parse_ip%map = -1

  elseif(parse_in_ip .and. name == 'GAP_data') then

     call QUIP_FoX_get_value(attributes, 'e0', value, status)
     if(status == 0) then
        read (value, *) parse_ip%e0(1)
        parse_ip%e0 = parse_ip%e0(1)
     endif

     parse_in_vdW_data = .true.

  elseif(parse_in_ip .and. parse_in_vdW_data .and. name == 'e0') then
     call QUIP_FoX_get_value(attributes, 'Z', value, status)
     if(status == 0) then
        read (value, *) Z
     else
        call system_abort('IPModel_vdW_read_params_xml cannot find Z')
     endif
     if( Z > size(parse_ip%e0) ) call system_abort('IPModel_vdW_read_params_xml: attribute Z = '//Z//' > '//size(parse_ip%e0))

     call QUIP_FoX_get_value(attributes, 'value', value, status)
     if(status == 0) then
        read (value, *) parse_ip%e0(Z)
     else
        call system_abort('IPModel_vdW_read_params_xml cannot find value in e0')
     endif

  elseif(parse_in_ip .and. name == 'per_type_data') then
     parse_ip%n_types = parse_ip%n_types + 1
     call reallocate(parse_ip%Z,parse_ip%n_types,copy=.true.)
     call reallocate(parse_ip%r0,parse_ip%n_types,copy=.true.)
     call reallocate(parse_ip%alpha0,parse_ip%n_types,copy=.true.)
     call reallocate(parse_ip%c6,parse_ip%n_types,copy=.true.)

     call QUIP_FoX_get_value(attributes, 'Z', value, status)
     if(status == 0) then
        read (value, *) parse_ip%Z(parse_ip%n_types)
     else
        call system_abort('IPModel_vdW_read_params_xml cannot find Z')
     endif

     call QUIP_FoX_get_value(attributes, 'r0', value, status)
     if(status == 0) then
        read (value, *) parse_ip%r0(parse_ip%n_types)
     else
        call system_abort('IPModel_vdW_read_params_xml cannot find r0')
     endif

     call QUIP_FoX_get_value(attributes, 'alpha0', value, status)
     if(status == 0) then
        read (value, *) parse_ip%alpha0(parse_ip%n_types)
     else
        call system_abort('IPModel_vdW_read_params_xml cannot find alpha0')
     endif

     call QUIP_FoX_get_value(attributes, 'c6', value, status)
     if(status == 0) then
        read (value, *) parse_ip%c6(parse_ip%n_types)
     else
        call system_abort('IPModel_vdW_read_params_xml cannot find c6')
     endif
  endif

end subroutine IPModel_startElement_handler

subroutine IPModel_endElement_handler(URI, localname, name)
  character(len=*), intent(in)   :: URI  
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name 

  if (parse_in_ip) then
    if(name == 'vdW_params') then
       parse_in_ip = .false.
       parse_in_ip_done = .true.
    elseif(name == 'vdW_data') then
       parse_in_vdW_data = .false.

    elseif(name == 'per_type_data') then

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

subroutine IPModel_vdW_read_params_xml(this, param_str)
  type(IPModel_vdW), intent(inout), target :: this
  character(len=*), intent(in) :: param_str

  type(xml_t) :: fxml

  if (len(trim(param_str)) <= 0) then
     call system_abort('IPModel_vdW_read_params_xml: invalid param_str length '//len(trim(param_str)) )
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
     call  system_abort('IPModel_vdW_read_params_xml: could not initialise vdW potential. No vdW_params present?')
     this%initialised = .true.
  endif

end subroutine IPModel_vdW_read_params_xml

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!% Printing of vdW parameters: number of different types, cutoff radius, atomic numbers, etc.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

subroutine IPModel_vdW_Print (this, file, dict)
  type(IPModel_vdW), intent(inout) :: this
  type(Inoutput), intent(inout),optional :: file
  type(Dictionary), intent(inout), optional :: dict
  integer :: i

#ifdef HAVE_GAP
  call Print("IPModel_vdW : Gaussian Approximation Potential", file=file)
  call Print("IPModel_vdW : label = "//this%label, file=file)
  call Print("IPModel_vdW : cutoff = "//this%cutoff, file=file)
  call Print("IPModel_vdW : E_scale = "//this%E_scale, file=file)
  call Print("IPModel_vdW : command_line = "//string(this%command_line),file=file)

#else
#endif

end subroutine IPModel_vdW_Print

end module IPModel_vdW_module
