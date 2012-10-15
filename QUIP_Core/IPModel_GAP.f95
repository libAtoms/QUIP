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

  real(dp) :: E_scale = 0.0_dp                                 !% scale factor for the potential 

  ! bispectrum parameters
  integer :: j_max = 0
  real(dp) :: z0 = 0.0_dp
  integer :: n_species = 0                                       !% Number of atomic types.
  integer, dimension(:), allocatable :: Z
  real(dp), dimension(116) :: z_eff = 0.0_dp
  real(dp), dimension(116) :: w_Z = 1.0_dp
  real(dp) :: e0 = 0.0_dp

  ! qw parameters
  integer :: qw_l_max = 0
  integer :: qw_f_n = 0
  logical :: qw_do_q = .false.
  logical :: qw_do_w = .false.
  real(dp), allocatable :: qw_cutoff(:)
  integer, allocatable :: qw_cutoff_f(:)
  real(dp), allocatable :: qw_cutoff_r1(:)

  integer :: cosnx_l_max, cosnx_n_max

  real(dp), dimension(:), allocatable :: pca_mean, NormFunction
  real(dp), dimension(:,:), allocatable :: pca_matrix, RadialTransform

  logical :: do_pca = .false.

  logical :: do_lammps = .false., do_fast_gap = .false., do_fast_gap_dgemv = .false.

  character(len=256) :: coordinates             !% Coordinate system used in GAP database

  character(len=STRING_LENGTH) :: label

#ifdef HAVE_GP_PREDICT
  type(gpSparse) :: my_gp
  type(descriptor), dimension(:), allocatable :: my_descriptor
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

  integer :: i_coordinate

  call Finalise(this)

  ! now initialise the potential
#ifndef HAVE_GP_PREDICT
  call system_abort('IPModel_GAP_Initialise_str: compiled without HAVE_GP_PREDICT')
#else

  call initialise(params)
  this%label=''

  call param_register(params, 'label', '', this%label, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'lammps', 'false', this%do_lammps, help_string="The way forces on local/ghost atoms are calculated. " // &
   "By default, cross terms are accumulated on the central atom. If >>lammps<< is set, cross-terms are accumulated on the neighbours.")
  call param_register(params, 'fast_gap', 'false', this%do_fast_gap, help_string="If true, activate various speedups to gp_predict.  " // &
   "Switches >>lammps<< to true, which may change amount of buffer needed around atoms in atom_mask_name")
  call param_register(params, 'fast_gap_dgemv', 'false', this%do_fast_gap_dgemv, help_string="If true, used dgemv instead of many " // &
   "dot_product() calls in GAP forces. Introduces some overhead in sorting of sparse points by Z.")
  call param_register(params, 'E_scale', '1.0', this%E_scale, help_string="rescaling factor for the potential")

  if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_SW_Initialise_str args_str')) &
  call system_abort("IPModel_GAP_Initialise_str failed to parse label from args_str="//trim(args_str))
  call finalise(params)

  call IPModel_GAP_read_params_xml(this, param_str)
  call gp_readXML(this%my_gp, param_str,label=trim(this%label))
  allocate(this%my_descriptor(this%my_gp%n_coordinate))

  this%cutoff = 0.0_dp
  do i_coordinate = 1, this%my_gp%n_coordinate
     call initialise(this%my_descriptor(i_coordinate),string(this%my_gp%coordinate(i_coordinate)%descriptor_str))
     this%cutoff = max(this%cutoff,cutoff(this%my_descriptor(i_coordinate)))
  enddo

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
  type(IPModel_GAP), intent(inout) :: this
  type(Atoms), intent(in) :: at
  real(dp), intent(out), optional :: e, local_e(:) !% \texttt{e} = System total energy, \texttt{local_e} = energy of each atom, vector dimensioned as \texttt{at%N}.  
  real(dp), intent(out), optional :: f(:,:), local_virial(:,:)   !% Forces, dimensioned as \texttt{f(3,at%N)}, local virials, dimensioned as \texttt{local_virial(9,at%N)} 
  real(dp), intent(out), optional :: virial(3,3)   !% Virial
  character(len=*), intent(in), optional :: args_str 
  type(MPI_Context), intent(in), optional :: mpi
  integer, intent(out), optional :: error

#ifdef HAVE_GP_PREDICT
  real(dp), pointer :: w_e(:)
  real(dp) :: e_i, f_gp_k, water_monomer_energy
  real(dp), dimension(:), allocatable   :: local_e_in, w
  real(dp), dimension(:,:,:), allocatable   :: virial_in
  real(dp), dimension(:,:,:,:), allocatable   :: dvec
  real(dp), dimension(:,:), allocatable   :: vec
  real(dp), dimension(:,:,:), allocatable   :: jack
  integer, dimension(:,:), allocatable :: water_monomer_index
  integer :: d, i, j, k, l, m, n, nei_max, iAo, iBo, n_water_pair

  real(dp), dimension(:,:), allocatable :: covariance, pca_matrix, f_in

  integer, dimension(3) :: shift
  real(dp), dimension(3) :: pos, f_gp
  type(Dictionary) :: params
  logical, dimension(:), pointer :: atom_mask_pointer
  logical :: has_atom_mask_name, do_atom_mask_lookup, do_lammps, do_fast_gap, do_fast_gap_dgemv, force_calc_new_x_star, new_x_star
  character(STRING_LENGTH) :: atom_mask_name
  real(dp) :: r_scale, E_scale

  real(dp), dimension(:), allocatable :: sparseScore
  logical :: do_rescale_r, do_rescale_E, do_sparseScore

  integer :: n_atoms_eff, ii, ji, jj, i_coordinate
  integer, allocatable, dimension(:) :: atom_mask_lookup, atom_mask_reverse_lookup

  type(descriptor_data) :: my_descriptor_data
  real(dp), dimension(:), allocatable :: gradPredict

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

  ! Has to be allocated as it's in the reduction clause.
  allocate(local_e_in(at%N))
  local_e_in = this%e0 !0.0_dp

  allocate(f_in(3,at%N))
  f_in = 0.0_dp

  allocate(virial_in(3,3,at%N))
  virial_in = 0.0_dp

  if (.not. assign_pointer(at, "weight", w_e)) nullify(w_e)

  atom_mask_pointer => null()

  do_lammps = .false.
  do_sparseScore = .false.
  if(present(args_str)) then
     call initialise(params)
     
     call param_register(params, 'atom_mask_name', 'NONE',atom_mask_name,has_value_target=has_atom_mask_name, &
     help_string="Name of a logical property in the atoms object. For atoms where this property is true, energies, forces, virials etc. are " // &
      "calculated.  Cross-terms are attributed to atoms depending on the >>lammps<< flag.")
     call param_register(params, 'lammps', ""//this%do_lammps, do_lammps, help_string="The way forces on local/ghost atoms are calculated. " // &
      "By default, cross terms are accumulated on the central atom. If >>lammps<< is set, cross-terms are accumulated on the neighbours.")
     call param_register(params, 'fast_gap', ""//this%do_fast_gap, do_fast_gap, help_string="If true, activate various speedups to gp_predict.  " // &
      "Switches >>lammps<< to true, which may change amount of buffer needed around atoms in atom_mask_name")
     call param_register(params, 'fast_gap_dgemv', ""//this%do_fast_gap_dgemv, do_fast_gap_dgemv, help_string="If true, used dgemv instead of many " // &
      "dot_product() calls in GAP forces. Introduces some overhead in sorting of sparse points by Z.")
     call param_register(params, 'r_scale', '1.0',r_scale, has_value_target=do_rescale_r, help_string="Recaling factor for distances. Default 1.0.")
     call param_register(params, 'E_scale', '1.0',E_scale, has_value_target=do_rescale_E, help_string="Recaling factor for energy. Default 1.0.")
     call param_register(params, 'sparseScore', 'F', do_sparseScore, help_string="Compute score for each test point.")

     if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='IPModel_GAP_Calc args_str')) &
     call system_abort("IPModel_GAP_Calc failed to parse args_str='"//trim(args_str)//"'")
     call finalise(params)


     if( has_atom_mask_name ) then
        if (.not. assign_pointer(at, trim(atom_mask_name) , atom_mask_pointer)) &
        call system_abort("IPModel_GAP_Calc did not find "//trim(atom_mask_name)//" propery in the atoms object.")
     else
        atom_mask_pointer => null()
     endif
     if (do_rescale_r .or. do_rescale_E) then
        RAISE_ERROR("IPModel_GAP_Calc: rescaling of potential at the calc() stage with r_scale and E_scale not yet implemented!", error)
     end if

  endif

  do i_coordinate = 1, this%my_gp%n_coordinate

     d = descriptor_dimensions(this%my_descriptor(i_coordinate))

     if(do_sparseScore) then
        call gpCoordinates_initialise_SparseScore(this%my_gp%coordinate(i_coordinate))
     endif

     if(present(f) .or. present(virial) .or. present(local_virial)) allocate(gradPredict(d))

     call calc(this%my_descriptor(i_coordinate),at,my_descriptor_data, &
        do_descriptor=.true.,do_grad_descriptor=present(f) .or. present(virial) .or. present(local_virial))

     if(do_sparseScore) allocate(sparseScore(size(my_descriptor_data%x)))

!$omp parallel default(none) private(i,gradPredict, e_i,n,j,pos,f_gp) shared(this,at,i_coordinate,my_descriptor_data,e,virial,local_virial,local_e,do_sparseScore,sparseScore,f) reduction(+:local_e_in,f_in,virial_in)

!$omp do
     do i = 1, size(my_descriptor_data%x)
        if( .not. my_descriptor_data%x(i)%has_data ) cycle

        call system_timer('IPModel_GAP_Calc_gp_predict')

        if(present(f) .or. present(virial) .or. present(local_virial)) then
           call reallocate(gradPredict,size(my_descriptor_data%x(i)%data(:)),zero=.true.)
           e_i =  gp_predict(this%my_gp%coordinate(i_coordinate) , xStar=my_descriptor_data%x(i)%data(:), gradPredict =  gradPredict, sparseScore=sparseScore(i), do_sparseScore=do_sparseScore)
        else
           e_i =  gp_predict(this%my_gp%coordinate(i_coordinate) , xStar=my_descriptor_data%x(i)%data(:), sparseScore=sparseScore(i), do_sparseScore=do_sparseScore)
        endif
        call system_timer('IPModel_GAP_Calc_gp_predict')
        if(present(e) .or. present(local_e)) then
           local_e_in( min(i,size(local_e_in)) ) = local_e_in( min(i,size(local_e_in)) ) + &
                e_i * my_descriptor_data%x(i)%covariance_cutoff
        endif
        if(present(f) .or. present(virial) .or. present(local_virial)) then
           do n = lbound(my_descriptor_data%x(i)%ii,1), ubound(my_descriptor_data%x(i)%ii,1)
              if( .not. my_descriptor_data%x(i)%has_grad_data(n) ) cycle
              j = my_descriptor_data%x(i)%ii(n)
              pos = my_descriptor_data%x(i)%pos(:,n)
              f_gp = matmul( gradPredict,my_descriptor_data%x(i)%grad_data(:,:,n)) * my_descriptor_data%x(i)%covariance_cutoff + &
              e_i * my_descriptor_data%x(i)%grad_covariance_cutoff(:,n)
              if( present(f) ) then
                 f_in(:,j) = f_in(:,j) - f_gp
                 if(my_descriptor_data%atom_based) f_in(:,i) = f_in(:,i) + f_gp
              endif
              if( present(virial) .or. present(local_virial) ) then
                 virial_in(:,:,j) = virial_in(:,:,j) - (pos .outer. f_gp)
                 if(my_descriptor_data%atom_based) virial_in(:,:,i) = virial_in(:,:,i) + (at%pos(:,i) .outer. f_gp)
              endif
           enddo
        endif
     enddo
!$omp end do
     if(allocated(gradPredict)) deallocate(gradPredict)
!$omp end parallel

     if(do_sparseScore) then
        do i = 1, size(my_descriptor_data%x)
           call print('DESCRIPTOR '//i//' SPARSE_SCORE = '//sparseScore(i))
        enddo
        deallocate(sparseScore)
     endif

     call finalise(my_descriptor_data)

  enddo

  if(present(f)) f = f_in

  if (present(mpi)) then
     if(present(f)) call sum_in_place(mpi,f)
     if(present(virial) .or. present(local_virial)) call sum_in_place(mpi,virial_in)
     if(present(e) .or. present(local_e) ) call sum_in_place(mpi,local_e_in)
  endif

  if(present(f)) f = this%E_scale*f
  if(present(e)) e = this%E_scale*sum(local_e_in)
  if(present(local_e)) local_e = this%E_scale*local_e_in
  if(present(virial)) virial = this%E_scale*sum(virial_in,dim=3)

  if(present(local_virial)) then
     do i = 1, at%N
        local_virial(:,i) = this%E_scale*reshape(virial_in(:,:,i),(/9/))
     enddo
  endif

  if(allocated(local_e_in)) deallocate(local_e_in)
  if(allocated(f_in)) deallocate(f_in)
  if(allocated(virial_in)) deallocate(virial_in)
  if(allocated(w)) deallocate(w)
  if(allocated(covariance)) deallocate(covariance)
  if(allocated(dvec)) deallocate(dvec)
  if(allocated(vec)) deallocate(vec)
  if(allocated(jack)) deallocate(jack)
  if(allocated(water_monomer_index)) deallocate(water_monomer_index)
  atom_mask_pointer => null()
  if (allocated(atom_mask_lookup)) deallocate(atom_mask_lookup)
  if (allocated(atom_mask_reverse_lookup)) deallocate(atom_mask_reverse_lookup)

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
        parse_ip%label = trim(value) ! if we found a label, AND didn't have one originally, pass it back to the object.
     endif

     if(parse_in_ip) then
        if(parse_ip%initialised) call finalise(parse_ip)
     endif

     call QUIP_FoX_get_value(attributes, 'svn_version', value, status)
     if( (status == 0) .and. (string_to_int(value) > current_version()) ) then
	call system_abort( &
	   'Database was created with a different version of the code.' // &
	   'Version of code used to generate the database is '//trim(value) //'.'// &
	   'Version of current code is '//current_version() // &
	   'Please update your code.')
     endif

  elseif(parse_in_ip .and. name == 'GAP_data') then

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

     allocate( parse_ip%Z(parse_ip%n_species) )

  elseif(parse_in_ip .and. name == 'water_monomer_params') then

     call QUIP_FoX_get_value(attributes, 'cutoff', value, status)
     if(status == 0) then
        read (value, *) parse_ip%cutoff
     else
        call system_abort('IPModel_GAP_read_params_xml cannot find cutoff')
     endif

  elseif(parse_in_ip .and. name == 'hf_dimer_params') then

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

  elseif(parse_in_ip .and. name == 'cosnx_params') then
  
     call QUIP_FoX_get_value(attributes, 'l_max', value, status)
     if(status == 0) then
        read (value, *) parse_ip%cosnx_l_max
     else
        call system_abort('IPModel_GAP_read_params_xml cannot find l_max')
     endif
  
     call QUIP_FoX_get_value(attributes, 'n_max', value, status)
     if(status == 0) then
        read (value, *) parse_ip%cosnx_n_max
     else
        call system_abort('IPModel_GAP_read_params_xml cannot find n_max')
     endif
  
     call QUIP_FoX_get_value(attributes, 'cutoff', value, status)
     if(status == 0) then
        read (value, *) parse_ip%cutoff
     else
        call system_abort('IPModel_GAP_read_params_xml cannot find cutoff')
     endif
  
     allocate(parse_ip%NormFunction(parse_ip%cosnx_n_max), parse_ip%RadialTransform(parse_ip%cosnx_n_max,parse_ip%cosnx_n_max))

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

  elseif(parse_in_ip .and. name == 'NormFunction') then
  
     parse_n_row = parse_ip%cosnx_n_max
     call zero(parse_cur_data)
  
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

    elseif(name == 'hf_dimer_params') then

    elseif(name == 'water_monomer_params') then

    elseif(name == 'water_dimer_params') then

    elseif(name == 'qw_so3_params') then

    elseif(name == 'radial_function') then

    elseif(name == 'per_type_data') then

    elseif(name == 'PCA_mean') then
       
       val = string(parse_cur_data)
       read(val,*) parse_ip%pca_mean

    elseif(name == 'row') then

       val = string(parse_cur_data)
       read(val,*) parse_ip%pca_matrix(:,parse_cur_row)

    elseif(name == 'NormFunction') then
       
       val = string(parse_cur_data)
       read(val,*) parse_ip%NormFunction
    
    elseif(name == 'RadialTransform_row') then
    
       val = string(parse_cur_data)
       read(val,*) parse_ip%RadialTransform(:,parse_cur_row)

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
  call Print("IPModel_GAP : E_scale = "//this%E_scale, file=file)
  call Print("IPModel_GAP : n_species = "//this%n_species, file=file)

  do i = 1, this%n_species
     call Print("IPModel_GAP : Z = "//this%Z(i), file=file)
     call Print("IPModel_GAP : z_eff = "//this%z_eff(this%Z(i)), file=file)
!     call Print("IPModel_GAP : delta = "//this%my_gp%delta(i), file=file)
!     call Print("IPModel_GAP : theta = "//this%my_gp%theta(:,i), file=file)
  enddo

  call Print("IPModel_GAP : command_line = "//string(this%command_line))
#endif

end subroutine IPModel_GAP_Print

end module IPModel_GAP_module
