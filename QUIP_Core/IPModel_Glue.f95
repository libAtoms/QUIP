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
!X IPModel_Glue
!X
!% Generic implementation of Glue potentials. Glue potentials are described
!% in the XML params file.
!X
!X Copyright Tim Green 2010. Licensed under the Gnu Public License version 3.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#include "error.inc"

module IPModel_Glue_module

use error_module
use system_module, only : dp, inoutput, print, verbosity_push_decrement, verbosity_pop, operator(//)
use dictionary_module
use paramreader_module
use linearalgebra_module
use spline_module
use atoms_types_module
use atoms_module

use mpi_context_module
use QUIP_Common_module

implicit none
private

include 'IPModel_interface.h'

type SplineDataContainer
    real(dp), allocatable, dimension(:,:) :: data
    real(dp) :: y1, yn
end type SplineDataContainer

public :: IPModel_Glue
type IPModel_Glue
  integer :: n_types = 0 ! Number of different species we have

  ! Maps from atomic number to type and vice versa
  integer, allocatable :: atomic_num(:), type_of_atomic_num(:)

  real(dp) :: cutoff = 0.0_dp
  real(dp), allocatable, dimension(:) :: poly

  type(SplineDataContainer), allocatable, dimension(:) :: spline_data_potential, spline_data_density
  type(SplineDataContainer), allocatable, dimension(:,:) :: spline_data_pair
  type(Spline), allocatable :: potential(:)
  type(Spline), allocatable :: pair(:,:)
  type(Spline), allocatable :: density(:)

  character(len=STRING_LENGTH) :: label

end type IPModel_Glue

logical, private :: parse_in_ip, parse_matched_label, parse_in_density_potential, parse_in_neighbours_potential, parse_in_density, parse_in_potential_pair
type(IPModel_Glue), private, pointer :: parse_ip
integer :: parse_curr_type_i = 0, parse_curr_type_j = 0, parse_curr_point, parse_n_neighbours

interface Initialise
  module procedure IPModel_Glue_Initialise_str
end interface Initialise

interface Finalise
  module procedure IPModel_Glue_Finalise
end interface Finalise

interface Print
  module procedure IPModel_Glue_Print
end interface Print

interface Calc
  module procedure IPModel_Glue_Calc
end interface Calc

contains

subroutine IPModel_Glue_Initialise_str(this, args_str, param_str)
  type(IPModel_Glue), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str

  type(Dictionary) :: params

  call Finalise(this)

  call initialise(params)
  this%label=''
  call param_register(params, 'label', '', this%label, help_string="No help yet.  This source file was $LastChangedBy$")
  if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_Glue_Initialise_str args_str')) then
    call system_abort("IPModel_Glue_Initialise_str failed to parse label from args_str="//trim(args_str))
  endif
  call finalise(params)

  call IPModel_Glue_read_params_xml(this, param_str)

end subroutine IPModel_Glue_Initialise_str

subroutine IPModel_Glue_Finalise(this)
  type(IPModel_Glue), intent(inout) :: this
  integer :: i, j

  ! Add finalisation code here

  if (allocated(this%atomic_num)) deallocate(this%atomic_num)
  if (allocated(this%type_of_atomic_num)) deallocate(this%type_of_atomic_num)
  
  if (allocated(this%spline_data_density)) then
     do i=1,this%n_types
        if(allocated(this%spline_data_density(i)%data)) deallocate(this%spline_data_density(i)%data)
     enddo
     deallocate(this%spline_data_density)
  endif

  if (allocated(this%spline_data_potential)) then
     do i=1,this%n_types
        if(allocated(this%spline_data_potential(i)%data)) deallocate(this%spline_data_potential(i)%data)
     enddo
     deallocate(this%spline_data_potential)
  endif

  if (allocated(this%spline_data_pair)) then
     do i=1,this%n_types
        do j=1,this%n_types
           if(allocated(this%spline_data_pair(j,i)%data)) deallocate(this%spline_data_pair(j,i)%data)
        enddo
     enddo
     deallocate(this%spline_data_pair)
  endif

  if(allocated(this%poly)) deallocate(this%poly)

  if (allocated(this%potential)) then
     do i=1,this%n_types
        call finalise(this%potential(i))
     enddo
     deallocate(this%potential)
  endif

  if (allocated(this%density)) then
     do i=1,this%n_types
        call finalise(this%density(i))
     enddo
     deallocate(this%density)
  endif

  if (allocated(this%pair)) then
     do i=1,this%n_types
        do j = 1, this%n_types
            call finalise(this%pair(j,i))
        enddo
     enddo
     deallocate(this%pair)
  endif

  this%n_types = 0
  this%label = ''
end subroutine IPModel_Glue_Finalise

subroutine IPModel_Glue_Calc(this, at, e, local_e, f, virial, local_virial, args_str, mpi, error)
  type(IPModel_Glue), intent(inout):: this
  type(Atoms), intent(inout)      :: at
  real(dp), intent(out), optional :: e, local_e(:)
  real(dp), intent(out), optional :: f(:,:), local_virial(:,:)   !% Forces, dimensioned as \texttt{f(3,at%N)}, local virials, dimensioned as \texttt{local_virial(9,at%N)} 
  real(dp), intent(out), optional :: virial(3,3)
  character(len=*), optional      :: args_str
  type(MPI_Context), intent(in), optional :: mpi
  integer, intent(out), optional :: error

  type(Dictionary) :: params
  logical, dimension(:), pointer :: atom_mask_pointer
  real(dp), dimension(:), allocatable :: local_e_in, rho_local 
  real(dp), dimension(:,:), allocatable :: f_in
  real(dp), dimension(:,:,:), allocatable :: local_virial_in
  logical :: has_atom_mask_name
  character(STRING_LENGTH) :: atom_mask_name
  real(dp) :: r_scale, E_scale
  logical :: do_rescale_r, do_rescale_E
  ! Loop variables
  integer :: i, ji, j, ti, tj  

  real(dp) :: r_ij_mag, r_ij_hat(3), pair_e_ij, dpair_e_ij ! Neighbour vector info

  ! For calculating forces and the virial tensor
  real(dp) :: dpotential_drho_drho_i_drij(3), dpotential_drho
    
  INIT_ERROR(error)


  if (present(e)) e = 0.0_dp
  if (present(local_e)) then
     call check_size('Local_E',local_e,(/at%N/),'IPModel_Glue_Calc', error)
  endif

  if(present(e) .or. present(local_e)) then
     allocate(local_e_in(at%N))
     local_e_in = 0.0_dp
  else
     allocate(local_e_in(1))
  endif

  if (present(f)) then
     call check_size('Force',f,(/3,at%N/),'IPModel_Glue_Calc', error)
     allocate(f_in(3,at%N))
     f_in = 0.0_dp
  else
     allocate(f_in(1,1))
  end if

  if (present(local_virial)) then
     call check_size('Local_virial',local_virial,(/9,at%N/),'IPModel_Glue_Calc', error)
     local_virial = 0.0_dp
  endif

  if(present(virial) .or. present(local_virial)) then
     allocate(local_virial_in(3,3,at%N))
     local_virial_in = 0.0_dp
  else
     allocate(local_virial_in(1,1,1))
  endif

  atom_mask_pointer => null()

  if(present(args_str)) then
     call initialise(params)
     call param_register(params, 'atom_mask_name', 'NONE',atom_mask_name,has_value_target=has_atom_mask_name, help_string="No help yet.  This source file was $LastChangedBy$")
     call param_register(params, 'r_scale', '1.0',r_scale, has_value_target=do_rescale_r, help_string="Recaling factor for distances. Default 1.0.")
     call param_register(params, 'E_scale', '1.0',E_scale, has_value_target=do_rescale_E, help_string="Recaling factor for energy. Default 1.0.")

     if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='IPModel_Glue_Calc args_str')) then
        RAISE_ERROR("IPModel_Glue_Calc failed to parse args_str='"//trim(args_str)//"'", error)
     endif
     call finalise(params)


     if( has_atom_mask_name ) then
        if (.not. assign_pointer(at, trim(atom_mask_name) , atom_mask_pointer)) then
           RAISE_ERROR("IPModel_Glue_Calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
        endif
     else
        atom_mask_pointer => null()
     endif
     if (do_rescale_r .or. do_rescale_E) then
        RAISE_ERROR("IPModel_Glue_Calc: rescaling of potential with r_scale and E_scale not yet implemented!", error)
     end if
  endif

  ! Iterate over atoms
  allocate(rho_local(at%N))
  rho_local = 0.0_dp

!$omp parallel do default(none) shared(this,at,atom_mask_pointer,rho_local) private(i,j,ji,ti,tj,r_ij_mag)
  do i = 1, at%N
     if(associated(atom_mask_pointer)) then
        if(.not. atom_mask_pointer(i)) cycle
     endif

     ! Get the type of this species from its atomic number
     ti = get_type(this%type_of_atomic_num, at%Z(i))
 
     ! Iterate over our nighbours
     do ji = 1, n_neighbours(at, i)
        j = neighbour(at, i, ji, distance = r_ij_mag)
        tj = get_type(this%type_of_atomic_num, at%Z(j))
 
        if (r_ij_mag < glue_cutoff(this, tj)) then ! Skip atoms beyond the cutoff
           rho_local(i) = rho_local(i) + eam_density(this, tj, r_ij_mag)
        endif
     enddo
  enddo
!$omp end parallel do

  ! Iterate over atoms
!$omp parallel do default(none) shared(this,at,atom_mask_pointer,rho_local,local_e_in,e,f,virial,local_e,local_virial) &
!$omp private(i,ji,j,ti,tj,r_ij_mag,r_ij_hat,dpotential_drho, dpotential_drho_drho_i_drij, pair_e_ij, dpair_e_ij) reduction(+:f_in,local_virial_in)
  do i = 1, at%N
     if(associated(atom_mask_pointer)) then
        if(.not. atom_mask_pointer(i)) cycle
     endif

    ! Get the type of this species from its atomic number
    ti = get_type(this%type_of_atomic_num, at%Z(i))

    dpotential_drho = eam_spline_potential_deriv(this, ti, rho_local(i))

    ! Iterate over our nighbours
    do ji = 1, n_neighbours(at, i)
       j = neighbour(at, i, ji, r_ij_mag, cosines=r_ij_hat)
       tj = get_type(this%type_of_atomic_num, at%Z(j))

       if (r_ij_mag < glue_cutoff(this, tj)) then ! Skip atoms beyond the cutoff
          if( present(f) .or. present(virial) .or. present(local_virial) ) then
             dpotential_drho_drho_i_drij = dpotential_drho * eam_density_deriv(this, tj, r_ij_mag) * r_ij_hat
          endif

          if(present(f)) then
             f_in(:,j) = f_in(:,j) - dpotential_drho_drho_i_drij
             f_in(:,i) = f_in(:,i) + dpotential_drho_drho_i_drij
          endif

          if(present(virial) .or. present(local_virial)) local_virial_in(:,:,j) = local_virial_in(:,:,j) - ( dpotential_drho_drho_i_drij .outer. r_ij_hat ) * r_ij_mag
       endif

       if( r_ij_mag < pair_cutoff(this,ti,tj) ) then
          pair_e_ij = eam_spline_pair(this,ti,tj,r_ij_mag)
          if( present(local_e) .or. present(e) ) local_e_in(i) = local_e_in(i) + 0.5_dp * pair_e_ij
          if( present(f) .or. present(virial) .or. present(local_virial) ) dpair_e_ij = eam_spline_pair_deriv(this,ti, tj, r_ij_mag)
          if( present(f) ) f_in(:,i) = f_in(:,i) + dpair_e_ij * r_ij_hat
          if( present(virial) .or. present(local_virial) ) local_virial_in(:,:,j) = local_virial_in(:,:,j) - 0.5_dp * dpair_e_ij * (r_ij_hat .outer. r_ij_hat) * r_ij_mag
       endif

    enddo ! ji

    if(present(local_e) .or. present(e)) local_e_in(i) = local_e_in(i) + eam_spline_potential(this, ti, rho_local(i))
  end do ! i
!$omp end parallel do  


  if(present(e)) e = sum(local_e_in)
  if(present(f)) f = f_in
  if(present(local_e)) local_e = local_e_in
  if(present(virial)) virial = sum(local_virial_in,dim=3)
  if(present(local_virial)) local_virial = reshape(local_virial_in,(/9,at%N/))

  if(allocated(rho_local)) deallocate(rho_local)
  if(allocated(local_e_in)) deallocate(local_e_in)
  if(allocated(f_in)) deallocate(f_in)
  if(allocated(local_virial_in)) deallocate(local_virial_in)

end subroutine IPModel_Glue_Calc

function glue_cutoff(this, ti, error)
  type(IPModel_Glue), intent(in) :: this
  integer, intent(in) :: ti
  integer, intent(out), optional :: error
  real(dp) :: glue_cutoff
  
  INIT_ERROR(error)

  if( .not. this%density(ti)%initialised ) then
      RAISE_ERROR("glue_cutoff: spline not initialised", error)
  endif

  glue_cutoff = max_knot(this%density(ti))

end function glue_cutoff

function pair_cutoff(this, ti, tj, error)
  type(IPModel_Glue), intent(in) :: this
  integer, intent(in) :: ti, tj
  integer, intent(out), optional :: error
  real(dp) :: pair_cutoff
  
  INIT_ERROR(error)

  if( .not. this%pair(ti,tj)%initialised ) then
     pair_cutoff = 0.0_dp 
  else
     pair_cutoff = max_knot(this%pair(ti,tj))
  endif

end function pair_cutoff

function eam_density(this, ti, r, error)
  type(IPModel_Glue), intent(in) :: this
  integer, intent(in) :: ti
  real(dp), intent(in) :: r
  integer, intent(out), optional :: error
  real(dp) :: eam_density
  
  INIT_ERROR(error)

  if( .not. this%density(ti)%initialised ) then
      RAISE_ERROR("eam_density: spline not initialised", error)
  endif

  if (r < glue_cutoff(this, ti)) then
    !eam_density = spline_value(this%density(ti),r) 
    eam_density = this%poly(ti)*(glue_cutoff(this, ti)-r)**3
  else
    eam_density = 0.0_dp
  endif
end function eam_density

function eam_density_deriv(this, ti, r, error)
  type(IPModel_Glue), intent(in) :: this
  integer, intent(in) :: ti
  real(dp), intent(in) :: r
  integer, intent(out), optional :: error
  real(dp) :: eam_density_deriv

  INIT_ERROR(error)

  if( .not. this%density(ti)%initialised ) then
      RAISE_ERROR("eam_density_deriv: spline not initialised", error)
  endif

  if (r < glue_cutoff(this, ti)) then
    !eam_density_deriv = spline_deriv(this%density(ti),r)
    eam_density_deriv = -3.0_dp * this%poly(ti)*(glue_cutoff(this, ti)-r)**2
  else
    eam_density_deriv = 0.0_dp
  endif
end function eam_density_deriv

function eam_spline_potential(this, ti, rho, error)
  ! Evaluate the potential spline at r for the ti^th species
  type(IPModel_Glue), intent(in) :: this
  integer, intent(in) :: ti
  real(dp), intent(in) :: rho
  integer, intent(out), optional :: error
  real(dp) :: eam_spline_potential

  INIT_ERROR(error)

  if( .not. this%potential(ti)%initialised ) then
      RAISE_ERROR("eam_spline_potential: spline not initialised", error)
  endif

  eam_spline_potential = spline_value(this%potential(ti), rho)

end function eam_spline_potential

function eam_spline_potential_deriv(this, ti, rho, error)
  ! Evaluate the potential spline derivitive at r for the ti^th species
  type(IPModel_Glue), intent(in) :: this
  integer, intent(in) :: ti
  real(dp), intent(in) :: rho
  integer, intent(out), optional :: error
  real(dp) :: eam_spline_potential_deriv

  INIT_ERROR(error)

  if( .not. this%potential(ti)%initialised ) then
      RAISE_ERROR("eam_spline_potential_deriv: spline not initialised", error)
  endif

  eam_spline_potential_deriv = spline_deriv(this%potential(ti), rho)

end function eam_spline_potential_deriv

function eam_spline_pair(this, ti, tj, rho, error)
  ! Evaluate the potential spline at r for the ti^th species
  type(IPModel_Glue), intent(in) :: this
  integer, intent(in) :: ti, tj
  real(dp), intent(in) :: rho
  integer, intent(out), optional :: error
  real(dp) :: eam_spline_pair

  INIT_ERROR(error)

  if( .not. this%pair(ti,tj)%initialised ) then
     eam_spline_pair = 0.0_dp
  else
     eam_spline_pair = spline_value(this%pair(ti,tj), rho)
  endif


end function eam_spline_pair

function eam_spline_pair_deriv(this, ti, tj, rho, error)
  ! Evaluate the potential spline derivitive at r for the ti^th species
  type(IPModel_Glue), intent(in) :: this
  integer, intent(in) :: ti, tj
  real(dp), intent(in) :: rho
  integer, intent(out), optional :: error
  real(dp) :: eam_spline_pair_deriv

  INIT_ERROR(error)

  if( .not. this%pair(ti,tj)%initialised ) then
     eam_spline_pair_deriv = 0.0_dp
  else
     eam_spline_pair_deriv = spline_deriv(this%pair(ti, tj), rho)
  endif


end function eam_spline_pair_deriv

subroutine IPModel_Glue_Print(this, file)
  type(IPModel_Glue), intent(in) :: this
  type(Inoutput), intent(inout),optional :: file

  integer :: ti

  call Print("IPModel_Glue : Glue Potential", file=file)
  call Print("IPModel_Glue : n_types = " // this%n_types // " cutoff = " // this%cutoff, file=file)

  do ti=1, this%n_types
    call Print ("IPModel_Glue : type " // ti // " atomic_num " // this%atomic_num(ti), file=file)
    call verbosity_push_decrement()
    call Print ("IPModel_Glue : " // &
        "cutoff " // this%cutoff , &
         file=file)
    call verbosity_pop()
  end do

end subroutine IPModel_Glue_Print

subroutine IPModel_Glue_read_params_xml(this, param_str)
  type(IPModel_Glue), intent(inout), target :: this
  character(len=*), intent(in) :: param_str
  integer :: ti, tj
  real(dp) :: max_cutoff
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
    call system_abort("IPModel_Glue_read_params_xml parsed file, but n_types = 0")
  endif

  ! Initialise the spline data structures with the spline data
  do ti = 1, this%n_types
  
     if(allocated(this%spline_data_density(ti)%data)) then
        call initialise(this%density(ti), this%spline_data_density(ti)%data(:,1), this%spline_data_density(ti)%data(:,2), &
           this%spline_data_density(ti)%y1, this%spline_data_density(ti)%yn)
     endif
     if(allocated(this%spline_data_potential(ti)%data)) then
        call initialise(this%potential(ti), this%spline_data_potential(ti)%data(:,1), this%spline_data_potential(ti)%data(:,2), 2.0e30_dp,2.0e30_dp)
     endif

     do tj = 1, this%n_types
        if(allocated(this%spline_data_pair(tj,ti)%data)) then
           call initialise(this%pair(tj,ti), this%spline_data_pair(tj,ti)%data(:,1), this%spline_data_pair(tj,ti)%data(:,2), 0.0_dp, 0.0_dp)
           call initialise(this%pair(ti,tj), this%spline_data_pair(tj,ti)%data(:,1), this%spline_data_pair(tj,ti)%data(:,2), 0.0_dp, 0.0_dp)
        endif
     enddo

     this%poly(ti) = dot_product(this%spline_data_density(ti)%data(:,2), (glue_cutoff(this,ti) - &
        this%spline_data_density(ti)%data(:,1))**3) / &
        dot_product((glue_cutoff(this,ti)-this%spline_data_density(ti)%data(:,1))**3, &
        (glue_cutoff(this,ti)-this%spline_data_density(ti)%data(:,1))**3)
  end do

  max_cutoff = 0.0_dp
  ! Find the largest cutoff
  do ti = 1, this%n_types
     max_cutoff = max(glue_cutoff(this,ti),max_cutoff)
     do tj = 1, this%n_types
        max_cutoff = max(pair_cutoff(this,tj,ti),max_cutoff)
     enddo
  enddo
  this%cutoff = max_cutoff

end subroutine IPModel_Glue_read_params_xml

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
  character(len=STRING_LENGTH) :: value

  integer :: ti, tj, num_points
  real(dp) :: v

  if (name == 'Glue_params') then ! new Template stanza

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
        call system_abort("Can't find n_types in Glue_params")
      endif

      allocate(parse_ip%atomic_num(parse_ip%n_types))
      parse_ip%atomic_num = 0
 
      ! Allocate n_types spline data objects for further allocation
      allocate(parse_ip%spline_data_density(parse_ip%n_types))
      allocate(parse_ip%spline_data_potential(parse_ip%n_types))
      allocate(parse_ip%spline_data_pair(parse_ip%n_types,parse_ip%n_types))

      ! Allocate the splines
      allocate(parse_ip%pair(parse_ip%n_types,parse_ip%n_types))
      allocate(parse_ip%potential(parse_ip%n_types))
      allocate(parse_ip%density(parse_ip%n_types))
      allocate(parse_ip%poly(parse_ip%n_types))
    endif


  elseif (parse_in_ip .and. name == 'per_type_data') then

    call QUIP_FoX_get_value(attributes, "type", value, status)
    if (status /= 0) call system_abort ("IPModel_Glue_read_params_xml cannot find type")
    read (value, *) ti ! This is the index of this particular species

    parse_curr_type_i = ti

    call QUIP_FoX_get_value(attributes, "atomic_num", value, status)
    if (status /= 0) call system_abort ("IPModel_Glue_read_params_xml cannot find atomic_num")
    read (value, *) parse_ip%atomic_num(ti) ! This is the atomic number of this particular species

    ! Re-do the atomic number -> species index map
    if (allocated(parse_ip%type_of_atomic_num)) deallocate(parse_ip%type_of_atomic_num)
    allocate(parse_ip%type_of_atomic_num(maxval(parse_ip%atomic_num)))
    parse_ip%type_of_atomic_num = 0
    do ti=1, parse_ip%n_types
      if (parse_ip%atomic_num(ti) > 0) &
        parse_ip%type_of_atomic_num(parse_ip%atomic_num(ti)) = ti
    end do

  elseif (parse_in_ip .and. name == 'per_pair_data') then

    call QUIP_FoX_get_value(attributes, "type1", value, status)
    if (status /= 0) call system_abort ("IPModel_Glue_read_params_xml cannot find type")
    read (value, *) ti ! This is the index of this particular species

    parse_curr_type_i = ti

    call QUIP_FoX_get_value(attributes, "type2", value, status)
    if (status /= 0) call system_abort ("IPModel_Glue_read_params_xml cannot find type")
    read (value, *) tj ! This is the index of this particular species

    parse_curr_type_j = tj

  elseif (parse_in_ip .and. parse_curr_type_i /= 0 .and. name == 'density') then
    ! Density is modelled as scale*exp(-extent*x). Could probably determine scale from normalization, but I won't.

    parse_in_density = .true.
    parse_curr_point = 1

    call QUIP_FoX_get_value(attributes, "num_points", value, status)
    if (status /= 0) call system_abort ("IPModel_Glue_read_params_xml density but no num_points")
    read (value, *) num_points ! scale of the electron density for this element

    call QUIP_FoX_get_value(attributes, "density_y1", value, status)
    if (status /= 0) call system_abort ("IPModel_Glue_read_params_xml density but no density_y1")
    read (value, *) parse_ip%spline_data_density(parse_curr_type_i)%y1 ! scale of the electron density for this element

    call QUIP_FoX_get_value(attributes, "density_yn", value, status)
    if (status /= 0) call system_abort ("IPModel_Glue_read_params_xml density but no density_yn")
    read (value, *) parse_ip%spline_data_density(parse_curr_type_i)%yn ! scale of the electron density for this element

    ! Allocate num_points in the current type
    allocate(parse_ip%spline_data_density(parse_curr_type_i)%data(num_points,2))
  
  elseif (parse_in_ip .and. parse_curr_type_i /= 0 .and. name == 'potential_density') then
    
    ! Potential supplied as density vs. energy
    parse_in_density_potential = .true.
    parse_curr_point = 1

    call QUIP_FoX_get_value(attributes, "num_points", value, status)
    if (status /= 0) call system_abort ("IPModel_Glue_read_params_xml potential_density but no num_points")
    read (value, *) num_points ! scale of the electron density for this element

    ! Allocate num_points in the current type
    allocate(parse_ip%spline_data_potential(parse_curr_type_i)%data(num_points,2))
  
  elseif (parse_in_ip .and. parse_curr_type_i /= 0 .and. name == 'potential_neighbours') then
    
    ! Potential supplied as distance to nearest neighbour vs. energy, along with number of neighbours
    parse_in_neighbours_potential = .true.
    parse_curr_point = 1

    call QUIP_FoX_get_value(attributes, "num_neighbours", value, status)
    if (status /= 0) call system_abort ("IPModel_Glue_read_params_xml missing num_neighbours on potential_neighbouring")
    read (value, *) parse_n_neighbours

    call QUIP_FoX_get_value(attributes, "num_points", value, status)
    if (status /= 0) call system_abort ("IPModel_Glue_read_params_xml potential_neighbours but no num_points")
    read (value, *) num_points ! scale of the electron density for this element

    ! Allocate num_points in the current type
    allocate(parse_ip%spline_data_potential(parse_curr_type_i)%data(num_points,2))
    
  elseif (parse_in_ip .and. parse_curr_type_i /= 0 .and. parse_curr_type_j /= 0 .and. name == 'potential_pair') then
    
    ! Pair potential
    parse_in_potential_pair = .true.
    parse_curr_point = 1

    call QUIP_FoX_get_value(attributes, "num_points", value, status)
    if (status /= 0) call system_abort ("IPModel_Glue_read_params_xml potential_pair but no num_points")
    read (value, *) num_points ! scale of the electron density for this element

    ! Allocate num_points in the current type
    allocate(parse_ip%spline_data_pair(parse_curr_type_i,parse_curr_type_j)%data(num_points,2))
  
  elseif (parse_in_ip .and. name == 'point') then
    
    if (parse_in_density_potential) then

      !
      !  Potential defined as a rho->E map
      !
      
      if (parse_curr_point > size(parse_ip%spline_data_density(parse_curr_type_i)%data,1)) call system_abort ("IPModel_Glue got too " // &
	"many points " // parse_curr_point // " type " // parse_curr_type_i // " in potential_density")

      call QUIP_FoX_get_value(attributes, "rho", value, status)
      if (status /= 0) call system_abort ("IPModel_Glue_read_params_xml cannot find rho")
      read (value, *) v
      parse_ip%spline_data_potential(parse_curr_type_i)%data(parse_curr_point,1) = v

      call QUIP_FoX_get_value(attributes, "E", value, status)
      if (status /= 0) call system_abort ("IPModel_Glue_read_params_xml cannot find E")
      read (value, *) v
      parse_ip%spline_data_potential(parse_curr_type_i)%data(parse_curr_point,2) = v

   elseif (parse_in_potential_pair) then

      if (parse_curr_point > size(parse_ip%spline_data_pair(parse_curr_type_i,parse_curr_type_j)%data,1)) call system_abort ("IPModel_Glue got too " // &
	"many points " // parse_curr_point // " type " // parse_curr_type_i // " in potential_pair")

      call QUIP_FoX_get_value(attributes, "r", value, status)
      if (status /= 0) call system_abort ("IPModel_Glue_read_params_xml cannot find r")
      read (value, *) v
      parse_ip%spline_data_pair(parse_curr_type_i,parse_curr_type_j)%data(parse_curr_point,1) = v

      call QUIP_FoX_get_value(attributes, "E", value, status)
      if (status /= 0) call system_abort ("IPModel_Glue_read_params_xml cannot find E")
      read (value, *) v
      parse_ip%spline_data_pair(parse_curr_type_i,parse_curr_type_j)%data(parse_curr_point,2) = v

    elseif (parse_in_neighbours_potential) then

      !
      !  Potential defined as a a->E map for N nearest neighbours
      !

      if (parse_curr_point > size(parse_ip%spline_data_potential(parse_curr_type_i)%data,1)) call system_abort ("IPModel_Glue got too " // &
"many points " // parse_curr_point // " type " // parse_curr_type_i // " in potential_neighbours")

      call QUIP_FoX_get_value(attributes, "a", value, status)
      if (status /= 0) call system_abort ("IPModel_Glue_read_params_xml cannot find a")
      read (value, *) v
      
      ! This implies that the denstiy must be defined above the potential
      parse_ip%spline_data_potential(parse_curr_type_i)%data(parse_curr_point,1) = eam_density(parse_ip, parse_curr_type_i, v) * parse_n_neighbours

      call QUIP_FoX_get_value(attributes, "E", value, status)
      if (status /= 0) call system_abort ("IPModel_Glue_read_params_xml cannot find E")
      read (value, *) v
      parse_ip%spline_data_potential(parse_curr_type_i)%data(parse_curr_point,2) = v

      ! call print(parse_ip%spline_data(parse_curr_type_i)%spline_potential(1, parse_curr_point) // " " // parse_ip%spline_data(parse_curr_type_i)%spline_potential(2, parse_curr_point) // " " )

    elseif (parse_in_density) then

      !
      !  Potential defined as a a->E map for N nearest neighbours
      !

      if (parse_curr_point > size(parse_ip%spline_data_density(parse_curr_type_i)%data,1)) call system_abort ("IPModel_Glue got too " // &
"many points " // parse_curr_point // " type " // parse_curr_type_i // " in density")

      call QUIP_FoX_get_value(attributes, "a", value, status)
      if (status /= 0) call system_abort ("IPModel_Glue_read_params_xml cannot find a")
      read (value, *) v
      
      ! This implies that the denstiy must be defined above the potential
      parse_ip%spline_data_density(parse_curr_type_i)%data(parse_curr_point,1) = v

      call QUIP_FoX_get_value(attributes, "rho", value, status)
      if (status /= 0) call system_abort ("IPModel_Glue_read_params_xml cannot find rho")
      read (value, *) v
      parse_ip%spline_data_density(parse_curr_type_i)%data(parse_curr_point,2) = v

      ! call print(parse_ip%spline_data(parse_curr_type_i)%spline_potential(1, parse_curr_point) // " " // parse_ip%spline_data(parse_curr_type_i)%spline_potential(2, parse_curr_point) // " " )

    endif

    parse_curr_point = parse_curr_point + 1

  endif

end subroutine IPModel_startElement_handler

subroutine IPModel_endElement_handler(URI, localname, name)
  character(len=*), intent(in)   :: URI
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name

  if (parse_in_ip) then
    if (name == 'Glue_params') then
      parse_in_ip = .false.
    elseif (name == 'potential_density') then
      parse_in_density_potential = .false.
    elseif (name == 'potential_pair') then
      parse_in_potential_pair = .false.
    elseif (name == 'potential_neighbours') then
      parse_in_neighbours_potential = .false.
    elseif (name == 'density') then
      parse_in_density = .false.
    elseif (name == 'per_type_data') then
      parse_curr_type_i = 0
    elseif (name == 'per_pair_data') then
      parse_curr_type_i = 0
      parse_curr_type_j = 0
    endif
  endif

end subroutine IPModel_endElement_handler

end module IPModel_Glue_module
