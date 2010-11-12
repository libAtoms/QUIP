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

use libatoms_module

use mpi_context_module
use QUIP_Common_module

implicit none
private

include 'IPModel_interface.h'

public :: SplineDataContainer
type SplineDataContainer
    real(dp), allocatable :: spline_potential(:,:)
    real(dp), allocatable :: spline_density(:,:)
    real(dp) :: density_y1, density_yn
end type SplineDataContainer

public :: IPModel_Glue
type IPModel_Glue
  integer :: n_types = 0 ! Number of different species we have

  ! Maps from atomic number to type and vice versa
  integer, allocatable :: atomic_num(:), type_of_atomic_num(:)

  real(dp) :: cutoff = 0.0_dp
  real(dp), allocatable, dimension(:) :: poly

  real(dp), allocatable :: density_extent(:), density_scale(:)

  type(SplineDataContainer), allocatable :: spline_data(:)
  type(Spline), allocatable :: potential(:)
  type(Spline), allocatable :: density(:)

  character(len=FIELD_LENGTH) :: label

end type IPModel_Glue

logical, private :: parse_in_ip, parse_matched_label, parse_in_density_potential, parse_in_neighbours_potential, parse_in_density
type(IPModel_Glue), private, pointer :: parse_ip
integer :: parse_curr_type = 0, parse_curr_point, parse_n_neighbours

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
  character(len=FIELD_LENGTH) label

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
  integer :: i

  ! Add finalisation code here

  if (allocated(this%atomic_num)) deallocate(this%atomic_num)
  if (allocated(this%type_of_atomic_num)) deallocate(this%type_of_atomic_num)
  
  if (allocated(this%density_extent)) deallocate(this%density_extent)
  if (allocated(this%density_scale)) deallocate(this%density_scale)

  if (allocated(this%spline_data)) then
    do i=1,this%n_types
      if(allocated(this%spline_data(i)%spline_potential)) deallocate(this%spline_data(i)%spline_potential)
      if(allocated(this%spline_data(i)%spline_density)) deallocate(this%spline_data(i)%spline_density)
    end do
    deallocate(this%spline_data)
    deallocate(this%poly)
  endif

  if (allocated(this%potential)) then
    do i=1,this%n_types
      call finalise(this%potential(i))
      call finalise(this%density(i))
    end do
    deallocate(this%potential)
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
  logical :: has_atom_mask_name
  character(FIELD_LENGTH) :: atom_mask_name
  ! Loop variables
  integer :: i, ji, j, ti, tj  

  ! If only Fortran allow you to declare variables somewhere else
  real(dp) :: r_ij_mag, r_ij_hat(3) ! Neighbour vector info
  real(dp) :: rho_local ! Accumulator for local electron density

  ! For calculating forces and the virial tensor
  real(dp) :: drho_i_drij
  real(dp) :: drho_i_dri(3), potential_deriv
  real(dp), dimension(3,3) :: drho_i_drij_outer_rij, virial_i
    
  INIT_ERROR(error)

  if (present(e)) e = 0.0_dp
  if (present(local_e)) then
     call check_size('Local_E',local_e,(/at%N/),'IPModel_Glue_Calc', error)
     local_e = 0.0_dp
  endif
  if (present(f)) then
     call check_size('Force',f,(/3,at%Nbuffer/),'IPModel_Glue_Calc', error)
     f = 0.0_dp
  end if
  if (present(virial)) virial = 0.0_dp
  if (present(local_virial)) then
     call check_size('Local_virial',local_virial,(/9,at%Nbuffer/),'IPModel_Glue_Calc', error)
     local_virial = 0.0_dp
  endif

  atom_mask_pointer => null()
  if(present(args_str)) then
     call initialise(params)
     call param_register(params, 'atom_mask_name', 'NONE',atom_mask_name,has_value_target=has_atom_mask_name, help_string="No help yet.  This source file was $LastChangedBy$")
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
  endif

  ! Iterate over atoms
  do i = 1, at%N
     if(associated(atom_mask_pointer)) then
        if(.not. atom_mask_pointer(i)) cycle
     endif

    rho_local = 0.0_dp ! Local density from neighbours
    drho_i_dri = 0.0_dp
    drho_i_drij = 0.0_dp
    drho_i_drij_outer_rij = 0.0_dp
    
    ! Get the type of this species from its atomic number
    ti = get_type(this%type_of_atomic_num, at%Z(i))

    ! Iterate over our nighbours
    do ji=1, atoms_n_neighbours(at, i)
      j = atoms_neighbour(at, i, ji, r_ij_mag, cosines=r_ij_hat)
      tj = get_type(this%type_of_atomic_num, at%Z(j))
      if (r_ij_mag < glue_cutoff(this, tj)) then ! Skip atoms beyond the cutoff
          rho_local = rho_local + eam_density(this, tj, r_ij_mag)
          drho_i_drij = eam_density_deriv(this, tj, r_ij_mag)
          drho_i_dri = drho_i_dri + drho_i_drij * r_ij_hat
          drho_i_drij_outer_rij = drho_i_drij_outer_rij + drho_i_drij*(r_ij_hat .outer. r_ij_hat) * r_ij_mag
      endif
    end do ! ji
    potential_deriv = eam_spline_potential_deriv(this, ti, rho_local)
    if(present(f)) f(:,i) = f(:,i) + drho_i_dri * potential_deriv

    ! Iterate over neighbours again to sum forces
    do ji=1, atoms_n_neighbours(at, i)
      j = atoms_neighbour(at, i, ji, r_ij_mag, cosines=r_ij_hat)
      tj = get_type(this%type_of_atomic_num, at%Z(j))
     
      if(present(f) .and. r_ij_mag < glue_cutoff(this, tj)) then         
          f(:,j) = f(:,j) - potential_deriv * eam_density_deriv(this, tj, r_ij_mag) * r_ij_hat
      endif
    end do ! ji
    if(present(local_e)) local_e(i) = eam_spline_potential(this, ti, rho_local)
    if(present(e)) e = e + eam_spline_potential(this, ti, rho_local)
    if(present(virial) .or. present(local_virial)) virial_i = potential_deriv * drho_i_drij_outer_rij
    if(present(virial))  virial = virial - virial_i
    if(present(local_virial)) local_virial(:,i) = local_virial(:,i) - reshape(virial_i,(/9/))
  end do ! i
end subroutine IPModel_Glue_Calc

function glue_cutoff(this, ti)
  type(IPModel_Glue), intent(in) :: this
  integer, intent(in) :: ti
  real(dp) :: glue_cutoff
  
  !glue_cutoff = this%density_extent(ti)*10
  glue_cutoff = max_knot(this%density(ti))

end function glue_cutoff

function eam_density(this, ti, r)
  type(IPModel_Glue), intent(in) :: this
  integer, intent(in) :: ti
  real(dp), intent(in) :: r
  real(dp) :: eam_density
  
  if (r < glue_cutoff(this, ti)) then
    !eam_density = this%density_scale(ti) * exp(-r/this%density_extent(ti))
    !eam_density = spline_value(this%density(ti),r) 
    eam_density = this%poly(ti)*(glue_cutoff(this, ti)-r)**3
  else
    eam_density = 0.0_dp
  endif
end function eam_density

function eam_density_deriv(this, ti, r)
  type(IPModel_Glue), intent(in) :: this
  integer, intent(in) :: ti
  real(dp), intent(in) :: r
  real(dp) :: eam_density_deriv

  if (r < glue_cutoff(this, ti)) then
    !eam_density_deriv = - this%density_scale(ti)/this%density_extent(ti) * exp(-r/this%density_extent(ti))
    !eam_density_deriv = spline_deriv(this%density(ti),r)
    eam_density_deriv = -3.0_dp * this%poly(ti)*(glue_cutoff(this, ti)-r)**2
  else
    eam_density_deriv = 0.0_dp
  endif
end function eam_density_deriv

function eam_spline_potential(this, ti, rho)
  ! Evaluate the potential spline at r for the ti^th species
  type(IPModel_Glue), intent(in) :: this
  integer, intent(in) :: ti
  integer :: last
  real(dp), intent(in) :: rho
  real(dp) :: eam_spline_potential

  real(dp) :: min_rho
  real(dp) :: max_rho

  !min_rho = min_knot(this%potential(ti))
  !max_rho = max_knot(this%potential(ti))

  !if (rho < min_rho) then
  !  !eam_spline_potential = this%potential(ti)%y(1) + (rho - this%potential(ti)%x(1)) * spline_deriv(this%potential(ti), min_rho)
  !  eam_spline_potential = this%potential(ti)%y(1) + (rho - this%potential(ti)%x(1)) * this%potential(ti)%yp1
  !  eam_spline_potential = this%potential(ti)%y(1) + (rho - this%potential(ti)%x(1)) * this%potential(ti)%yp1
  !elseif (rho >= max_rho) then
  !  last = size(this%potential(ti)%x)
  !  !eam_spline_potential = this%potential(ti)%y(last) + (rho - this%potential(ti)%x(last)) * spline_deriv(this%potential(ti), max_rho)
  !  eam_spline_potential = this%potential(ti)%y(last) + (rho - this%potential(ti)%x(last)) * this%potential(ti)%ypn
  !else
  !  eam_spline_potential = spline_value(this%potential(ti), rho)
  !endif

  eam_spline_potential = spline_value(this%potential(ti), rho)

end function eam_spline_potential

function eam_spline_potential_deriv(this, ti, rho)
  ! Evaluate the potential spline derivitive at r for the ti^th species
  type(IPModel_Glue), intent(in) :: this
  integer, intent(in) :: ti
  real(dp), intent(in) :: rho
  real(dp) :: eam_spline_potential_deriv

  real(dp) :: min_rho
  real(dp) :: max_rho

  !min_rho = min_knot(this%potential(ti))
  !max_rho = max_knot(this%potential(ti))

  !if (rho < min_rho) then
  !  !eam_spline_potential_deriv = spline_deriv(this%potential(ti), min_rho)
  !  eam_spline_potential_deriv = this%potential(ti)%yp1
  !elseif (rho >= max_rho) then
  !  !eam_spline_potential_deriv = spline_deriv(this%potential(ti), max_rho)
  !  eam_spline_potential_deriv = this%potential(ti)%ypn
  !else
  !  eam_spline_potential_deriv = spline_deriv(this%potential(ti), rho)
  !endif
  eam_spline_potential_deriv = spline_deriv(this%potential(ti), rho)

end function eam_spline_potential_deriv


subroutine IPModel_Glue_Print(this, file)
  type(IPModel_Glue), intent(in) :: this
  type(Inoutput), intent(inout),optional :: file

  integer :: ti, tj

  call Print("IPModel_Glue : Template Potential", file=file)
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

subroutine sort_spline(this, ti)
  type(IPModel_Glue), intent(inout), target :: this
  integer :: m, ti, num_potential_points
  real(dp) :: temp1, temp2
  logical :: finished

  ! Knuth forgive me, for I have written bubble sort. Reimplement this as quicksort or something that isn't dreadful
  num_potential_points = size(this%spline_data(ti)%spline_potential(1,:))
  finished = .false.

  do while (finished .neqv. .true.)
      finished = .true.
      do m=2, num_potential_points
          if (this%spline_data(ti)%spline_potential(1,m-1) > this%spline_data(ti)%spline_potential(1,m)) then
              temp1 = this%spline_data(ti)%spline_potential(1,m-1)
              temp2 = this%spline_data(ti)%spline_potential(2,m-1)
              this%spline_data(ti)%spline_potential(1,m-1) = this%spline_data(ti)%spline_potential(1,m)
              this%spline_data(ti)%spline_potential(2,m-1) = this%spline_data(ti)%spline_potential(2,m)
              this%spline_data(ti)%spline_potential(1,m) = temp1
              this%spline_data(ti)%spline_potential(2,m) = temp2
              finished = .false.
          end if
      end do
  end do
end subroutine sort_spline

subroutine IPModel_Glue_read_params_xml(this, param_str)
  type(IPModel_Glue), intent(inout), target :: this
  character(len=*), intent(in) :: param_str
  integer :: ti, num_potential_points, m
  real(dp) :: grad_upper, grad_lower
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

  do ti=1, this%n_types
    call sort_spline(this, ti)
    !do m=1, size(this%spline_data(ti)%spline_potential(1,:))
    !  call print (ti // " " // this%spline_data(ti)%spline_potential(1,m) // " " // this%spline_data(ti)%spline_potential(2,m))
    !end do
  end do

  ! Initialise the spline data structures with the spline data
  do ti=1, this%n_types
    num_potential_points = size(this%spline_data(ti)%spline_potential(1,:))
  
    grad_upper = (this%spline_data(ti)%spline_potential(2,num_potential_points) - this%spline_data(ti)%spline_potential(2,num_potential_points-1)) / &
                 (this%spline_data(ti)%spline_potential(1,num_potential_points) - this%spline_data(ti)%spline_potential(1,num_potential_points-1))

    grad_lower = (this%spline_data(ti)%spline_potential(2,2) - this%spline_data(ti)%spline_potential(2,1)) / &
                 (this%spline_data(ti)%spline_potential(1,2) - this%spline_data(ti)%spline_potential(1,1))

    !call initialise(this%potential(ti), this%spline_data(ti)%spline_potential(1,:), this%spline_data(ti)%spline_potential(2,:), grad_lower, grad_upper)
    call initialise(this%potential(ti), this%spline_data(ti)%spline_potential(1,:), this%spline_data(ti)%spline_potential(2,:),2.0e30_dp,2.0e30_dp)
    call initialise(this%density(ti), this%spline_data(ti)%spline_density(1,:), this%spline_data(ti)%spline_density(2,:), this%spline_data(ti)%density_y1, this%spline_data(ti)%density_yn)

    this%poly(ti) = dot_product(this%spline_data(ti)%spline_density(2,:), (glue_cutoff(this,ti) - &
      this%spline_data(ti)%spline_density(1,:))**3) / &
      dot_product((glue_cutoff(this,ti)-this%spline_data(ti)%spline_density(1,:))**3, &
      (glue_cutoff(this,ti)-this%spline_data(ti)%spline_density(1,:))**3)
  end do

  ! Find the largest cutoff
  do ti=1, this%n_types
    if(glue_cutoff(this,ti) > this%cutoff) then
      this%cutoff = glue_cutoff(this,ti)
    endif
  end do

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
  character(len=FIELD_LENGTH) :: value

  logical shifted
  integer ti, tj
  real(dp) :: v
  integer num_points

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
        call system_abort("Can't find n_types in Template_params")
      endif

      allocate(parse_ip%atomic_num(parse_ip%n_types))
      parse_ip%atomic_num = 0
 
      allocate(parse_ip%density_extent(parse_ip%n_types))
      parse_ip%density_extent = 0

      allocate(parse_ip%density_scale(parse_ip%n_types))
      parse_ip%density_scale = 0

      ! Allocate n_types spline data objects for further allocation
      allocate(parse_ip%spline_data(parse_ip%n_types))

      ! Allocate the splines
      allocate(parse_ip%potential(parse_ip%n_types))
      allocate(parse_ip%density(parse_ip%n_types))
      allocate(parse_ip%poly(parse_ip%n_types))
    endif


  elseif (parse_in_ip .and. name == 'per_type_data') then

    call QUIP_FoX_get_value(attributes, "type", value, status)
    if (status /= 0) call system_abort ("IPModel_Glue_read_params_xml cannot find type")
    read (value, *) ti ! This is the index of this particular species

    parse_curr_type = ti

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

  elseif (parse_in_ip .and. parse_curr_type /= 0 .and. name == 'density') then
    ! Density is modelled as scale*exp(-extent*x). Could probably determine scale from normalization, but I won't.

    !call QUIP_FoX_get_value(attributes, "a", value, status)
    !if (status /= 0) call system_abort ("IPModel_Glue_read_params_xml density but no a")
    !read (value, *) parse_ip%density_extent(parse_curr_type) ! characteristic extent of the electron density for this element in angstrom

    !call QUIP_FoX_get_value(attributes, "scale", value, status)
    !if (status /= 0) call system_abort ("IPModel_Glue_read_params_xml density but no scale")
    !read (value, *) parse_ip%density_scale(parse_curr_type) ! scale of the electron density for this element

    parse_in_density = .true.
    parse_curr_point = 1

    call QUIP_FoX_get_value(attributes, "num_points", value, status)
    if (status /= 0) call system_abort ("IPModel_Glue_read_params_xml density but no num_points")
    read (value, *) num_points ! scale of the electron density for this element

    call QUIP_FoX_get_value(attributes, "density_y1", value, status)
    if (status /= 0) call system_abort ("IPModel_Glue_read_params_xml density but no density_y1")
    read (value, *) parse_ip%spline_data(parse_curr_type)%density_y1 ! scale of the electron density for this element

    call QUIP_FoX_get_value(attributes, "density_yn", value, status)
    if (status /= 0) call system_abort ("IPModel_Glue_read_params_xml density but no density_yn")
    read (value, *) parse_ip%spline_data(parse_curr_type)%density_yn ! scale of the electron density for this element

    ! Allocate num_points in the current type
    allocate(parse_ip%spline_data(parse_curr_type)%spline_density(2, num_points))
  
  elseif (parse_in_ip .and. parse_curr_type /= 0 .and. name == 'potential_density') then
    
    ! Potential supplied as density vs. energy
    parse_in_density_potential = .true.
    parse_curr_point = 1

    call QUIP_FoX_get_value(attributes, "num_points", value, status)
    if (status /= 0) call system_abort ("IPModel_Glue_read_params_xml potential_density but no num_points")
    read (value, *) num_points ! scale of the electron density for this element

    ! Allocate num_points in the current type
    allocate(parse_ip%spline_data(parse_curr_type)%spline_potential(2, num_points))
  
  elseif (parse_in_ip .and. parse_curr_type /= 0 .and. name == 'potential_neighbours') then
    
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
    allocate(parse_ip%spline_data(parse_curr_type)%spline_potential(2, num_points))
    
  elseif (parse_in_ip .and. name == 'point') then
    
    if (parse_in_density_potential) then

      !
      !  Potential defined as a rho->E map
      !
      
      if (parse_curr_point > size(parse_ip%spline_data(parse_curr_type)%spline_potential(1,:))) call system_abort ("IPModel_Glue got too " // &
	"many points " // parse_curr_point // " type " // parse_curr_type // " in potential_density")

      call QUIP_FoX_get_value(attributes, "rho", value, status)
      if (status /= 0) call system_abort ("IPModel_Glue_read_params_xml cannot find rho")
      read (value, *) v
      parse_ip%spline_data(parse_curr_type)%spline_potential(1, parse_curr_point) = v

      call QUIP_FoX_get_value(attributes, "E", value, status)
      if (status /= 0) call system_abort ("IPModel_Glue_read_params_xml cannot find E")
      read (value, *) v
      parse_ip%spline_data(parse_curr_type)%spline_potential(2, parse_curr_point) = v

    elseif (parse_in_neighbours_potential) then

      !
      !  Potential defined as a a->E map for N nearest neighbours
      !

      if (parse_curr_point > size(parse_ip%spline_data(parse_curr_type)%spline_potential(1,:))) call system_abort ("IPModel_Glue got too " // &
"many points " // parse_curr_point // " type " // parse_curr_type // " in potential_neighbours")

      call QUIP_FoX_get_value(attributes, "a", value, status)
      if (status /= 0) call system_abort ("IPModel_Glue_read_params_xml cannot find a")
      read (value, *) v
      
      ! This implies that the denstiy must be defined above the potential
      parse_ip%spline_data(parse_curr_type)%spline_potential(1, parse_curr_point) = eam_density(parse_ip, parse_curr_type, v) * parse_n_neighbours

      call QUIP_FoX_get_value(attributes, "E", value, status)
      if (status /= 0) call system_abort ("IPModel_Glue_read_params_xml cannot find E")
      read (value, *) v
      parse_ip%spline_data(parse_curr_type)%spline_potential(2, parse_curr_point) = v

      ! call print(parse_ip%spline_data(parse_curr_type)%spline_potential(1, parse_curr_point) // " " // parse_ip%spline_data(parse_curr_type)%spline_potential(2, parse_curr_point) // " " )

    elseif (parse_in_density) then

      !
      !  Potential defined as a a->E map for N nearest neighbours
      !

      if (parse_curr_point > size(parse_ip%spline_data(parse_curr_type)%spline_density(1,:))) call system_abort ("IPModel_Glue got too " // &
"many points " // parse_curr_point // " type " // parse_curr_type // " in density")

      call QUIP_FoX_get_value(attributes, "a", value, status)
      if (status /= 0) call system_abort ("IPModel_Glue_read_params_xml cannot find a")
      read (value, *) v
      
      ! This implies that the denstiy must be defined above the potential
      parse_ip%spline_data(parse_curr_type)%spline_density(1, parse_curr_point) = v

      call QUIP_FoX_get_value(attributes, "rho", value, status)
      if (status /= 0) call system_abort ("IPModel_Glue_read_params_xml cannot find rho")
      read (value, *) v
      parse_ip%spline_data(parse_curr_type)%spline_density(2, parse_curr_point) = v

      ! call print(parse_ip%spline_data(parse_curr_type)%spline_potential(1, parse_curr_point) // " " // parse_ip%spline_data(parse_curr_type)%spline_potential(2, parse_curr_point) // " " )

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
    elseif (name == 'potential_neighbours') then
      parse_in_neighbours_potential = .false.
    elseif (name == 'density') then
      parse_in_density = .false.
    elseif (name == 'per_type_data') then
      parse_curr_type = 0
    endif
  endif

end subroutine IPModel_endElement_handler

end module IPModel_Glue_module
