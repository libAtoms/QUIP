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
!X IPModel_GAP module  
!X
!% Module for Gaussian Approximation Potential.
!%
!% The IPModel_GAP object contains all the parameters read from a
!% 'GAP_params' XML stanza.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
module IPModel_GAP_module

use libatoms_module

use mpi_context_module
use QUIP_Common_module

!use gp_module
!use bispectrum_module

implicit none

private 

include 'IPModel_interface.h'

! this stuff is here for now, but it should live somewhere else eventually
! lower down in the GP
integer, parameter :: N_max = 11  
integer, parameter :: L_max = 3

public :: IPModel_GAP
type IPModel_GAP
  integer :: n_types = 0         !% Number of atomic types. 
  integer, allocatable :: atomic_num(:), type_of_atomic_num(:)  !% Atomic number dimensioned as \texttt{n_types}. 
  real(dp) :: cutoff = 0.0_dp    !% Cutoff for computing connection.

  character(len=256) datafile !% File name containing the GAP database

  character(len=FIELD_LENGTH) label
  type(mpi_context) :: mpi

  !type(gp) :: my_gp

  ! this stuff is here for now, but it should live somewhere else eventually
  ! lower down in the GP
  real(dp), dimension(N_max) :: gr0 
  real(dp), dimension(N_max) :: gsigma 


end type IPModel_GAP

logical :: parse_in_ip, parse_matched_label
type(IPModel_GAP), pointer :: parse_ip

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

subroutine IPModel_GAP_Initialise_str(this, args_str, param_str, mpi)
  type(IPModel_GAP), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str
  type(mpi_context), intent(in), optional :: mpi

  type(Dictionary) :: params

  call Finalise(this)

  call initialise(params)
  this%label = ''
  call param_register(params, 'label', '', this%label)
  if (.not. param_read_line(params, args_str, ignore_unknown=.true.)) then
    call system_abort("IPModel_GAP_Initialise_str failed to parse label from args_str="//trim(args_str))
  endif
  call finalise(params)

  call IPModel_GAP_read_params_xml(this, param_str)

  if (present(mpi)) this%mpi = mpi

  ! now initialise the potential
  this%gsigma = 5.0_dp
  !call uni_div(this%gr0,2.0_dp,3.0_dp)
  !call initialise(this%gr0,this%gsigma)

  !call gp_read(this%my_gp, this%datafile)

  this%cutoff = 3.0_dp

end subroutine IPModel_GAP_Initialise_str

subroutine IPModel_GAP_Finalise(this)
  type(IPModel_GAP), intent(inout) :: this

  if (allocated(this%atomic_num)) deallocate(this%atomic_num)
  if (allocated(this%type_of_atomic_num)) deallocate(this%type_of_atomic_num)

  this%n_types = 0
  this%label = ''
end subroutine IPModel_GAP_Finalise

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!% The potential calculator: this routine computes energy, forces and the virial.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

subroutine IPModel_GAP_Calc(this, at, e, local_e, f, virial)
  type(IPModel_GAP), intent(inout) :: this
  type(Atoms), intent(in) :: at
  real(dp), intent(out), optional :: e, local_e(:) !% \texttt{e} = System total energy, \texttt{local_e} = energy of each atom, vector dimensioned as \texttt{at%N}.  
  real(dp), intent(out), optional :: f(:,:)        !% Forces, dimensioned as \texttt{f(3,at%N)} 
  real(dp), intent(out), optional :: virial(3,3)   !% Virial

  real(dp), pointer :: w_e(:)
  real(dp) :: e_i, error_e
  real(dp), dimension(:,:), allocatable :: atoms_pos
  real(dp), dimension(:), allocatable   :: vec
  integer :: i, j, n

  !type(fourier_coefficient) :: f_hat
  !type(bispectrum_coefficient) :: bispectrum


  if (present(e)) e = 0.0_dp
  if (present(local_e)) local_e = 0.0_dp
  if (present(virial)) virial = 0.0_dp
  if (present(f)) then 
     if(size(f,1) .ne. 3 .or. size(f,2) .ne. at%N) call system_abort('IPModel_GAP_Calc: f is the wrong size')
     f = 0.0_dp
  end if

  ! no forces or virials for now
  if(present(f) .or. present(virial)) &
       call system_abort('IPModel_GAP_Cals: no forces or virials yet!')

  if (.not. assign_pointer(at, "weight", w_e)) nullify(w_e)

  ! temporary variable for bispectrum
  !allocate( vec(this%my_gp%d) )

  do i=1, at%N
    if (this%mpi%active) then
      if (mod(i-1, this%mpi%n_procs) /= this%mpi%my_proc) cycle
    endif

    ! pre-libatoms legacy: coordinates in array
    ! 1st: central atom
    ! other atoms relative to the central atom
    allocate( atoms_pos(3,atoms_n_neighbours(at,i)+1) )

    atoms_pos(:,1) = 0.0_dp
    do n = 1, atoms_n_neighbours(at,i)
       j = atoms_neighbour(at,i,n)
       atoms_pos(:,n+1) = diff_min_image(at,j,i)
    enddo
    
    ! f_hat: fourier transformed coordinates
    !call fourier_transform(atoms_pos,f_hat,N_max=N_max,L_max=L_max)
    ! bispectrum: bispectrum of f_hat
    !call calc_bispectrum(f_hat,bispectrum)
    ! transforming from bispectrum format to array
    !call bispectrum2vec( bispectrum,vec )
    
    ! predicting function value:
    ! estimated value: e_i
    ! estimated error: error_e
    ! variable: vec
    ! GP data: my_gp
    !call gp_predict(e_i, error_e, vec, this%my_gp)
    !e_i = gp_mean(vec, this%my_gp)

    if(present(local_e)) local_e(i) = e_i
    if (associated(w_e)) e_i = w_e(i)*e_i
    if(present(e)) e = e + e_i

    deallocate( atoms_pos )
  end do

  deallocate( vec )
  !call finalise(f_hat)
  !call finalise(bispectrum)

  if (present(e)) e = sum(this%mpi, e)
  if (present(local_e)) call sum_in_place(this%mpi, local_e)
  if (present(virial)) call sum_in_place(this%mpi, virial)
  if (present(f)) call sum_in_place(this%mpi, f)

end subroutine IPModel_GAP_Calc



!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!% XML param reader functions.
!% An example for XML stanza is given below, please notice that
!% they are simply dummy parameters for testing purposes, with no physical meaning.
!%
!%> <GAP_params datafile="file" n_types="1" label="default">
!%> <per_type_data type="1" atomic_num="14" />
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

  logical shifted
  integer ti, tj

  if (name == 'GAP_params') then ! new GAP stanza
     
    if (parse_in_ip) &
          call system_abort("IPModel_startElement_handler entered GAP_params with parse_in true. Probably a bug in FoX (4.0.1, e.g.)")
     
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
	call system_abort("Can't find n_types in GAP_params")
      endif

      allocate(parse_ip%atomic_num(parse_ip%n_types))
      parse_ip%atomic_num = 0

      parse_ip%datafile = ''
      call QUIP_FoX_get_value(attributes, 'datafile', value, status)
      if (status == 0) then
         parse_ip%datafile = trim(value)
      else
	call system_abort("Can't find datafile in GAP_params")
      endif

    endif

  elseif (parse_in_ip .and. name == 'per_type_data') then

    call QUIP_FoX_get_value(attributes, "type", value, status)
    if (status /= 0) call system_abort ("IPModel_GAP_read_params_xml cannot find type")
    read (value, *) ti

    call QUIP_FoX_get_value(attributes, "atomic_num", value, status)
    if (status /= 0) call system_abort ("IPModel_GAP_read_params_xml cannot find atomic_num")
    read (value, *) parse_ip%atomic_num(ti)

    ! build reverse lookup table to atomic types
    if (allocated(parse_ip%type_of_atomic_num)) deallocate(parse_ip%type_of_atomic_num)
    allocate(parse_ip%type_of_atomic_num(maxval(parse_ip%atomic_num)))
    parse_ip%type_of_atomic_num = 0
    do ti=1, parse_ip%n_types
      if (parse_ip%atomic_num(ti) > 0) &
	parse_ip%type_of_atomic_num(parse_ip%atomic_num(ti)) = ti
    end do

  endif

end subroutine IPModel_startElement_handler

subroutine IPModel_endElement_handler(URI, localname, name)
  character(len=*), intent(in)   :: URI  
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name 

  if (parse_in_ip) then
    if (name == 'GAP_params') then
      parse_in_ip = .false.
    end if
  endif

end subroutine IPModel_endElement_handler

subroutine IPModel_GAP_read_params_xml(this, param_str)
  type(IPModel_GAP), intent(inout), target :: this
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
    call system_abort("IPModel_GAP_read_params_xml parsed file, but n_types = 0")
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

  integer :: ti, tj

  call Print("IPModel_GAP : Gaussian Approximation Potential", file=file)
  call Print("IPModel_GAP : datafile = "//trim(this%datafile)//" n_types = " // this%n_types // " cutoff = " // this%cutoff, file=file)

  do ti=1, this%n_types
    call Print ("IPModel_GAP : type " // ti // " atomic_num " // this%atomic_num(ti), file=file)
  end do

end subroutine IPModel_GAP_Print

end module IPModel_GAP_module
