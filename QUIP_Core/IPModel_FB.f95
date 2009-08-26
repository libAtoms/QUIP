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
!X IPModel_FB module  
!X
!% 
!% 
!% 
!% 
!% 
!% 
!%
!% 
!% 
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

module IPModel_FB_module

use libatoms_module
use IPEwald_module

use mpi_context_module
use QUIP_Common_module

implicit none

private 

include 'IPModel_interface.h'

public :: IPModel_FB
type IPModel_FB
  integer :: n_types = 0         !% Number of atomic types. 
  integer, allocatable :: atomic_num(:), type_of_atomic_num(:)  !% Atomic number dimensioned as \texttt{n_types}. 

  real(dp) :: cutoff = 0.0_dp    !% Cutoff for computing connection.

  real(dp), dimension(:), allocatable :: charge
  real(dp), dimension(:,:), allocatable :: A, B, C, r_cut

  character(len=FIELD_LENGTH) label
  type(mpi_context) :: mpi

end type IPModel_FB

logical :: parse_in_ip, parse_matched_label
integer :: parse_cur_type_i, parse_cur_type_j
type(IPModel_FB), pointer :: parse_ip

interface Initialise
  module procedure IPModel_FB_Initialise_str
end interface Initialise

interface Finalise
  module procedure IPModel_FB_Finalise
end interface Finalise

interface Print
  module procedure IPModel_FB_Print
end interface Print

interface Calc
  module procedure IPModel_FB_Calc
end interface Calc

contains

subroutine IPModel_FB_Initialise_str(this, args_str, param_str, mpi)
  type(IPModel_FB), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str
  type(mpi_context), intent(in), optional :: mpi

  type(Dictionary) :: params

  call Finalise(this)

  call initialise(params)
  this%label = ''
  call param_register(params, 'label', '', this%label)
  if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_FB_Initialise_str args_str')) then
    call system_abort("IPModel_FB_Initialise_str failed to parse label from args_str="//trim(args_str))
  endif
  call finalise(params)

  call IPModel_FB_read_params_xml(this, param_str)

  this%cutoff = maxval(this%r_cut)

  if (present(mpi)) this%mpi = mpi

end subroutine IPModel_FB_Initialise_str

subroutine IPModel_FB_Finalise(this)
  type(IPModel_FB), intent(inout) :: this

  if (allocated(this%atomic_num)) deallocate(this%atomic_num)
  if (allocated(this%type_of_atomic_num)) deallocate(this%type_of_atomic_num)

  if (allocated(this%A)) deallocate(this%A)
  if (allocated(this%B)) deallocate(this%B)
  if (allocated(this%C)) deallocate(this%C)
  if (allocated(this%charge)) deallocate(this%charge)
  if (allocated(this%r_cut)) deallocate(this%r_cut)

  this%n_types = 0
  this%label = ''
end subroutine IPModel_FB_Finalise

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!% The potential calculator: this routine computes energy, forces and the virial.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

subroutine IPModel_FB_Calc(this, at, e, local_e, f, virial)
  type(IPModel_FB), intent(in) :: this
  type(Atoms), intent(in) :: at
  real(dp), intent(out), optional :: e !% \texttt{e} = System total energy
  real(dp), dimension(:), intent(out), optional :: local_e !% \texttt{local_e} = energy of each atom, vector dimensioned as \texttt{at%N}.  
  real(dp), dimension(:,:), intent(out), optional :: f        !% Forces, dimensioned as \texttt{f(3,at%N)} 
  real(dp), dimension(3,3), intent(out), optional :: virial   !% Virial

  integer :: i, j, n, ti, tj
  real(dp) :: r_ij, de, de_dr
  real(dp), dimension(3) :: u_ij
  real(dp), dimension(:), pointer :: charge

  type(Atoms) :: at_ew

  if (present(e)) e = 0.0_dp
  if (present(local_e)) local_e = 0.0_dp
  if (present(virial)) virial = 0.0_dp
  if (present(f)) f = 0.0_dp 

  at_ew = at
  call add_property(at_ew, 'charge', 0.0_dp)
  if( .not. assign_pointer(at_ew, 'charge', charge) ) call system_abort('')
  charge = 0.0_dp

  do i = 1, at_ew%N
    ti = get_type(this%type_of_atomic_num, at_ew%Z(i))
    charge(i) = this%charge(ti)
  enddo
  call ewald_calc(at_ew,e=e,f=f,virial=virial)
  call finalise(at_ew)

  do i = 1, at%N
    if (this%mpi%active) then
      if (mod(i-1, this%mpi%n_procs) /= this%mpi%my_proc) cycle
    endif
    ti = get_type(this%type_of_atomic_num, at%Z(i))
    do n = 1, atoms_n_neighbours(at, i)
      j = atoms_neighbour(at, i, n, distance=r_ij, cosines = u_ij)
      if (i < j) cycle
      tj = get_type(this%type_of_atomic_num, at%Z(j))
      if (present(e) .or. present(local_e)) then
	de = this%A(ti,tj) * exp( - r_ij / this%B(ti,tj) ) - this%C(ti,tj) / r_ij**6
	if (present(local_e)) then
	  local_e(i) = local_e(i) + 0.5_dp*de
	  local_e(j) = local_e(j) + 0.5_dp*de
	endif
	if (present(e)) e = e + de
      endif
      if (present(f) .or. present(virial)) then
	de_dr = -this%A(ti,tj) * exp( - r_ij / this%B(ti,tj) ) / this%B(ti,tj) + &
              & 6.0_dp * this%C(ti,tj) / r_ij**7
	if (present(f)) then
	  f(:,i) = f(:,i) + de_dr*u_ij
	  f(:,j) = f(:,j) - de_dr*u_ij
	endif
	if (present(virial)) virial = virial - de_dr*(u_ij .outer. u_ij)*r_ij
      endif
    end do
  end do

  if (present(e)) e = sum(this%mpi, e)
  if (present(local_e)) call sum_in_place(this%mpi, local_e)
  if (present(virial)) call sum_in_place(this%mpi, virial)
  if (present(f)) call sum_in_place(this%mpi, f)

end subroutine IPModel_FB_Calc

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!% XML param reader functions.
!% An example for XML stanza is given below, please notice that
!% they are simply dummy parameters for testing purposes, with no physical meaning.
!%
!%> <LJ_params n_types="2" label="default">
!%> <per_type_data type="1" atomic_num="29" />
!%> <per_type_data type="2" atomic_num="79" />
!%> <per_pair_data type1="1" type2="1" sigma="4.0" eps6="1.0" 
!%>       eps12="1.0" cutoff="6.0" shifted="T" />
!%> <per_pair_data type1="2" type2="2" sigma="5.0" eps6="2.0" 
!%>       eps12="2.0" cutoff="7.5" shifted="T" />
!%> <per_pair_data type1="1" type2="2" sigma="4.5" eps6="1.5" 
!%>       eps12="1.5" cutoff="6.75" shifted="T" />
!%> </LJ_params>
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
subroutine IPModel_startElement_handler(URI, localname, name, attributes)
  character(len=*), intent(in)   :: URI  
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name 
  type(dictionary_t), intent(in) :: attributes

  integer :: status
  character(len=1024) :: value

  integer :: ti, Zi, Zj

  if (name == 'FB_params') then ! new FB stanza

    if (parse_in_ip) &
      call system_abort("IPModel_startElement_handler entered FB_params with parse_in true. Probably a bug in FoX (4.0.1, e.g.)")

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
	call system_abort("Can't find n_types in FB_params")
      endif

      allocate(parse_ip%atomic_num(parse_ip%n_types))
      parse_ip%atomic_num = 0
      allocate(parse_ip%A(parse_ip%n_types,parse_ip%n_types))
      parse_ip%A = 0.0_dp
      allocate(parse_ip%B(parse_ip%n_types,parse_ip%n_types))
      parse_ip%B = 0.0_dp
      allocate(parse_ip%C(parse_ip%n_types,parse_ip%n_types))
      parse_ip%C = 0.0_dp
      allocate(parse_ip%charge(parse_ip%n_types))
      parse_ip%charge = 0.0_dp
      allocate(parse_ip%r_cut(parse_ip%n_types,parse_ip%n_types))
      parse_ip%r_cut = 0.0_dp
    endif

  elseif (parse_in_ip .and. name == 'per_type_data') then

    call QUIP_FoX_get_value(attributes, "type", value, status)
    if (status /= 0) call system_abort ("IPModel_FB_read_params_xml cannot find type")
    read (value, *) ti

    call QUIP_FoX_get_value(attributes, "atomic_num", value, status)
    if (status /= 0) call system_abort ("IPModel_FB_read_params_xml cannot find atomic_num")
    read (value, *) parse_ip%atomic_num(ti)

    call QUIP_FoX_get_value(attributes, "charge", value, status)
    if (status /= 0) call system_abort ("IPModel_FB_read_params_xml cannot find charge")
    read (value, *) parse_ip%charge(ti)

    if (allocated(parse_ip%type_of_atomic_num)) deallocate(parse_ip%type_of_atomic_num)
    allocate(parse_ip%type_of_atomic_num(maxval(parse_ip%atomic_num)))
    parse_ip%type_of_atomic_num = 0
    do ti=1, parse_ip%n_types
      if (parse_ip%atomic_num(ti) > 0) &
	parse_ip%type_of_atomic_num(parse_ip%atomic_num(ti)) = ti
    end do

  elseif (parse_in_ip .and. name == 'per_pair_data') then

    call QUIP_FoX_get_value(attributes, "atomic_num_i", value, status)
    if (status /= 0) call system_abort ("IPModel_FB_read_params_xml cannot find atomic_num_i")
    read (value, *) Zi

    call QUIP_FoX_get_value(attributes, "atomic_num_j", value, status)
    if (status /= 0) call system_abort ("IPModel_FB_read_params_xml cannot find atomic_num_i")
    read (value, *) Zj

    parse_cur_type_i = get_type(parse_ip%type_of_atomic_num,Zi)
    parse_cur_type_j = get_type(parse_ip%type_of_atomic_num,Zj)

    call QUIP_FoX_get_value(attributes, "A", value, status)
    if (status /= 0) call system_abort ("IPModel_FB_read_params_xml cannot find A")
    read (value, *) parse_ip%A(parse_cur_type_i,parse_cur_type_j)

    call QUIP_FoX_get_value(attributes, "B", value, status)
    if (status /= 0) call system_abort ("IPModel_FB_read_params_xml cannot find B")
    read (value, *) parse_ip%B(parse_cur_type_i,parse_cur_type_j)

    call QUIP_FoX_get_value(attributes, "C", value, status)
    if (status /= 0) call system_abort ("IPModel_FB_read_params_xml cannot find C")
    read (value, *) parse_ip%C(parse_cur_type_i,parse_cur_type_j)

    call QUIP_FoX_get_value(attributes, "r_cut", value, status)
    if (status /= 0) call system_abort ("IPModel_FB_read_params_xml cannot find r_cut")
    read (value, *) parse_ip%r_cut(parse_cur_type_i,parse_cur_type_j)

    if (parse_cur_type_i /= parse_cur_type_j) then
        parse_ip%A(parse_cur_type_j,parse_cur_type_i) = parse_ip%A(parse_cur_type_i,parse_cur_type_j)
        parse_ip%B(parse_cur_type_j,parse_cur_type_i) = parse_ip%B(parse_cur_type_i,parse_cur_type_j)
        parse_ip%C(parse_cur_type_j,parse_cur_type_i) = parse_ip%C(parse_cur_type_i,parse_cur_type_j)
        parse_ip%r_cut(parse_cur_type_j,parse_cur_type_i) = parse_ip%r_cut(parse_cur_type_i,parse_cur_type_j)
    end if

  endif

end subroutine IPModel_startElement_handler

subroutine IPModel_endElement_handler(URI, localname, name)
  character(len=*), intent(in)   :: URI  
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name 

  if (parse_in_ip) then
    if (name == 'FB_params') then
      parse_in_ip = .false.
    end if
  endif

end subroutine IPModel_endElement_handler

subroutine IPModel_FB_read_params_xml(this, param_str)
  type(IPModel_FB), intent(inout), target :: this
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
    call system_abort("IPModel_FB_read_params_xml parsed file, but n_types = 0")
  endif

end subroutine IPModel_FB_read_params_xml


!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!% Printing of FB parameters: number of different types, cutoff radius, atomic numbers, etc.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

subroutine IPModel_FB_Print (this, file)
  type(IPModel_FB), intent(in) :: this
  type(Inoutput), intent(inout),optional :: file

  integer :: ti, tj

  call Print("IPModel_FB : Flikkema Bromley", file=file)
  call Print("IPModel_FB : n_types = " // this%n_types // " cutoff = " // this%cutoff, file=file)

  do ti=1, this%n_types
    call Print ("IPModel_FB : type " // ti // " atomic_num " // this%atomic_num(ti) // " charge " // this%charge(ti), file=file)
    call verbosity_push_decrement()
    do tj=1, this%n_types
      call Print ("IPModel_FB : interaction " // ti // " " // tj // " A " // this%A(ti,tj) // " B " // &
      & this%B(ti,tj) // " C " // this%C(ti,tj) // &
      & " r_cut " // this%r_cut(ti,tj), file=file)
    end do
    call verbosity_pop()
  end do

end subroutine IPModel_FB_Print

end module IPModel_FB_module

!<FB_params n_types="2" label="default">
!<!-- Flikkema & Bromley  -->
!<per_type_data type="1" atomic_num="8" charge="-1.2"/>
!<per_type_data type="2" atomic_num="14" charge="2.4"/>
!<per_pair_data atomic_num_i="8" atomic_num_j="8" A="1428.406" B="0.358" C="41.374" r_cut="6.0" />
!<per_pair_data atomic_num_i="8" atomic_num_j="14" A="10454.202" B="0.208" C="63.047" r_cut="6.0" />
!<per_pair_data atomic_num_i="14" atomic_num_j="14" A="79502.113" B="0.201" C="446.780" r_cut="6.0" />
!</FB_params>
