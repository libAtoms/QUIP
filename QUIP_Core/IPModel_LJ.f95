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
!X IPModel_LJ module  
!X
!% Module for Lennard-Jones pair potential.
!% \begin{equation} 
!%   \nonumber
!%     V(r) = 4 \epsilon \left[  \left( \frac{\sigma}{r} \right)^{12} - \left( \frac{\sigma}{r} \right)^6 \right] 
!% \end{equation} 
!% For parameters see Ashcroft and Mermin, {\it Solid State Physics}. 
!%
!% The IPModel_LJ object contains all the parameters read from a
!% 'LJ_params' XML stanza.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
module IPModel_LJ_module

use libatoms_module

use mpi_context_module
use QUIP_Common_module

implicit none

private 

include 'IPModel_interface.h'

public :: IPModel_LJ
type IPModel_LJ
  integer :: n_types = 0         !% Number of atomic types. 
  integer, allocatable :: atomic_num(:), type_of_atomic_num(:)  !% Atomic number dimensioned as \texttt{n_types}. 

  real(dp) :: cutoff = 0.0_dp    !% Cutoff for computing connection.

  real(dp), allocatable :: sigma(:,:), eps6(:,:), eps12(:,:), cutoff_a(:,:), shift(:,:) !% IP parameters.

  character(len=FIELD_LENGTH) label
  type(mpi_context) :: mpi

end type IPModel_LJ

logical :: parse_in_ip, parse_matched_label
type(IPModel_LJ), pointer :: parse_ip

interface Initialise
  module procedure IPModel_LJ_Initialise_str
end interface Initialise

interface Finalise
  module procedure IPModel_LJ_Finalise
end interface Finalise

interface Print
  module procedure IPModel_LJ_Print
end interface Print

interface Calc
  module procedure IPModel_LJ_Calc
end interface Calc

contains

subroutine IPModel_LJ_Initialise_str(this, args_str, param_str, mpi)
  type(IPModel_LJ), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str
  type(mpi_context), intent(in), optional :: mpi

  type(Dictionary) :: params

  call Finalise(this)

  call initialise(params)
  this%label = ''
  call param_register(params, 'label', '', this%label)
  if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_LJ_Initialise_str args_str')) then
    call system_abort("IPModel_LJ_Initialise_str failed to parse label from args_str="//trim(args_str))
  endif
  call finalise(params)

  call IPModel_LJ_read_params_xml(this, param_str)

  this%cutoff = maxval(this%cutoff_a)

  if (present(mpi)) this%mpi = mpi

end subroutine IPModel_LJ_Initialise_str

subroutine IPModel_LJ_Finalise(this)
  type(IPModel_LJ), intent(inout) :: this

  if (allocated(this%atomic_num)) deallocate(this%atomic_num)
  if (allocated(this%type_of_atomic_num)) deallocate(this%type_of_atomic_num)

  if (allocated(this%sigma)) deallocate(this%sigma)
  if (allocated(this%eps6)) deallocate(this%eps6)
  if (allocated(this%eps12)) deallocate(this%eps12)
  if (allocated(this%cutoff_a)) deallocate(this%cutoff_a)
  if (allocated(this%shift)) deallocate(this%shift)

  this%n_types = 0
  this%label = ''
end subroutine IPModel_LJ_Finalise

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!% The potential calculator: this routine computes energy, forces and the virial.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

subroutine IPModel_LJ_Calc(this, at, e, local_e, f, virial, args_str)
  type(IPModel_LJ), intent(inout) :: this
  type(Atoms), intent(inout) :: at
  real(dp), intent(out), optional :: e, local_e(:) !% \texttt{e} = System total energy, \texttt{local_e} = energy of each atom, vector dimensioned as \texttt{at%N}.  
  real(dp), intent(out), optional :: f(:,:)        !% Forces, dimensioned as \texttt{f(3,at%N)} 
  real(dp), intent(out), optional :: virial(3,3)   !% Virial
  character(len=*), intent(in), optional      :: args_str

  real(dp), pointer :: w_e(:)
  integer i, ji, j, ti, tj
  real(dp) :: dr(3), dr_mag
  real(dp) :: de, de_dr
  logical :: i_is_min_image

  integer :: i_calc, n_extra_calcs
  character(len=20) :: extra_calcs_list(10)

  logical :: do_flux = .false.
  real(dp), pointer :: velo(:,:)
  real(dp) :: flux(3)

  if (present(e)) e = 0.0_dp
  if (present(local_e)) local_e = 0.0_dp
  if (present(virial)) virial = 0.0_dp
  if (present(f)) then 
     if(size(f,1) .ne. 3 .or. size(f,2) .ne. at%N) call system_abort('IPMOdel_LJ_Calc: f is the wrong size')
     f = 0.0_dp
  end if

  if (present(args_str)) then
    if (len_trim(args_str) > 0) then
      n_extra_calcs = parse_extra_calcs(args_str, extra_calcs_list)
      if (n_extra_calcs > 0) then
	do i_calc=1, n_extra_calcs
	  select case(trim(extra_calcs_list(i_calc)))
	    case("flux")
	      if (.not. assign_pointer(at, "velo", velo)) &
		call system_abort("IPModel_LJ_Calc Flux calculation requires velo field")
	      do_flux = .true.
	      flux = 0.0_dp
	    case default
	      call system_abort("Unsupported extra_calc '"//trim(extra_calcs_list(i_calc))//"'")
	  end select
	end do
      endif ! n_extra_calcs
    endif ! len_trim(args_str)
  endif ! present(args_str)

  if (.not. assign_pointer(at, "weight", w_e)) nullify(w_e)

  do i = 1, at%N
    i_is_min_image = is_min_image(at,i)

    if (this%mpi%active) then
      if (mod(i-1, this%mpi%n_procs) /= this%mpi%my_proc) cycle
    endif

    do ji = 1, atoms_n_neighbours(at, i)
      j = atoms_neighbour(at, i, ji, dr_mag, cosines = dr)
      
      if (dr_mag .feq. 0.0_dp) cycle
      if ((i < j) .and. i_is_min_image) cycle

      ti = get_type(this%type_of_atomic_num, at%Z(i))
      tj = get_type(this%type_of_atomic_num, at%Z(j))

      if (present(e) .or. present(local_e)) then
	de = IPModel_LJ_pairenergy(this, ti, tj, dr_mag)
	if (present(local_e)) then
	  local_e(i) = local_e(i) + 0.5_dp*de
          if(i_is_min_image) local_e(j) = local_e(j) + 0.5_dp*de
	endif
	if (present(e)) then
	  if (associated(w_e)) then
	    de = de*0.5_dp*(w_e(i)+w_e(j))
	  endif
          if(i_is_min_image) then
             e = e + de
          else
             e = e + 0.5_dp*de
          endif
	endif
      endif
      if (present(f) .or. present(virial) .or. do_flux) then
	de_dr = IPModel_LJ_pairenergy_deriv(this, ti, tj, dr_mag)
	if (associated(w_e)) then
	  de_dr = de_dr*0.5_dp*(w_e(i)+w_e(j))
	endif
	if (present(f)) then
	  f(:,i) = f(:,i) + de_dr*dr
	  if(i_is_min_image) f(:,j) = f(:,j) - de_dr*dr
	endif
	if (do_flux) then
	  ! -0.5 (v_i + v_j) . F_ij * dr_ij
	  flux = flux - 0.5_dp*sum((velo(:,i)+velo(:,j))*(de_dr*dr))*(dr*dr_mag)
	endif
	if (present(virial)) then
	  if(i_is_min_image) then
             virial = virial - de_dr*(dr .outer. dr)*dr_mag
          else
             virial = virial - 0.5_dp*de_dr*(dr .outer. dr)*dr_mag
          endif
	endif
      endif
    end do
  end do

  if (present(e)) e = sum(this%mpi, e)
  if (present(local_e)) call sum_in_place(this%mpi, local_e)
  if (present(virial)) call sum_in_place(this%mpi, virial)
  if (present(f)) call sum_in_place(this%mpi, f)
  if (do_flux) then
    flux = flux / cell_volume(at)
    call sum_in_place(this%mpi, flux)
    call set_value(at%params, "Flux", flux)
  endif

end subroutine IPModel_LJ_Calc

!% This routine computes the two-body term for a pair of atoms  separated by a distance r.
function IPModel_LJ_pairenergy(this, ti, tj, r)
  type(IPModel_LJ), intent(in) :: this
  integer, intent(in) :: ti, tj   !% Atomic types.
  real(dp), intent(in) :: r       !% Distance.
  real(dp) :: IPModel_LJ_pairenergy

  real(dp) :: tpow

  if ((r .feq. 0.0_dp) .or. (r > this%cutoff_a(ti,tj))) then
    IPModel_LJ_pairenergy = 0.0
    return
  endif

  tpow = (this%sigma(ti,tj)/r)**6

  IPModel_LJ_pairenergy = (this%eps12(ti,tj)*tpow*tpow - this%eps6(ti,tj)*tpow) + this%shift(ti,tj)
end function IPModel_LJ_pairenergy

!% Derivative of the two-body term.
function IPModel_LJ_pairenergy_deriv(this, ti, tj, r)
  type(IPModel_LJ), intent(in) :: this
  integer, intent(in) :: ti, tj   !% Atomic types.
  real(dp), intent(in) :: r       !% Distance.
  real(dp) :: IPModel_LJ_pairenergy_deriv

  real(dp) :: tpow

  if ((r .feq. 0.0_dp) .or. (r > this%cutoff_a(ti,tj))) then
    IPModel_LJ_pairenergy_deriv = 0.0
    return
  endif

  tpow = (this%sigma(ti,tj)/r)**6

  IPModel_LJ_pairenergy_deriv = (-12.0_dp*this%eps12(ti,tj)*tpow*tpow + 6.0_dp*this%eps6(ti,tj)*tpow)/r
end function IPModel_LJ_pairenergy_deriv

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

  logical shifted
  integer ti, tj

  if (name == 'LJ_params') then ! new LJ stanza

    if (parse_in_ip) &
      call system_abort("IPModel_startElement_handler entered LJ_params with parse_in true. Probably a bug in FoX (4.0.1, e.g.)")

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
	call system_abort("Can't find n_types in LJ_params")
      endif

      allocate(parse_ip%atomic_num(parse_ip%n_types))
      parse_ip%atomic_num = 0
      allocate(parse_ip%sigma(parse_ip%n_types,parse_ip%n_types))
      parse_ip%sigma = 1.0_dp
      allocate(parse_ip%eps6(parse_ip%n_types,parse_ip%n_types))
      parse_ip%eps6 = 0.0_dp
      allocate(parse_ip%eps12(parse_ip%n_types,parse_ip%n_types))
      parse_ip%eps12 = 0.0_dp
      allocate(parse_ip%cutoff_a(parse_ip%n_types,parse_ip%n_types))
      parse_ip%cutoff_a = 0.0_dp
      allocate(parse_ip%shift(parse_ip%n_types,parse_ip%n_types))
      parse_ip%shift = 0.0_dp
    endif

  elseif (parse_in_ip .and. name == 'per_type_data') then

    call QUIP_FoX_get_value(attributes, "type", value, status)
    if (status /= 0) call system_abort ("IPModel_LJ_read_params_xml cannot find type")
    read (value, *) ti

    call QUIP_FoX_get_value(attributes, "atomic_num", value, status)
    if (status /= 0) call system_abort ("IPModel_LJ_read_params_xml cannot find atomic_num")
    read (value, *) parse_ip%atomic_num(ti)

    if (allocated(parse_ip%type_of_atomic_num)) deallocate(parse_ip%type_of_atomic_num)
    allocate(parse_ip%type_of_atomic_num(maxval(parse_ip%atomic_num)))
    parse_ip%type_of_atomic_num = 0
    do ti=1, parse_ip%n_types
      if (parse_ip%atomic_num(ti) > 0) &
	parse_ip%type_of_atomic_num(parse_ip%atomic_num(ti)) = ti
    end do

    parse_ip%shift = 0.0_dp

  elseif (parse_in_ip .and. name == 'per_pair_data') then

    call QUIP_FoX_get_value(attributes, "type1", value, status)
    if (status /= 0) call system_abort ("IPModel_LJ_read_params_xml cannot find type1")
    read (value, *) ti
    call QUIP_FoX_get_value(attributes, "type2", value, status)
    if (status /= 0) call system_abort ("IPModel_LJ_read_params_xml cannot find type2")
    read (value, *) tj

    call QUIP_FoX_get_value(attributes, "sigma", value, status)
    if (status /= 0) call system_abort ("IPModel_LJ_read_params_xml cannot find sigma")
    read (value, *) parse_ip%sigma(ti,tj)
    call QUIP_FoX_get_value(attributes, "eps6", value, status)
    if (status /= 0) call system_abort ("IPModel_LJ_read_params_xml cannot find eps6")
    read (value, *) parse_ip%eps6(ti,tj)
    call QUIP_FoX_get_value(attributes, "eps12", value, status)
    if (status /= 0) call system_abort ("IPModel_LJ_read_params_xml cannot find eps12")
    read (value, *) parse_ip%eps12(ti,tj)
    call QUIP_FoX_get_value(attributes, "cutoff", value, status)
    if (status /= 0) call system_abort ("IPModel_LJ_read_params_xml cannot find cutoff")
    read (value, *) parse_ip%cutoff_a(ti,tj)
    call QUIP_FoX_get_value(attributes, "shifted", value, status)
    if (status /= 0) call system_abort ("IPModel_LJ_read_params_xml cannot find shifted")
    read (value, *) shifted

    if (shifted) then
      parse_ip%shift(ti,tj) = - IPModel_LJ_pairenergy(parse_ip, ti, tj, parse_ip%cutoff_a(ti,tj))
    endif

    if (ti /= tj) then
      parse_ip%eps6(tj,ti) = parse_ip%eps6(ti,tj)
      parse_ip%eps12(tj,ti) = parse_ip%eps12(ti,tj)
      parse_ip%sigma(tj,ti) = parse_ip%sigma(ti,tj)
      parse_ip%cutoff_a(tj,ti) = parse_ip%cutoff_a(ti,tj)
      parse_ip%shift(tj,ti) = parse_ip%shift(ti,tj)
    end if

  endif

end subroutine IPModel_startElement_handler

subroutine IPModel_endElement_handler(URI, localname, name)
  character(len=*), intent(in)   :: URI  
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name 

  if (parse_in_ip) then
    if (name == 'LJ_params') then
      parse_in_ip = .false.
    end if
  endif

end subroutine IPModel_endElement_handler

subroutine IPModel_LJ_read_params_xml(this, param_str)
  type(IPModel_LJ), intent(inout), target :: this
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
    call system_abort("IPModel_LJ_read_params_xml parsed file, but n_types = 0")
  endif

end subroutine IPModel_LJ_read_params_xml


!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!% Printing of LJ parameters: number of different types, cutoff radius, atomic numbers, etc.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

subroutine IPModel_LJ_Print (this, file)
  type(IPModel_LJ), intent(in) :: this
  type(Inoutput), intent(inout),optional :: file

  integer :: ti, tj

  call Print("IPModel_LJ : Lennard-Jones", file=file)
  call Print("IPModel_LJ : n_types = " // this%n_types // " cutoff = " // this%cutoff, file=file)

  do ti=1, this%n_types
    call Print ("IPModel_LJ : type " // ti // " atomic_num " // this%atomic_num(ti), file=file)
    call verbosity_push_decrement()
    do tj=1, this%n_types
      call Print ("IPModel_LJ : interaction " // ti // " " // tj // " sigma " // this%sigma(ti,tj) // " eps6,12 " // &
	this%eps6(ti,tj) // " " // this%eps12(ti,tj) // " cutoff_a " // this%cutoff_a(ti,tj) // " shift " // &
	this%shift(ti,tj), file=file)
    end do
    call verbosity_pop()
  end do

end subroutine IPModel_LJ_Print

function parse_extra_calcs(args_str, extra_calcs_list) result(n_extra_calcs)
  character(len=*), intent(in) :: args_str
  character(len=*), intent(out) :: extra_calcs_list(:)
  integer :: n_extra_calcs

  character(len=FIELD_LENGTH) :: extra_calcs_str
  type(Dictionary) :: params

  n_extra_calcs = 0
  call initialise(params)
  call param_register(params, "extra_calcs", "", extra_calcs_str)
  if (param_read_line(params, args_str, ignore_unknown=.true.,task='parse_extra_calcs')) then
    if (len_trim(extra_calcs_str) > 0) then
      call split_string_simple(extra_calcs_str, extra_calcs_list, n_extra_calcs, ":")
    end if
  end if
  call finalise(params)

end function parse_extra_calcs

end module IPModel_LJ_module
