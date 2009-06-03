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
!X IPModel_ASAP
!X
!% Interface to ASAP potential.
!% P. Tangney and S. Scandolo,
!% An ab initio parametrized interatomic force field for silica
!% J. Chem. Phys, 117, 8898 (2002). 
!%
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
module IPModel_ASAP_module

use libatoms_module

use mpi_context_module
use QUIP_Common_module

implicit none
private

include 'IPModel_interface.h'

logical, private :: asap_initialised = .false.

integer, private, parameter :: PARAM_LINE_LENGTH = 255, N_PARAM_LINE = 94

public :: IPModel_ASAP
type IPModel_ASAP
  integer :: n_types = 0, n_atoms = 0
  integer, allocatable :: atomic_num(:), type_of_atomic_num(:)

  real(dp) :: cutoff = 0.0_dp

  character(len=FIELD_LENGTH) :: label
  character(len=PARAM_LINE_LENGTH) param_lines(N_PARAM_LINE)
  integer :: n_param_lines = 0, param_col = 1
  type(mpi_context) :: mpi
  logical :: initialised

end type IPModel_ASAP

logical :: parse_in_ip, parse_matched_label, parse_in_ip_params
type(IPModel_ASAP), pointer :: parse_ip

interface Initialise
  module procedure IPModel_ASAP_Initialise_str
end interface Initialise

interface Finalise
  module procedure IPModel_ASAP_Finalise
end interface Finalise

interface Print
  module procedure IPModel_ASAP_Print
end interface Print

interface Calc
  module procedure IPModel_ASAP_Calc
end interface Calc

contains

subroutine IPModel_ASAP_Initialise_str(this, args_str, param_str, mpi)
  type(IPModel_ASAP), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str
  type(mpi_context), intent(in), optional :: mpi

  type(Dictionary) :: params
  character(len=FIELD_LENGTH) label

  this%initialised = .false.

  call Finalise(this)

  call initialise(params)
  this%label=''
  call param_register(params, 'label', '', this%label)
  if (.not. param_read_line(params, args_str, ignore_unknown=.true.)) then
    call system_abort("IPModel_ASAP_Initialise_str failed to parse label from args_str="//trim(args_str))
  endif
  call finalise(params)

  call IPModel_ASAP_read_params_xml(this, param_str)
  this%n_atoms = 0
  this%initialised = .true.
  
  if (present(mpi)) this%mpi = mpi

end subroutine IPModel_ASAP_Initialise_str

subroutine IPModel_ASAP_Finalise(this)
  type(IPModel_ASAP), intent(inout) :: this

#ifdef HAVE_ASAP
  if (asap_initialised) then
     call asap_singlepoint_finalise()
     asap_initialised = .false.
  end if
  this%initialised = .false.
#else
  call system_abort('ASAP potential is not compiled in. Recompile with HAVE_ASAP=1')
#endif

  if (allocated(this%atomic_num)) deallocate(this%atomic_num)
  if (allocated(this%type_of_atomic_num)) deallocate(this%type_of_atomic_num)

  this%n_types = 0
  this%label = ''
end subroutine IPModel_ASAP_Finalise


subroutine IPModel_ASAP_Calc(this, at, e, local_e, f, virial, args_str)
   type(IPModel_ASAP), intent(inout):: this
   type(Atoms), intent(inout)      :: at
   real(dp), intent(out), optional :: e, local_e(:)
   real(dp), intent(out), optional :: f(:,:)
   real(dp), intent(out), optional :: virial(3,3)
   character(len=*), optional      :: args_str

   type(Dictionary) :: params
   integer, allocatable :: spind(:)
   real(dp) :: asap_e, asap_stress(3,3)
   real(dp), allocatable :: asap_f(:,:)
   integer :: i
   logical :: restart

   allocate(asap_f(3,at%N))
   allocate(spind(at%N))

#ifdef HAVE_ASAP

   call initialise(params)
   this%label=''
   call param_register(params, 'restart', 'F', restart)
   if (.not. param_read_line(params, args_str, ignore_unknown=.true.)) then
      call system_abort("IPModel_ASAP_Initialise_str failed to parse args_str="//trim(args_str))
   endif
   call finalise(params)

   if (.not. asap_initialised .or. this%n_atoms /= at%n) then
      
      if (asap_initialised) then
         call asap_singlepoint_finalise()
         asap_initialised = .false.
      end if

      this%n_atoms = at%n
      call asap_singlepoint_init(this%n_atoms, this%n_types, this%param_lines)
      asap_initialised = .true.
   end if

   spind = 0
   do i=1,this%n_types
      where(at%Z == this%atomic_num(i)) spind = i
   end do

   ! ASAP uses atomic units - lengths are in Bohr, energies in Hartree,
   ! forces in Hartree/Bohr and stress in Hartree/Bohr**3
   call asap_singlepoint_calc(at%pos/BOHR, spind, at%lattice/BOHR, &
        asap_e, asap_f, asap_stress, restart)

   ! Convert to {eV, A, fs} units
   if (present(e)) e = asap_e*HARTREE
   if (present(f)) f = asap_f*(HARTREE/BOHR)
   if (present(virial)) virial = asap_stress*(HARTREE/(BOHR**3))
   
   deallocate(asap_f)
   deallocate(spind)

#else
  call system_abort('ASAP potential is not compiled in. Recompile with HAVE_ASAP=1')
#endif

end subroutine IPModel_ASAP_Calc


subroutine IPModel_ASAP_Print(this, file)
  type(IPModel_ASAP), intent(in) :: this
  type(Inoutput), intent(inout),optional :: file

  integer :: ti, tj, i

  call Print("IPModel_ASAP : ASAP Potential", file=file)
  call Print("IPModel_ASAP : n_types = " // this%n_types //" n_atoms = "//this%n_atoms// " cutoff = " // this%cutoff, file=file)
  call Print("IPModel_ASAP : n_param_lines = " // this%n_param_lines)

  do i=1,this%n_param_lines
     call print(this%param_lines(i))
  end do

  do ti=1, this%n_types
    call Print ("IPModel_ASAP : type " // ti // " atomic_num " // this%atomic_num(ti), file=file)
    call verbosity_push_decrement()
    call Print ("IPModel_ASAP : " // &
        "cutoff " // this%cutoff , &
         file=file)
    call verbosity_pop()
  end do

end subroutine IPModel_ASAP_Print

subroutine IPModel_ASAP_read_params_xml(this, param_str)
  type(IPModel_ASAP), intent(inout), target :: this
  character(len=*), intent(in) :: param_str

  type(xml_t) :: fxml
  integer i

  if (len(trim(param_str)) <= 0) return

  parse_in_ip = .false. 
  parse_matched_label = .false.
  parse_ip => this

  this%n_param_lines = 1
  this%param_col = 1
  do i=1,N_PARAM_LINE
     this%param_lines(i) = repeat(' ',PARAM_LINE_LENGTH)
  end do

  call open_xml_string(fxml, param_str)
  call parse(fxml,  &
    startElement_handler = IPModel_startElement_handler, &
    endElement_handler = IPModel_endElement_handler, &
    characters_handler = IPModel_character_handler)
  call close_xml_t(fxml)

  if (this%n_types == 0) then
    call system_abort("IPModel_ASAP_read_params_xml parsed file, but n_types = 0")
  endif

end subroutine IPModel_ASAP_read_params_xml

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
  integer ti, tj, i, num_fields

  if (name == 'ASAP_params') then ! new ASAP stanza

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
        call system_abort("Can't find n_types in ASAP_params")
      endif

      allocate(parse_ip%atomic_num(parse_ip%n_types))
      parse_ip%atomic_num = 0

      call QUIP_FoX_get_value(attributes, "cutoff", value, status)
      if (status /= 0) call system_abort ("IPModel_ASAP_read_params_xml cannot find cutoff")
      read (value, *) parse_ip%cutoff
    endif

  elseif (parse_in_ip .and. name == 'per_type_data') then

    call QUIP_FoX_get_value(attributes, "type", value, status)
    if (status /= 0) call system_abort ("IPModel_ASAP_read_params_xml cannot find type")
    read (value, *) ti

    call QUIP_FoX_get_value(attributes, "atomic_num", value, status)
    if (status /= 0) call system_abort ("IPModel_ASAP_read_params_xml cannot find atomic_num")
    read (value, *) parse_ip%atomic_num(ti)

    if (allocated(parse_ip%type_of_atomic_num)) deallocate(parse_ip%type_of_atomic_num)
    allocate(parse_ip%type_of_atomic_num(maxval(parse_ip%atomic_num)))
    parse_ip%type_of_atomic_num = 0
    do ti=1, parse_ip%n_types
      if (parse_ip%atomic_num(ti) > 0) &
        parse_ip%type_of_atomic_num(parse_ip%atomic_num(ti)) = ti
    end do

  else if (parse_in_ip .and. name == 'params') then

     call Print('entering <params>')
     parse_in_ip_params = .true.

  endif

end subroutine IPModel_startElement_handler

subroutine IPModel_character_handler(chars)
  character(len=*), intent(in) :: chars

  integer :: i, l

  if (parse_in_ip_params) then
     do i=1,len(chars)
        if (chars(i:i) == char(10)) then
           if (i /= 1 .and. i /= len(chars)) then
              parse_ip%n_param_lines = parse_ip%n_param_lines + 1
              if (parse_ip%n_param_lines > N_PARAM_LINE) then
                 call system_abort('IPModel_ASAP: too many param lines '//parse_ip%n_param_lines)
              end if
              parse_ip%param_col = 1
           end if
        else
           parse_ip%param_lines(parse_ip%n_param_lines)(parse_ip%param_col:parse_ip%param_col) = chars(i:i)
           parse_ip%param_col = parse_ip%param_col + 1
           if (parse_ip%param_col > PARAM_LINE_LENGTH) call system_abort('IPModel_ASAPL: param line too long')
        end if
     end do
  end if

end subroutine IPModel_character_handler

subroutine IPModel_endElement_handler(URI, localname, name)
  character(len=*), intent(in)   :: URI
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name

  if (parse_in_ip) then
    if (name == 'ASAP_params') then
      parse_in_ip = .false.
    end if
    
    if (name == 'params') then
       call Print('leaving <params>')
       parse_in_ip_params = .false.
    end if

  endif

end subroutine IPModel_endElement_handler

end module IPModel_ASAP_module
