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
!X IPModel_LMTO_TBE
!X
!% Interface to empirical Tight Binding code within LM suite
!% (https://bitbucket.org/lmto/lm/branch/libtbe)
!%
!% Requires linking with LMTO code (HAVE_LMTO_TBE=1 in Makefile.inc).
!%
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#include "error.inc"

module IPModel_LMTO_TBE_module

use error_module
use units_module, only: RYDBERG, BOHR
use system_module, only : dp, inoutput, print, verbosity_push_decrement, verbosity_pop, operator(//)
use dictionary_module
use paramreader_module
use linearalgebra_module
use atoms_types_module
use atoms_module
#ifdef HAVE_LMTO_TBE
use libtbe, only: tbinit, tbsc
#endif

use mpi_context_module
use QUIP_Common_module

implicit none
private

include 'IPModel_interface.h'

public :: IPModel_LMTO_TBE
type IPModel_LMTO_TBE
  integer :: n_types = 0
  integer, allocatable :: atomic_num(:), type_of_atomic_num(:)
  real(dp) :: cutoff = 0.0_dp !% cutoff is always zero as LMTO code does its own neighbour calculation
  character(len=STRING_LENGTH) :: label, control_file
  type(Atoms) :: oldat

end type IPModel_LMTO_TBE

logical, private :: parse_in_ip, parse_matched_label
type(IPModel_LMTO_TBE), private, pointer :: parse_ip

interface Initialise
  module procedure IPModel_LMTO_TBE_Initialise_str
end interface Initialise

interface Finalise
  module procedure IPModel_LMTO_TBE_Finalise
end interface Finalise

interface Print
  module procedure IPModel_LMTO_TBE_Print
end interface Print

interface Calc
  module procedure IPModel_LMTO_TBE_Calc
end interface Calc

contains

subroutine IPModel_LMTO_TBE_Initialise_str(this, args_str, param_str)
  type(IPModel_LMTO_TBE), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str

  type(Dictionary) :: params

  call Finalise(this)

  call initialise(params)
  this%label=''
  call param_register(params, 'label', '', this%label, help_string="No help yet.  This source file was $LastChangedBy$")
  if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_LMTO_TBE_Initialise_str args_str')) then
    call system_abort("IPModel_LMTO_TBE_Initialise_str failed to parse label from args_str="//trim(args_str))
  endif
  call finalise(params)

  call IPModel_LMTO_TBE_read_params_xml(this, param_str)

#ifdef HAVE_LMTO_TBE
  call tbinit(this%control_file)
#else
  call system_abort("IPModel_LMTO_TBE can only be used when compiled with HAVE_LMTO_TBE enabled")
#endif


end subroutine IPModel_LMTO_TBE_Initialise_str

subroutine IPModel_LMTO_TBE_Finalise(this)
  type(IPModel_LMTO_TBE), intent(inout) :: this

  if (allocated(this%atomic_num)) deallocate(this%atomic_num)
  if (allocated(this%type_of_atomic_num)) deallocate(this%type_of_atomic_num)

  this%n_types = 0
  this%label = ''
end subroutine IPModel_LMTO_TBE_Finalise


subroutine IPModel_LMTO_TBE_Calc(this, at, e, local_e, f, virial, local_virial, args_str, mpi, error)
   type(IPModel_LMTO_TBE), intent(inout):: this
   type(Atoms), intent(inout)      :: at
   real(dp), intent(out), optional :: e, local_e(:)
   real(dp), intent(out), optional :: f(:,:), local_virial(:,:)   !% Forces, dimensioned as \texttt{f(3,at%N)}, local virials, dimensioned as \texttt{local_virial(9,at%N)} 
   real(dp), intent(out), optional :: virial(3,3)
   character(len=*), optional      :: args_str
   type(MPI_Context), intent(in), optional :: mpi
   integer, intent(out), optional :: error

   integer i
   logical :: reorder = .false.
   integer, allocatable :: atypes(:)
   real(dp), allocatable :: force(:,:), atpos(:,:)
   real(dp) :: energy, atlattice(3,3)
   integer, pointer :: index(:), oldindex(:)

   INIT_ERROR(error)

#ifdef HAVE_LMTO_TBE
   if (present(e)) e = 0.0_dp
   if (present(local_e)) then
      call check_size('Local_E',local_e,(/at%N/),'IPModel_LMTO_TBE_Calc', error)
      local_e = 0.0_dp
   endif
   if (present(f)) then
      call check_size('Force',f,(/3,at%Nbuffer/),'IPModel_LMTO_TBE_Calc', error)
      f = 0.0_dp
   end if
   if (present(virial)) virial = 0.0_dp
   if (present(local_virial)) then
      call check_size('Local_virial',local_virial,(/9,at%Nbuffer/),'IPModel_LMTO_TBE_Calc', error)
      local_virial = 0.0_dp
   endif

   ! check if Atoms is "equivalent" (isoatomic) to previous Atoms
   ! if not we need to re-initialise the TB code
   if (at%n /= this%oldat%n) then
      call print('LMTO_TBE_Calc: natoms changed, reordering')
      reorder = .true.
   else if (maxval(abs(at%lattice - this%oldat%lattice)) > 1e-7_dp) then
      call print('LMTO_TBE_Calc: lattice components changed, reordering')
      reorder = .true.
   else if (any(at%z /= this%oldat%z)) then
      call print('LMTO_TBE_Calc: atomic numbers changed, reordering')
      reorder = .true.
   else if (has_property(at, 'index')) then
      ! clusters have an 'index' property that points back to original atom
      if (has_property(this%oldat, 'index')) then
         call assign_property_pointer(at, 'index', index, error)
         PASS_ERROR(error)
         call assign_property_pointer(this%oldat, 'index', oldindex, error)
         PASS_ERROR(error)
         if (any(index /= oldindex)) then
            call print('LMTO_TBE_Calc: atomic numbers changed, reordering')
            reorder = .true.
         end if
      else
         call print('LMTO_TBE_Calc: new atoms has "index" property, old does not, reordering')
         reorder = .true.
      end if
   end if

   allocate(atypes(at%n), force(3, at%n), atpos(3, at%n))
   do i=1, at%N
      atypes(i) = get_type(this%type_of_atomic_num, at%Z(i))
   end do

   call print('LMTO_TBE_Calc: invoking tbsc() on '//at%n//' atoms with reorder='//reorder)

   atlattice = at%lattice/BOHR
   atpos = at%pos/BOHR
   ! libtbe expects the force to be initialised as it will accumulate forces
   force = 0.0_dp

   call tbsc(atlattice, at%n, atpos, atypes, energy, force, &
        reorder=reorder, incomm=mpi%communicator)
   if (present(e)) then
      e = energy*RYDBERG
   end if
   if (present(f)) then
      f = force*(RYDBERG/BOHR)
   end if

   deallocate(atypes, force, atpos)

   call atoms_copy_without_connect(this%oldat, at)
#else
   RAISE_ERROR("IPModel_LMTO_TBE_Calc() called but LMTO_TBE support not compiled in", error)
#endif

end subroutine IPModel_LMTO_TBE_Calc


subroutine IPModel_LMTO_TBE_Print(this, file)
  type(IPModel_LMTO_TBE), intent(in) :: this
  type(Inoutput), intent(inout),optional :: file

  integer :: ti

  call Print("IPModel_LMTO_TBE : LMTO_TBE Potential", file=file)
  call Print("IPModel_LMTO_TBE : n_types = " // this%n_types, file=file)

  do ti=1, this%n_types
    call Print ("IPModel_LMTO_TBE : type " // ti // " atomic_num " // this%atomic_num(ti), file=file)
    call verbosity_push_decrement()
    call Print ("IPModel_LMTO_TBE", file=file)
    call verbosity_pop()
  end do

end subroutine IPModel_LMTO_TBE_Print

subroutine IPModel_LMTO_TBE_read_params_xml(this, param_str)
  type(IPModel_LMTO_TBE), intent(inout), target :: this
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
    call system_abort("IPModel_LMTO_TBE_read_params_xml parsed file, but n_types = 0")
  endif

end subroutine IPModel_LMTO_TBE_read_params_xml

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

  integer ti

  if (name == 'LMTO_TBE_params') then ! new LMTO_TBE stanza

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
        call system_abort("Can't find n_types in LMTO_TBE_params")
      endif

      allocate(parse_ip%atomic_num(parse_ip%n_types))
      parse_ip%atomic_num = 0

      call QUIP_FoX_get_value(attributes, "control_file", value, status)
      if (status /= 0) call system_abort ("IPModel_LMTO_TBE_read_params_xml cannot find control_file")
      read (value, *) parse_ip%control_file
      
    endif

  elseif (parse_in_ip .and. name == 'per_type_data') then

    call QUIP_FoX_get_value(attributes, "type", value, status)
    if (status /= 0) call system_abort ("IPModel_LMTO_TBE_read_params_xml cannot find type")
    read (value, *) ti

    call QUIP_FoX_get_value(attributes, "atomic_num", value, status)
    if (status /= 0) call system_abort ("IPModel_LMTO_TBE_read_params_xml cannot find atomic_num")
    read (value, *) parse_ip%atomic_num(ti)

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
    if (name == 'LMTO_TBE_params') then
      parse_in_ip = .false.
    end if
  endif

end subroutine IPModel_endElement_handler

end module IPModel_LMTO_TBE_module
