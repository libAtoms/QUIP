! H0 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! H0 X
! H0 X   libAtoms+QUIP: atomistic simulation library
! H0 X
! H0 X   Portions of this code were written by
! H0 X     Albert Bartok-Partay, Silvia Cereda, Gabor Csanyi, James Kermode,
! H0 X     Ivan Solt, Wojciech Szlachta, Csilla Varnai, Steven Winfield.
! H0 X
! H0 X   Copyright 2006-2017.
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
!X IPModel_DispTS
!X
!% Pairwise dispersion correction from Tkatchenko and Scheffler,
!% PRL 102(7) 073005 (2009).
!%
!% Free-atom reference data must be supplied in the XML file; Hirshfeld volumes
!% must be supplied in the Atoms object.
!%
!% Free-atom reference data is available e.g. from Chu and Dalgarno,
!% JCP 121, 4083 (2004) (used in original implementation) or from
!% Gould and Bučko, JCTC 12(8), 3603 (2016) (better coverage)
!%
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#include "error.inc"

module IPModel_DispTS_module

use error_module
use system_module, only : dp, inoutput, print, verbosity_push_decrement, verbosity_pop, operator(//)
use dictionary_module
use paramreader_module
use linearalgebra_module
use atoms_types_module
use atoms_module

use mpi_context_module
use QUIP_Common_module

implicit none
private

include 'IPModel_interface.h'

public :: IPModel_DispTS
type IPModel_DispTS
  integer :: n_types = 0
  integer, allocatable :: atomic_num(:), type_of_atomic_num(:)
  real(dp), allocatable :: c6_free(:), alpha_free(:), r_vdW_free(:)

  real(dp) :: cutoff = 0.0_dp

  character(len=STRING_LENGTH) :: label

end type IPModel_DispTS

logical, private :: parse_in_ip, parse_matched_label
type(IPModel_DispTS), private, pointer :: parse_ip

interface Initialise
  module procedure IPModel_DispTS_Initialise_str
end interface Initialise

interface Finalise
  module procedure IPModel_DispTS_Finalise
end interface Finalise

interface Print
  module procedure IPModel_DispTS_Print
end interface Print

interface Calc
  module procedure IPModel_DispTS_Calc
end interface Calc

contains

subroutine IPModel_DispTS_Initialise_str(this, args_str, param_str)
  type(IPModel_DispTS), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str

  type(Dictionary) :: params

  call Finalise(this)

  call initialise(params)
  this%label=''
  call param_register(params, 'label', '', this%label, help_string="No help yet.  This source file was $LastChangedBy$")
  if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_DispTS_Initialise_str args_str')) then
    call system_abort("IPModel_DispTS_Initialise_str failed to parse label from args_str="//trim(args_str))
  endif
  call finalise(params)

  call IPModel_DispTS_read_params_xml(this, param_str)

  !  Add initialisation code here

end subroutine IPModel_DispTS_Initialise_str

subroutine IPModel_DispTS_Finalise(this)
  type(IPModel_DispTS), intent(inout) :: this

  ! Add finalisation code here

  if (allocated(this%atomic_num)) deallocate(this%atomic_num)
  if (allocated(this%type_of_atomic_num)) deallocate(this%type_of_atomic_num)
  if (allocated(this%c6_free)) deallocate(this%c6_free)
  if (allocated(this%alpha_free)) deallocate(this%alpha_free)
  if (allocated(this%r_vdW_free)) deallocate(this%r_vdW_free)

  this%n_types = 0
  this%label = ''
end subroutine IPModel_DispTS_Finalise


subroutine IPModel_DispTS_Calc(this, at, e, local_e, f, virial, local_virial, args_str, mpi, error)
   type(IPModel_DispTS), intent(inout):: this
   type(Atoms), intent(inout)      :: at
   real(dp), intent(out), optional :: e, local_e(:)
   real(dp), intent(out), optional :: f(:,:), local_virial(:,:)   !% Forces, dimensioned as \texttt{f(3,at%N)}, local virials, dimensioned as \texttt{local_virial(9,at%N)} 
   real(dp), intent(out), optional :: virial(3,3)
   character(len=*), optional      :: args_str
   type(MPI_Context), intent(in), optional :: mpi
   integer, intent(out), optional :: error

   real(dp), pointer, dimension(:) :: my_hirshfeld_volume
   character(STRING_LENGTH)        :: hirshfeld_vol_name

   INIT_ERROR(error)

   if (present(e)) e = 0.0_dp
   if (present(local_e)) then
      call check_size('Local_E',local_e,(/at%N/),'IPModel_DispTS_Calc', error)
      local_e = 0.0_dp
   endif
   if (present(f)) then
      call check_size('Force',f,(/3,at%Nbuffer/),'IPModel_DispTS_Calc', error)
      f = 0.0_dp
   end if
   if (present(virial)) virial = 0.0_dp
   if (present(local_virial)) then
      call check_size('Local_virial',local_virial,(/9,at%Nbuffer/),'IPModel_DispTS_Calc', error)
      local_virial = 0.0_dp
   endif

   if (present(args_str)) then
       if (len_trim(args_str) > 0) then
           call initialise(params)
           call param_register(params, 'hirshfeld_vol_name', 'hirshfeld_rel_volume', hirshfeld_vol_name, &
                               help_string='Name of the Atoms property containing relative Hirshfeld volumes $v/v_{free}$')

           if (.not. param_read_line(params, args_str, ignore_unknown=.true., task='IPModel_MBD_Calc args_str')) then
               RAISE_ERROR("IPModel_DispTS_Calc failed to parse args_str '"//trim(args_str)//"'", error)
           endif
           call finalise(params)
           call assign_property_pointer(at, trim(hirshfeld_vol_name), my_hirshfeld_volume, error)
           PASS_ERROR_WITH_INFO("IPModel_DispTS_Calc could not find '"//trim(hirshfeld_vol_name)//"' property in the Atoms object", error)
        endif
   else
       call assign_property_pointer(at, 'hirshfeld_rel_volume', my_hirshfeld_volume, error)
   endif

   !TODO smooth cutoff
   !TODO tie in with Ewald code for long range

   RAISE_ERROR('IPModel_Calc - not implemented',error)

end subroutine IPModel_DispTS_Calc


subroutine IPModel_DispTS_Print(this, file)
  type(IPModel_DispTS), intent(in) :: this
  type(Inoutput), intent(inout),optional :: file

  integer :: ti

  call Print("IPModel_DispTS : T-S dispersion correction potential", file=file)
  call Print("IPModel_DispTS : n_types = " // this%n_types // " cutoff = " // this%cutoff, file=file)

  do ti=1, this%n_types
    call Print("IPModel_DispTS : type " // ti // " atomic_num " // this%atomic_num(ti), file=file)
    call Print("IPModel_DispTS : type " // ti // " free-atom: C6 " // this%c6_free(ti) &
               // " polarizability " // this%alpha_free(ti) &
               // " vdW radius " // this%r_vdW_free(ti), file=file)
  end do

end subroutine IPModel_DispTS_Print

subroutine IPModel_DispTS_read_params_xml(this, param_str)
  type(IPModel_DispTS), intent(inout), target :: this
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
    call system_abort("IPModel_DispTS_read_params_xml parsed file, but n_types = 0")
  endif

end subroutine IPModel_DispTS_read_params_xml

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

  if (name == 'DispTS_params') then ! new DispTS stanza

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
        call system_abort("Can't find n_types in DispTS_params")
      endif

      allocate(parse_ip%atomic_num(parse_ip%n_types))
      allocate(parse_ip%c6_free(parse_ip%n_types))
      allocate(parse_ip%alpha_free(parse_ip%n_types))
      allocate(parse_ip%r_vdW_free(parse_ip%n_types))
      parse_ip%atomic_num = 0
      parse_ip%c6_free = 0.0_dp
      parse_ip%alpha_free = 0.0_dp
      parse_ip%r_vdW_free = 0.0_dp

      call QUIP_FoX_get_value(attributes, "cutoff", value, status)
      if (status /= 0) call system_abort ("IPModel_DispTS_read_params_xml cannot find cutoff")
      read (value, *) parse_ip%cutoff
    endif


  elseif (parse_in_ip .and. name == 'per_type_data') then

    call QUIP_FoX_get_value(attributes, "type", value, status)
    if (status /= 0) call system_abort ("IPModel_DispTS_read_params_xml cannot find type")
    read (value, *) ti

    call QUIP_FoX_get_value(attributes, "atomic_num", value, status)
    if (status /= 0) call system_abort ("IPModel_DispTS_read_params_xml cannot find atomic_num")
    read (value, *) parse_ip%atomic_num(ti)

    call QUIP_FoX_get_value(attributes, "c6_free", value, status)
    if (status /= 0) call system_abort ("IPModel_DispTS_read_params_xml cannot find c6_free")
    read (value, *) parse_ip%c6_free(ti)

    call QUIP_FoX_get_value(attributes, "alpha_free", value, status)
    if (status /= 0) call system_abort ("IPModel_DispTS_read_params_xml cannot find alpha_free")
    read (value, *) parse_ip%alpha_free(ti)

    call QUIP_FoX_get_value(attributes, "r_vdW_free", value, status)
    if (status /= 0) call system_abort ("IPModel_DispTS_read_params_xml cannot find r_vdW_free")
    read (value, *) parse_ip%r_vdW_free(ti)

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
    if (name == 'DispTS_params') then
      parse_in_ip = .false.
    end if
  endif

end subroutine IPModel_endElement_handler

end module IPModel_DispTS_module
