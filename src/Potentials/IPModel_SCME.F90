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
!X IPModel_SCME
!X
!% SCME module for the single-centre multipole expansion water module, developed by Thor Wikfeldt
!%
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#include "error.inc"

module IPModel_SCME_module

use error_module
use system_module, only : dp, inoutput, print, verbosity_push_decrement, verbosity_pop, operator(//)
use dictionary_module
use linearalgebra_module, only: is_diagonal
use topology_module, only: find_water_monomer
use paramreader_module
use linearalgebra_module
use atoms_types_module
use atoms_module

use mpi_context_module
use QUIP_Common_module

#ifdef HAVE_SCME
use scme, only : scme_calculate
#endif

implicit none
private

include 'IPModel_interface.h'

public :: IPModel_SCME
type IPModel_SCME
  !integer :: n_types = 0
  !integer, allocatable :: atomic_num(:), type_of_atomic_num(:)

  real(dp) :: cutoff = 0.0_dp
  logical :: full_interaction_order, use_repulsion, use_PS_PES

  character(len=STRING_LENGTH) :: label

end type IPModel_SCME

logical, private :: parse_in_ip, parse_matched_label
type(IPModel_SCME), private, pointer :: parse_ip

interface Initialise
  module procedure IPModel_SCME_Initialise_str
end interface Initialise

interface Finalise
  module procedure IPModel_SCME_Finalise
end interface Finalise

interface Print
  module procedure IPModel_SCME_Print
end interface Print

interface Calc
  module procedure IPModel_SCME_Calc
end interface Calc

contains

subroutine IPModel_SCME_Initialise_str(this, args_str, param_str)
  type(IPModel_SCME), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str

  type(Dictionary) :: params

  call Finalise(this)

  call initialise(params)
  this%label=''
  call param_register(params, 'label', '', this%label, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'full_interaction_order', 'T', this%full_interaction_order, help_string="Whether to truncate the interaction order calculation at 5th order")
  call param_register(params, 'use_repulsion', 'F', this%use_repulsion, help_string="Whether to use repulsion in SCME")
  call param_register(params, 'use_PS_PES', 'T', this%use_PS_PES, help_string="Whether to use the PS potential energy surface")
  if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_SCME_Initialise_str args_str')) then
    call system_abort("IPModel_SCME_Initialise_str failed to parse label from args_str="//trim(args_str))
  endif
  call finalise(params)

  if( trim(this%label) /= "version_20170802" ) then
     call system_abort("IPModel_SCME_Initialise_str: SCME with updated parameters/damping. Make sure your potential is compatible. Proceed with caution, email Albert for instructions if in doubt.")
  endif
  !call IPModel_SCME_read_params_xml(this, param_str)
  this%cutoff = 2.0_dp

  !  Add initialisation code here

end subroutine IPModel_SCME_Initialise_str

subroutine IPModel_SCME_Finalise(this)
  type(IPModel_SCME), intent(inout) :: this

  ! Add finalisation code here
  this%cutoff = 0.0_dp

  this%label = ''
end subroutine IPModel_SCME_Finalise







!/////////////////////////////////
subroutine IPModel_SCME_Calc(this, at, e, local_e, f, virial, local_virial, args_str, mpi, error)
   type(IPModel_SCME), intent(inout):: this
   type(Atoms), intent(inout)      :: at
   real(dp), intent(out), optional :: e, local_e(:)
   real(dp), intent(out), optional :: f(:,:)            !j
   real(dp), intent(out), optional :: local_virial(:,:) !j
   !j real(dp), intent(out), optional :: f(:,:), local_virial(:,:)   !% Forces, dimensioned as \texttt{f(3,at%N)}, local virials, dimensioned as \texttt{local_virial(9,at%N)} 
   real(dp), intent(out), optional :: virial(3,3)
   character(len=*), optional      :: args_str
   type(MPI_Context), intent(in), optional :: mpi
   integer, intent(out), optional :: error

   integer :: nWater, i, a, m
   integer, dimension(:,:), allocatable :: water_monomer_index
   real(dp) :: uTot
   real(dp), dimension(3) :: lattice
   real(dp), dimension(:), allocatable :: raOri
   real(dp), dimension(:,:), allocatable :: coords !j
   real(dp), dimension(:,:), allocatable :: fa     !j

   INIT_ERROR(error)
#ifndef HAVE_SCME
   RAISE_ERROR('IPModel_SCME_Calc - not linked to the scme library',error)
#endif

   if (present(e)) e = 0.0_dp
   
   if (present(local_e)) then
      RAISE_ERROR('IPModel_SCME_Calc - local energies not implemented',error)
      call check_size('Local_E',local_e,(/at%N/),'IPModel_SCME_Calc', error)
      local_e = 0.0_dp
   endif
   
   if (present(f)) then
      RAISE_ERROR('IPModel_SCME_Calc - forces not implemented',error)
      call check_size('Force',f,(/3,at%Nbuffer/),'IPModel_SCME_Calc', error)
      f = 0.0_dp
   end if
   
   if (present(virial)) then
      RAISE_ERROR('IPModel_SCME_Calc - virials not implemented',error)
      virial = 0.0_dp
   endif

   if (present(local_virial)) then
      RAISE_ERROR('IPModel_SCME_Calc - local virials not implemented',error)
      call check_size('Local_virial',local_virial,(/9,at%Nbuffer/),'IPModel_SCME_Calc', error)
      local_virial = 0.0_dp
   endif

   nWater = count(at%Z==8)
   allocate(water_monomer_index(3,nWater))
   call find_water_monomer(at,water_monomer_index,error=error)

   !j allocate(raOri(3*at%N),fa(3*at%N)) 
   allocate(raOri(3*at%N),fa(3,at%N)) !j
   allocate( coords(3,at%N) ) !j
   
   if( is_diagonal(at%lattice) ) then
      do i = 1, 3
         lattice(i) = at%lattice(i,i)
      enddo
   else
      RAISE_ERROR('IPModel_SCME_Calc - lattice must be orthorhombic',error)
   endif

   do i = 1, nWater !j molecules.  I should change raOrig to at(3,n_atoms) or at(3,3,n_water)
      
      ! This is the new 3 by N_atoms coordinate matrix
      coords(:,(i-1)*3+1) = at%pos(:,water_monomer_index(2,i)) !j coords is hho ordered in scme
      coords(:,(i-1)*3+2) = at%pos(:,water_monomer_index(3,i))
      coords(:,(i-1)*3+3) = at%pos(:,water_monomer_index(1,i))
      
      ! This is the old version:
      do a = 1, 3 !j xyz
         raOri( (i-1)*6 + a ) = at%pos(a,water_monomer_index(2,i))            !j h
         raOri( (i-1)*6 + a + 3 ) = at%pos(a,water_monomer_index(3,i))        !j h
         raOri( (i-1)*3 + nWater*6 + a ) = at%pos(a,water_monomer_index(1,i)) !j o
      enddo
   enddo
 
#ifdef HAVE_SCME
   !call scme_calculate(at%N,raOri, lattice, fa, uTot)
   call scme_calculate(at%N,coords, lattice, fa, uTot,in_FULL=this%full_interaction_order, in_USE_REP=this%use_repulsion, in_USE_PS_PES=this%use_PS_PES)
#endif
   
   
   if(present(f))then 
    do m = 1,nWater !jö Rearange from HHO (in fa in SCME) to OHH ( in f in QUIP)
      f(:,(m-1)*3+1)=fa(:,(m-1)*3+3) !O
      f(:,(m-1)*3+2)=fa(:,(m-1)*3+1) !H1
      f(:,(m-1)*3+3)=fa(:,(m-1)*3+2) !H1
    enddo
   endif 
   
   !if(present(f)) f = fa !j copy over the forces to the quip-defined array f from scme-specified fa ! How come the order gets correct? 
   
   if (present(e)) e = uTot !j if (present): Refers to optional arguments

   if(allocated(water_monomer_index)) deallocate(water_monomer_index)
   if(allocated(raOri)) deallocate(raOri)
   if(allocated(fa)) deallocate(fa)

end subroutine IPModel_SCME_Calc
!///////////////////////////////////////////////////////////






subroutine IPModel_SCME_Print(this, file)
  type(IPModel_SCME), intent(in) :: this
  type(Inoutput), intent(inout),optional :: file

  integer :: ti

  call Print("IPModel_SCME : SCME Potential", file=file)

end subroutine IPModel_SCME_Print

subroutine IPModel_SCME_read_params_xml(this, param_str)
  type(IPModel_SCME), intent(inout), target :: this
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

end subroutine IPModel_SCME_read_params_xml

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

  if (name == 'SCME_params') then ! new SCME stanza

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

    if(parse_in_ip) then
       call finalise(parse_ip)
    endif
  endif

end subroutine IPModel_startElement_handler

subroutine IPModel_endElement_handler(URI, localname, name)
  character(len=*), intent(in)   :: URI
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name

  if (parse_in_ip) then
    if (name == 'SCME_params') then
      parse_in_ip = .false.
    end if
  endif

end subroutine IPModel_endElement_handler

end module IPModel_SCME_module
