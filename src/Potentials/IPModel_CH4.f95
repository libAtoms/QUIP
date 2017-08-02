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
!X IPModel_CH4
!X
!% CH4 module for use when implementing a new Interatomic Potential
!%
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#include "error.inc"

module IPModel_CH4_module

use error_module
use system_module, only : dp, inoutput, print, verbosity_push_decrement, verbosity_pop, operator(//)
use dictionary_module
use paramreader_module
use linearalgebra_module
use atoms_types_module
use units_module
use atoms_module
use topology_module

use mpi_context_module
use QUIP_Common_module

implicit none
private

include 'IPModel_interface.h'

public :: IPModel_CH4
type IPModel_CH4
  real(dp) :: cutoff = 0.0_dp

  character(len=STRING_LENGTH) :: label

end type IPModel_CH4

logical, private :: parse_in_ip, parse_matched_label
type(IPModel_CH4), private, pointer :: parse_ip

interface Initialise
  module procedure IPModel_CH4_Initialise_str
end interface Initialise

interface Finalise
  module procedure IPModel_CH4_Finalise
end interface Finalise

interface Print
  module procedure IPModel_CH4_Print
end interface Print

interface Calc
  module procedure IPModel_CH4_Calc
end interface Calc

interface
      function ch4pot_func_10(r)
      real(8), intent(in) :: r(10)
      real(8) :: ch4pot_func_10
      endfunction ch4pot_func_10
endinterface

contains

subroutine IPModel_CH4_Initialise_str(this, args_str, param_str)
  type(IPModel_CH4), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str

  type(Dictionary) :: params

  call Finalise(this)

  call initialise(params)
  this%label=''
  call param_register(params, 'label', '', this%label, help_string="No help yet.  This source file was $LastChangedBy$")
  if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_CH4_Initialise_str args_str')) then
    call system_abort("IPModel_CH4_Initialise_str failed to parse label from args_str="//trim(args_str))
  endif
  call finalise(params)

  this%cutoff = 2.0

end subroutine IPModel_CH4_Initialise_str

subroutine IPModel_CH4_Finalise(this)
  type(IPModel_CH4), intent(inout) :: this

  ! Add finalisation code here

  this%label = ''
  this%cutoff = 0.0_dp
end subroutine IPModel_CH4_Finalise


subroutine IPModel_CH4_Calc(this, at, e, local_e, f, virial, local_virial, args_str, mpi, error)
   type(IPModel_CH4), intent(inout):: this
   type(Atoms), intent(inout)      :: at
   real(dp), intent(out), optional :: e, local_e(:)
   real(dp), intent(out), optional :: f(:,:), local_virial(:,:)   !% Forces, dimensioned as \texttt{f(3,at%N)}, local virials, dimensioned as \texttt{local_virial(9,at%N)} 
   real(dp), intent(out), optional :: virial(3,3)
   character(len=*), optional      :: args_str
   type(MPI_Context), intent(in), optional :: mpi
   integer, intent(out), optional :: error

   integer, dimension(5) :: signature
   integer, dimension(:,:), allocatable :: monomer_index
   logical, dimension(:), allocatable :: is_associated
   real(dp), dimension(3,4) :: pos
   real(dp), dimension(4,4) :: bmat, bmati
   real(dp), dimension(10) :: rten
   integer :: i, j, iC, iH, n_monomers

   INIT_ERROR(error)

   if (present(e)) e = 0.0_dp
   if (present(local_e)) then
      call check_size('Local_E',local_e,(/at%N/),'IPModel_CH4_Calc', error)
      local_e = 0.0_dp
   endif
   if (present(f)) then
      call check_size('Force',f,(/3,at%Nbuffer/),'IPModel_CH4_Calc', error)
      f = 0.0_dp
   end if
   if (present(virial)) virial = 0.0_dp
   if (present(local_virial)) then
      call check_size('Local_virial',local_virial,(/9,at%Nbuffer/),'IPModel_CH4_Calc', error)
      local_virial = 0.0_dp
   endif

   signature = (/6,1,1,1,1/)
   allocate(is_associated(at%N))

   
   is_associated = .false.
   call find_general_monomer(at,monomer_index,signature,is_associated,this%cutoff,error=error)

   n_monomers = size(monomer_index,2)

   call ch4_initialise(bmat,bmati)

   do i = 1, n_monomers
      iC = monomer_index(1,i)
      do j = 1, 4
         iH = monomer_index(j+1,i)
         pos(:,j) = diff_min_image(at,iC,iH) / BOHR
      enddo
      call bond_cart2radau_int10(pos,bmati,rten)
      if(present(e)) e = e + ch4pot_func_10(rten) * INVERSE_CM
   enddo
      
   if(allocated(monomer_index)) deallocate(monomer_index)
   if(allocated(is_associated)) deallocate(is_associated)

end subroutine IPModel_CH4_Calc


subroutine IPModel_CH4_Print(this, file)
  type(IPModel_CH4), intent(in) :: this
  type(Inoutput), intent(inout),optional :: file

  integer :: ti

  call Print("IPModel_CH4 : CH4 Potential", file=file)
  call Print("IPModel_CH4 : cutoff = " // this%cutoff, file=file)

end subroutine IPModel_CH4_Print

end module IPModel_CH4_module
