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
!X IPModel_MTP
!X
!% MTP interatomic potential: interface to Alex Shapeev's C++ code. 
!%
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#include "error.inc"

module IPModel_MTP_module

use iso_c_binding, only : C_NULL_CHAR
use error_module
use system_module, only : dp, inoutput, print, operator(//)
use dictionary_module
use paramreader_module
use linearalgebra_module
use atoms_types_module
use atoms_module
use CInOutput_module

use mpi_context_module
use QUIP_Common_module

implicit none
private

include 'IPModel_interface.h'

public :: IPModel_MTP
type IPModel_MTP
  real(dp) :: cutoff = 0.0_dp
end type IPModel_MTP

logical, private :: parse_in_ip, parse_matched_label
type(IPModel_MTP), private, pointer :: parse_ip

interface Initialise
  module procedure IPModel_MTP_Initialise_str
end interface Initialise

interface Finalise
  module procedure IPModel_MTP_Finalise
end interface Finalise

interface Print
  module procedure IPModel_MTP_Print
end interface Print

interface Calc
  module procedure IPModel_MTP_Calc
end interface Calc

#ifdef HAVE_MTP
interface
   subroutine alex_initialise(filename, cutoff, error) bind(C)
      use iso_c_binding, only: c_double, c_char, c_int
      character(kind=c_char) :: filename(*)
      real(kind=c_double) :: cutoff
      integer(kind=c_int) :: error
   endsubroutine alex_initialise
endinterface

interface
   subroutine alex_compute(neigh_count,neigh_len,input_pos,out_force,out_site_en, error) bind(C)
      use iso_c_binding, only: c_ptr, c_int
      integer(kind=c_int) :: neigh_count
      type(c_ptr), value :: neigh_len
      type(c_ptr), value :: input_pos
      type(c_ptr), value :: out_force
      type(c_ptr), value :: out_site_en
      integer(kind=c_int) :: error
   endsubroutine alex_compute
endinterface
#endif

contains

subroutine IPModel_MTP_Initialise_str(this, args_str, param_str, error)
  type(IPModel_MTP), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str
  type(Dictionary) :: params
  integer, optional, intent(out):: error

  character(len=STRING_LENGTH) :: alex_filename
  integer :: mtp_error

  INIT_ERROR(error)
  call Finalise(this)

  call initialise(params)
  call param_register(params, 'alex_filename', '', alex_filename, help_string='filename for coefficiets')
  if(.not. param_read_line(params, args_str, ignore_unknown=.true., task='IPModel_MTP_Initialise args_str')) then
     RAISE_ERROR("IPModel_MTP_Init failed to parse args_str='"//trim(args_str)//"'", error)
  end if
  call finalise(params)

#ifdef HAVE_MTP
  call alex_initialise(trim(alex_filename)//C_NULL_CHAR, this%cutoff,mtp_error)
  if( mtp_error /= 0 ) then
     RAISE_ERROR("IPModel_MTP_Init received error "//mtp_error//" from alex_initialise()",error)
  endif
#else
  RAISE_ERROR("support for MTP not compiled in. Obtain MTP code first, and HAVE_MTP = 1 in your Makefile.inc, and recompile QUIP",error)
#endif

end subroutine IPModel_MTP_Initialise_str

subroutine IPModel_MTP_Finalise(this)
  type(IPModel_MTP), intent(inout) :: this

  ! Add finalisation code here

end subroutine IPModel_MTP_Finalise


subroutine IPModel_MTP_Calc(this, at, e, local_e, f, virial, local_virial, args_str, mpi, error)
   use iso_c_binding, only : c_int, c_double, c_loc
   type(IPModel_MTP), intent(inout):: this
   type(Atoms), intent(inout)      :: at
   real(dp), intent(out), optional :: e, local_e(:)
   real(dp), intent(out), optional :: f(:,:), local_virial(:,:)   !% Forces, dimensioned as \texttt{f(3,at%N)}, local virials, dimensioned as \texttt{local_virial(9,at%N)} 
   real(dp), intent(out), optional :: virial(3,3)
   character(len=*), optional      :: args_str
   type(MPI_Context), intent(in), optional :: mpi
   integer, intent(out), optional :: error

   ! Add calc() code here


   integer :: i, j, n, neighbour_index
   real(dp) :: r_ij
   real(dp), dimension(3) :: diff_ij, f_ij
   integer(kind=c_int), dimension(:), allocatable, target :: alex_n_neighbours
   real(kind=c_double), dimension(:), allocatable, target :: alex_energy, alex_r_neighbours, alex_f_neighbours
   integer :: mtp_error

   INIT_ERROR(error)

   allocate(alex_n_neighbours(at%N), alex_energy(at%N))

   do i = 1, at%N
      alex_n_neighbours(i) = n_neighbours(at,i,max_dist=this%cutoff)
   enddo

   allocate(alex_r_neighbours(3*sum(alex_n_neighbours)), alex_f_neighbours(3*sum(alex_n_neighbours)))

   neighbour_index = 0

   do i = 1, at%N
      do n = 1, n_neighbours(at,i)
         j = neighbour(at, i, n, distance=r_ij, diff=diff_ij)
         if( r_ij > this%cutoff ) cycle

         neighbour_index = neighbour_index + 1
         alex_r_neighbours(3*(neighbour_index-1)+1:3*neighbour_index) = diff_ij
      enddo
   enddo

#ifdef HAVE_MTP
   call alex_compute(at%N,c_loc(alex_n_neighbours), c_loc(alex_r_neighbours), c_loc(alex_f_neighbours), c_loc(alex_energy), mtp_error)
   if( mtp_error /= 0 ) then
      call write(at, 'stdout', prefix='MTP_ERROR')
      RAISE_ERROR("IPModel_MTP_Calc received error "//mtp_error//" from alex_compute()",error)
   endif
#endif

   if (present(e)) e = sum(alex_energy)
   if (present(local_e)) then
      call check_size('Local_E',local_e,(/at%N/),'IPModel_MTP_Calc', error)
      local_e = alex_energy
   endif

   if (present(f).or.present(virial)) then
      if(present(f)) then
         call check_size('Force',f,(/3,at%N/),'IPModel_MTP_Calc', error)
         f = 0.0_dp
      endif
      if(present(virial)) virial = 0.0_dp
      
      neighbour_index = 0
      do i = 1, at%N
         do n = 1, n_neighbours(at,i)
            j = neighbour(at, i, n, distance=r_ij, diff=diff_ij)
            if( r_ij > this%cutoff ) cycle

            neighbour_index = neighbour_index + 1
            f_ij = alex_f_neighbours(3*(neighbour_index-1)+1:3*neighbour_index)
            if(present(f)) then
               f(:,j) = f(:,j) - f_ij
               f(:,i) = f(:,i) + f_ij
            endif
            if(present(virial)) then
               virial = virial - ( diff_ij .outer. f_ij )
            endif
         enddo
      enddo

   end if

   if (present(local_virial)) then
      call check_size('Local_virial',local_virial,(/9,at%Nbuffer/),'IPModel_MTP_Calc', error)
      local_virial = 0.0_dp
   endif

   if(allocated(alex_n_neighbours)) deallocate(alex_n_neighbours)
   if(allocated(alex_energy)) deallocate(alex_energy)
   if(allocated(alex_r_neighbours)) deallocate(alex_r_neighbours)
   if(allocated(alex_f_neighbours)) deallocate(alex_f_neighbours)

end subroutine IPModel_MTP_Calc


subroutine IPModel_MTP_Print(this, file)
  type(IPModel_MTP), intent(in) :: this
  type(Inoutput), intent(inout),optional :: file

  call Print("IPModel_MTP : MTP Potential", file=file)
  call Print("IPModel_MTP : cutoff = " // this%cutoff, file=file)

end subroutine IPModel_MTP_Print


end module IPModel_MTP_module
