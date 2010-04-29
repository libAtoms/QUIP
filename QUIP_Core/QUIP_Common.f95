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
!X QUIP_Common_module
!X
!% The QUIP_Common module contains functions to get the atom type from 
!% the atomic number or
!% a token from a string.  
!%  
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

module QUIP_Common_module

use FoX_sax, only: xml_t, dictionary_t, haskey, getvalue, open_xml_string, close_xml_t, parse
use System_module

implicit none

private

public :: get_type
public :: xml_t, dictionary_t, open_xml_string, close_xml_t, parse
public :: QUIP_FoX_get_value

public :: pauli_sigma
complex(dp), parameter :: pauli_sigma(2,2,3) = reshape( (/ &
    cmplx(0.0_dp,0.0_dp,dp), cmplx(1.0_dp,0.0_dp,dp), cmplx(1.0_dp,0.0_dp,dp), cmplx(0.0_dp,0.0_dp,dp), &
    cmplx(0.0_dp,0.0_dp,dp), cmplx(0.0_dp,1.0_dp,dp), -cmplx(0.0_dp,1.0_dp,dp), cmplx(0.0_dp,0.0_dp,dp), &
    cmplx(1.0_dp,0.0_dp,dp), cmplx(0.0_dp,0.0_dp,dp), cmplx(0.0_dp,0.0_dp,dp), cmplx(-1.0_dp,0.0_dp,dp) /) , &
    (/ 2, 2, 3 /) )
    

contains

!% Get the type from the atomic number, with error checking.
function get_type(type_of_atomic_num, Z)
  integer, intent(in) :: type_of_atomic_num(:)
  integer, intent(in) :: Z
  integer get_type

  if (Z < 1 .or. Z > size(type_of_atomic_num)) &
    call system_abort ('get_type: Atomic number '//Z //' out of range')
  if (type_of_atomic_num(Z) == 0) &
    call system_abort ('get_type: Atomic number '//Z//' does not correspond to a defined type')

  get_type = type_of_atomic_num(Z)

end function get_type

subroutine QUIP_FoX_get_value(attributes, key, val, status)
  type(dictionary_t), intent(in) :: attributes
  character(len=*), intent(in) :: key
  character(len=*), intent(inout) :: val
  integer, intent(out), optional :: status

  if (HasKey(attributes,key)) then
    val = GetValue(attributes, trim(key))
    if (present(status)) status = 0
  else
    val = ""
    if (present(status)) status = 1
  endif
end subroutine

end module QUIP_Common_module
