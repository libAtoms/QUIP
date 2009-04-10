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
