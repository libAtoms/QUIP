!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X     libAtoms: atomistic simulation library
!X     
!X     Copyright 2006-2007.
!X
!X     Authors: Gabor Csanyi, Steven Winfield, James Kermode
!X     Contributors: Noam Bernstein, Alessio Comisso
!X
!X     The source code is released under the GNU General Public License,
!X     version 2, http://www.gnu.org/copyleft/gpl.html
!X
!X     If you would like to license the source code under different terms,
!X     please contact Gabor Csanyi, gabor@csanyi.net
!X
!X     When using this software, please cite the following reference:
!X
!X     http://www.libatoms.org
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X  test_force_model module
!X  
!%  Given at atoms structure, 'test_force_model' does a numerical test
!%  of the force against the finite difference of the energies
!%
!%  The name of the force model to be tested is hard-wired here.
!X 
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

! $Id: test_force_model.f95,v 1.8 2007-04-18 10:28:47 jrk33 Exp $

! $Log: not supported by cvs2svn $
! Revision 1.7  2007/04/17 17:38:01  jrk33
! Standardised subroutine and function references and printing argument order.
!
! Revision 1.6  2007/04/17 09:57:19  gc121
! put copyright statement in each file
!
! Revision 1.5  2007/04/11 15:44:17  saw44
! Updated/Added comments for the documentation generator
!
! Revision 1.4  2007/03/12 17:08:55  jrk33
! Reformatted docs
!
! Revision 1.3  2007/03/01 13:51:46  jrk33
! Documentation comments reformatted and edited throughout. Anything starting "!(no space)%"
!  is picked up by the documentation generation script
!
! Revision 1.2  2007/01/03 16:58:13  jrk33
! Changed #ifdef MPI to #ifdef _MPI to reflect changes in libAtoms
!
! Revision 1.1  2006/12/12 15:20:17  gc121
! added a generic gradient model test
!
! Revision 1.8  2006/10/16 09:26:31  gc121
! added calconnectfast call
!
! Revision 1.7  2006/06/29 17:57:01  gc121
! sync
!
! Revision 1.6  2006/06/29 11:10:43  jrk33
! Removed unused variables
!
! Revision 1.5  2006/06/23 10:51:24  jrk33
! f = FORCE_ROUTINE(at) -> call FORCE_ROUTINE(at, f)
!
! Revision 1.4  2006/02/16 23:08:46  gc121
! typo
!
! Revision 1.3  2006/02/15 18:59:35  gc121
! added BKS potential call
!
! Revision 1.2  2006/02/15 14:17:57  saw44
! Added save attribute to "at" variable after email from Alessio saying xlf compiler requires this
!
! Revision 1.1  2006/01/30 10:43:42  gc121
! finite diff. test on a force model
!
!

! Because FORTRAN lacks usable function pointers, we have to 
! hard-wire the model to be tested. Put the actual function
! names here. The corresponding module also needs to be "use"-d

module test_force_model_module
  use system_module
  use atoms_module
  use minimization_module
  type(atoms), private, save::at

  real(dp), allocatable, dimension(:) :: test_force_weights

#define ENERGY_ROUTINE dummy_energy
#define ENERGY_ROUTINE_NAME "dummy_energy"
#define FORCE_ROUTINE dummy_force
#define FORCE_ROUTINE_NAME "dummy_force"

contains

  !% OMIT
  function e(p)
    real(dp)::p(:)
    real(dp)::e

    at%pos = reshape(p, (/3,at%N/))
    call calc_dists(at)
    e = ENERGY_ROUTINE(at) ! energy function

  end function e
  
  !% OMIT
  function eprime(p)
    real(dp)::p(:)
    real(dp)::eprime(size(p))
    real(dp)::f(3,at%N)

    at%pos = reshape(p, (/3,at%N/))
    call calc_dists(at)
    call FORCE_ROUTINE(at, f) ! force function

    eprime = -reshape(f, (/at%N*3/)) 
  end function eprime
  
  subroutine test_force_model(a)
    logical::t
    type(atoms)::a
    at = a
    call calc_connect_fast(at)
    t = test_gradient(reshape(at%pos, (/at%N*3/)), e, eprime, NERD)
  end subroutine test_force_model
  

  !% OMIT
  function dummy_energy(at) result(e)
    type(Atoms), intent(in) :: at
    real(dp) :: e

    e = 0.0_dp

  end function dummy_energy


  !% OMIT
  subroutine dummy_force(at, f)
    type(Atoms), intent(in) :: at
    real(dp), dimension(3,at%N) :: f

    f = 0.0_dp
  end subroutine dummy_force


end module Test_Force_Model_module

