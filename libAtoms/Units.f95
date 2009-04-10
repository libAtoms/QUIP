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
!X  Units module
!X  
!%  This module holds a collection of our units and conversion factors.
!%  We use Angstroms (\AA), eVs, and femtoseconds (fs).
!%
!% \begin{description}
!%  \item[Length:]   Angstroms    $=    \mathrm{a.u.} \times \mathtt{Bohr}$
!%  \item[Energy:]   eV           $=    \mathrm{a.u.} \times \mathtt{Hartree}$
!%  \item[Time:]     fs
!%  \item[Mass:]     $ E T^2 / L^2 = \mathrm{a.u.} \times \mathtt{Hartree} \times \mathtt{AU\_fs}^2 / \mathtt{Bohr}^2 $
!% \end{description}
!%
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

!X $Id: Units.f95,v 1.11 2008-06-09 16:45:48 ab686 Exp $

!X $Log: not supported by cvs2svn $
!X Revision 1.10  2007/04/18 12:43:14  jrk33
!X Updated doc comment
!X
!X Revision 1.9  2007/04/17 17:26:18  jrk33
!X Changed module name to lower case
!X
!X Revision 1.8  2007/04/17 09:57:19  gc121
!X put copyright statement in each file
!X
!X Revision 1.7  2007/04/11 15:44:17  saw44
!X Updated/Added comments for the documentation generator
!X
!X Revision 1.6  2007/03/12 17:04:25  jrk33
!X Added GPA conversion constant
!X
!X Revision 1.5  2007/03/01 13:51:46  jrk33
!X Documentation comments reformatted and edited throughout. Anything starting "!(no space)%"
!  is picked up by the documentation generation script
!X
!X Revision 1.4  2007/02/28 15:47:22  saw44
!X Added RADIANS_PER_DEGREE
!X
!X Revision 1.3  2007/02/19 11:04:34  saw44
!X Added DEGREES_PER_RADIAN
!X
!X Revision 1.2  2006/12/04 11:56:39  nb326
!X Add Rydberg
!X
!X Revision 1.1.1.1  2006/12/04 11:11:30  gc121
!X Imported sources
!X
!X Revision 1.7  2006/06/20 17:23:19  gc121
!X added new copyright notice to include James, Gian, Mike and Alessandro
!X
!X Revision 1.6  2006/06/08 14:04:33  saw44
!X Added Avogadros number and conversion from kcal/mol to eV
!X
!X Revision 1.5  2006/01/27 16:13:17  gc121
!X fixed units to be in eV,A,fs
!X
!X Revision 1.4  2006/01/19 15:02:12  gc121
!X added copyright headers and cvs magic tags to files that were missing them
!X

module units_module

use system_module ! for definition of dp
implicit none

real(dp), parameter :: ELECTRONMASS_GPERMOL =  5.48579903e-4_dp !% grams/mol
real(dp), parameter :: ELEM_CHARGE = 1.60217653e-19_dp !% coulombs
real(dp), parameter :: HARTREE = 27.2113961_dp !% eV
real(dp), parameter :: RYDBERG = 0.5_dp*HARTREE !% eV
real(dp), parameter :: BOHR = 0.529177249_dp !% Angstrom
real(dp), parameter :: HBAR_EVSEC = 6.5821220e-16_dp !% hbar in eV seconds
real(dp), parameter :: HBAR_AU = 1.0_dp              !% hbar in a.u.
real(dp), parameter :: HBAR = (HBAR_EVSEC*1e+15_dp)    !% hbar in eV fs
real(dp), parameter :: ONESECOND = 1e15_dp           !% 1 second in fs
real(dp), parameter :: ONESECOND_AU = (1.0_dp/(HBAR_EVSEC/(HBAR_AU*HARTREE))) !% 1 second in a.u.
real(dp), parameter :: AU_FS = (1.0_dp/ONESECOND_AU*ONESECOND) !% a.u. time in fs
real(dp), parameter :: MASSCONVERT = (1.0_dp/ELECTRONMASS_GPERMOL*HARTREE*AU_FS*AU_FS/(BOHR*BOHR)) !% = 1e7 / (N_A * ELEM_CHARGE)
real(dp), parameter :: BOLTZMANN_K = 8.617385e-5_dp !% eV/Kelvin
real(dp), parameter :: PI = 3.14159265358979323846264338327950288_dp
real(dp), parameter :: N_A = 6.0221479e23_dp !% Avogadro's number
real(dp), parameter :: KCAL_MOL = 4.3383e-2_dp !% eV
real(dp), parameter :: DEGREES_PER_RADIAN = 180.0_dp / PI
real(dp), parameter :: RADIANS_PER_DEGREE = PI / 180.0_dp
real(dp), parameter :: GPA = 1.6022e-19_dp*1.0e30_dp/1.0e9_dp !% Convert from \textsc{libAtoms} units to Gigapascals
real(dp), parameter :: EPSILON_0 = 8.854187817e-12_dp / ELEM_CHARGE * 1.0e-10_dp !% epsilon_0 in e / V Angstrom

complex(dp), parameter :: CPLX_ZERO = (0.0_dp,0.0_dp)
complex(dp), parameter :: CPLX_IMAG = (0.0_dp,1.0_dp)
complex(dp), parameter :: CPLX_ONE  = (1.0_dp,0.0_dp)

end module units_module
