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
!X  Units module
!X  
!%  This module holds a collection of our units and conversion factors.
!%  We use Angstroms ($\mathrm{\AA}$), eV, and femtoseconds (fs).
!%
!%  - Length:   Angstroms, ``a.u. * BOHR``
!%  - Energy:   eV,        ``a.u. * HARTREE``
!%  - Time:     fs
!%  - Mass:     $E T^2/L^2$ = ``a.u. * HARTREE * AU_FS**2 / BOHR**2``
!%
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

module units_module

use system_module, only : dp ! for definition of dp
implicit none
public

#ifdef LEGACY_UNITS
! Values of units pre-2017.
! Since results will vary, keep these around for a while in case
! someone needs consistency with these.
! Tests will fail if these values are used.
real(dp), parameter :: ELECTRONMASS_GPERMOL =  5.48579903e-4_dp !% grams/mol
real(dp), parameter :: ELEM_CHARGE = 1.60217653e-19_dp !% coulombs
real(dp), parameter :: HARTREE = 27.2113961_dp !% eV
real(dp), parameter :: RYDBERG = 0.5_dp*HARTREE !% eV
real(dp), parameter :: BOHR = 0.529177249_dp !% Angstrom
real(dp), parameter :: CUBIC_BOHR = BOHR*BOHR/HARTREE !% QUIP Units of polarisability
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
real(dp), parameter :: DEBYE = 1.0e-21_dp/299792458.0_dp/ELEM_CHARGE*1e10_dp !% 1D $= 10^{-18}$ statcoulomb-centrimetre in e-A
real(dp), parameter :: SQRT_TWO = sqrt(2.0_dp)
! not included originally, make sure they're not undefined
real(dp), parameter :: GRAM = 1e-3_dp  !% in kg
real(dp), parameter :: ELECTRONMASS = ELECTRONMASS_GPERMOL*GRAM/N_A
real(dp), parameter :: EPSILON_0_AU = EPSILON_0*ELEM_CHARGE/1.0e-10_dp
real(dp), parameter :: VACUUM_C = 299792458.0_dp  !% exact in m/s
#else
! CODATA 2014 taken from
! http://arxiv.org/pdf/1507.07956.pdf
! Fundamentals
real(dp), parameter :: ELECTRONMASS = 9.10938356e-31_dp  !% electron mass / kg
real(dp), parameter :: N_A = 6.022140857e23_dp  !% Avogadro constant
real(dp), parameter :: ELEM_CHARGE = 1.6021766208e-19_dp  !% e / Coulombs
real(dp), parameter :: HARTREE = 27.21138602_dp  !% EH / eV
real(dp), parameter :: BOHR = 0.52917721067_dp  !% in Angstrom
real(dp), parameter :: HBAR_EVSEC = 6.582119514e-16_dp  !% hbar in eV seconds
real(dp), parameter :: BOLTZMANN_K = 8.6173303e-5_dp  !% eV/kelvin
real(dp), parameter :: EPSILON_0_AU = 8.854187817e-12_dp  !% exact in F/m
real(dp), parameter :: VACUUM_C = 299792458.0_dp  !% exact in m/s

! Others
real(dp), parameter :: GRAM = 1e-3_dp  !% in kg
real(dp), parameter :: PI = 3.14159265358979323846264338327950288_dp
real(dp), parameter :: SQRT_TWO = sqrt(2.0_dp)

! derived units for QUIP
real(dp), parameter :: ELECTRONMASS_GPERMOL = ELECTRONMASS*N_A/GRAM  !% grams/mol
real(dp), parameter :: RYDBERG = 0.5_dp*HARTREE !% eV
real(dp), parameter :: CUBIC_BOHR = BOHR*BOHR/HARTREE !% QUIP Units of polarisability
real(dp), parameter :: HBAR_AU = 1.0_dp              !% hbar in a.u.
real(dp), parameter :: HBAR = (HBAR_EVSEC*1e+15_dp)    !% hbar in eV fs
real(dp), parameter :: ONESECOND = 1e15_dp           !% 1 second in fs
real(dp), parameter :: ONESECOND_AU = (1.0_dp/(HBAR_EVSEC/(HBAR_AU*HARTREE))) !% 1 second in a.u.
real(dp), parameter :: AU_FS = (1.0_dp/ONESECOND_AU*ONESECOND) !% a.u. time in fs
real(dp), parameter :: MASSCONVERT = (1.0_dp/ELECTRONMASS_GPERMOL*HARTREE*AU_FS*AU_FS/(BOHR*BOHR)) !% = 1e7 / (N_A * ELEM_CHARGE)
real(dp), parameter :: KCAL_MOL = 4.184*1000/(ELEM_CHARGE*N_A) !% Thermochemical definition in eV
real(dp), parameter :: DEGREES_PER_RADIAN = 180.0_dp / PI
real(dp), parameter :: RADIANS_PER_DEGREE = PI / 180.0_dp
real(dp), parameter :: GPA = ELEM_CHARGE*1.0e30_dp/1.0e9_dp !% Convert from \textsc{libAtoms} units to Gigapascals
real(dp), parameter :: EPSILON_0 = EPSILON_0_AU/ELEM_CHARGE*1.0e-10_dp !% epsilon_0 in e / V Angstrom
real(dp), parameter :: DEBYE = 1.0e-21_dp/VACUUM_C/ELEM_CHARGE*1e10_dp !% 1D $= 10^{-18}$ statcoulomb-centrimetre in e-A
#endif

complex(dp), parameter :: CPLX_ZERO = (0.0_dp,0.0_dp)
complex(dp), parameter :: CPLX_IMAG = (0.0_dp,1.0_dp)
complex(dp), parameter :: CPLX_ONE  = (1.0_dp,0.0_dp)

! Make the values of some compiler macros available at runtime
#ifdef HAVE_LOTF
integer, parameter :: have_lotf = 1
#else
integer, parameter :: have_lotf = 0
#endif

#ifdef HAVE_CP2K
integer, parameter :: have_cp2k = 1
#else
integer, parameter :: have_cp2k = 0
#endif

#ifdef HAVE_NETCDF
integer, parameter :: have_netcdf = 1
#else
integer, parameter :: have_netcdf = 0
#endif

end module units_module
