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
!X  Periodic Table module
!X  
!%  This module contains a list of elements, their masses and covalent radii.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

module periodictable_module

use system_module      ! for definition of real(dp)
use units_module
!
! The Periodic Table
!

implicit none
private

public :: ElementName, ElementMass, ElementValence, ElementCovRad, atomic_number, atomic_number_from_symbol, atomic_number_from_mass

character(3),parameter,dimension(0:116) :: ElementName =   (/"xx ",                                    &
   "H  ","He ","Li ","Be ","B  ","C  ","N  ","O  ","F  ","Ne ","Na ","Mg ","Al ","Si ","P  ","S  ",    &
   "Cl ","Ar ","K  ","Ca ","Sc ","Ti ","V  ","Cr ","Mn ","Fe ","Co ","Ni ","Cu ","Zn ","Ga ","Ge ",    &
   "As ","Se ","Br ","Kr ","Rb ","Sr ","Y  ","Zr ","Nb ","Mo ","Tc ","Ru ","Rh ","Pd ","Ag ","Cd ",    &
   "In ","Sn ","Sb ","Te ","I  ","Xe ","Cs ","Ba ","La ","Ce ","Pr ","Nd ","Pm ","Sm ","Eu ","Gd ",    &
   "Tb ","Dy ","Ho ","Er ","Tm ","Yb ","Lu ","Hf ","Ta ","W  ","Re ","Os ","Ir ","Pt ","Au ","Hg ",    &
   "Tl ","Pb ","Bi ","Po ","At ","Rn ","Fr ","Ra ","Ac ","Th ","Pa ","U  ","Np ","Pu ","Am ","Cm ",    &
   "Bk ","Cf ","Es ","Fm ","Md ","No ","Lr ","Rf ","Db ","Sg ","Bh ","Hs ","Mt ","Ds ","Rg ","Uub",    &
   "Uut","Uuq","Uup","Uuh" /) !% Mapping of atomic number to element name


! NOTE: constants used in array initializers below are SINGLE
! PRECISION, so values of ElementMass and ElementCovRad are slightly
! incorrect, e.g. masses differ by ~1.0e-7 from results of double
! precision multiplication of atomic masses by MASSCONVERT. Adding
! "_dp" after each constant fixes this but causes a number of
! regression tests to fail.

! Units: grams per Mole * MASSCONVERT (conforming to eV,A,fs system)

real(dp),parameter,dimension(116) :: ElementMass =                                                     &
(/1.00794, 4.00260, 6.941, 9.012187, 10.811, 12.0107, 14.00674, 15.9994, 18.99840, 20.1797, 22.98977,  &
24.3050, 26.98154, 28.0855, 30.97376, 32.066, 35.4527, 39.948, 39.0983, 40.078, 44.95591, 47.867,      &
50.9415, 51.9961, 54.93805, 55.845, 58.93320, 58.6934, 63.546, 65.39, 69.723, 72.61, 74.92160, 78.96,  &
79.904, 83.80, 85.4678, 87.62, 88.90585, 91.224, 92.90638, 95.94, 98.0, 101.07, 102.90550, 106.42,     &
107.8682, 112.411, 114.818, 118.710, 121.760, 127.60, 126.90447, 131.29, 132.90545, 137.327, 138.9055, &
140.116, 140.90765, 144.24, 145.0, 150.36, 151.964, 157.25, 158.92534, 162.50, 164.93032, 167.26,      &
168.93421, 173.04, 174.967, 178.49, 180.9479, 183.84, 186.207, 190.23, 192.217, 195.078, 196.96655,    &
200.59, 204.3833, 207.2, 208.98038, 209.0, 210.0, 222.0, 223.0, 226.0, 227.0, 232.0381, 231.03588,     &
238.0289, 237.0, 244.0, 243.0, 247.0, 247.0, 251.0, 252.0, 257.0, 258.0, 259.0, 262.0, 261.0, 262.0,   &
263.0, 264.0, 265.0, 268.0, 271.0, 272.0, 285.0, 284.0, 289.0, 288.0, 292.0/)*MASSCONVERT 
!% Element mass in grams per Mole $\times$ 'MASSCONVERT' (conforming to eV,\AA,fs unit system).

! Units: Angstroms

real(dp),parameter,dimension(116) :: ElementCovRad =                                                   &
(/0.320,0.310,1.630,0.900,0.820,0.770,0.750,0.730,0.720,0.710,1.540,1.360,1.180,1.110,1.060,1.020,     &
0.990,0.980,2.030,1.740,1.440,1.320,1.220,1.180,1.170,1.170,1.160,1.150,1.170,1.250,1.260,1.220,1.200, &
1.160,1.140,1.120,2.160,1.910,1.620,1.450,1.340,1.300,1.270,1.250,1.250,1.280,1.340,1.480,1.440,1.410, &
1.400,1.360,1.330,1.310,2.350,1.980,1.690,1.650,1.650,1.840,1.630,1.620,1.850,1.610,1.590,1.590,1.580, &
1.570,1.560,2.000,1.560,1.440,1.340,1.300,1.280,1.260,1.270,1.300,1.340,1.490,1.480,1.470,1.460,1.460, &
2.000,2.000,2.000,2.000,2.000,1.650,2.000,1.420,2.000,2.000,2.000,2.000,2.000,2.000,2.000,2.000,2.000, &
2.000,2.000,2.000,2.000,2.000,2.000,2.000,2.000,2.000,2.000,2.000,2.000,2.000,2.000,2.000/)
!% Covalent radii in \AA.

integer,parameter,dimension(116) :: ElementValence =  &
(/1,-1, 1, 2, 3, 4, 3, 2, 1,-1, 1, 2, 3, 4, 3, 2,     &
  1,-1, 1, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,     &
 -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,     &
 -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,     &
 -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,     &
 -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,     &
 -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,     &
 -1,-1,-1,-1/)


interface atomic_number
   !% Do a reverse lookup of atomic number from either symbol or mass
   module procedure atomic_number_from_symbol, atomic_number_from_mass
end interface atomic_number

contains

  !Look up the atomic number for a given atomic symbol
  function atomic_number_from_symbol(atomic_symbol)
    character(*), intent(in) :: atomic_symbol
    integer                  :: atomic_number_from_symbol
    integer                  :: i

    if (verify(trim(adjustl(atomic_symbol)),"0123456789") == 0) then ! an integer
       read (atomic_symbol, *) atomic_number_from_symbol
       if (atomic_number_from_symbol < 1 .or. atomic_number_from_symbol > size(ElementName)) then
	  atomic_number_from_symbol = 0
       endif
       return
    else ! not an integer, hopefully an element abbreviation
       do i = 1, 116
	  if (trim(lower_case(adjustl(atomic_symbol)))==trim(lower_case(ElementName(i)))) then
	     atomic_number_from_symbol = i
	     return
	  end if
       end do
    end if

    !If unsuccessful, return 0
    atomic_number_from_symbol = 0

  end function atomic_number_from_symbol

  !Look up the atomic number for a given atomic mass (IN GRAMS PER MOLE)
  !Note: this may fail for some of the transuranic elements... so put those ununpentium simulations on hold for a while ;-)
  function atomic_number_from_mass(atomic_mass)
    real(dp), intent(in) :: atomic_mass
    integer              :: atomic_number_from_mass
    integer              :: i
    real(dp), parameter  :: TOL = 0.01_dp

    do i = 1, 116
       if (abs(atomic_mass - ElementMass(i)/MASSCONVERT) < TOL) then
          atomic_number_from_mass = i
          return
       end if
    end do

    !If unsuccessful, return 0
    atomic_number_from_mass = 0

  end function atomic_number_from_mass

  !ElementName formatting, used by atoms_read_xyz
  !First leter uppercase, others lowercase
  function ElementFormat(lower_UPPER) result(UPPER_lower)
    character(*), intent(in)             :: lower_UPPER
    character(len=len_trim(lower_UPPER)) :: UPPER_lower
    character(len=*), parameter          :: lc = 'abcdefghijklmnopqrstuvwxyz', &
                                            UC = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    integer                              :: i,j

    UPPER_lower = lower_UPPER
    j = index(lc,lower_UPPER(1:1))
    if (j>0) UPPER_lower(1:1) = UC(j:j)
    do i = 2, len_trim(lower_UPPER)
       j = index(UC,lower_UPPER(i:i))
       if (j>0) UPPER_lower(i:i) = lc(j:j)
    enddo

  end function ElementFormat

end module periodictable_module
