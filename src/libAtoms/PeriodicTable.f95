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

! IUPAC 2017 element symbols
character(3),parameter,dimension(0:118) :: ElementName =   (/"xx ",                                             &
   "H  ","He ",                                                                                                 &
   "Li ","Be ","B  ","C  ","N  ","O  ","F  ","Ne ",                                                             &
   "Na ","Mg ","Al ","Si ","P  ","S  ","Cl ","Ar ",                                                             &
   "K  ","Ca ","Sc ","Ti ","V  ","Cr ","Mn ","Fe ","Co ","Ni ","Cu ","Zn ","Ga ","Ge ","As ","Se ","Br ","Kr ", &
   "Rb ","Sr ","Y  ","Zr ","Nb ","Mo ","Tc ","Ru ","Rh ","Pd ","Ag ","Cd ","In ","Sn ","Sb ","Te ","I  ","Xe ", &
   "Cs ","Ba ","La ","Ce ","Pr ","Nd ","Pm ","Sm ","Eu ","Gd ","Tb ","Dy ","Ho ","Er ","Tm ","Yb ","Lu ","Hf ", &
   "Ta ","W  ","Re ","Os ","Ir ","Pt ","Au ","Hg ","Tl ","Pb ","Bi ","Po ","At ","Rn ",                         &
   "Fr ","Ra ","Ac ","Th ","Pa ","U  ","Np ","Pu ","Am ","Cm ","Bk ","Cf ","Es ","Fm ","Md ","No ","Lr ","Rf ", &
   "Db ","Sg ","Bh ","Hs ","Mt ","Ds ","Rg ","Cn ","Nh ","Fl ","Mc ","Lv ","Ts ","Og "                          &
   /) !% Mapping of atomic number to element name


! NOTE: constants used in array initializers below are SINGLE
! PRECISION, so values of ElementMass and ElementCovRad are slightly
! incorrect, e.g. masses differ by ~1.0e-7 from results of double
! precision multiplication of atomic masses by MASSCONVERT. Adding
! "_dp" after each constant fixes this but causes a number of
! regression tests to fail.

! Units: grams per Mole * MASSCONVERT (conforming to eV,A,fs system)

real(dp),parameter,dimension(118) :: ElementMass = (/                                                 &
1.008_dp,4.002602_dp,6.94_dp,9.0121831_dp,10.81_dp,12.011_dp,14.007_dp,15.999_dp,18.998403163_dp,     &
20.1797_dp,22.98976928_dp,24.305_dp,26.9815385_dp,28.085_dp,30.973761998_dp,32.06_dp,35.45_dp,        &
39.948_dp,39.0983_dp,40.078_dp,44.955908_dp,47.867_dp,50.9415_dp,51.9961_dp,54.938044_dp,55.845_dp,   &
58.933194_dp,58.6934_dp,63.546_dp,65.38_dp,69.723_dp,72.63_dp,74.921595_dp,78.971_dp,79.904_dp,       &
83.798_dp,85.4678_dp,87.62_dp,88.90584_dp,91.224_dp,92.90637_dp,95.95_dp,97.90721_dp,101.07_dp,       &
102.9055_dp,106.42_dp,107.8682_dp,112.414_dp,114.818_dp,118.71_dp,121.76_dp,127.6_dp,126.90447_dp,    &
131.293_dp,132.90545196_dp,137.327_dp,138.90547_dp,140.116_dp,140.90766_dp,144.242_dp,144.91276_dp,   &
150.36_dp,151.964_dp,157.25_dp,158.92535_dp,162.5_dp,164.93033_dp,167.259_dp,168.93422_dp,173.054_dp, &
174.9668_dp,178.49_dp,180.94788_dp,183.84_dp,186.207_dp,190.23_dp,192.217_dp,195.084_dp,196.966569_dp,&
200.592_dp,204.38_dp,207.2_dp,208.9804_dp,208.98243_dp,209.98715_dp,222.01758_dp,223.01974_dp,        &
226.02541_dp,227.02775_dp,232.0377_dp,231.03588_dp,238.02891_dp,237.04817_dp,244.06421_dp,            &
243.06138_dp,247.07035_dp,247.07031_dp,251.07959_dp,252.083_dp,257.09511_dp,258.09843_dp,259.101_dp,  &
262.11_dp,267.122_dp,268.126_dp,271.134_dp,270.133_dp,269.1338_dp,278.156_dp,281.165_dp,281.166_dp,   &
285.177_dp,286.182_dp,289.19_dp,289.194_dp,293.204_dp,293.208_dp,294.214_dp /)*MASSCONVERT
!% Element mass in grams per Mole $\times$ 'MASSCONVERT' (conforming to eV,\AA,fs unit system).

! Units: Angstroms

real(dp),parameter,dimension(118) :: ElementCovRad =                                                   &
(/0.320,0.310,1.630,0.900,0.820,0.770,0.750,0.730,0.720,0.710,1.540,1.360,1.180,1.110,1.060,1.020,     &
0.990,0.980,2.030,1.740,1.440,1.320,1.220,1.180,1.170,1.170,1.160,1.150,1.170,1.250,1.260,1.220,1.200, &
1.160,1.140,1.120,2.160,1.910,1.620,1.450,1.340,1.300,1.270,1.250,1.250,1.280,1.340,1.480,1.440,1.410, &
1.400,1.360,1.330,1.310,2.350,1.980,1.690,1.650,1.650,1.840,1.630,1.620,1.850,1.610,1.590,1.590,1.580, &
1.570,1.560,2.000,1.560,1.440,1.340,1.300,1.280,1.260,1.270,1.300,1.340,1.490,1.480,1.470,1.460,1.460, &
2.000,2.000,2.000,2.000,2.000,1.650,2.000,1.420,2.000,2.000,2.000,2.000,2.000,2.000,2.000,2.000,2.000, &
2.000,2.000,2.000,2.000,2.000,2.000,2.000,2.000,2.000,2.000,2.000,2.000,2.000,2.000,2.000,2.000,2.000/)
!% Covalent radii in \AA.

integer,parameter,dimension(118) :: ElementValence =  &
(/1,-1, 1, 2, 3, 4, 3, 2, 1,-1, 1, 2, 3, 4, 3, 2,     &
  1,-1, 1, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,     &
 -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,     &
 -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,     &
 -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,     &
 -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,     &
 -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,     &
 -1,-1,-1,-1,-1,-1/)


interface atomic_number
   !% Do a reverse lookup of atomic number from either symbol or mass
   module procedure atomic_number_from_symbol, atomic_number_from_mass
end interface atomic_number

contains

  !Look up the atomic number for a given atomic symbol
  elemental function atomic_number_from_symbol(atomic_symbol)
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
       do i = 1, 118
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
  !Note: this may fail for some of the transuranic elements...
  !so put those ununpentium simulations on hold for a while ;-)
  function atomic_number_from_mass(atomic_mass)
    real(dp), intent(in) :: atomic_mass
    integer              :: atomic_number_from_mass
    integer              :: i
    real(dp), parameter  :: TOL = 0.01_dp

    do i = 1, 118
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
