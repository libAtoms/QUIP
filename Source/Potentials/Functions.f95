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

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X Functions Module 
!X
!% The Functions Module contains useful subroutines and functions
!% to compute the Fermi occupation $n(E)$ of given energy levels,
!% and error functions. 
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
module functions_module

use System_module

implicit none
private

public erf, erfc
public :: f_Fermi, f_Fermi_deriv, smooth_cutoff

contains

!% Calculate the Fermi occupation for en energy level ('E'),
!% at a given temperature ('T'),
!% and chemical potential ('mu').
elemental function f_Fermi(mu, T, E)
  real(dp), intent(in) :: mu, T, E
  real(dp) :: f_Fermi

  if (abs(T) < 1.0e-10_dp) then
    if (E-mu < 0.0_dp) then
      f_Fermi = 1.0_dp
    else if (E-mu > 0.0_dp) then
      f_Fermi = 0.0_dp
    else
      f_Fermi = 0.5_dp
    endif
  else
    f_Fermi = 1.0_dp/(1.0_dp+exp((E-mu)/T))
  endif
end function f_Fermi

elemental function f_Fermi_deriv(mu, T, E)
  real(dp), intent(in) :: mu, T, E
  real(dp) :: f_Fermi_deriv

  real(dp) :: arg

  if (abs(T) < 1.0e-10_dp) then
    f_Fermi_deriv = 0.0_dp
  else
    arg = abs((E-mu)/T)
    f_Fermi_deriv = -1.0_dp/((1.0_dp+exp(-arg))**2 * exp(arg) * T)
  endif
end function f_Fermi_deriv


!| double precision function erf
!| calculate error function of argument
!| returns error function of argument
!|
!|I double precision x: argument
!|
!| Noam Bernstein 9/12/2001

double precision function erf(x)
implicit none
    double precision x

    if (x .lt. 0.0D0) then
	erf = -gammp_half(x*x)
    else
	erf = gammp_half(x*x)
    end if

end function

!| double precision function erfc
!| calculate error function complement of argument
!| returns error function complement of argument
!|
!|I double precision x: argument
!|
!| Noam Bernstein 9/12/2001

double precision function erfc(x)
implicit none
    double precision x

    if (x .lt. 0.0D0) then
	erfc = 1.0D0 + gammp_half(x*x)
    else
	erfc = gammq_half(x*x)
    end if
end function

!| double precision function gammp_half
!| calculate incomplete gamma p function of (0.5,x)
!| returns incomplete gamma p function of (0.5,x)
!|
!|I double precision x: argument
!| 
!| Noam Bernstein 9/12/2001

double precision function gammp_half(x)
implicit none
    double precision x

    double precision Gamma_half ! = sqrt(pi)
    parameter ( Gamma_half = 1.77245385090552D0 )

    gammp_half = 2.0D0/Gamma_half * sqrt(x) * p_F0(x)
    
end function

!| double precision function gammq_half
!| calculate incomplete gamma q function (complement of gamma p) of (0.5,x)
!| returns incomplete gamma q function (complement of gamma p) of (0.5,x)
!|
!|I double precision x: argument
!| 
!| Noam Bernstein 9/12/2001

double precision function gammq_half(x)
implicit none
    double precision x

    double precision Gamma_half ! = sqrt(pi)
    parameter ( Gamma_half = 1.77245385090552D0 )

    if (x .ne. 0) then
	gammq_half = 2.0D0/Gamma_half * sqrt(x) * p_F0C(x)
    else
	gammq_half = 1.0D0
    end if
    
end function

!| double precision function p_F0(x)
!| calculate F_m(x) for m=0 using FFFMT3
!| returns value of F_0(x)
!|
!|I double precision x: argument
!|
!| Noam Bernstein 9/13/2001

double precision function p_F0(x)
implicit none
    double precision x

    double precision xx(1)
    double precision s1(1,5)

    xx(1) = x
    call FFFMT3_F90(1, 1, xx, s1)

    p_F0 = s1(1,1)

end function

!| double precision function p_F0C(x)
!| calculate sqrt(pi/x)/2-F_m(x) for m=0 using FFFMT3C
!| returns value of sqrt(pi/x)/2-F_0(x)
!|
!|I double precision x: argument
!|
!| Noam Bernstein 9/13/2001

double precision function p_F0C(x)
implicit none
    double precision x

    double precision xx(1)
    double precision s1(1,5)

    xx(1) = x
    call FFFMT3C_F90(1, 1, xx, s1)

    p_F0C = s1(1,1)

end function

!| subroutine FFFMT3_F90(x)
!% Calculate $F_m(x)$ for $m=0$ using 'FFFMT3C'
!|
!|I integer MX, NPTS: maximum, actual number of pts to be evaluated
!|I double precision xx(NPTS): arguments
!|O double precision s1(NPTS,5): returned values
!|
!| Noam Bernstein 9/13/2001
!| From Mark Pederson's code

subroutine FFFMT3_F90(MX,NPTS,XX,S1)
! MODIFIED BY MARK R PEDERSON (1990)
! To FORTRAN 90 by Noam Bernstein
implicit none
    integer MX, NPTS
    double precision XX(MX),S1(MX,5)

    logical, save:: FIRST = .true.
    double precision, save:: F(10,1601)

    double precision ONEOVRI(5)
    double precision G(20)

    double precision:: ROOTPI = 1.77245385090551602728D0
    double precision:: ONEOVER3 = 1.0D0/3.0D0

    integer I, L, K, KK, IT, IPTS
    double precision T, X, DEN, TERM, SUM, FACT, Q, EPS, EX
    double precision RT, POL1, RECT1, DT

    IF (NPTS.GT.MX) THEN
	PRINT *,'FFFMT3: MX MUST BE AT LEAST: ',NPTS
	stop
    END IF
    IF (FIRST) THEN
	! ONEOVER3 = 1.0D0/3.0D0
	FACT = 1.0D0
	DO I = 1,5
	    ONEOVRI(I) = 1.0D0/I
	    FACT = FACT*ONEOVRI(I)
	END DO
	FIRST = .FALSE.
	EPS = 1.0D-12
	DO I = 1,10
	    F(I,1) = 1.0D0/(2*I-1)
	END DO
	T = 0.0D0
	DO L = 2,1601
	    T = (1.0D-2)*(L-1)
	    X = 2*T
	    DEN = 39.0D0
	    TERM = 1.0D0/DEN
	    SUM = TERM
	    !DO 130 I = 2,100
	    DO I = 2,100
		DEN = DEN+2.0D0
		TERM = TERM*X/DEN
		SUM = SUM+TERM
		Q = TERM/SUM
		! IF (Q-EPS)  140,140,130
		if (Q .le. EPS) exit
	    END DO
	    ! 130   CONTINUE
	    ! 140   EX = EXP(-T)
	    EX = EXP(-T)
	    G(20) = EX*SUM
	    !
	    !  USE DOWNWARD RECURSION
	    !
	    DO I = 1,19
		K = 21-I
		KK = K-1
		G(KK) = (X*G(K)+EX)/(2*KK-1)
	    END DO
	    DO I = 1,10
		F(I,L) = G(I)
	    END DO
	END DO
	DO L = 1,1601
	    F(5,L) = F(5,L)/2.0D0
	    F(6,L) = F(6,L)/6.0D0
	END DO
    END IF
    !
    DO IPTS = 1,NPTS
	IF (XX(IPTS) .GT. 16.0D0) THEN
	    POL1 =  -EXP(-XX(IPTS))
	    RT = SQRT(XX(IPTS))*ROOTPI
	    RECT1 = 0.5D0/XX(IPTS)
	    S1(IPTS,1) = RECT1*(RT+POL1)
	    S1(IPTS,2) = RECT1*(S1(IPTS,1)+POL1)
	    S1(IPTS,3) = RECT1*(3*S1(IPTS,2)+POL1)
	ELSE
	    IT = INT(100*XX(IPTS)+1.5D0)
	    POL1 = EXP(-XX(IPTS))
	    RECT1 = 2*XX(IPTS)
	    DT = (1.0D-2)*(IT-1)-XX(IPTS)
	    S1(IPTS,3) = F(3,IT)+DT*(F(4,IT)+DT*(F(5,IT)+DT*F(6,IT)))
	    S1(IPTS,2) = (S1(IPTS,3)*RECT1+POL1)*ONEOVER3
	    S1(IPTS,1) = (S1(IPTS,2)*RECT1+POL1)
	END IF
    END DO
    RETURN
end subroutine

!| subroutine FFFMT3C_F90(x)
!% Calculate $\frac{1}{2}\sqrt{\frac{\pi}{x}}-F_m(x)$ for $m=0$ using FFFMT3C
!|
!|I integer MX, NPTS: maximum, actual number of pts to be evaluated
!|I double precision xx(NPTS): arguments
!|O double precision s1(NPTS,5): returned values
!|
!| Noam Bernstein 9/13/2001
!| Modified for F_0(x) complement from Mark Pederson's code

subroutine FFFMT3C_F90(MX,NPTS,XX,S1)
! MODIFIED BY MARK R PEDERSON (1990)
! To FORTRAN 90 by Noam Bernstein
implicit none
    integer MX, NPTS
    double precision XX(MX),S1(MX,5)

    logical, save:: FIRST = .true.
    double precision, save:: F(10,1601)

    double precision ONEOVRI(5)
    double precision G(20)

    double precision:: ROOTPI = 1.77245385090551602728D0
    double precision:: ONEOVER3 = 1.0D0/3.0D0

    integer I, L, K, KK, IT, IPTS
    double precision T, X, DEN, TERM, SUM, FACT, Q, EPS, EX
    double precision RT, POL1, RECT1, DT

    IF (NPTS.GT.MX) THEN
	PRINT *,'FFFMT3: MX MUST BE AT LEAST: ',NPTS
	stop
    END IF
    IF (FIRST) THEN
	! ONEOVER3 = 1.0D0/3.0D0
	FACT = 1.0D0
	DO I = 1,5
	    ONEOVRI(I) = 1.0D0/I
	    FACT = FACT*ONEOVRI(I)
	END DO
	FIRST = .FALSE.
	EPS = 1.0D-12
	DO I = 1,10
	    F(I,1) = 1.0D0/(2*I-1)
	END DO
	T = 0.0D0
	DO L = 2,1601
	    T = (1.0D-2)*(L-1)
	    X = 2*T
	    DEN = 39.0D0
	    TERM = 1.0D0/DEN
	    SUM = TERM
	    !DO 130 I = 2,100
	    DO I = 2,100
		DEN = DEN+2.0D0
		TERM = TERM*X/DEN
		SUM = SUM+TERM
		Q = TERM/SUM
		! IF (Q-EPS)  140,140,130
		if (Q .le. EPS) exit
	    END DO
	    ! 130   CONTINUE
	    ! 140   EX = EXP(-T)
	    EX = EXP(-T)
	    G(20) = EX*SUM
	    !
	    !  USE DOWNWARD RECURSION
	    !
	    DO I = 1,19
		K = 21-I
		KK = K-1
		G(KK) = (X*G(K)+EX)/(2*KK-1)
	    END DO
	    DO I = 1,10
		F(I,L) = G(I)
	    END DO
	END DO
	DO L = 1,1601
	    F(5,L) = F(5,L)/2.0D0
	    F(6,L) = F(6,L)/6.0D0
	END DO
    END IF
    !
    DO IPTS = 1,NPTS
	IF (XX(IPTS) .GT. 16.0D0) THEN
	    POL1 =  -EXP(-XX(IPTS))
	    RT = SQRT(XX(IPTS))*ROOTPI
	    RECT1 = 0.5D0/XX(IPTS)
	    !! S1(IPTS,1) = RECT1*(RT+POL1)
	    ! for computing the complement of F0
	    S1(IPTS,1) = -RECT1*POL1
	    ! S1(IPTS,2) = RECT1*(S1(IPTS,1)+POL1)
	    ! S1(IPTS,3) = RECT1*(3*S1(IPTS,2)+POL1)
	ELSE
	    IT = INT(100*XX(IPTS)+1.5D0)
	    POL1 = EXP(-XX(IPTS))
	    RECT1 = 2*XX(IPTS)
	    DT = (1.0D-2)*(IT-1)-XX(IPTS)
	    S1(IPTS,3) = F(3,IT)+DT*(F(4,IT)+DT*(F(5,IT)+DT*F(6,IT)))
	    S1(IPTS,2) = (S1(IPTS,3)*RECT1+POL1)*ONEOVER3
	    S1(IPTS,1) = (S1(IPTS,2)*RECT1+POL1)
	    ! for computing the complement of F0
	    S1(IPTS,1) = ROOTPI/(2.0D0*sqrt(XX(IPTS))) - S1(IPTS,1)
	END IF
    END DO
    RETURN
end subroutine


!% Smooth cutoff function from D.J. Cole \emph{et al.}, J. Chem. Phys. {\bf 127}, 204704 (2007).
subroutine smooth_cutoff(x,R,D,fc,dfc_dx)

  real(dp), intent(in) :: x,R,D
  real(dp), intent(out) :: fc,dfc_dx
  real(dp), parameter :: pi = dacos(-1.0d0)

  if (x .lt. (R-D)) then
     fc = 1.0d0
     dfc_dx = 0.0d0
     return
  else if (x .gt. (R+D)) then
     fc = 0.0d0
     dfc_dx = 0.0d0
     return
  else
     fc = 1.0d0 - (x-R+D)/(2.0d0*D) + 1.0d0/(2.0d0*pi)*dsin((pi/D)*(x-R+D))
     dfc_dx = 1.0d0/(2.0d0*D)* (dcos((pi/D)*(x-R+D)) - 1.0d0)
     return
  end if
end subroutine smooth_cutoff


end module functions_module
