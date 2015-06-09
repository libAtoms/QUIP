!X   Portions of this code were written by Mark Pederson as part of
!X   the NRLMOL program, developed as part of his employment for 
!X   the U.S. Government, and are not subject to copyright in the USA.

      ! compute <psi_i | x,y,z,1 | psi_j>
      ! where <psi_i | r> is a Gaussian times a spherical harmonic
      ! returns 10x10 matrix for different spherical harmonics
      ! 1 - s
      ! 2-4 px, py, pz
      ! 5-10 xx yy zz xy xz yz
      SUBROUTINE GINTED(A1,A2,A,B,W)
      IMPLICIT  REAL*8 (A-H,O-Z)
      DIMENSION A(3),B(3),W(10,10,4),L3(3,4)
      DIMENSION V(3,3,5,3)
      DIMENSION ML(10,3)
      SAVE
      DATA  L3/2,1,1,  1,2,1,  1,1,2,  1,1,1/
      DATA ML/1,2,1,1,3,1,1,2,2,1,
     1        1,1,2,1,1,3,1,2,1,2,
     2        1,1,1,2,1,1,3,1,2,2/
      AT=1.0D0/(A1+A2)
      Q1=0.5D0*AT
      Q2=Q1*Q1
      XS=0.0D0
      DO 10 IFT=1,3
      XS=XS+(A1*A2*((A(IFT)-B(IFT))**2))*AT
      D=(A1*A(IFT)+A2*B(IFT))*AT
      F1=D-A(IFT)
      F2=D-B(IFT)
      F3=D
      F001=F3
      F010=F2
      F011=F010*F3
      F020=F010*F2
      F021=F020*F3
      F100=F1
      F101=F100*F3
      F110=F100*F2
      F111=F110*F3
      F120=F110*F2
      F121=F120*F3
      F200=F100*F1
      F201=F200*F3
      F210=F200*F2
      F211=F210*F3
      F220=F210*F2
      F221=F220*F3
      V(1,1,1,IFT)=1
      V(1,2,1,IFT)=+F010
      V(1,3,1,IFT)=+F020+Q1
      V(2,1,1,IFT)=+F100
      V(2,2,1,IFT)=+F110+Q1
      V(2,3,1,IFT)=+F120+Q1*F100+2*Q1*F010
      V(3,1,1,IFT)=+F200+Q1
      V(3,2,1,IFT)=+F210+2*Q1*F100+Q1*F010
      V(3,3,1,IFT)=+F220+Q1*F200+4*Q1*F110+Q1*F020+3*Q2
      V(1,1,2,IFT)=+F001
      V(1,2,2,IFT)=+F011+Q1
      V(1,3,2,IFT)=+F021+2*Q1*F010+Q1*F001
      V(2,1,2,IFT)=+F101+Q1
      V(2,2,2,IFT)=+F111+Q1*F100+Q1*F010+Q1*F001
      V(2,3,2,IFT)=+F121+2*Q1*F110+Q1*F101+Q1*F020+2*Q1*F011+3*Q2
      V(3,1,2,IFT)=+F201+2*Q1*F100+Q1*F001
      V(3,2,2,IFT)=+F211+Q1*F200+2*Q1*F110+2*Q1*F101+Q1*F011+3*Q2
      V(3,3,2,IFT)=+F221+2*Q1*F210+Q1*F201+2*Q1*F120+4*Q1*F111
     &             +6*Q2*F100+Q1*F021+6*Q2*F010+3*Q2*F001
 10   CONTINUE
      D=EXP(-XS)*(SQRT(3.14159265358979324D0*AT)**3)
      DO 40 K=1,4 
       DO 35 J=1,10
        DO 30 I=1,10
         W(I,J,K)=V(ML(I,1),ML(J,1),L3(1,K),1)*
     &            V(ML(I,2),ML(J,2),L3(2,K),2)*
     &            V(ML(I,3),ML(J,3),L3(3,K),3)
   30   CONTINUE
   35  CONTINUE
   40 CONTINUE
      DO 50 K=1,4
       DO 48 J=1,10
        DO 46 I=1,10
         W(I,J,K)=W(I,J,K)*D
  46    CONTINUE
  48   CONTINUE
  50  CONTINUE
      END
