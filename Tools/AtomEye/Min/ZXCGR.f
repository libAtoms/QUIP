C-----------------------------------------------------------------------ZXCG0030
C                                                                       ZXCG0040
C   COMPUTER            - HP9000/DOUBLE                                 ZXCG0050
C                                                                       ZXCG0060
C   LATEST REVISION     - NOVEMBER 1, 1979                              ZXCG0070
C                                                                       ZXCG0080
C   PURPOSE             - A CONJUGATE GRADIENT ALGORITHM FOR FINDING    ZXCG0090
C                           THE MINIMUM OF A FUNCTION OF N VARIABLES    ZXCG0100
C                                                                       ZXCG0110
C   USAGE               - CALL ZXCGR (FUNCT,N,ACC,MAXFN,DFPRED,X,G,F,W, ZXCG0120
C                           IER)                                        ZXCG0130
C                                                                       ZXCG0140
C   ARGUMENTS    FUNCT  - A USER SUPPLIED SUBROUTINE WHICH CALCULATES   ZXCG0150
C                           THE OBJECTIVE FUNCTION AND ITS GRADIENT     ZXCG0160
C                           FOR GIVEN PARAMETER VALUES                  ZXCG0170
C                           X(1),X(2),...,X(N).                         ZXCG0180
C                           THE CALLING SEQUENCE HAS THE FOLLOWING FORM ZXCG0190
C                           CALL FUNCT (N,X,F,G)                        ZXCG0200
C                           WHERE X AND G ARE VECTORS OF LENGTH N.      ZXCG0210
C                           THE SCALAR F IS FOR THE OBJECTIVE FUNCTION. ZXCG0220
C                           G(1), G(2), ..., G(N) ARE FOR THE COMPONENTSZXCG0230
C                           OF THE GRADIENT OF F.                       ZXCG0240
C                           FUNCT MUST APPEAR IN AN EXTERNAL STATEMENT  ZXCG0250
C                           IN THE CALLING PROGRAM. FUNCT MUST NOT      ZXCG0260
C                           ALTER THE VALUES OF X(I),I=1,...,N OR N.    ZXCG0270
C                N      - THE NUMBER OF PARAMETERS OF THE OBJECTIVE     ZXCG0280
C                           FUNCTION. (INPUT) (I.E.,THE LENGTH OF X)    ZXCG0290
C                ACC    - CONVERGENCE CRITERION. (INPUT)                ZXCG0300
C                           THE CALCULATION ENDS WHEN THE SUM OF SQUARESZXCG0310
C                           OF THE COMPONENTS OF G IS LESS THAN ACC.    ZXCG0320
C                MAXFN  - MAXIMUM NUMBER OF FUNCTION EVALUATIONS (I.E., ZXCG0330
C                           CALLS TO SUBROUTINE FUNCT) ALLOWED. (INPUT) ZXCG0340
C                           IF MAXFN IS SET TO ZERO, THEN THERE IS      ZXCG0350
C                           NO RESTRICTION ON THE NUMBER OF FUNCTION    ZXCG0360
C                           EVALUATIONS.                                ZXCG0370
C                DFPRED - A ROUGH ESTIMATE OF THE EXPECTED REDUCTION    ZXCG0380
C                           IN F, WHICH IS USED TO DETERMINE THE SIZE   ZXCG0390
C                           OF THE INITIAL CHANGE TO X. (INPUT)         ZXCG0400
C                           NOTE THAT DFPRED IS THE EXPECTED REDUCTION  ZXCG0410
C                           ITSELF, AND DOES NOT DEPEND ON ANY RATIOS.  ZXCG0420
C                           A BAD VALUE OF DFPRED CAUSES AN ERROR       ZXCG0430
C                           MESSAGE, WITH IER=129, AND A RETURN ON THE  ZXCG0440
C                           FIRST ITERATION. (SEE THE DESCRIPTION OF    ZXCG0450
C                           IER BELOW)                                  ZXCG0460
C                X      - VECTOR OF LENGTH N CONTAINING PARAMETER       ZXCG0470
C                           VALUES.                                     ZXCG0480
C                         ON INPUT, X MUST CONTAIN THE INITIAL          ZXCG0490
C                           PARAMETER ESTIMATES.                        ZXCG0500
C                         ON OUTPUT, X CONTAINS THE FINAL PARAMETER     ZXCG0510
C                           ESTIMATES AS DETERMINED BY ZXCGR.           ZXCG0520
C                G      - A VECTOR OF LENGTH N CONTAINING THE           ZXCG0530
C                           COMPONENTS OF THE GRADIENT OF F AT THE      ZXCG0540
C                           FINAL PARAMETER ESTIMATES. (OUTPUT)         ZXCG0550
C                F      - A SCALAR CONTAINING THE VALUE OF THE FUNCTION ZXCG0560
C                           AT THE FINAL PARAMETER ESTIMATES. (OUTPUT)  ZXCG0570
C                W      - WORK VECTOR OF LENGTH 6*N.                    ZXCG0580
C                IER    - ERROR PARAMETER. (OUTPUT)                     ZXCG0590
C                         IER = 0 IMPLIES THAT CONVERGENCE WAS          ZXCG0600
C                           ACHIEVED AND NO ERRORS OCCURRED.            ZXCG0610
C                         TERMINAL ERROR                                ZXCG0620
C                           IER = 129 IMPLIES THAT THE LINE SEARCH OF   ZXCG0630
C                             AN INTEGRATION WAS ABANDONED. THIS        ZXCG0640
C                             ERROR MAY BE CAUSED BY AN ERROR IN THE    ZXCG0650
C                             GRADIENT.                                 ZXCG0660
C                           IER = 130 IMPLIES THAT THE CALCULATION      ZXCG0670
C                             CANNOT CONTINUE BECAUSE THE SEARCH        ZXCG0680
C                             DIRECTION IS UPHILL.                      ZXCG0690
C                           IER = 131 IMPLIES THAT THE ITERATION WAS    ZXCG0700
C                             TERMINATED BECAUSE MAXFN WAS EXCEEDED.    ZXCG0710
C                           IER = 132 IMPLIES THAT THE CALCULATION      ZXCG0720
C                             WAS TERMINATED BECAUSE TWO CONSECUTIVE    ZXCG0730
C                             ITERATIONS FAILED TO REDUCE F.            ZXCG0740
C                                                                       ZXCG0750
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32                         ZXCG0760
C                       - SINGLE/H36,H48,H60                            ZXCG0770
C                                                                       ZXCG0780
C   REQD. IMSL ROUTINES - UERTST,UGETIO                                 ZXCG0790
C                                                                       ZXCG0800
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND           ZXCG0810
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL      ZXCG0820
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP  ZXCG0830
C                                                                       ZXCG0840
C   REMARKS  1.  THE ROUTINE INCLUDES NO THOROUGH CHECKS ON THE PART    ZXCG0850
C                OF THE USER PROGRAM THAT CALCULATES THE DERIVATIVES    ZXCG0860
C                OF THE OBJECTIVE FUNCTION. THEREFORE, BECAUSE          ZXCG0870
C                DERIVATIVE CALCULATION IS A FREQUENT SOURCE OF         ZXCG0880
C                ERROR, THE USER SHOULD VERIFY INDEPENDENTLY THE        ZXCG0890
C                CORRECTNESS OF THE DERIVATIVES THAT ARE GIVEN TO       ZXCG0900
C                THE ROUTINE.                                           ZXCG0910
C            2.  BECAUSE OF THE CLOSE RELATION BETWEEN THE CONJUGATE    ZXCG0920
C                GRADIENT METHOD AND THE METHOD OF STEEPEST DESCENTS,   ZXCG0930
C                IT IS VERY HELPFUL TO CHOOSE THE SCALE OF THE          ZXCG0940
C                VARIABLES IN A WAY THAT BALANCES THE MAGNITUDES OF     ZXCG0950
C                THE COMPONENTS OF A TYPICAL DERIVATE VECTOR. IT        ZXCG0960
C                CAN BE PARTICULARLY INEFFICIENT IF A FEW COMPONENTS    ZXCG0970
C                OF THE GRADIENT ARE MUCH LARGER THAN THE REST.         ZXCG0980
C            3.  IF THE VALUE OF THE PARAMETER ACC IN THE ARGUMENT      ZXCG0990
C                LIST OF THE ROUTINE IS SET TO ZERO, THEN THE           ZXCG1000
C                SUBROUTINE WILL CONTINUE ITS CALCULATION UNTIL IT      ZXCG1010
C                STOPS REDUCING THE OBJECTIVE FUNCTION. IN THIS CASE    ZXCG1020
C                THE USUAL BEHAVIOUR IS THAT CHANGES IN THE             ZXCG1030
C                OBJECTIVE FUNCTION BECOME DOMINATED BY COMPUTER        ZXCG1040
C                ROUNDING ERRORS BEFORE PRECISION IS LOST IN THE        ZXCG1050
C                GRADIENT VECTOR. THEREFORE, BECAUSE THE POINT OF       ZXCG1060
C                VIEW HAS BEEN TAKEN THAT THE USER REQUIRES THE         ZXCG1070
C                LEAST POSSIBLE VALUE OF THE FUNCTION, A VALUE OF       ZXCG1080
C                THE OBJECTIVE FUNCTION THAT IS SMALL DUE TO            ZXCG1090
C                COMPUTER ROUNDING ERRORS CAN PREVENT FURTHER           ZXCG1100
C                PROGRESS. HENCE THE PRECISION IN THE FINAL VALUES      ZXCG1110
C                OF THE VARIABLES MAY BE ONLY ABOUT HALF THE            ZXCG1120
C                NUMBER OF SIGNIFICANT DIGITS IN THE COMPUTER           ZXCG1130
C                ARITHMETIC, BUT THE LEAST VALUE OF F IS USUALLY        ZXCG1140
C                FOUND TO QUITE HIGH ACCURACY.                          ZXCG1150
C                                                                       ZXCG1160
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.       ZXCG1170
C                                                                       ZXCG1180
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN ZXCG1190
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,    ZXCG1200
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.        ZXCG1210
C                                                                       ZXCG1220
C-----------------------------------------------------------------------ZXCG1230
C                                                                       ZXCG1240
      SUBROUTINE ZXCGRF77  (FUNCT,N,ACC,MAXFN,DFPRED,X,G,F,W,IER)          ZXCG1250
      external FUNCT
C                                  SPECIFICATIONS FOR ARGUMENTS         ZXCG1260
      INTEGER            N,MAXFN,IER                                    ZXCG1270
      DOUBLE PRECISION   ACC,DFPRED,X(N),G(N),F,W(1)                    ZXCG1280
C                                  SPECIFICATIONS FOR LOCAL VARIABLES   ZXCG1290
      INTEGER            MAXLIN,MXFCON,I,IGINIT,IGOPT,IRETRY,IRSDG,     ZXCG1300
     1                   IRSDX,ITERC,ITERFM,ITERRS,IXOPT,NCALLS,NFBEG,  ZXCG1310
     2                   NFOPT                                          ZXCG1320
      DOUBLE PRECISION   BETA,DDSPLN,DFPR,FCH,FINIT,FMIN,GAMDEN,GAMA,   ZXCG1330
     1                   GINIT,GMIN,GNEW,GSPLN,GSQRD,SBOUND,STEP,STEPCH,ZXCG1340
     2                   STMIN,SUM,WORK                                 ZXCG1350
      DATA               MAXLIN/5/,MXFCON/2/                            ZXCG1360
C                                  FIRST EXECUTABLE STATEMENT           ZXCG1370
      IER = 0                                                           ZXCG1380
C                                  THE WORKING SPACE ARRAY IS SPLIT     ZXCG1390
C                                    INTO SIX VECTORS OF LENGTH N. THE  ZXCG1400
C                                    FIRST PART IS USED FOR THE SEARCH  ZXCG1410
C                                    DIRECTION OF AN ITERATION. THE     ZXCG1420
C                                    SECOND AND THIRD PARTS CONTAIN THE ZXCG1430
C                                    INFORMATION THAT IS REQUIRED BY    ZXCG1440
C                                    THE CONJUGACY CONDITIONS OF THE    ZXCG1450
C                                    RESTART PROCEDURE. THE FOURTH PART ZXCG1460
C                                    CONTAINS THE GRADIENT AT THE START ZXCG1470
C                                    OF AN ITERATION. THE FIFTH PART    ZXCG1480
C                                    CONTAINS THE PARAMETERS THAT GIVE  ZXCG1490
C                                    THE LEAST CALCULATED VALUE OF F.   ZXCG1500
C                                    THE SIXTH PART CONTAINS THE        ZXCG1510
C                                    GRADIENT VECTOR WHERE F IS LEAST.  ZXCG1520
      IRSDX = N                                                         ZXCG1530
      IRSDG = IRSDX+N                                                   ZXCG1540
      IGINIT = IRSDG+N                                                  ZXCG1550
      IXOPT = IGINIT+N                                                  ZXCG1560
      IGOPT = IXOPT+N                                                   ZXCG1570
C                                  SET SOME PARAMETERS TO BEGIN THE     ZXCG1580
C                                    CALCULATION. ITERC AND             ZXCG1590
C                                    NCALLS COUNT THE NUMBER OF         ZXCG1600
C                                    ITERATIONS AND CALLS OF FUNCT.     ZXCG1610
C                                    ITERFM IS THE NUMBER OF THE MOST   ZXCG1620
C                                    RECENT ITERATION THAT DECREASES F. ZXCG1630
      ITERC = 0                                                         ZXCG1640
      NCALLS = 0                                                        ZXCG1650
      ITERFM = ITERC                                                    ZXCG1660
C                                  CALL SUBROUTINE FUNCT. LET THE       ZXCG1670
C                                    INITIAL SEARCH DIRECTION BE MINUS  ZXCG1680
C                                    THE GRADIENT VECTOR. USUALLY THE   ZXCG1690
C                                    PARAMETER ITERRS GIVES THE         ZXCG1700
C                                    ITERATION NUMBER OF THE MOST       ZXCG1710
C                                    RECENT RESTART, BUT IT IS SET TO   ZXCG1720
C                                    ZERO WHEN THE STEEPEST DESCENT     ZXCG1730
C                                    DIRECTION IS USED.                 ZXCG1740
    5 NCALLS = NCALLS+1                                                 ZXCG1750
      CALL FUNCT(N,X,F,G)                                                    ZXCG1760
      IF (NCALLS.GE.2) GO TO 20                                         ZXCG1770
   10 DO 15 I=1,N                                                       ZXCG1780
   15 W(I) = -G(I)                                                      ZXCG1790
      ITERRS = 0                                                        ZXCG1800
      IF (ITERC.GT.0) GO TO 80                                          ZXCG1810
C                                  SET SUM TO G SQUARED. GMIN AND GNEW  ZXCG1820
C                                    ARE THE OLD AND THE NEW            ZXCG1830
C                                    DIRECTIONAL DERIVATIVES ALONG THE  ZXCG1840
C                                    CURRENT SEARCH DIRECTION. LET FCH  ZXCG1850
C                                    BE THE DIFFERENCE BETWEEN F AND    ZXCG1860
C                                    THE PREVIOUS BEST VALUE OF THE     ZXCG1870
C                                    OBJECTIVE FUNCTION.                ZXCG1880
   20 GNEW = 0.0D0                                                      ZXCG1890
      SUM = 0.0D0                                                       ZXCG1900
      DO 25 I=1,N                                                       ZXCG1910
         GNEW = GNEW+W(I)*G(I)                                          ZXCG1920
   25 SUM = SUM+G(I)**2                                                 ZXCG1930
      IF (NCALLS.EQ.1) GO TO 35                                         ZXCG1940
      FCH = F-FMIN                                                      ZXCG1950
C                                  STORE THE VALUES OF X, F AND G, IF   ZXCG1960
C                                    THEY ARE THE BEST THAT HAVE BEEN   ZXCG1970
C                                    CALCULATED SO FAR, AND NOTE G      ZXCG1980
C                                    SQUARED AND THE VALUE OF NCALLS.   ZXCG1990
C                                    TEST FOR CONVERGENCE.              ZXCG2000
      IF (FCH) 35,30,50                                                 ZXCG2010
   30 IF (GNEW/GMIN.LT.-1.0D0) GO TO 45                                 ZXCG2020
   35 FMIN = F                                                          ZXCG2030
      GSQRD = SUM                                                       ZXCG2040
      NFOPT = NCALLS                                                    ZXCG2050
      DO 40 I=1,N                                                       ZXCG2060
         W(IXOPT+I) = X(I)                                              ZXCG2070
   40 W(IGOPT+I) = G(I)                                                 ZXCG2080
   45 IF (SUM.LE.ACC) GO TO 9005                                        ZXCG2090
C                                  TEST IF THE VALUE OF MAXFN ALLOWS    ZXCG2100
C                                    ANOTHER CALL OF FUNCT.             ZXCG2110
   50 IF (NCALLS.NE.MAXFN) GO TO 55                                     ZXCG2120
      IER = 131                                                         ZXCG2130
      GO TO 9000                                                        ZXCG2140
   55 IF (NCALLS.GT.1) GO TO 100                                        ZXCG2150
C                                  SET DFPR TO THE ESTIMATE OF THE      ZXCG2160
C                                    REDUCTION IN F GIVEN IN THE        ZXCG2170
C                                    ARGUMENT LIST, IN ORDER THAT THE   ZXCG2180
C                                    INITIAL CHANGE TO THE PARAMETERS   ZXCG2190
C                                    IS OF A SUITABLE SIZE. THE VALUE   ZXCG2200
C                                    OF STMIN IS USUALLY THE            ZXCG2210
C                                    STEP-LENGTH OF THE MOST RECENT     ZXCG2220
C                                    LINE SEARCH THAT GIVES THE LEAST   ZXCG2230
C                                    CALCULATED VALUE OF F.             ZXCG2240
      DFPR = DFPRED                                                     ZXCG2250
      STMIN = DFPRED/GSQRD                                              ZXCG2260
C                                  BEGIN THE ITERATION                  ZXCG2270
   80 ITERC = ITERC+1                                                   ZXCG2280
C                                  STORE THE INITIAL FUNCTION VALUE AND ZXCG2290
C                                    GRADIENT, CALCULATE THE INITIAL    ZXCG2300
C                                    DIRECTIONAL DERIVATIVE, AND BRANCH ZXCG2310
C                                    IF ITS VALUE IS NOT NEGATIVE. SET  ZXCG2320
C                                    SBOUND TO MINUS ONE TO INDICATE    ZXCG2330
C                                    THAT A BOUND ON THE STEP IS NOT    ZXCG2340
C                                    KNOWN YET, AND SET NFBEG TO THE    ZXCG2350
C                                    CURRENT VALUE OF NCALLS. THE       ZXCG2360
C                                    PARAMETER IRETRY SHOWS THE NUMBER  ZXCG2370
C                                    OF ATTEMPTS AT SATISFYING THE BETA ZXCG2380
C                                    CONDITION.                         ZXCG2390
      FINIT = F                                                         ZXCG2400
      GINIT = 0.0D0                                                     ZXCG2410
      DO 85 I=1,N                                                       ZXCG2420
         W(IGINIT+I) = G(I)                                             ZXCG2430
   85 GINIT = GINIT+W(I)*G(I)                                           ZXCG2440
      IF (GINIT.GE.0.0D0) GO TO 165                                     ZXCG2450
      GMIN = GINIT                                                      ZXCG2460
      SBOUND = -1.0D0                                                   ZXCG2470
      NFBEG = NCALLS                                                    ZXCG2480
      IRETRY = -1                                                       ZXCG2490
C                                  SET STEPCH SO THAT THE INITIAL       ZXCG2500
C                                    STEP-LENGTH IS CONSISTENT WITH THE ZXCG2510
C                                    PREDICTED REDUCTION IN F, SUBJECT  ZXCG2520
C                                    TO THE CONDITION THAT IT DOES NOT  ZXCG2530
C                                    EXCEED THE STEP-LENGTH OF THE      ZXCG2540
C                                    PREVIOUS ITERATION. LET STMIN BE   ZXCG2550
C                                    THE STEP TO THE LEAST CALCULATED   ZXCG2560
C                                    VALUE OF F.                        ZXCG2570
      STEPCH = DMIN1(STMIN,DABS(DFPR/GINIT))                            ZXCG2580
      STMIN = 0.0D0                                                     ZXCG2590
C                                  CALL SUBROUTINE FUNCT AT THE VALUE   ZXCG2600
C                                    OF X THAT IS DEFINED BY THE NEW    ZXCG2610
C                                    CHANGE TO THE STEP-LENGTH, AND LET ZXCG2620
C                                    THE NEW STEP-LENGTH BE STEP. THE   ZXCG2630
C                                    VARIABLE WORK IS USED AS WORK      ZXCG2640
C                                    SPACE.                             ZXCG2650
   90 STEP = STMIN+STEPCH                                               ZXCG2660
      WORK = 0.0D0                                                      ZXCG2670
      DO 95 I=1,N                                                       ZXCG2680
         X(I) = W(IXOPT+I)+STEPCH*W(I)                                  ZXCG2690
   95 WORK = DMAX1(WORK,DABS(X(I)-W(IXOPT+I)))                          ZXCG2700
      IF (WORK.GT.0.0D0) GO TO 5                                        ZXCG2710
C                                  TERMINATE THE LINE SEARCH IF STEPCH  ZXCG2720
C                                    IS EFFECTIVELY ZERO.               ZXCG2730
      IF (NCALLS.GT.NFBEG+1) GO TO 115                                  ZXCG2740
      IF (DABS(GMIN/GINIT)-0.2D0) 170,170,115                           ZXCG2750
C                                  LET SPLN BE THE QUADRATIC SPLINE     ZXCG2760
C                                    THAT INTERPOLATES THE CALCULATED   ZXCG2770
C                                    FUNCTION VALUES AND DIRECTIONAL    ZXCG2780
C                                    DERIVATIVES AT THE POINTS STMIN    ZXCG2790
C                                    AND STEP OF THE LINE SEARCH, WHERE ZXCG2800
C                                    THE KNOT OF THE SPLINE IS AT       ZXCG2810
C                                    0.5*(STMIN+STEP). REVISE STMIN,    ZXCG2820
C                                    GMIN AND SBOUND, AND SET DDSPLN TO ZXCG2830
C                                    THE SECOND DERIVATIVE OF SPLN AT   ZXCG2840
C                                    THE NEW STMIN. HOWEVER, IF FCH IS  ZXCG2850
C                                    ZERO, IT IS ASSUMED THAT THE       ZXCG2860
C                                    MAXIMUM ACCURACY IS ALMOST         ZXCG2870
C                                    ACHIEVED, SO DDSPLN IS CALCULATED  ZXCG2880
C                                    USING ONLY THE CHANGE IN THE       ZXCG2890
C                                    GRADIENT.                          ZXCG2900
  100 WORK = (FCH+FCH)/STEPCH-GNEW-GMIN                                 ZXCG2910
      DDSPLN = (GNEW-GMIN)/STEPCH                                       ZXCG2920
      IF (NCALLS.GT.NFOPT) SBOUND = STEP                                ZXCG2930
      IF (NCALLS.GT.NFOPT) GO TO 105                                    ZXCG2940
      IF (GMIN*GNEW.LE.0.0D0) SBOUND = STMIN                            ZXCG2950
      STMIN = STEP                                                      ZXCG2960
      GMIN = GNEW                                                       ZXCG2970
      STEPCH = -STEPCH                                                  ZXCG2980
  105 IF (FCH.NE.0.0D0) DDSPLN = DDSPLN+(WORK+WORK)/STEPCH              ZXCG2990
C                                                                       ZXCG3000
C                                  TEST FOR CONVERGENCE OF THE LINE     ZXCG3010
C                                    SEARCH, BUT FORCE AT LEAST TWO     ZXCG3020
C                                    STEPS TO BE TAKEN IN ORDER NOT TO  ZXCG3030
C                                    LOSE QUADRATIC TERMINATION.        ZXCG3040
      IF (GMIN.EQ.0.0D0) GO TO 170                                      ZXCG3050
      IF (NCALLS.LE.NFBEG+1) GO TO 120                                  ZXCG3060
      IF (DABS(GMIN/GINIT).LE.0.2D0) GO TO 170                          ZXCG3070
C                                  APPLY THE TEST THAT DEPENDS ON THE   ZXCG3080
C                                    PARAMETER MAXLIN.                  ZXCG3090
  110 IF (NCALLS.LT.NFOPT+MAXLIN) GO TO 120                             ZXCG3100
  115 IER = 129                                                         ZXCG3110
      GO TO 170                                                         ZXCG3120
C                                  SET STEPCH TO THE GREATEST CHANGE TO ZXCG3130
C                                    THE CURRENT VALUE OF STMIN THAT IS ZXCG3140
C                                    ALLOWED BY THE BOUND ON THE LINE   ZXCG3150
C                                    SEARCH. SET GSPLN TO THE GRADIENT  ZXCG3160
C                                    OF THE QUADRATIC SPLINE AT         ZXCG3170
C                                    (STMIN+STEPCH). HENCE CALCULATE    ZXCG3180
C                                    THE VALUE OF STEPCH THAT MINIMIZES ZXCG3190
C                                    THE SPLINE FUNCTION, AND THEN      ZXCG3200
C                                    OBTAIN THE NEW FUNCTION AND        ZXCG3210
C                                    GRADIENT VECTOR, FOR THIS VALUE OF ZXCG3220
C                                    THE CHANGE TO THE STEP-LENGTH.     ZXCG3230
  120 STEPCH = 0.5D0*(SBOUND-STMIN)                                     ZXCG3240
      IF (SBOUND.LT.-0.5D0) STEPCH = 9.0D0*STMIN                        ZXCG3250
      GSPLN = GMIN+STEPCH*DDSPLN                                        ZXCG3260
      IF (GMIN*GSPLN.LT.0.0D0) STEPCH = STEPCH*GMIN/(GMIN-GSPLN)        ZXCG3270
      GO TO 90                                                          ZXCG3280
C                                  CALCULATE THE VALUE OF BETA THAT     ZXCG3290
C                                    OCCURS IN THE NEW SEARCH           ZXCG3300
C                                    DIRECTION.                         ZXCG3310
  125 SUM = 0.0D0                                                       ZXCG3320
      DO 130 I=1,N                                                      ZXCG3330
  130 SUM = SUM+G(I)*W(IGINIT+I)                                        ZXCG3340
      BETA = (GSQRD-SUM)/(GMIN-GINIT)                                   ZXCG3350
C                                  TEST THAT THE NEW SEARCH DIRECTION   ZXCG3360
C                                    CAN BE MADE DOWNHILL. IF IT        ZXCG3370
C                                    CANNOT, THEN MAKE ONE ATTEMPT TO   ZXCG3380
C                                    IMPROVE THE ACCURACY OF THE LINE   ZXCG3390
C                                    SEARCH.                            ZXCG3400
      IF (DABS(BETA*GMIN).LE.0.2D0*GSQRD) GO TO 135                     ZXCG3410
      IRETRY = IRETRY+1                                                 ZXCG3420
      IF (IRETRY.LE.0) GO TO 110                                        ZXCG3430
C                                  APPLY THE TEST THAT DEPENDS ON THE   ZXCG3440
C                                    PARAMETER MXFCON.                  ZXCG3450
C                                    SET DFPR TO THE PREDICTED          ZXCG3460
C                                    REDUCTION IN F ON THE NEXT         ZXCG3470
C                                    ITERATION.                         ZXCG3480
  135 IF (F.LT.FINIT) ITERFM = ITERC                                    ZXCG3490
      IF (ITERC.LT.ITERFM+MXFCON) GO TO 140                             ZXCG3500
      IER = 132                                                         ZXCG3510
      GO TO 9000                                                        ZXCG3520
  140 DFPR = STMIN*GINIT                                                ZXCG3530
C                                  BRANCH IF A RESTART PROCEDURE IS     ZXCG3540
C                                    REQUIRED DUE TO THE ITERATION      ZXCG3550
C                                    NUMBER OR DUE TO THE SCALAR        ZXCG3560
C                                    PRODUCT OF CONSECUTIVE GRADIENTS.  ZXCG3570
      IF (IRETRY.GT.0) GO TO 10                                         ZXCG3580
      IF (ITERRS.EQ.0) GO TO 155                                        ZXCG3590
      IF (ITERC-ITERRS.GE.N) GO TO 155                                  ZXCG3600
      IF (DABS(SUM).GE.0.2D0*GSQRD) GO TO 155                           ZXCG3610
C                                  CALCULATE THE VALUE OF GAMA THAT     ZXCG3620
C                                    OCCURS IN THE NEW SEARCH           ZXCG3630
C                                    DIRECTION, AND SET SUM TO A SCALAR ZXCG3640
C                                    PRODUCT FOR THE TEST BELOW. THE    ZXCG3650
C                                    VALUE OF GAMDEN IS SET BY THE      ZXCG3660
C                                    RESTART PROCEDURE.                 ZXCG3670
      GAMA = 0.0D0                                                      ZXCG3680
      SUM = 0.0D0                                                       ZXCG3690
      DO 145 I=1,N                                                      ZXCG3700
         GAMA = GAMA+G(I)*W(IRSDG+I)                                    ZXCG3710
  145 SUM = SUM+G(I)*W(IRSDX+I)                                         ZXCG3720
      GAMA = GAMA/GAMDEN                                                ZXCG3730
C                                  RESTART IF THE NEW SEARCH DIRECTION  ZXCG3740
C                                    IS NOT SUFFICIENTLY DOWNHILL.      ZXCG3750
C                                                                       ZXCG3760
      IF (DABS(BETA*GMIN+GAMA*SUM).GE.0.2D0*GSQRD) GO TO 155            ZXCG3770
C                                                                       ZXCG3780
C                                  CALCULATE THE NEW SEARCH DIRECTION.  ZXCG3790
      DO 150 I=1,N                                                      ZXCG3800
  150 W(I) = -G(I)+BETA*W(I)+GAMA*W(IRSDX+I)                            ZXCG3810
      GO TO 80                                                          ZXCG3820
C                                  APPLY THE RESTART PROCEDURE.         ZXCG3830
  155 GAMDEN = GMIN-GINIT                                               ZXCG3840
      DO 160 I=1,N                                                      ZXCG3850
         W(IRSDX+I) = W(I)                                              ZXCG3860
         W(IRSDG+I) = G(I)-W(IGINIT+I)                                  ZXCG3870
  160 W(I) = -G(I)+BETA*W(I)                                            ZXCG3880
      ITERRS = ITERC                                                    ZXCG3890
      GO TO 80                                                          ZXCG3900
C                                  SET IER TO INDICATE THAT THE SEARCH  ZXCG3910
C                                    DIRECTION IS UPHILL.               ZXCG3920
  165 IER = 130                                                         ZXCG3930
C                                  ENSURE THAT F, X AND G ARE OPTIMAL.  ZXCG3940
  170 IF (NCALLS.EQ.NFOPT) GO TO 180                                    ZXCG3950
      F = FMIN                                                          ZXCG3960
      DO 175 I=1,N                                                      ZXCG3970
         X(I) = W(IXOPT+I)                                              ZXCG3980
  175 G(I) = W(IGOPT+I)                                                 ZXCG3990
  180 IF (IER.EQ.0) GO TO 125                                           ZXCG4000
 9000 CONTINUE                                                          ZXCG4010
c      CALL UERTST (IER,'ZXCGRF77 ')                                    ZXCG4020
 9005 RETURN                                                            ZXCG4030
      END                                                               ZXCG4040
