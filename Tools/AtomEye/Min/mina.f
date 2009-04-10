      DOUBLE PRECISION FUNCTION MINAF77 (FE,N,NDIV,DEL,L,U,GUESS,X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXN=1024)
C     ------------------------------------------------------------------
C     MINA FINDS AN APPROXIMATE MINIMUM OF A REAL FUNCTION OF      
C     N VARIABLES, GIVEN AN INITIAL ESTIMATE OF THE POSITION OF   
C     THE MINIMUM AND RANGES FOR EACH OF THE VARIABLES.            
C     MINA USES A SELECTIVE DIRECTED SEARCH OF A SURROUNDING       
C     N-DIMENSIONAL GRID OF POINTS TO FIND A DIRECTION IN WHICH   
C     THE FUNCTION DECREASES.  IT THEN PROCEEDS IN THIS DIRECTION  
C     AS FAR AS THE FUNCTION DECREASES, THEN DETERMINES A NEW      
C     DIRECTION TO TRAVEL.  WHEN NO SUCH DIRECTION IS FOUND THE    
C     SEARCH INCREMENT FACTOR IS DECREASED AND THE PROCESS         
C     IS REPEATED.                                                 
C                                                                     
C     DESCRIPTION OF ARGUMENTS                                        
C     THE USER MUST DIMENSION ALL ARRAYS APPEARING IN THE CALL LIST
C     L(N), U(N), GUESS(N), X(N)                              
C     
C     INPUT--                                                      
C      FE  - NAME OF FUNCTION OF N VARIABLES TO BE MINIMIZED.     
C            (THIS NAME MUST APPEAR IN AN EXTERNAL STATEMENT.)     
C            FORM OF THE CALLING SEQUENCE MUST BE FUNCTION FE(X),  
C            WHERE X IS AN ARRAY OF N VARIABLE VALUES. THE        
C            ORDERING OF THE VARIABLES IS ARBITRARY, EXCEPT        
C            THAT IT MUST AGREE WITH THE ORDERING USED IN          
C            ARRAYS A AND GUESS.                                   
C     N    - NUMBER OF VARIABLES.  (N .GE. 1)                     
C     NDIV - NUMBER OF REFINEMENTS OF THE SEARCH INCREMENTS TO USE.
C            AT EACH REFINEMENT, THE INCREMENT IN EACH DIMENSION   
C            IS DIVIDED BY 10.  (USUALLY NDIV IS ABOUT 3 OR 4.)    
C     DEL  - FRACTION OF VARIABLE RANGE (IN EACH DIMENSION) TO USE 
C            AS THE INITIAL INCREMENT (IN THAT DIMENSION)
C     L,U  - ARRAYS OF SEARCH BOUND, DIMENSIONED N
C            L(I) SHOULD BE THE LOWER BOUND OF THE I-TH VARIABLE.
C            U(I) SHOULD BE THE UPPER BOUND OF THE I-TH VARIABLE.
C     GUESS - ARRAY OF N INITIAL VALUES.  GUESS(I) SHOULD BE THE   
C             INITIAL VALUE TO USE FOR THE I-TH VARIABLE.           
C     
C     OUTPUT--                                                     
C     X    - ARRAY (DIMENSIONED N) GIVING THE VALUES OF THE       
C            VARIABLES AT THE MINIMUM.  X(I) WILL BE THE VALUE     
C            OF THE I-TH VARIABLE.                                 
C     RETURNS FUNCTION VALUE AT THE MINIMUM                         
C     ------------------------------------------------------------------
      DOUBLE PRECISION L(MAXN),U(MAXN),GUESS(MAXN),X(MAXN)
      DOUBLE PRECISION XNOW(MAXN),XNEW(MAXN),R(MAXN)                         
      IF (N.LE.MAXN) GO TO 2                                         
      WRITE (6,*) 'error: minaf77: N is greater than MAXN = ', MAXN
      CALL EXIT(1)
    2 NX = N                                                         
      IDIV = 0                                                        
      DO 5 I=1,NX                                                     
      XNOW(I) = GUESS(I)                                              
      IF (XNOW(I).LT.L(I)) XNOW(I) = L(I)                         
      IF (XNOW(I).GT.U(I)) XNOW(I) = U(I)                         
      IF (L(I)-U(I)) 5,5,4                                        
    4 WRITE (6,*) 'error: minaf77: range min of x', I, ' =', L(I),
     A            ' > than max', U(I) 
      CALL EXIT(1)
    5 R(I) = U(I)-L(I)                                            
      DELTA = DEL                                                     
      IF (DELTA.LE.0.0) DELTA = 0.1                                   
      FEOW = FE(XNOW)                                                 
C     FIND NEW DIRECTION                                              
    7 DO 8 I=1,NX                                                     
    8 XNEW(I) = XNOW(I)                                               
      FOLD = FEOW                                                     
   10 DO 40 I=1,NX                                                    
      IF (XNOW(I).GE.U(I)) GO TO 20                                 
      XNEW(I) = MIN(XNOW(I)+DELTA*R(I),U(I))
      FEEW = FE(XNEW)                                                 
      IF (FEEW.LT.FEOW) GO TO 30                                      
   20 IF (XNOW(I) .LE. L(I)) GO TO 25                               
      XNEW(I) = MAX(XNOW(I)-DELTA*R(I),L(I))
      FEEW = FE(XNEW)                                                 
      IF (FEEW.LT.FEOW) GO TO 30                                      
   25 XNEW(I) = XNOW(I)                                               
      GO TO 40                                                        
   30 FEOW = FEEW                                                     
   40 CONTINUE                                                        
      ISTEP = 1                                                       
C     REFINE IF NEEDED                                                
      IF (FEOW.LT.FOLD) GO TO 50                                      
      IF (IDIV.GE.NDIV) GO TO 100                                     
      DELTA = DELTA*0.1                                               
      IDIV = IDIV+1                                                   
      GO TO 10                                                        
C     TRY TO CONTINUE IN CHOSEN DIRECTION                             
   50 ICHNG = 0                                                       
      FAC = 1.0                                                       
      IF ((ISTEP/10)*10.EQ.ISTEP) FAC = 2.0                           
      DO 60 I=1,NX                                                    
      DX = (XNEW(I)-XNOW(I))*FAC                                      
      XNOW(I) = XNEW(I)                                               
      IF (DX) 52,54,56                                                
   52 XNEW(I) = MAX(XNOW(I)+DX,L(I))
      IF (XNEW(I).LT.XNOW(I)) ICHNG = 1                               
      GO TO 60                                                        
   54 XNEW(I) = XNOW(I)                                               
      GO TO 60                                                        
   56 XNEW(I) = MIN(XNOW(I)+DX,U(I))
      IF (XNEW(I).GT.XNOW(I)) ICHNG = 1                               
   60 CONTINUE                                                        
      IF (ICHNG.EQ.0) GO TO 7                                         
      FEEW = FE(XNEW)                                                 
      IF (FEEW.GE.FEOW) GO TO 7                                       
      FEOW = FEEW                                                     
      ISTEP = ISTEP+1                                                 
      GO TO 50                                                        
C     RETURN ANSWERS                                                  
  100 FOFX = FOLD                                                     
      DO 110 I=1,NX                                                   
  110 X(I) = XNOW(I)                                                  
      MINAF77 = FOFX
      RETURN
      END                                                             
C     ----------------------- MINAF77 ends --------------------------------
      
      
      SUBROUTINE SIMIN (F,K,EPS,ANS,S,NEV,ICONT,Y)                    
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ------------------------------------------------------------------
C     DESCRIPTION OF PARAMETERS                                       
C      --INPUT--                                                      
C        F  - NAME OF FUNCTION OF K VARIABLES TO BE MINIMIZED.        
C             (THIS NAME MUST APPEAR IN AN EXTERNAL STATEMENT.)       
C             FORM OF THE CALLING SEQUENCE MUST BE FUNCTION F(X),     
C             WHERE X IS AN ARRAY OF K VARIABLES.                     
C        K  - THE NUMBER OF VARIABLES.  K MUST BE AT LEAST 2.         
C             NORMALLY K SHOULD BE LESS THAN ABOUT 10, AS SIMIN       
C             BECOMES LESS EFFECTIVE FOR LARGER VALUES OF K.          
C        EPS- THE CONVERGENCE CRITERION.  LET YAVG BE THE AVERAGE     
C             VALUE OF THE FUNCTION F AT THE K+1 POINTS OF THE        
C             SIMPLEX, AND LET R BE THEIR STANDARD ERROR.  (THAT IS,  
C             THE ROOT-MEAN-SQUARE OF THE SET OF VALUES (Y(I)-YAVG),  
C             WHERE Y(I) IS THE FUNCTION VALUE AT THE I-TH POINT OF   
C             THE SIMPLEX.)  THEN--                                   
C             IF EPS.GT.0, CONVERGENCE IS OBTAINED IF  R.LE.EPS.      
C             IF EPS.LT.0, CONVERGENCE IS IF  R.LE.ABS(EPS*YAVG).     
C             IF EPS=0, THE PROCESS WILL NOT CONVERGE BUT INSTEAD WILL
C             QUIT WHEN NEV FUNCTION EVALUATIONS HAVE BEEN USED.      
C        ANS- AN ARRAY OF LENGTH K CONTAINING A GUESS FOR THE LOCATION
C             OF A MINIMUM OF F.                                      
C        S  - A SCALE PARAMETER, WHICH MAY BE A SIMPLE VARIABLE OR AN 
C             ARRAY OF LENGTH K.  USE OF AN ARRAY IS SIGNALLED BY     
C             SETTING S(1) NEGATIVE.                                  
C             -SIMPLE VARIABLE CASE.  HERE S IS THE LENGTH OF EACH    
C             SIDE OF THE INITIAL SIMPLEX.  THUS, THE INITIAL SEARCH  
C             RANGE IS THE SAME FOR ALL THE VARIABLES.                
C             -ARRAY CASE.  HERE THE LENGTH OF SIDE I OF THE INITIAL  
C             SIMPLEX IS ABS(S(I)).  THUS, THE INITIAL SEARCH RANGE   
C             MAY BE DIFFERENT FOR DIFFERENT VARIABLES.               
C             NOTE-- THE VALUE(S) USED FOR S ARE NOT VERY CRITICAL.   
C             ANY REASONABLE GUESS SHOULD DO O.K.                     
C        NEV- THE MAXIMUM NUMBER OF FUNCTION EVALUATIONS TO BE USED.  
C             (THE ACTUAL NUMBER USED MAY EXCEED THIS SLIGHTLY SO THE 
C             LAST SEARCH ITERATION MAY BE COMPLETED.)                
C        ICONT - ICONT SHOULD BE ZERO ON ANY CALL TO SIMIN WHICH      
C             IS NOT A CONTINUATION OF A PREVIOUS CALL.               
C             IF ICONT=1 THE PROBLEM WILL BE CONTINUED.  IN THIS      
C             CASE THE WORK ARRAY Y MUST BE THE SAME ARRAY THAT WAS   
C             USED IN THE CALL THAT IS BEING CONTINUED (AND THE VALUES
C             IN IT MUST BE UNCHANGED).  THE REASON FOR THIS IS THAT  
C             IF ICONT=1 THEN THE ARGUMENT S IS IGNORED AND THE SIMPLE
C             AND RELATED FUNCTION VALUES THAT WERE STORED IN ARRAY Y 
C             DURING A PREVIOUS EXECUTION ARE USED TO CONTINUE THAT   
C             PREVIOUS PROBLEM.                                       
C        Y  - A WORK ARRAY CONTAINING AT LEAST K*K + 5*K + 1 WORDS.   
C             IF ICONT=1 THIS MUST BE THE SAME ARRAY USED IN THE CALL 
C             THAT IS BEING CONTINUED.                                
C      --OUTPUT--                                                     
C        ANS- ANS WILL CONTAIN THE LOCATION OF THE POINT WITH THE     
C             SMALLEST VALUE OF THE FUNCTION THAT WAS FOUND.          
C        S  - IN THE SIMPLE VARIABLE CASE S WILL BE RETURNED AS THE   
C             AVERAGE DISTANCE FROM THE VERTICES TO THE CENTROID OF   
C             THE SIMPLEX.                                            
C             IN THE ARRAY CASE S(I) WILL BE RETURNED AS THE AVERAGE  
C             DISTANCE IN THE I-TH DIMENSION OF VERTICES FROM         
C             THE CENTROID.  (S(1) WILL BE NEGATED.)                  
C             NOTE-- THE VALUE(S) RETURNED IN S ARE USEFUL FOR        
C             ASSESSING THE FLATNESS OF THE FUNCTION NEAR THE         
C             MINIMUM.  THE LARGER THE VALUE OF S (FOR A GIVEN        
C             VALUE OF EPS), THE FLATTER THE FUNCTION.                
C        NEV- NEV WILL BE THE COUNT OF THE ACTUAL NUMBER OF FUNCTION  
C             EVALUATIONS USED.                                       
C        Y  - WILL CONTAIN ALL DATA NEEDED TO CONTINUE THE MINIMIZATIO
C             SEARCH EFFICIENTLY IN A SUBSEQUENT CALL.                
C             NOTE -- THE FIRST K+1 ELEMENTS OF Y WILL CONTAIN THE    
C             FUNCTION VALUES AT THE K+1 POINTS OF THE LATEST SIMPLEX.
C             THE NEXT K*(K+1) ELEMENTS OF Y WILL BE THE K+1 POINTS   
C             OF THE SIMPLEX (IN EXACT CORRESPONDENSE TO THE ARRAY    
C             P DISCUSSED IN REFERENCE 1 ABOVE).  THE REMAINING 3*K   
C             WORDS ARE TEMPORARY WORKING STORAGE ONLY.               
C     ------------------------------------------------------------------
      DIMENSION ANS(K),S(K),Y(1)                                      
      EXTERNAL F                                                      
      IF(K.GE.2 .AND. S(1).NE.0.) GO TO 10                            
      WRITE (6,*) 'S(1)=0 OR K IS LESS THAN 2'
      RETURN                                                          
   10 WRITE (6,*) 'K IS LARGE THAN 100'
      IP = K+2                                                        
      IC = IP+K*(K+1)                                                 
      IR = IC+K                                                       
      IRR = IR+K                                                      
      CALL SIMINA(F,K,EPS,ANS,S,NEV,ICONT,Y,Y(IP),Y(IC),Y(IR),Y(IRR)) 
      RETURN                                                          
      END                                                             
C     ----------------------- SIMIN ends -------------------------------
      
      
      SUBROUTINE SIMINA (F,K,EPS,ANS,S,NEV,ICONT,Y,P,PC,PR,PRR)        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ------------------------------------------------------------------
C     SIMINA IS A CORE MINIMIZATION ROUTINE CALLED ONLY BY SIMIN.     
C     ------------------------------------------------------------------
      DIMENSION ANS(K),S(K),Y(3),P(K,3),PC(K),PR(K),PRR(K)            
      DATA ALPHA,BETA,GAMMA/1.0,0.5,2.0/                              
      KOUNT = 0                                                       
      KK = K                                                          
      IF (KK.LE.1) GO TO 99                                           
      ONEK = 1.0/FLOAT(KK)                                            
      KP1 = KK+1                                                      
      ONEKP1 = 1.0/FLOAT(KP1)                                         
      TOL = FLOAT(KP1)*EPS**2                                         
C     INITIAL SIMPLEX                                                 
      IF (ICONT.GE.1) GO TO 10                                        
      IF (S(1)) 4,99,1                                                
    1 SKP1 = S(1)*ONEKP1                                              
      DO 2 I=1,KP1                                                    
      DO 2 J=1,KK                                                     
    2 P(J,I) = ANS(J) - SKP1                                          
      DO 3 J=1,KK                                                     
    3 P(J,J+1) = P(J,J+1) + S(1)                                      
      GO TO 7                                                         
    4 DO 5 I=1,KP1                                                    
      DO 5 J=1,KK                                                     
    5 P(J,I) = ANS(J) - ABS(S(J))*ONEKP1                              
      DO 6 J=1,KK                                                     
    6 P(J,J+1) = P(J,J+1) + ABS(S(J))                                 
C     FUNCTION VALUES FOR INITIAL SIMPLEX                             
    7 I1 = 1                                                          
      DO 8 I=1,KP1                                                    
      Y(I) = F(P(1,I))                                                
      IF (Y(I).GT.Y(I1)) I1 = I                                       
    8 CONTINUE                                                        
      YANS = F(ANS)                                                   
      KOUNT = KP1+1                                                   
      IF (YANS.GE.Y(I1)) GO TO 10                                     
      Y(I1) = YANS                                                    
      DO 9 J=1,KK                                                     
    9 P(J,I1) = ANS(J)                                                
C     RE-START / NEXT ITERATION                                       
C     IF K.LT.0 VALUES IN THE P AND Y ARRAYS (AND ONLY THESE VALUES)  
C     WILL NOT HAVE BEEN DEFINED IN THIS CALL.  THIS IS NON-ANSI USAGE
C     FIRST FIND LARGEST, SECOND LARGEST, AND SMALLEST FUNCTION VALUES
   10 I1 = 1                                                          
      IL = 1                                                          
      DO 12 I=2,KP1                                                   
      IF (Y(I).LT.Y(IL)) IL = I                                       
      IF (Y(I).GT.Y(I1)) I1 = I                                       
   12 CONTINUE                                                        
      I2 = IL                                                         
      DO 13 I=1,KP1                                                   
      IF (I.EQ.I1) GO TO 13                                           
      IF (Y(I).GT.Y(I2)) I2 = I                                       
   13 CONTINUE                                                        
C     COMPUTE CENTROID, LEAVING OUT P(*,I1)                           
      DO 15 J=1,KK                                                    
      SUM = 0.0                                                       
      DO 14 I=1,KP1                                                   
      IF (I.EQ.I1) GO TO 14                                           
      SUM = SUM + P(J,I)                                              
   14 CONTINUE                                                        
   15 PC(J) = SUM*ONEK                                                
C     FORM REFLECTED POINT AND TEST                                   
      DO 20 J=1,KK                                                    
   20 PR(J) = PC(J) + ALPHA*(PC(J)-P(J,I1))                           
      YR = F(PR)                                                      
      KOUNT = KOUNT+1                                                 
      IF (YR.LT.Y(IL)) GO TO 30                                       
      IF (YR.GE.Y(I2)) GO TO 40                                       
C     ACCEPT REFLECTED POINT                                          
   21 Y(I1) = YR                                                      
      DO 22 J=1,KK                                                    
   22 P(J,I1) = PR(J)                                                 
      GO TO 60                                                        
C     EXPAND IN FAVORABLE DIRECTION AND TEST                          
   30 DO 31 J=1,KK                                                    
   31 PRR(J) = PR(J) + GAMMA*(PR(J)-PC(J))                            
      YRR = F(PRR)                                                    
      KOUNT = KOUNT+1                                                 
      IF (YRR.GE.YR) GO TO 21                                         
C     ACCEPT EXPANDED POINT                                           
      Y(I1) = YRR                                                     
      DO 32 J=1,KK                                                    
   32 P(J,I1) = PRR(J)                                                
      GO TO 60                                                        
C     DECIDE WHETHER TO ACCEPT REFLECTED POINT.                       
   40 IF (YR.GE.Y(I1)) GO TO 42                                       
      Y(I1) = YR                                                      
      DO 41 J=1,KK                                                    
   41 P(J,I1) = PR(J)                                                 
C     TRY CONTRACTION.                                                
   42 DO 43 J=1,KK                                                    
   43 PR(J) = PC(J) + BETA*(P(J,I1)-PC(J))                            
      YCT = F(PR)                                                     
      KOUNT = KOUNT+1                                                 
      IF (YCT.GT.Y(I1)) GO TO 50                                      
      Y(I1) = YCT                                                     
      DO 44 J=1,KK                                                    
   44 P(J,I1) = PR(J)                                                 
      GO TO 60                                                        
C     ALL EFFORTS FAILED.  SHRINK THE SIMPLEX ABOUT BEST POINT.       
   50 DO 52 I=1,KP1                                                   
      IF (I.EQ.IL) GO TO 52                                           
      DO 51 J=1,KK                                                    
   51 P(J,I) = 0.5*(P(J,I)+P(J,IL))                                   
      Y(I) = F(P(1,I))                                                
   52 CONTINUE                                                        
      KOUNT = KOUNT+KP1                                               
C     CHECK FOR CONVERGENCE                                           
   60 IF (KOUNT.GE.NEV) GO TO 65                                      
      IF (EPS.EQ.0.0) GO TO 10                                        
      SUM = 0.0                                                       
      DO 61 I=1,KP1                                                   
   61 SUM = SUM + Y(I)                                                
      YAVG = SUM*ONEKP1                                               
      SUM = 0.0                                                       
      DO 62 I=1,KP1                                                   
   62 SUM = SUM + (Y(I)-YAVG)**2                                      
      IF (EPS) 64,63,63                                               
   63 IF (SUM-TOL) 65,65,10                                           
   64 IF (SUM-TOL*ABS(YAVG)) 65,65,10                                 
C     CONVERGENCE OBTAINED.                                           
C     COMPUTE CENTROID                                                
   65 DO 68 J=1,KK                                                    
      SUM = 0.0                                                       
      DO 67 I=1,KP1                                                   
   67 SUM = SUM+P(J,I)                                                
   68 PC(J) = SUM*ONEKP1                                              
      IF (S(1)) 73,69,69                                              
C     COMPUTE S(1) AS AVERAGE DISTANCE OF VERTICES FROM CENTROID.     
   69 DIST = 0.0                                                      
      DO 71 I=1,KP1                                                   
      SUM = 0.0                                                       
      DO 70 J=1,KK                                                    
   70 SUM = SUM + (P(J,I)-PC(J))**2                                   
   71 DIST = DIST + SQRT(SUM)                                         
      S(1) = DIST*ONEKP1                                              
      GO TO 80                                                        
C     COMPUTE S(J) AS AVERAGE DISTANCE IN J-TH DIMENSION OF           
C     VERTICES FROM THE CENTROID.                                     
   73 DO 75 J=1,KK                                                    
      SUM = 0.0                                                       
      DO 74 I=1,KP1                                                   
   74 SUM = SUM + ABS(P(J,I)-PC(J))                                   
   75 S(J) = SUM*ONEKP1                                               
      S(1) = -S(1)                                                    
C     RETURN P(*,IL) AS ANSWER                                        
   80 IL = 1                                                          
      DO 82 I=2,KP1                                                   
      IF (Y(I).LT.Y(IL)) IL = I                                       
   82 CONTINUE                                                        
      DO 84 J=1,KK                                                    
   84 ANS(J) = P(J,IL)                                                
      NEV = KOUNT                                                     
      RETURN                                                          
C     ERROR MESSAGE                                                   
   99 WRITE (6,*) 'S(1)=0. OR K IS LESS THAN 2.'
      RETURN                                                          
      END                                                             
C     ---------------------- SIMINA ends -------------------------------
