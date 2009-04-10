/***************************************/
/* libMin: Minimization Library        */
/*         -lm -lScalar -lVecMat       */
/*                                     */
/* Dec. 6, 1999 Ju Li <liju99@mit.edu> */
/***************************************/

#ifndef _Min_H
#define _Min_H

#include <stdio.h>
#include <Scalar.h>
#include <VecMat.h>

/* Minimization of single-variable f(x) on interval [a,b]: do parabolic */
/* interpolation when applicable; returns min; xmin is stored is *xmin  */
/* which is guaranteed to be within tol from the real minimum on [a,b]. */
double onemin (double (*f)(double), double xa, double xb, double tol,
               double *xmin);

double FORTRAN_SYMBOL(minaf77)
    (double (*fe)(double []), int *n, int *ndiv, double *del,
     double L[], double U[], double xguess[], double x[]);

/*
  MINA FINDS AN APPROXIMATE MINIMUM OF A REAL FUNCTION OF      
  N VARIABLES, GIVEN AN INITIAL ESTIMATE OF THE POSITION OF   
  THE MINIMUM AND RANGES FOR EACH OF THE VARIABLES.            
  MINA USES A SELECTIVE DIRECTED SEARCH OF A SURROUNDING       
  N-DIMENSIONAL GRID OF POINTS TO FIND A DIRECTION IN WHICH   
  THE FUNCTION DECREASES.  IT THEN PROCEEDS IN THIS DIRECTION  
  AS FAR AS THE FUNCTION DECREASES, THEN DETERMINES A NEW      
  DIRECTION TO TRAVEL.  WHEN NO SUCH DIRECTION IS FOUND THE    
  SEARCH INCREMENT FACTOR IS DECREASED AND THE PROCESS         
  IS REPEATED.                                                 

  DESCRIPTION OF ARGUMENTS                                        
  THE USER MUST DIMENSION ALL ARRAYS APPEARING IN THE CALL LIST
  L(N), U(N), GUESS(N), X(N)                              

  INPUT--                                                      
  FE  - NAME OF FUNCTION OF N VARIABLES TO BE MINIMIZED.     
        (THIS NAME MUST APPEAR IN AN EXTERNAL STATEMENT.)     
	FORM OF THE CALLING SEQUENCE MUST BE FUNCTION FE(X),  
	WHERE X IS AN ARRAY OF N VARIABLE VALUES. THE        
	ORDERING OF THE VARIABLES IS ARBITRARY, EXCEPT        
	THAT IT MUST AGREE WITH THE ORDERING USED IN          
	ARRAYS A AND GUESS.                                   
 N    - NUMBER OF VARIABLES.  (N .GE. 1)                     
 NDIV - NUMBER OF REFINEMENTS OF THE SEARCH INCREMENTS TO USE.
        AT EACH REFINEMENT, THE INCREMENT IN EACH DIMENSION   
	IS DIVIDED BY 10.  (USUALLY NDIV IS ABOUT 3 OR 4.)    
 DEL  - FRACTION OF VARIABLE RANGE (IN EACH DIMENSION) TO USE 
        AS THE INITIAL INCREMENT (IN THAT DIMENSION)
 L,U  - ARRAYS OF SEARCH BOUND, DIMENSIONED N
        L(I) SHOULD BE THE LOWER BOUND OF THE I-TH VARIABLE.
        U(I) SHOULD BE THE UPPER BOUND OF THE I-TH VARIABLE.
 GUESS - ARRAY OF N INITIAL VALUES.  GUESS(I) SHOULD BE THE   
        INITIAL VALUE TO USE FOR THE I-TH VARIABLE.
       
  OUTPUT--                                                     
  X   - ARRAY (DIMENSIONED N) GIVING THE VALUES OF THE       
        VARIABLES AT THE MINIMUM.  X(I) WILL BE THE VALUE     
        OF THE I-TH VARIABLE.                                 
  RETURNS FUNCTION VALUE AT THE MINIMUM                         
*/

double FORTRAN_SYMBOL(simannf77)
    (double (*fe)(double []), int *n,
     double *T, double *RT, int *ISEED,
     double L[], double U[], double xguess[], double x[]);

/* cg.c: */

/*********************************/
/* Conjugate Gradient Minimizers */
/*********************************/

/****************************************/
/* Notation:                            */
/* eV to denote unit of potential value */
/* A  to denote unit of x[]             */
/* *1 to denote specific quantity       */
/* *N to denote extensive quantity      */
/****************************************/

/* http://131.111.48.24/mackay/c/macopt.html; mackay@mrao.cam.ac.uk */
#define MACOPT_END_IF_SMALL_STEP      0
#define MACOPT_END_IF_SMALL_GRADIENT  1
typedef struct
{
    int convergence_condition; /* use step_square_sum_tolerance criterion */
    bool rich;               /* can afford more computation? */
    int linmin_maxits;     /* in maclinmin() */
    double lastx_default;  /* lastx is set to this if maclinmin() is reset */
    double linmin_g1;    /* factors for growing and shrinking the interval */
    double linmin_g2;
    double linmin_g3;
} Macopt_Option;

typedef struct
{
    Macopt_Option macopt;
    /* zxcgr: */
    int max_potential_evaluations;
    /* zxcgr,macopt: (eV/A)^2*N */
    double gradient_square_sum_tolerance;
    /* macopt: (A^2*N) */
    double step_square_sum_tolerance;
    /* zxcgr: eV*N */
    double rough_potential_decrease_expectation;
} CG_Option;
extern CG_Option CG_option;

typedef struct
{
    int N;        /* dimension of parameter space */
    double lastx; 
    double gtyp ; /* stores the rms gradient for linmin */
    double *g, *h, *xi, *pt, *gx, *gy, *gunused ;
    int restart ; /* whether to restart macopt - fresh cg directions */
} Macopt_State;

typedef struct
{
    double (*potential) (double *x, double *gradient);
    int number_of_evaluations;
    int error;
    char *error_mesg;
    double potential_min;   /* eV*N */
    double *gradient;       /* eV/A */
    Macopt_State macopt;
} CG_State;
extern CG_State CG_state;

void FORTRAN_SYMBOL(zxcgrf77)
    (void (*potential) (int *N, double *x, double *value, double *gradient),
     int *N, double *ACC, int *MAXFN, double *DFPRED, double *X,
     double *G, double *F, double *W, int *IER);

/* Driver of IMSL conjugate gradient minimizer zxcgr_() */
double zxcgr (double (*potential) (double *x, double *gradient),
              int N, double *x);

/* http://131.111.48.24/mackay/c/macopt.html; mackay@mrao.cam.ac.uk */
double macopt (double (*potential) (double *x, double *gradient),
             int N, double *x);

#endif
