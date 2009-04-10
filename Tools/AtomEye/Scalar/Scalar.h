/*******************************************/
/* libScalar:  -lm                         */
/*             -lIO                        */
/*                                         */
/* Scalar Arithmetic Operations and Macros */
/*                                         */
/* Nov 11 1999 Ju Li <liju99@mit.edu>      */
/*******************************************/

#ifndef _Scalar_h
#define _Scalar_h

#include <math.h>
#include <time.h>
#include <IO.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

/** Mathematical Constants **/

#define HALF      (0.5)
#define THIRD     (1./3.)
#define FOURTH    (0.25)
#define FIFTH     (0.2)
#define SIXTH     (1./6.)
#define SEVENTH   (1./7.)
#define EIGHTH    (0.125)
#define NINTH     (1./9.)
#define TENTH     (1./10.)
#define ELEVENTH  (1./11.)
#define TWELFTH   (1./12.)

/* extended precision to 60 digits: see Notes/mathconst.mws */
#define SQRT2 (1.41421356237309504880168872420969807856967187537694807317668)
#define SQRT3 (1.73205080756887729352744634150587236694280525381038062805581)
#define SQRT5 (2.23606797749978969640917366873127623544061835961152572427090)
#define SQRT6 (2.44948974278317809819728407470589139196594748065667012843269)
#define SQRT7 (2.64575131106459059050161575363926042571025918308245018036833)
#define SQRT8 (2.82842712474619009760337744841939615713934375075389614635336)
#define CBRT2 (1.25992104989487316476721060727822835057025146470150798008198)
#define CBRT3 (1.44224957030740838232163831078010958839186925349935057754642)
#define CBRT4 (1.58740105196819947475170563927230826039149332789985300980829)
#define CBRT5 (1.70997594667669698935310887254386010986805511054305492438286)
#define CBRT6 (1.81712059283213965889121175632726050242821046314121967148133)
#define CBRT7 (1.91293118277238910119911683954876028286243905034587576621065)
#define CBRT9 (2.08008382305190411453005682435788538633780534037326210969759)
/* (-1+sqrt(5.))/2 */
#define GOLDEN_RATIO \
  (.618033988749894848204586834365638117720309179805762862135450)
#define GOLDEN_RATIO_RECIPROCAL (GOLDEN_RATIO + 1)

#ifndef PI
#define PI    (3.14159265358979323846264338327950288419716939937510582097494)
#endif
#define RADIAN_TO_DEGREE(r) ((r)/PI*180.)
#define DEGREE_TO_RADIAN(d) ((d)/180.*PI)

#define SIN(x) sin(DEGREE_TO_RADIAN(x))
#define COS(x) cos(DEGREE_TO_RADIAN(x))
#define TAN(x) tan(DEGREE_TO_RADIAN(x))
#define ASIN(x) RADIAN_TO_DEGREE(asin(x))
#define ACOS(x) RADIAN_TO_DEGREE(acos(x))
#define ATAN(x) RADIAN_TO_DEGREE(atan(x))

#ifndef EPS
/* IEEE 754 floating-point arithmetic relative error bound */
#define EPS (2.220446e-16)
#endif

/* Deemed reasonably large number of floating point operations */
#ifndef REASONABLE_FLOPS
#define REASONABLE_FLOPS (1000)
#endif
/* on variables to estimate their likely numerical relative error */

/* When the condition number gets bad, for instance if subtracting */
/* two close numbers, we should override the above to larger value */

#ifndef TINY /* as compared to unity */
#define TINY (REASONABLE_FLOPS*EPS)
#endif

#ifndef HUGE /* as compared to unity */
#define HUGE (1./TINY)
#endif

#ifndef HUGE_INT
#define HUGE_INT (INT_MAX/2)
#endif

/* Close to the floating-point overflow */
/* IEEE 754 32-bit floating-point standard: 1.175494e-38 to 3.402823e38 */
#define SINGLE_PRECISION_INFINITY (1e37)
/* IEEE 754 64-bit floating-point standard: 2.225074e-308 to 1.797693e308 */
#define DOUBLE_PRECISION_INFINITY (1e307)

#ifndef SMALL /* as compared to unity */
#define SMALL (sqrt(EPS))
#endif

#ifndef LARGE /* as compared to unity */
 #define LARGE (1./SMALL)
#endif

/* a wise choice for estimating quadratic curvature */
#ifndef DELTA /* as compared to unity */
#define DELTA (sqrt(sqrt(EPS)))
#endif

/* very large exponent for later purposes of defining inf-norm, etc. */
#ifndef EXPONENT_INFINITY
#define EXPONENT_INFINITY (16.)
#endif

#define EQ(x,y)  ( (x) == (y) )
#define GE(x,y)  ( (x) >= (y) )
#define GT(x,y)  ( (x) >  (y) )
#define LE(x,y)  ( (x) <= (y) )
#define LT(x,y)  ( (x) <  (y) )
/* x = [lower_bound, upper_bound] */
#define XIN(x,lower_bound,upper_bound) (GE(x,lower_bound) && LE(x,upper_bound))
#define XINW(x,width) XIN(x,0.,width)
/* more likely to overflow than underflow, I think */
#define XOU(x,lower_bound,upper_bound) (GT(x,upper_bound) || LT(x,lower_bound))
#define XOUW(x,width) XOU(x,0.,width)
/* with no assumption about which is bigger, bound1 or bound2 */
#define INSIDE(x,bound1,bound2) (XIN(x,bound1,bound2) || XIN(x,bound2,bound1))
#define OUSIDE(x,bound1,bound2) (XOU(x,bound1,bound2) && XOU(x,bound2,bound1))
/* i = [lower_bound, upper_bound) */
#define IN(i,lower_bound,upper_bound) (GE(i,lower_bound) && LT(i,upper_bound))
#define INW(i,n) IN(i,0,n)  /* i = 0..n-1, C convention */
#define OU(i,lower_bound,upper_bound) (GE(i,upper_bound) || LT(i,lower_bound))
#define OUW(i,n) OU(i,0,n)

#define ABS(x)  ((x)>0?(x):-(x))
#define ABS_LE(x,positive_bound) XIN(x,-(positive_bound),(positive_bound))
#define ABS_GT(x,positive_bound) XOU(x,-(positive_bound),(positive_bound))
#define ABS_LT(x,positive_bound) \
  (GT(x,-(positive_bound)) && LT(x,(positive_bound)))
#define ABS_GE(x,positive_bound) \
  (GE(x,(positive_bound)) || LE(x,-(positive_bound)))

#ifndef MIN 
#define MIN(x,y)      ((x)<(y)?(x):(y))
#endif
#ifndef MAX
#define MAX(x,y)      ((x)>(y)?(x):(y))
#endif
#define SWAP(x,y,tmp) ((tmp)=(x), (x)=(y), (y)=(tmp))
#define SWITCH(p,q,state0,state1) CAN ( if ((p) == (state0)) \
  { (p)=(state1); (q)=(state0); } else { (p)=(state0); (q)=(state1); } )
#define SWITCHSTATE(p,q,state) SWITCH(p,q,state,(state)+1)
#define ENSURE(min,max,tmp) { if (min>max) SWAP(min,max,tmp); }
#define FIND3(x0,x1,x2,op,xm) CAN( xm = x0; \
  if (xm op x1) xm = x1; if (xm op x2) xm = x2 )
#define FIND4(x0,x1,x2,x3,op,xm) CAN( xm = x0; \
  if (xm op x1) xm = x1; if (xm op x2) xm = x2; if (xm op x3) xm = x3 )
#define BEPOSITIVE(x)       CAN( if ((x) < 0) x = -(x) )
#define SQUARE(i)     ((i)*(i))
#define CUBE(i)       ((i)*(i)*(i))
#define QUAD(i)       ((i)*(i)*(i)*(i))
#define CEIL(i,j)     (((i)<=0)?0:((i)-1)/(j)+1)
#define SEIL(i,j)     ((size_t)CEIL(i,j))
#define ROUNDUP_TO(i,j)         (CEIL(i,j)*(j))
#define SEPARATION(x,y)         (fabs((x)-(y)))
#define SMALLSEPARATION(x,y,z)  (fabs((x)-(y))<z)

#define SPHERE_VOLUME(radius) (4.*PI/3.*CUBE(radius))
#define safe_avg(sum,count) ((double)(sum)/(((count)!=0)?(count):1))

#define NONEED_TRIM(x)  ( ((x) >= 0.) && ((x) <  1.) )
#define NEED_TRIM(x)    ( ((x) <  0.) || ((x) >= 1.) )
#define Trim(x) { if ((x)<0) do (x)++; while((x)<0); else \
  if ((x)>=1) do (x)--; while ((x)>=1); }
/* faster than TRIM(x) when x is around 1 */
#define TRIM(x) ((x)-floor((double)(x)))  /* x's image on [0,1) */
/* faster than Trim(x) when x is very large */

#define NONEED_IMAGE(x)  ( ((x) >= -0.5) && ((x) <  0.5) )
#define NEED_IMAGE(x)    ( ((x) <  -0.5) || ((x) >= 0.5) )
#define Image(x) { if ((x)<-0.5) do (x)++; while((x)<-0.5); else \
  if ((x)>=0.5) do (x)--; while ((x)>=0.5); }
/* faster than IMAGE(x) when x is around 1 */
#define IMAGE(x)  (TRIM(x+0.5)-0.5) /* x's image on [-0.5,0.5) */
/* faster than Image(x) when x is very large */

#define DISTANCE2(x,y,z) (SQUARE(x)+SQUARE(y)+SQUARE(z))
#define DISTANCE(x,y,z)   sqrt(DISTANCE2(x,y,z))
#define DISTANCED2(x,y)  (SQUARE(x)+SQUARE(y))
#define DISTANCE2D(x,y)   sqrt(DISTANCED2(x,y))

/* b^2 = c^2 - a^2:  the answer b will be nonnegative */
#define SafePythagorasComplement(c,a)                   \
    ( (c)>fabs(a) ? sqrt(SQUARE(c)-SQUARE(a)) : 0 )

/* x is floating-point */
#define Floor(x) (((x)==(int)(x))?(int)(x):(((x)>0)?((int)(x)):((int)(x)-1)))
#define FLOOR(x,i) { (i)=(int)(x); if ((i)>(x)) (i)--; }

#if defined(_alpha)
/* Compaq C warns about negative shifts which cannot actually happen */
#define SAFE_RSHIFT(a,n)  ( ((n)>=0) ? ((a)>>ABS(n)) : ((a)<<ABS(n)) )
#define SAFE_LSHIFT(a,n)  ( ((n)>=0) ? ((a)<<ABS(n)) : ((a)>>ABS(n)) )
#else
#define SAFE_RSHIFT(a,n)  ( ((n)>=0) ? ((a)>>(n)) : ((a)<<(-(n))) )
#define SAFE_LSHIFT(a,n)  ( ((n)>=0) ? ((a)<<(n)) : ((a)>>(-(n))) )
#endif

/* x, target can now be of any numeric type */
#define ERR(x,target) Ferr((double)(x), (double)(target))

#define ISTINY(a) (fabs(a)<TINY)
#define NOTTINY(a) (fabs(a)>=TINY)
#define ISSMALL(a) (fabs(a)<SMALL)
#define NOTSMALL(a) (fabs(a)>=SMALL)
#define ISDELTA(a) (fabs(a)<DELTA)
#define NOTDELTA(a) (fabs(a)>=DELTA)
#define NEAR(x,y)    ISTINY((x)-(y))
#define NOTNEAR(x,y) NOTTINY((x)-(y))

/* whether x and y are equal to floating point accuracy */
#define EQUAL(x,y) ((fabs((double)(x-y))<TINY) || (fabs(ERR(x,y))<TINY))

/* If a number is tiny, set it to zero. USE WITH CARE! */
#define SET_ZERO_IF_TINY(a) { if ISTINY(a) (a)=0.; }

/* Genesis of all randomness */
#define DEFAULT_RANDOM_SEED (1975425) /* for marking and reproducing results */

/* randomizing subsequent Frandom() sequence with seed */
#define Randomize(seed) srandom((unsigned int)(seed))
/* if never called, Frandom() automatically uses seed 1 */
#define TimeRandomize() Randomize(time((time_t *)NULL))

#ifdef _Linux
#define Frandom() ((double)(random()+0.5)/(1.+RAND_MAX))
#define FRANDOM() ((double)(random()+0.5)/(1.+RAND_MAX)-0.5)
#else
#define Frandom() ((double)(random()+0.5)/pow(2.,31.))
#define FRANDOM() ((double)(random()+0.5)/pow(2.,31.)-0.5)
#endif
/* Genesis of all randomness */

#define HALF_DECISION() (Frandom() > 0.5)

/* randomly pick an integer in min..max (set {min,min+1... max}) */
#define Fran(min,max) (INT(min)+(int)floor((INT(max)-INT(min)+1)*Frandom()))
/* randomly pick an integer in min..max-1 (set {min,min+1... max-1}) */
#define fran(min,max) (INT(min)+(int)floor((INT(max)-INT(min))*Frandom()))

/* return normally distributed random number with <x>=E, <(x-E)^2>=sigma2 */
#define Frandnorm(E,sigma2) \
(sqrt(-2.*log(Frandom()))*sin(2*acos(-1.)*Frandom())*sqrt((double)sigma2)+E)
/* The above is not most efficient. For array assignment, use Vfrandnorm() */

/* return normally distributed random number with <x>=0, <x^2>=1 */
#define FRANDNORM() (sqrt(-2.*log(Frandom()))*sin(2*acos(-1.)*Frandom()))

/* return exponentially distributed random number with <x>=tau */
#define Frandexp(tau) (-log(Frandom())*(tau))

/* return exponentially distributed random number with <x>=1 */
#define FRANDEXP() (-log(Frandom()))

/* if a is an integer factor of b */
#define isfactor(a,b) (((int)(b))%((int)(a))==0)
#define isnotfactor(a,b) (((int)(b))%((int)(a))!=0)

/* f(x) := x * Heaviside(x) */
#define ZERO_RAMP(x) (((x)>0)?(x):0)

/** Scalar.c: **/

/* greatest common divisor of a and b */
long gcd (long a, long b);
#define GCD(a,b)  gcd((long)(a),(long)(b))
/* least common multiple */
#define LCM(a,b) (long)(a)*(b)/gcd((long)(a),(long)(b))
/* check if a is a factor of b */
#define ISFACTOR(a,b) (((long)(b))%((long)(a))==0)
#define ISEVEN(a) (((unsigned int)(a)&1)==0)
#define ISODD(a)  (((unsigned int)(a)&1))

/* find k,l such that k * x + l * y = gcd(x,y) */
void diophantine2 (int x, int y, int *k, int *l);

/* set Bitmap "b", starting from the bit "offset", every "period" */
/* bits to the lowest "period" bits of "value", for "n" periods.  */
Bmap *Bperiod (Bmap *b, int offset, int n, int period, int value);

/* relative error of x compared to target */
#define rerr(x,target) ((x)/(target)-1)
double Ferr(double x, double target);

/* calculates |_ a/b _| with mathematically correct floor */
int floorDiv (int a, int b);

/* calculates |^ a/b ^| with mathamatically correct ceiling */
int ceilDiv (int a,int b);

/* return (n!) = GAMMA(n+1) in double precision */
double factorial (int n);

/* Take average on structures whose first element is the counter  */
#define avg_what(type,s)  (sizeof(type)/sizeof(double)),((double *)(s))
#define Avg_what(type,s)  (sizeof(type)/sizeof(double)),((double *)(&(s)))
void avg_clear (int n, double *s);
int avg_add (int n, double *s, ...);
int avg_done (int n, double *s);
int avg_recover (int n, double *s);
/* to avoid alignment trouble please define the counter as double */

/* Generate a string of n-1 random base64 characters ended by EOS. */
char *RandomBase64String (int n, char *a);

/* (-1)%3=-1, but positive_remainder(-1,3)=2 */
int positive_remainder(int a, int b);

/* quadrature.c: */

/* Integration of one-dimensional function: */

/* Gaussian quadrature with W(x)=1, N=10, NR pp.148 */
double qgaus (double (*func)(double), double a, double b);

/* Romberg quadrature: */
#define QROMB_JMAX               20
#define QROMB_K                  5
#define QROMB_DEFAULT_TOLERANCE  1e-14
extern double qromb_tolerance;
double qromb (double (*func)(double), double a, double b);


/* newtonian.c: */

/************************************************************/
/* K-evaluation Integration of Newtonian Dynamics:          */
/*                                                          */
/* Given x(0),dotx(0) and routine ddotx(x) which is derived */
/* from a potential, we want to get x(h),dotx(h) where h is */
/* called a full step, during which K calls to ddotx(x) are */
/* made. Let g = h / K. For a scheme to be competitive, it  */
/* should give better results at t_f = Lh >> h compared to  */
/* other schemes with equal g, for some g/accuracy choices. */
/*                                                          */
/* Algorithms with so-called FSAL property can combine two  */
/* consecutive h-steps to save one evaluation. We provide   */
/* _start, _fsal, _stop routines where _start and _fsal     */
/* need K evaluations and _stop may need one. For complete- */
/* ness, we also give _nofsal which needs K+1 evaluations.  */
/************************************************************/

/* k=0..K-1 (some _nofsal/_stop call k=K), n, x[n], acc[n]=workspace */
typedef void (*Newtonian_ddotx) (int, int, double *, double *);


/***************************/
/* Symplectic Algorithm 1: */
/***************************/
/* First x, then dotx evolution. No great need to use FSAL savings. */

/* This algorithm is not the best for heat current computation or */
/* total energy checksum because when ddotx() is called at k=0, x */
/* is already shifted. Ruth83, Schlier98_6a, Tsitouras99.         */

/****************************************************/
/* Computer "Experiments" on Classical Fluids. I.   */
/* Thermodynamical Properties of Lennard-Jones      */
/* Molecules, L. Verlet, Phys. Rev. 159, 98 (1967). */
/* A Canonical Integration Technique, R.D. Ruth,    */
/* IEEE Trans. Nucl. Sci. NS-30, 2669 (1983).       */
/****************************************************/
/* evaluations per h-step:  1   */
/* global truncation error: h^2 */

/* First x, then dotx evolution. No great need to use FSAL savings. */
void symplectic_Ruth83_nofsal
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx);

/* FSAL sequence: _start, _fsal, _stop */
void symplectic_Ruth83_start
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx);
void symplectic_Ruth83_fsal
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx);
void symplectic_Ruth83_stop
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx);

/**************************************************************/
/* Ch. Schlier and A. Seiter, J. Phys. Chem. 102 (1998) 9399. */
/**************************************************************/
/* evaluations per h-step:  9   */
/* global truncation error: h^6 */

/* First x, then dotx evolution. No great need to use FSAL savings. */
void symplectic_Schlier98_6a_nofsal
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx);

/* FSAL sequence: _start, _fsal, _stop */
void symplectic_Schlier98_6a_start
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx);
void symplectic_Schlier98_6a_fsal
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx);
void symplectic_Schlier98_6a_stop
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx);


/************************************************************/
/* Haruo Yoshida, Phys. Lett. A 150, 262 (1990).            */
/* Celestial Mechanics and Dynamical Astronomy 74 (1999)    */
/* 223-230, Ch. Tsitouras,                                  */
/* www.math.ntua.gr/people/tsitoura/publications.html       */
/************************************************************/
/* evaluations per h-step:  33   */
/* global truncation error: h^10 */

/* First x, then dotx evolution. No great need to use FSAL savings. */
void symplectic_Tsitouras99_nofsal
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx);

/* FSAL sequence: _start, _fsal, _stop */
void symplectic_Tsitouras99_start
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx);
void symplectic_Tsitouras99_fsal
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx);
void symplectic_Tsitouras99_stop
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx);


/***************************/
/* Symplectic Algorithm 2: */
/***************************/
/* First dotx, then x evolution. You do need FSAL savings. */

/* This algorithm is more appropriate for heat current computation */
/* or total energy checksum when ddotx() of _start is called at k  */
/* = 0, if _stop was called at the last step, because exact p(nh), */
/* q(nh) would then be passed to ddotx(). Calvo93, Schlier00_6b,   */
/* Schlier00_8b, Schlier00_8c. */


/************************************************************/
/* The Development of Variable-Step Symplectic Integrators, */
/* With Application to the Two-Body Problem, M.P. Calvo,    */
/* J.M. Sanz-Serna, SIAM J. Sci. Comput. 14, 936-952 (1993) */
/* Numerical Hamiltonian Problems, J.M. Sanz-Serna,         */
/* M.P. Calvo, Chapman & Hall, London (1994).               */
/************************************************************/
/* evaluations per h-step:  4   */
/* global truncation error: h^4 */

/* First dotx, then x evolution. No FSAL savings yet. */
void symplectic_Calvo93_nofsal
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx);

/* FSAL sequence: _start, _fsal, _stop */
void symplectic_Calvo93_start
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx);
void symplectic_Calvo93_fsal
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx);
void symplectic_Calvo93_stop
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx);


/*******************************************************/
/* High-Order Symplectic Integration: An Assessment,   */
/* Ch. Schlier, A. Seiter, Comp. Phys. Comm. (2000)    */
/* phya8.physik.uni-freiburg.de/abt/papers/papers.html */
/*******************************************************/
/* evaluations per h-step:  8   */
/* global truncation error: h^6 */

/* First dotx, then x evolution. No FSAL savings yet. */
void symplectic_Schlier00_6b_nofsal
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx);

/* FSAL sequence: _start, _fsal, _stop */
void symplectic_Schlier00_6b_start
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx);
void symplectic_Schlier00_6b_fsal
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx);
void symplectic_Schlier00_6b_stop
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx);


/*******************************************************/
/* Ch. Schlier and A. Seiter, Comp. Phys. Comm. (2000) */
/*******************************************************/
/* evaluations per h-step:  17  */
/* global truncation error: h^8 */

/* First dotx, then x evolution. No FSAL savings yet.  */
void symplectic_Schlier00_8b_nofsal
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx);

/* FSAL sequence: _start, _fsal, _stop */
void symplectic_Schlier00_8b_start
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx);
void symplectic_Schlier00_8b_fsal
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx);
void symplectic_Schlier00_8b_stop
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx);


/*******************************************************/
/* Ch. Schlier and A. Seiter, Comp. Phys. Comm. (2000) */
/*******************************************************/
/* evaluations per h-step:  17  */
/* global truncation error: h^8 */

/* First dotx, then x evolution. No FSAL savings yet.  */
void symplectic_Schlier00_8c_nofsal
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx);

/* FSAL sequence: _start, _fsal, _stop */
void symplectic_Schlier00_8c_start
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx);
void symplectic_Schlier00_8c_fsal
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx);
void symplectic_Schlier00_8c_stop
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx);


/**********************************************************************/
/* Routers to various routines by K, the number of evaluations per h. */
/* To make fair comparison, one should use equal g for different K's. */
/* It gives "exact" x(wh) and dotx(wh), w=1..L, after each full step. */
/**********************************************************************/
void symplectic_nofsal ( int K, double g, int n, double *x, double *dotx,
                         double *acc, Newtonian_ddotx ddotx );


/**********************************************************************/
/* Routers to various routines by K, the number of evaluations per h. */
/* To make fair comparison, one should use equal g for different K's. */
/* To integrate t from 0 to Lh, one should call 1 x symplectic_start, */
/* (L-1) x symplectic_fsal, and 1 x symplectic_stop. One cannot have  */
/* "exact" x(wh) and dotx(wh) except when w=0,L, but capturing k=0    */
/* calls to ddotx(k,n,x,acc) at w=1..L-1 gives approximate answers.   */
/**********************************************************************/
void symplectic_start ( int K, double g, int n, double *x, double *dotx,
                        double *acc, Newtonian_ddotx ddotx );
void symplectic_fsal ( int K, double g, int n, double *x, double *dotx,
                       double *acc, Newtonian_ddotx ddotx );
void symplectic_stop ( int K, double g, int n, double *x, double *dotx,
                       double *acc, Newtonian_ddotx ddotx );


/*******************************************************/
/* Well known non-symplectic integrators & Benchmarks: */
/*******************************************************/

/************************************************************/
/* classic 4th-order Runge-Kutta method: Numerical Recipes  */
/* in C: The Art of Scientific Computing, William H. Press, */
/* Saul A. Teukolsky, William T. Vetterling, Brian P.       */
/* Flannery, Cambridge University Press, Cambridge (1992).  */
/* http://mathworld.wolfram.com/Runge-KuttaMethod.html      */
/************************************************************/
/* evaluations per h-step:  4   */
/* global truncation error: h^4 */

/* suited for heat current computation & total energy checksum at k=0 */

/* allocate 8*n internal workspace. No integration is performed yet. */
void rk4_prep ( double h, int n, double *x, double *dotx, double *acc,
                Newtonian_ddotx ddotx );
/* "acc" would not be used, so it can be passed as NULL. */
void rk4 (double h, int n, double *x, double *dotx, double *acc,
          Newtonian_ddotx ddotx);
/* free 8*n internal workspace */
void rk4_free ( double h, int n, double *x, double *dotx, double *acc,
                Newtonian_ddotx ddotx );


/**********************************************************/
/* Gear predictor-corrector algorithm: Numerical Initial  */
/* Value Problems in Ordinary Differential Equation,      */
/* C.W. Gear, Englewood Cliffs N.J., Prentice-Hall, 1971. */
/* itp.nat.uni-magdeburg.de/~schinner/dip/node22.html     */
/**********************************************************/
/* k-value: keep k n-arrays: x^0,x^1,..,x^{k-1}:          */
/* initialized by "exact" time derivatives.               */
/* Before corrector, x0 has local truncation error h^k;   */
/* after corrector, x0 has local truncation error h^(k+1) */
/* therefore the global truncation error of x0 is h^k.    */
/* The global truncation error of dotx is h^{k-1}, etc.   */
/**********************************************************/

#define GEAR_TEST_H  50

/* 4-value Gear predictor-corrector algorithm */
/* evaluations per h-step:  1   */
/* local  truncation error: h^5 */
/* global truncation error: h^4 */

/* Initialize 4-value Gear integrator and allocate 2*n  */
/* internal workspace. No integration is performed yet. */
void gear4_prep ( double h, int n, double *x, double *dotx, double *acc,
                  Newtonian_ddotx ddotx );
void gear4 ( double h, int n, double *x, double *dotx, double *acc,
             Newtonian_ddotx ddotx );
/* free 2*n internal workspace */
void gear4_free ( double h, int n, double *x, double *dotx, double *acc,
                  Newtonian_ddotx ddotx );

/******************************************************/
/* 5-value Gear predictor-corrector algorithm:        */
/* Computer Simulation of Liquids, M.P. Allen and     */
/* D.J. Tildesley, Clarendon Press, (Oxford 1987).    */
/* ftp://ftp.dl.ac.uk/ccp5/ALLEN_TILDESLEY/F.02       */
/* itp.nat.uni-magdeburg.de/~schinner/dip/node22.html */
/******************************************************/
/* evaluations per h-step:  1   */
/* local  truncation error: h^6 */
/* global truncation error: h^5 */

/* Initialize 5-value Gear integrator and allocate 3*n  */
/* internal workspace. No integration is performed yet. */
void gear5_prep ( double h, int n, double *x, double *dotx, double *acc,
                  Newtonian_ddotx ddotx );
void gear5 ( double h, int n, double *x, double *dotx, double *acc,
             Newtonian_ddotx ddotx );
/* free 3*n internal workspace */
void gear5_free ( double h, int n, double *x, double *dotx, double *acc,
                  Newtonian_ddotx ddotx );


/* 6-value Gear predictor-corrector algorithm */
/* evaluations per h-step:  1   */
/* local  truncation error: h^7 */
/* global truncation error: h^6 */

/* Initialize 6-value Gear integrator and allocate 4*n  */
/* internal workspace. No integration is performed yet. */
void gear6_prep ( double h, int n, double *x, double *dotx, double *acc,
                  Newtonian_ddotx ddotx );
void gear6 ( double h, int n, double *x, double *dotx, double *acc,
             Newtonian_ddotx ddotx );
/* free 4*n internal workspace */
void gear6_free ( double h, int n, double *x, double *dotx, double *acc,
                  Newtonian_ddotx ddotx );

/* 7-value Gear predictor-corrector algorithm */
/* evaluations per h-step:  1   */
/* local  truncation error: h^8 */
/* global truncation error: h^7 */

/* Initialize 7-value Gear integrator and allocate 5*n  */
/* internal workspace. No integration is performed yet. */
void gear7_prep ( double h, int n, double *x, double *dotx, double *acc,
                  Newtonian_ddotx ddotx );
void gear7 ( double h, int n, double *x, double *dotx, double *acc,
             Newtonian_ddotx ddotx );
/* free 5*n internal workspace */
void gear7_free ( double h, int n, double *x, double *dotx, double *acc,
                  Newtonian_ddotx ddotx );

/* 8-value Gear predictor-corrector algorithm */
/* evaluations per h-step:  1   */
/* local  truncation error: h^9 */
/* global truncation error: h^8 */

/* Initialize 8-value Gear integrator and allocate 6*n  */
/* internal workspace. No integration is performed yet. */
void gear8_prep ( double h, int n, double *x, double *dotx, double *acc,
                  Newtonian_ddotx ddotx );
void gear8 ( double h, int n, double *x, double *dotx, double *acc,
             Newtonian_ddotx ddotx );
/* free 6*n internal workspace */
void gear8_free ( double h, int n, double *x, double *dotx, double *acc,
                  Newtonian_ddotx ddotx );

/* symplectic algorithm 1 */
#define NEWTONIAN_RUTH83          0
#define NEWTONIAN_SCHLIER98_6A   (NEWTONIAN_RUTH83 + 1)
#define NEWTONIAN_TSITOURAS99    (NEWTONIAN_SCHLIER98_6A + 1)
#define NEWTONIAN_SYMPLECTIC_ALGORITHM_1_MAX  (NEWTONIAN_TSITOURAS99 + 1)
/* symplectic algorithm 2 */
#define NEWTONIAN_CALVO93        NEWTONIAN_SYMPLECTIC_ALGORITHM_1_MAX
#define NEWTONIAN_SCHLIER00_6B   (NEWTONIAN_CALVO93 + 1)
#define NEWTONIAN_SCHLIER00_8B   (NEWTONIAN_SCHLIER00_6B + 1)
#define NEWTONIAN_SCHLIER00_8C   (NEWTONIAN_SCHLIER00_8B + 1)
#define NEWTONIAN_SYMPLECTIC_ALGORITHM_2_MAX  (NEWTONIAN_SCHLIER00_8C + 1)
/* non-symplectic algorithms */
#define NEWTONIAN_RK4            NEWTONIAN_SYMPLECTIC_ALGORITHM_2_MAX
#define NEWTONIAN_GEAR4          (NEWTONIAN_RK4 + 1)
#define NEWTONIAN_GEAR5          (NEWTONIAN_GEAR4 + 1)
#define NEWTONIAN_GEAR6          (NEWTONIAN_GEAR5 + 1)
#define NEWTONIAN_GEAR7          (NEWTONIAN_GEAR6 + 1)
#define NEWTONIAN_GEAR8          (NEWTONIAN_GEAR7 + 1)
#define NEWTONIAN_MAX            (NEWTONIAN_GEAR8 + 1)

typedef void (*Integrator) ( double, int, double *, double *, double *,
                             Newtonian_ddotx );
typedef struct
{
    char *name;
    bool is_symplectic;
    int K;
    int order;
    Integrator a,b,c;
} Newtonian;

extern Newtonian newtonian [NEWTONIAN_MAX];

/***************************************************************/
/* prep  +  [ firststep, L-1 nextstep, halt ] repeat  +  free  */
/*            |                                                */
/* k=0 property calculation for algorithm 2 and non-symplectic */
/***************************************************************/

/* Allocate and initialize states if necessary, but do NOT integrate. */
void newtonian_integrator_prep
( int id, double g, int n, double *x, double *dotx,
  double *acc, Newtonian_ddotx ddotx );

/* Get the integrator to run its first full step of a session. */
void newtonian_integrator_firststep
( int id, double g, int n, double *x, double *dotx,
  double *acc, Newtonian_ddotx ddotx );

/* Run other full steps of a session. */
void newtonian_integrator_nextstep
( int id, double g, int n, double *x, double *dotx,
  double *acc, Newtonian_ddotx ddotx );

/* Halt the session and do necessary mop-ups. To */
/* start a new session, call _firststep again.   */
void newtonian_integrator_halt
( int id, double g, int n, double *x, double *dotx,
  double *acc, Newtonian_ddotx ddotx );

/* Free workspace if any. No more sessions. To restart, call _prep. */
void newtonian_integrator_free
( int id, double g, int n, double *x, double *dotx,
  double *acc, Newtonian_ddotx ddotx );


/** spline.c: **/

#define SPLINE_CUBIC_NATURAL   0
#define SPLINE_CUBIC_PERIODIC  1
#define SPLINE_AKIMA_NATURAL   2
#define SPLINE_AKIMA_PERIODIC  3

typedef struct
{
    gsl_spline *spline;
    gsl_interp_accel *acc;
    double xmin, xmax, ymax, ymin;
} Spline;

/* Allocate and initialize Spline structure based on */
/* lookup table x[0..N-1],y[0..N-1] and spline type. */
Spline *spline_initialize (int N, double *x, double *y, int spline_type);

/* Evaluate spline value. */
double spline_value (Spline *sp, double x);

/* Evaluate first-order derivative. */
double spline_deriv (Spline *sp, double x);

/* Evaluate second-order derivative. */
double spline_deriv2 (Spline *sp, double x);

/* Free all data structure. */
void spline_finalize (Spline *sp);


/** histogram.c: **/

typedef struct
{
    double xmin, xmax;   /* domain of interest: {xmin,xmax} */
    int bins;            /* # of collection bins in the domain of interest */
    double xbin;         /* width of one bin */
    int totalcount;      /* counts of total submission */
    double domaincount;  /* counts that fall into the domain of interest */
    double *bincount;    /* count in each bin */
    double *x_sample;    /* center x of each bin */
    double *density_profile;  /* density profile which integrates to 1 */
} Histogram;

/* Given finite xmin < xmax and bins >= 1, reallocate */
/* data structure in h and reset all counters.        */
void histogram_reinitialize
(double xmin, double xmax, int bins, Histogram *h);

/* Update histogram by data array x[0..n-1]; return domaincount */
double histogram_intercept (int n, double *x, Histogram *h);

/* Return the interception rate = domaincount / totalcount */
double histogram_interception_rate (Histogram *h);

/* Update the density profile; return the interception rate */
double histogram_update_density_profile (Histogram *h);

/* Save histogram in matlab script "fname" as  */
/* symbol "token" and notify I/O stream "info" */
double histogram_save_as_matlab
(Histogram *h, FILE *info, char *token, char *fname);
#define histogram_save_as_MATLAB(h,info,token) histogram_save_as_matlab\
  (h,info,token,str2(token,"_h.m"))

/* Clear all counters as if it is just the beginning */
void histogram_reset_all_counters (Histogram *h);

/* Free the histogram data structure and reset all counters */
void histogram_free (Histogram *h);

#endif  /* _Scalar_h */
