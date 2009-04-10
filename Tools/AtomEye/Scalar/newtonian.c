/*******************************************/
/* libScalar:                              */
/*                                         */
/* Scalar Arithmetic Operations and Macros */
/*                                         */
/* Nov 11 1999 Ju Li <liju99@mit.edu>      */
/*******************************************/

#include "Scalar.h"

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

/************************************************************/
/* 1. Construction of Higher Order Symplectic Integrators,  */
/* Haruo Yoshida, Phys. Lett. A 150, 262 (1990).            */
/*                                                          */
/* 2. The Development of Variable-Step Symplectic           */
/* Integrators, With Application to the Two-Body Problem,   */
/* M.P. Calvo, J.M. Sanz-Serna, SIAM J. Sci. Comput. 14,    */
/* 936 (1993); Numerical Hamiltonian Problems, J.M.         */
/* Sanz-Serna, M.P. Calvo, Chapman & Hall, London (1994).   */
/*                                                          */
/* 3. Symplectic Integrator for Molecular Dynamics of a     */
/* Protein in Water, H. Ishida, Y. Nagai, A. Kidera,        */
/* Chemical Physics Letters 282, 115-120 (1998).            */
/*                                                          */
/* 4. Symplectic Integration of Classical Trajectories: A   */
/* Case Study, Ch. Schlier and A. Seiter, J. Phys. Chem.    */
/* 102 (1998) 9399; High-Order Symplectic Integration: An   */
/* Assessment, Ch. Schlier and A. Seiter, Comp. Phys. Comm. */
/* xxx (2000) p. yyy;                                       */
/* phya8.physik.uni-freiburg.de/abt/papers/papers.html      */
/*                                                          */
/* 5. A tenth order symplectic Runge-Kutta-Nystrom method,  */
/* Celestial Mechanics and Dynamical Astronomy 74 (1999)    */
/* 223-230, Ch. Tsitouras,                                  */
/* www.math.ntua.gr/people/tsitoura/publications.html       */
/************************************************************/


static void advance_x (double how_much, int n, double *dotx, double *x)
{
    register int i;
    for (i=n; i--;) x[i] += how_much * dotx[i];
    return;
} /* end advance_x() */


/* k-th evaluation of the acceleration */
static void advance_dotx
(int k, double how_much, int n, double *x, double *acc,
 Newtonian_ddotx ddotx, double *dotx)
{
    register int i;
    (*ddotx) (k, n, x, acc);
    for (i=n; i--;) dotx[i] += how_much * acc[i];
    return;
} /* end advance_x() */


/***************************/
/* Symplectic Algorithm 1: */
/***************************/

/* This algorithm is not the best for heat current computation or */
/* total energy checksum because when ddotx() is called at k=0, x */
/* is already shifted. Ruth83, Schlier98_6a, Tsitouras99.         */

/* First x, then dotx evolution. No great need to use FSAL savings. */
static void Algorithm1_nofsal
(int K, const double *a, double h, int n, double *x, double *dotx,
 double *acc, Newtonian_ddotx ddotx)
{
    int k;
    for (k=0; k<K; k++)
    {
        advance_x ( a[2*k] * h , n , dotx, x );
        /* k = 0..K-1 */
        advance_dotx ( k , a[2*k+1] * h , n , x , acc , ddotx , dotx );
    }
    advance_x ( a[2*K] * h , n , dotx, x );
    return;
} /* end Algorithm1_nofsal() */


/* FSAL sequence: _start, _fsal, _stop */
static void Algorithm1_start
(int K, const double *a, double h, int n, double *x, double *dotx,
 double *acc, Newtonian_ddotx ddotx)
{
    int k;
    for (k=0; k<K; k++)
    {
        advance_x ( a[2*k] * h , n , dotx, x );
        /* k = 0..K-1 */
        advance_dotx ( k , a[2*k+1] * h , n , x , acc , ddotx , dotx );
    }
    return;
} /* end Algorithm1_start() */

static void Algorithm1_fsal
(int K, const double *a, double h, int n, double *x, double *dotx,
 double *acc, Newtonian_ddotx ddotx)
{
    int k;
    advance_x ( (a[2*K]+a[0]) * h , n , dotx, x );
    /* k=0. _nofsal and _fsal makes no difference to ddotx */
    advance_dotx ( 0 , a[1] * h , n , x , acc , ddotx , dotx );
    for (k=1; k<K; k++)
    {
        advance_x ( a[2*k] * h , n , dotx, x );
        /* k = 1..K-1 */
        advance_dotx ( k , a[2*k+1] * h , n , x , acc , ddotx , dotx );
    }
    return;
} /* end Algorithm1_fsal() */

static void Algorithm1_stop
(int K, const double *a, double h, int n, double *x, double *dotx,
 double *acc, Newtonian_ddotx ddotx)
{ /* this _stop does NOT call k = K */
    advance_x ( a[2*K] * h , n , dotx, x );
    return;
} /* end Algorithm1_stop() */


/****************************************************/
/* Computer "Experiments" on Classical Fluids. I.   */
/* Thermodynamical Properties of Lennard-Jones      */
/* Molecules, L. Verlet, Phys. Rev. 159, 98 (1967). */
/* A Canonical Integration Technique, R.D. Ruth,    */
/* IEEE Trans. Nucl. Sci. NS-30, 2669 (1983).       */
/****************************************************/
#define Ruth83_K  1
/* global truncation error: h^2 */
static const double Ruth83 [2*Ruth83_K+1] =
{
    0.5,
    1,
    0.5
};

/* First x, then dotx evolution. No great need to use FSAL savings. */
void symplectic_Ruth83_nofsal
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx)
{
    Algorithm1_nofsal
        ( Ruth83_K, Ruth83, h, n, x, dotx, acc, ddotx );
    return;
} /* end symplectic_Ruth83_nofsal() */

/* FSAL sequence: _start, _fsal, _stop */
void symplectic_Ruth83_start
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx)
{
    Algorithm1_start
        ( Ruth83_K, Ruth83, h, n, x, dotx, acc, ddotx );
    return;
} /* end symplectic_Ruth83_start() */

void symplectic_Ruth83_fsal
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx)
{
    Algorithm1_fsal
        ( Ruth83_K, Ruth83, h, n, x, dotx, acc, ddotx );
    return;
} /* end symplectic_Ruth83_fsal() */

void symplectic_Ruth83_stop
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx)
{
    Algorithm1_stop
        ( Ruth83_K, Ruth83, h, n, x, dotx, acc, ddotx );
    return;
} /* end symplectic_Ruth83_stop() */


/* Ch. Schlier and A. Seiter, J. Phys. Chem. 102 (1998) 9399. */
#define Schlier98_6a_K  9
/* global truncation error: h^6 */
static const double Schlier98_6a[2*Schlier98_6a_K+1] =
{
    0.09517625454177405267746114335519342,
    0.66629689399770780134207498907168068,
    -0.12795028552368677941219191621429411,
    0.02461890095210508713078430308713062,
    0.10597295345325113143793587608716998,
    -0.41072553361795113231992873918199025,
    0.44822227660082748416851634186561201,
    0.65772926205091317768935130009339042,
    -0.02142119907216588887172144509368130,

    -0.87583904676554986768456370614042295,
    
    -0.02142119907216588887172144509368130,
    0.65772926205091317768935130009339042,
    0.44822227660082748416851634186561201,
    -0.41072553361795113231992873918199025,
    0.10597295345325113143793587608716998,
    0.02461890095210508713078430308713062,
    -0.12795028552368677941219191621429411,
    0.66629689399770780134207498907168068,
    0.09517625454177405267746114335519342
};


/* First x, then dotx evolution. No great need to use FSAL savings. */
void symplectic_Schlier98_6a_nofsal
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx)
{
    Algorithm1_nofsal
        ( Schlier98_6a_K, Schlier98_6a, h, n, x, dotx, acc, ddotx );
    return;
} /* end symplectic_Schlier98_6a_nofsal() */


/* FSAL sequence: _start, _fsal, _stop */
void symplectic_Schlier98_6a_start
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx)
{
    Algorithm1_start
        ( Schlier98_6a_K, Schlier98_6a, h, n, x, dotx, acc, ddotx );
    return;
} /* end symplectic_Schlier98_6a_start() */

void symplectic_Schlier98_6a_fsal
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx)
{
    Algorithm1_fsal
        ( Schlier98_6a_K, Schlier98_6a, h, n, x, dotx, acc, ddotx );
    return;
} /* end symplectic_Schlier98_6a_fsal() */

void symplectic_Schlier98_6a_stop
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx)
{
    Algorithm1_stop
        ( Schlier98_6a_K, Schlier98_6a, h, n, x, dotx, acc, ddotx );
    return;
} /* end symplectic_Schlier98_6a_stop() */
        

/************************************************************/
/* Haruo Yoshida, Phys. Lett. A 150, 262 (1990).            */
/* Celestial Mechanics and Dynamical Astronomy 74 (1999)    */
/* 223-230, Ch. Tsitouras,                                  */
/* www.math.ntua.gr/people/tsitoura/publications.html       */
/************************************************************/
/* w1  := 0.02690013604768968151437422144441685467297755661;  */
/* w2  := 0.939801567135683337900037741418097674632505563;    */
/* w3  := -0.00803583920385358749646880826318191806393661063; */
/* w4  := -0.866485197373761372803661767454208401679117010;   */
/* w5  := 0.1023112911193598731078563285067131541328142449;   */
/* w6  := -0.1970772151393080101376018465105491958660525085;  */
/* w7  := 0.617877713318069357335731125307019691019646679;    */
/* w8  := 0.1907272896000121001605903836891198270441436012;   */
/* w9  := 0.2072605028852482559382954630002620777969060377;   */
/* w10 := -0.395006197760920667393122535979679328161187572;   */
/* w11 := -0.582423447311644594710573905438945739845940956;   */
/* w12 := 0.742673314357319863476853599632017373530365297;    */
/* w13 := 0.1643375495204672910151440244080443210570501579;   */
/* w14 := -0.615116639060545182658778437156157368647972997;   */
/* w15 := 0.2017504140367640350582861633379013481712172488;   */
/* w16 := 0.45238717224346720617588658607423353932336395045;  */
/* see Notes/Tsitouras99.mws */

#define Tsitouras99_K  33
/* global truncation error: h^10 */
static const double Tsitouras99 [2*Tsitouras99_K+1] =
{
    .22619358612173360308794329303711676966168197522500,
    .45238717224346720617588658607423353932336395045,
    .32706879314011562061708637470606744374729059962500,
    .2017504140367640350582861633379013481712172488,
    -.20668311251189057380024613690912801023837787410000,
    -.615116639060545182658778437156157368647972997,
    -.22538954477003894582181720637405652379546141955000,
    .1643375495204672910151440244080443210570501579,
    .45350543193889357724599881202003084729370772745000,
    .742673314357319863476853599632017373530365297,
    .08012493352283763438313984709653581684221217050000,
    -.582423447311644594710573905438945739845940956,
    -.48871482253628263105184822070931253400356426400000,
    -.395006197760920667393122535979679328161187572,
    -.09387284743783620572741353648970862518214076715000,
    .2072605028852482559382954630002620777969060377,
    .19899389624263017804944292334469095242052481945000,
    .1907272896000121001605903836891198270441436012,
    .40430250145904072874816075449806975903189514010000,
    .617877713318069357335731125307019691019646679,
    .21040024908938067359906463939823524757679708525000,
    -.1970772151393080101376018465105491958660525085,
    -.047382962009974068514872759001918020866619131800000,
    .1023112911193598731078563285067131541328142449,
    -.38208695312720074984790271947374762377315138255000,
    -.866485197373761372803661767454208401679117010,
    -.43726051828880748015006528785869515987152681031500,
    -.00803583920385358749646880826318191806393661063,
    .46588286396591487520178446657745787828428447618500,
    .939801567135683337900037741418097674632505563,
    .48335085159168650970720598143125726465274155980500,
    .02690013604768968151437422144441685467297755661,
    -.46843234639020274572566122529289548178029390412500,

    -.96376482882809517296569667203020781823356536486,
    
    -.46843234639020274572566122529289548178029390412500,
    .02690013604768968151437422144441685467297755661,
    .48335085159168650970720598143125726465274155980500,
    .939801567135683337900037741418097674632505563,
    .46588286396591487520178446657745787828428447618500,
    -.00803583920385358749646880826318191806393661063,
    -.43726051828880748015006528785869515987152681031500,
    -.866485197373761372803661767454208401679117010,
    -.38208695312720074984790271947374762377315138255000,
    .1023112911193598731078563285067131541328142449,
    -.047382962009974068514872759001918020866619131800000,
    -.1970772151393080101376018465105491958660525085,
    .21040024908938067359906463939823524757679708525000,
    .617877713318069357335731125307019691019646679,
    .40430250145904072874816075449806975903189514010000,
    .1907272896000121001605903836891198270441436012,
    .19899389624263017804944292334469095242052481945000,
    .2072605028852482559382954630002620777969060377,
    -.09387284743783620572741353648970862518214076715000,
    -.395006197760920667393122535979679328161187572,
    -.48871482253628263105184822070931253400356426400000,
    -.582423447311644594710573905438945739845940956,
    .08012493352283763438313984709653581684221217050000,
    .742673314357319863476853599632017373530365297,
    .45350543193889357724599881202003084729370772745000,
    .1643375495204672910151440244080443210570501579,
    -.22538954477003894582181720637405652379546141955000,
    -.615116639060545182658778437156157368647972997,
    -.20668311251189057380024613690912801023837787410000,
    .2017504140367640350582861633379013481712172488,
    .32706879314011562061708637470606744374729059962500,
    .45238717224346720617588658607423353932336395045,
    .22619358612173360308794329303711676966168197522500
};

/* First x, then dotx evolution. No great need to use FSAL savings. */
void symplectic_Tsitouras99_nofsal
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx)
{
    Algorithm1_nofsal
        ( Tsitouras99_K, Tsitouras99, h, n, x, dotx, acc, ddotx );
    return;
} /* end symplectic_Tsitouras99_nofsal() */


/* FSAL sequence: _start, _fsal, _stop */
void symplectic_Tsitouras99_start
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx)
{
    Algorithm1_start
        ( Tsitouras99_K, Tsitouras99, h, n, x, dotx, acc, ddotx );
    return;
} /* end symplectic_Tsitouras99_start() */

void symplectic_Tsitouras99_fsal
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx)
{
    Algorithm1_fsal
        ( Tsitouras99_K, Tsitouras99, h, n, x, dotx, acc, ddotx );
    return;
} /* end symplectic_Tsitouras99_fsal() */

void symplectic_Tsitouras99_stop
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx)
{
    Algorithm1_stop
        ( Tsitouras99_K, Tsitouras99, h, n, x, dotx, acc, ddotx );
    return;
} /* end symplectic_Tsitouras99_stop() */


/***************************/
/* Symplectic Algorithm 2: */
/***************************/

/* This algorithm is more appropriate for heat current computation */
/* or total energy checksum when ddotx() of _start is called at k  */
/* = 0, if _stop was called at the last step, because exact p(nh), */
/* q(nh) would then be passed to ddotx(). Calvo93, Schlier00_6b,   */
/* Schlier00_8b, Schlier00_8c. */

/* First dotx, then x evolution. No FSAL savings yet. */
static void Algorithm2_nofsal
(int K, const double *a, double h, int n, double *x, double *dotx,
 double *acc, Newtonian_ddotx ddotx)
{
    int k;
    for (k=0; k<K; k++)
    { /* k = 0..K-1: at k = 0 it has both right p(nh) and q(nh) */
        advance_dotx ( k , a[2*k] * h , n , x , acc , ddotx , dotx );
        advance_x ( a[2*k+1] * h , n , dotx, x );
    }
    /* k = K: this has right q(nh), but not right p(nh) */
    advance_dotx ( K , a[2*K] * h , n , x , acc , ddotx , dotx );
    return;
} /* end Algorithm2_nofsal() */

/* FSAL sequence: _start, _fsal, _stop */
static void Algorithm2_start
(int K, const double *a, double h, int n, double *x, double *dotx,
 double *acc, Newtonian_ddotx ddotx)
{
    int k;
    for (k=0; k<K; k++)
    { /* k = 0..K-1: at k = 0 it has both right p(nh) and q(nh) */
        advance_dotx ( k , a[2*k] * h , n , x , acc , ddotx , dotx );
        advance_x ( a[2*k+1] * h , n , dotx, x );
    }
    /* finish with right q((n+1)h) but not p((n+1)h) */
    return;
} /* end Algorithm2_start() */

static void Algorithm2_fsal
(int K, const double *a, double h, int n, double *x, double *dotx,
 double *acc, Newtonian_ddotx ddotx)
{
    int k;
    /* call k = 0: q(nh) is right but not p(nh) */
    advance_dotx ( 0 , (a[2*K]+a[0]) * h , n , x , acc , ddotx , dotx );
    advance_x ( a[1] * h , n , dotx, x );
    for (k=1; k<K; k++)
    { /* k = 1..K-1 */
        advance_dotx ( k , a[2*k] * h , n , x , acc , ddotx , dotx );
        advance_x ( a[2*k+1] * h , n , dotx, x );
    }
    /* finish with right q((n+1)h) but not p((n+1)h) */
    return;
} /* end Algorithm2_fsal() */

static void Algorithm2_stop
(int K, const double *a, double h, int n, double *x, double *dotx,
 double *acc, Newtonian_ddotx ddotx)
{ /* call k = K: right q((n+1)h) but not p((n+1)h) */
    advance_dotx ( K , a[2*K] * h , n , x , acc , ddotx , dotx );
    /* right q((n+1)h) and p((n+1)h) */
    return;
} /* end Algorithm2_stop() */


/************************************************************/
/* The Development of Variable-Step Symplectic Integrators, */
/* With Application to the Two-Body Problem, M.P. Calvo,    */
/* J.M. Sanz-Serna, SIAM J. Sci. Comput. 14, 936-952 (1993) */
/* Numerical Hamiltonian Problems, J.M. Sanz-Serna,         */
/* M.P. Calvo, Chapman & Hall, London (1994).               */
/************************************************************/
#define Calvo93_K  4
/* global truncation error: h^4 */
static const double Calvo93 [2*Calvo93_K+1] =
{
    /* 0.0617588581356263250, */
    /* 0.2051776615422863900, */
    /* 0.3389780265536433551, */
    /* 0.4030212816042146300, */
    /* 0.6147913071755775662, */
    /* -0.1209208763389140700, */
    /* -0.1405480146593733802, */
    /* 0.5127219331924131000, */
    /* 0.1250198227945261338 */
    /* www.glue.umd.edu/~dsweet/ChaosClasses/IntegratorsDoc/crk4C_h.html */
    0.061758858135626325,
    0.205177661542286386,
    0.338978026553643355,
    0.403021281604214587,
    0.614791307175577566,
    -0.120920876338914008,
    -0.140548014659373380,
    0.512721933192413035,
    0.125019822794526133
    /* Symplectic Integrator for Molecular Dynamics of a */
    /* Protein in Water, H. Ishida, Y. Nagai, A. Kidera, */
    /* Chemical Physics Letters 282, 115-120 (1998).     */
};
/* note that Calvo's coefficients were optimized and is not symmetric */


/* First dotx, then x evolution. No FSAL savings yet. */
void symplectic_Calvo93_nofsal
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx)
{
    Algorithm2_nofsal
        ( Calvo93_K, Calvo93, h, n, x, dotx, acc, ddotx );
    return;
} /* end symplectic_Calvo93_nofsal() */

/* FSAL sequence: _start, _fsal, _stop */
void symplectic_Calvo93_start
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx)
{
    Algorithm2_start
        ( Calvo93_K, Calvo93, h, n, x, dotx, acc, ddotx );
    return;
} /* end symplectic_Calvo93_start() */

void symplectic_Calvo93_fsal
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx)
{
    Algorithm2_fsal
        ( Calvo93_K, Calvo93, h, n, x, dotx, acc, ddotx );
    return;
} /* end symplectic_Calvo93_fsal() */

void symplectic_Calvo93_stop
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx)
{
    Algorithm2_stop
        ( Calvo93_K, Calvo93, h, n, x, dotx, acc, ddotx );
    return;
} /* end symplectic_Calvo93_stop() */


/*******************************************************/
/* High-Order Symplectic Integration: An Assessment,   */
/* Ch. Schlier, A. Seiter, Comp. Phys. Comm. (2000)    */
/* phya8.physik.uni-freiburg.de/abt/papers/papers.html */
/*******************************************************/
#define Schlier00_6b_K  8
/* global truncation error: h^6 */
static const double Schlier00_6b[2*Schlier00_6b_K+1] =
{
    0.06942944346252987735848865824703402,
    0.28487837717280084052745346456657828,
    -0.13315519831598209409961309951373512,
    0.32783975759612945412054678367325547,
    0.00129038917981078974230481746443284,
    -0.38122104271932629475622784374211274,
    0.42243536567364142699881962380226825,
    0.26850290795039600010822759550227899,

    0.28,

    0.26850290795039600010822759550227899,
    0.42243536567364142699881962380226825,
    -0.38122104271932629475622784374211274,
    0.00129038917981078974230481746443284,
    0.32783975759612945412054678367325547,
    -0.13315519831598209409961309951373512,
    0.28487837717280084052745346456657828,
    0.06942944346252987735848865824703402
};


/* First dotx, then x evolution. No FSAL savings yet. */
void symplectic_Schlier00_6b_nofsal
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx)
{
    Algorithm2_nofsal
        ( Schlier00_6b_K, Schlier00_6b, h, n, x, dotx, acc, ddotx );
    return;
} /* end symplectic_Schlier00_6b_nofsal() */


/* FSAL sequence: _start, _fsal, _stop */
void symplectic_Schlier00_6b_start
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx)
{
    Algorithm2_start
        ( Schlier00_6b_K, Schlier00_6b, h, n, x, dotx, acc, ddotx );
    return;
} /* end symplectic_Schlier00_6b_start() */

void symplectic_Schlier00_6b_fsal
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx)
{
    Algorithm2_fsal
        ( Schlier00_6b_K, Schlier00_6b, h, n, x, dotx, acc, ddotx );
    return;
} /* end symplectic_Schlier00_6b_fsal() */

void symplectic_Schlier00_6b_stop
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx)
{
    Algorithm2_stop
        ( Schlier00_6b_K, Schlier00_6b, h, n, x, dotx, acc, ddotx );
    return;
} /* end symplectic_Schlier00_6b_stop() */
        

/* http://phya8.physik.uni-freiburg.de/abt/papers/papers.html */
#define Schlier00_8b_K   17
/* global truncation error: h^8 */
static const double Schlier00_8b [2*Schlier00_8b_K+1] =
{
    0.03676680389912337302666154929429291,
    0.11072655003739784175754797312279745,
    0.16040429374255560219395381214509780,
    0.61101267825171523627962718607785428,
    -0.00472877643941287918639412436088645,
    -0.19202809069032535396838334049379558,
    0.02983098489335056954884440558763334,
    -0.25979073929811660257162833544861286,
    0.19135844311091097984885756175207225,
    0.38384564066882093754274499421236298,
    -0.03781968145745128677723635761417376,
    0.32661664886778120135972921761872954,
    0.00351845996378093605518443870229385,
    -0.53463443374897025678663398242742174,
    0.13067013867271618676514580608303276,
    -0.39935632081078281354806842349635698,
    -0.01000066638557348147501709158936269,

    0.90721613344495961987012942166888585,

    -0.01000066638557348147501709158936269,
    -0.39935632081078281354806842349635698,
    0.13067013867271618676514580608303276,
    -0.53463443374897025678663398242742174,
    0.00351845996378093605518443870229385,
    0.32661664886778120135972921761872954,
    -0.03781968145745128677723635761417376,
    0.38384564066882093754274499421236298,
    0.19135844311091097984885756175207225,
    -0.25979073929811660257162833544861286,
    0.02983098489335056954884440558763334,
    -0.19202809069032535396838334049379558,
    -0.00472877643941287918639412436088645,
    0.61101267825171523627962718607785428,
    0.16040429374255560219395381214509780,
    0.11072655003739784175754797312279745,
    0.03676680389912337302666154929429291
};

/* Ch. Schlier and A. Seiter, Comp. Phys. Comm. (2000) */
/* First dotx, then x evolution. No FSAL savings yet.  */
void symplectic_Schlier00_8b_nofsal
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx)
{
    Algorithm2_nofsal
        ( Schlier00_8b_K, Schlier00_8b, h, n, x, dotx, acc, ddotx );
    return;
} /* end symplectic_Schlier00_8b_nofsal() */

/* FSAL sequence: _start, _fsal, _stop */
void symplectic_Schlier00_8b_start
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx)
{
    Algorithm2_start
        ( Schlier00_8b_K, Schlier00_8b, h, n, x, dotx, acc, ddotx );
    return;
} /* end symplectic_Schlier00_8b_start() */

void symplectic_Schlier00_8b_fsal
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx)
{
    Algorithm2_fsal
        ( Schlier00_8b_K, Schlier00_8b, h, n, x, dotx, acc, ddotx );
    return;
} /* end symplectic_Schlier00_8b_fsal() */

void symplectic_Schlier00_8b_stop
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx)
{
    Algorithm2_stop
        ( Schlier00_8b_K, Schlier00_8b, h, n, x, dotx, acc, ddotx );
    return;
} /* end symplectic_Schlier00_8b_stop() */


/* http://phya8.physik.uni-freiburg.de/abt/papers/papers.html */
#define Schlier00_8c_K   17
/* global truncation error: h^8 */
static const double Schlier00_8c [2*Schlier00_8c_K+1] =
{
    0.04463795052359022755913999625733590,
    0.13593258071690959145543264213495574,
    0.21988440427147072254445535069606167,
    0.13024946780523828601621193778196846,
    0.10250365693975069608261241007779814,
    0.43234521869358547487983257884877035,
    -0.00477482916916881658022489063962934,
    -0.58253476904040845493112837930861212,
    -0.03886264282111817697737420875189743,
    0.31548728537940479698273603797274199,
    0.18681583743297155471526153503972746,
    0.26500275499062083398346002963079872,
    -0.02405084735747361993573587982407554,
    -0.45040492499772251180922896712151891,
    -0.05897433015592386914575323926766330,
    -0.02168476171861335324934388684707580,
    0.07282080033590128173761892641234244,

    0.55121429634197067334405601381594315,

    0.07282080033590128173761892641234244,
    -0.02168476171861335324934388684707580,
    -0.05897433015592386914575323926766330,
    -0.45040492499772251180922896712151891,
    -0.02405084735747361993573587982407554,
    0.26500275499062083398346002963079872,
    0.18681583743297155471526153503972746,
    0.31548728537940479698273603797274199,
    -0.03886264282111817697737420875189743,
    -0.58253476904040845493112837930861212,
    -0.00477482916916881658022489063962934,
    0.43234521869358547487983257884877035,
    0.10250365693975069608261241007779814,
    0.13024946780523828601621193778196846,
    0.21988440427147072254445535069606167,
    0.13593258071690959145543264213495574,
    0.04463795052359022755913999625733590
};

/* Ch. Schlier and A. Seiter, Comp. Phys. Comm. (2000) */
/* First dotx, then x evolution. No FSAL savings yet.  */
void symplectic_Schlier00_8c_nofsal
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx)
{
    Algorithm2_nofsal
        ( Schlier00_8c_K, Schlier00_8c, h, n, x, dotx, acc, ddotx );
    return;
} /* end symplectic_Schlier00_8c_nofsal() */

/* FSAL sequence: _start, _fsal, _stop */
void symplectic_Schlier00_8c_start
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx)
{
    Algorithm2_start
        ( Schlier00_8c_K, Schlier00_8c, h, n, x, dotx, acc, ddotx );
    return;
} /* end symplectic_Schlier00_8c_start() */

void symplectic_Schlier00_8c_fsal
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx)
{
    Algorithm2_fsal
        ( Schlier00_8c_K, Schlier00_8c, h, n, x, dotx, acc, ddotx );
    return;
} /* end symplectic_Schlier00_8c_fsal() */

void symplectic_Schlier00_8c_stop
(double h, int n, double *x, double *dotx, double *acc,
 Newtonian_ddotx ddotx)
{
    Algorithm2_stop
        ( Schlier00_8c_K, Schlier00_8c, h, n, x, dotx, acc, ddotx );
    return;
} /* end symplectic_Schlier00_8c_stop() */


/**********************************************************************/
/* Routers to various routines by K, the number of evaluations per h. */
/* To make fair comparison, one should use equal g for different K's. */
/* It gives "exact" x(wh) and dotx(wh), w=1..L, after each full step. */
/**********************************************************************/
void symplectic_nofsal ( int K, double g, int n, double *x, double *dotx,
                         double *acc, Newtonian_ddotx ddotx )
{
    switch (K)
    {
        case Ruth83_K:  /* 1 */   /* not much more expensive */
            symplectic_Ruth83_nofsal (K*g, n, x, dotx, acc, ddotx);
            return;
        case Calvo93_K:  /* 4 */
            symplectic_Calvo93_nofsal (K*g, n, x, dotx, acc, ddotx);
            return;
        case Schlier00_6b_K:  /* 8 */
            symplectic_Schlier00_6b_nofsal (K*g, n, x, dotx, acc, ddotx);
            return;
        case Schlier98_6a_K:  /* 9 */   /* not much more expensive */
            symplectic_Schlier98_6a_nofsal (K*g, n, x, dotx, acc, ddotx);
            return;
        case Schlier00_8c_K:  /* 17 */
            symplectic_Schlier00_8c_nofsal (K*g, n, x, dotx, acc, ddotx);
            return;
        case Tsitouras99_K:  /* 33 */   /* not much more expensive */
            symplectic_Tsitouras99_nofsal (K*g, n, x, dotx, acc, ddotx);
            return;
        default:
            pe ("symplectic_nofsal: we have no %d-evaluation integrator.\n",
                K);
    }
    return;
} /* end symplectic_nofsal() */


/**********************************************************************/
/* Routers to various routines by K, the number of evaluations per h. */
/* To make fair comparison, one should use equal g for different K's. */
/* To integrate t from 0 to Lh, one should call 1 x symplectic_start, */
/* (L-1) x symplectic_fsal, and 1 x symplectic_stop. One cannot have  */
/* "exact" x(wh) and dotx(wh) except when w=0,L, but capturing k = 0  */
/* calls to ddotx(k,n,x,acc) at w=1..L-1 gives approximate answers.   */
/**********************************************************************/
void symplectic_start ( int K, double g, int n, double *x, double *dotx,
                        double *acc, Newtonian_ddotx ddotx )
{
    switch (K)
    {
        case Ruth83_K:  /* 1 */   /* not much more expensive */
            symplectic_Ruth83_start (K*g, n, x, dotx, acc, ddotx);
            return;
        case Calvo93_K:  /* 4 */
            symplectic_Calvo93_start (K*g, n, x, dotx, acc, ddotx);
            return;
        case Schlier00_6b_K:  /* 8 */
            symplectic_Schlier00_6b_start (K*g, n, x, dotx, acc, ddotx);
            return;
        case Schlier98_6a_K:  /* 9 : no great need */
            symplectic_Schlier98_6a_start (K*g, n, x, dotx, acc, ddotx);
            return;
        case Schlier00_8c_K:  /* 17 */
            symplectic_Schlier00_8c_start (K*g, n, x, dotx, acc, ddotx);
            return;
        case Tsitouras99_K:  /* 33 : no great need */
            symplectic_Tsitouras99_start (K*g, n, x, dotx, acc, ddotx);
            return;
        default:
            pe ("symplectic_start: we have no %d-evaluation integrator.\n",
                K);
    }
    return;
} /* end symplectic_start() */

void symplectic_fsal ( int K, double g, int n, double *x, double *dotx,
                       double *acc, Newtonian_ddotx ddotx )
{
    switch (K)
    {
        case Ruth83_K:  /* 1 */   /* not much more expensive */
            symplectic_Ruth83_fsal (K*g, n, x, dotx, acc, ddotx);
            return;
        case Calvo93_K:  /* 4 */
            symplectic_Calvo93_fsal (K*g, n, x, dotx, acc, ddotx);
            return;
        case Schlier00_6b_K:  /* 8 */
            symplectic_Schlier00_6b_fsal (K*g, n, x, dotx, acc, ddotx);
            return;
        case Schlier98_6a_K:  /* 9 : no great need */
            symplectic_Schlier98_6a_fsal (K*g, n, x, dotx, acc, ddotx);
            return;
        case Schlier00_8c_K:  /* 17 */
            symplectic_Schlier00_8c_fsal (K*g, n, x, dotx, acc, ddotx);
            return;
        case Tsitouras99_K:  /* 33 : no great need */
            symplectic_Tsitouras99_fsal (K*g, n, x, dotx, acc, ddotx);
            return;
        default:
            pe ("symplectic_fsal: we have no %d-evaluation integrator.\n",
                K);
    }
    return;
} /* end symplectic_fsal() */

void symplectic_stop ( int K, double g, int n, double *x, double *dotx,
                       double *acc, Newtonian_ddotx ddotx )
{
    switch (K)
    {
        case Ruth83_K:  /* 1 */   /* not much more expensive */
            symplectic_Ruth83_stop (K*g, n, x, dotx, acc, ddotx);
            return;
        case Calvo93_K:  /* 4 */
            symplectic_Calvo93_stop (K*g, n, x, dotx, acc, ddotx);
            return;
        case Schlier00_6b_K:  /* 8 */
            symplectic_Schlier00_6b_stop (K*g, n, x, dotx, acc, ddotx);
            return;
        case Schlier98_6a_K:  /* 9 : no great need */
            symplectic_Schlier98_6a_stop (K*g, n, x, dotx, acc, ddotx);
            return;
        case Schlier00_8c_K:  /* 17 */
            symplectic_Schlier00_8c_stop (K*g, n, x, dotx, acc, ddotx);
            return;
        case Tsitouras99_K:  /* 33 : no great need */
            symplectic_Tsitouras99_stop (K*g, n, x, dotx, acc, ddotx);
            return;
        default:
            pe ("symplectic_stop: we have no %d-evaluation integrator.\n",
                K);
    }
    return;
} /* end symplectic_stop() */


/*******************************************************/
/* Well known non-symplectic integrators & Benchmarks: */
/*******************************************************/

/* convert ddotx() in x[n]+dotx[n] to dotX() in X[2n] representation */
static void xdotx_to_X ( int n, double *x, double *dotx, double *X )
{
    register int i;
    for (i=n; i--;) X[i]   = x[i];
    for (i=n; i--;) X[n+i] = dotx[i];
    return;
} /* end xdotx_to_X() */

static void X_to_xdotx ( int n, double *X, double *x, double *dotx )
{
    register int i;
    for (i=n; i--;) x[i]    = X[i];
    for (i=n; i--;) dotx[i] = X[n+i];
    return;
} /* end X_to_xdotx() */

static void ddotx_to_dotX ( int k, int n, double *x, double *dotx,
                            Newtonian_ddotx ddotx, double *dotX )
{
    register int i;
    for (i=n; i--;) dotX[i] = dotx[i];
    (*ddotx) (k, n, x, dotX+n);
    return;
} /* end ddotx_to_dotX() */


/************************************************************/
/* classic 4th-order Runge-Kutta method: Numerical Recipes  */
/* in C: The Art of Scientific Computing, William H. Press, */
/* Saul A. Teukolsky, William T. Vetterling, Brian P.       */
/* Flannery, Cambridge University Press, Cambridge (1992).  */
/* http://mathworld.wolfram.com/Runge-KuttaMethod.html      */
/************************************************************/
static double *dydx = NULL;
static double *dym  = NULL;
static double *dyt  = NULL;
static double *yt   = NULL;

#define rk4_K  4
/* global truncation error: h^4 */

/* allocate 8*n internal workspace. No integration is performed yet. */
void rk4_prep ( double h, int n, double *x, double *dotx, double *acc,
                Newtonian_ddotx ddotx )
{
    REALLOC (rk4, dydx, 2*n, double);
    REALLOC (rk4, dym,  2*n, double);
    REALLOC (rk4, dyt,  2*n, double);
    REALLOC (rk4, yt,   2*n, double);
    return;
} /* end rk4_prep() */

/* "acc" would not be used, so it can be passed as NULL. */
void rk4 (double h, int n, double *x, double *dotx, double *acc,
          Newtonian_ddotx ddotx)
{
    register int i;
    register double hh;
    /* call k=0: suited for heat current computation & total energy checksum */
    ddotx_to_dotX ( 0, n, x, dotx, ddotx, dydx );  /* dydx := k1 / h */
    hh = h / 2;
    for (i=n; i--;) yt[i]   = x[i]    + hh * dydx[i];
    for (i=n; i--;) yt[n+i] = dotx[i] + hh * dydx[n+i];
    ddotx_to_dotX ( 1, n, yt, yt+n, ddotx, dyt );  /* dyt := k2 / h */
    for (i=n; i--;) yt[i]   = x[i]    + hh * dyt[i];
    for (i=n; i--;) yt[n+i] = dotx[i] + hh * dyt[n+i];
    ddotx_to_dotX ( 2, n, yt, yt+n, ddotx, dym );  /* dym := k3 / h */
    for (i=n; i--;) yt[i]   = x[i]    + h * dym[i];
    for (i=n; i--;) yt[n+i] = dotx[i] + h * dym[n+i];
    for (i=2*n; i--;) dym[i] += dyt[i];  /* dym := (k2 + k3) / h */
    ddotx_to_dotX ( 3, n, yt, yt+n, ddotx, dyt );  /* dyt := k4 / h */
    hh = h / 6.0;
    for (i=n; i--;) x[i]    += hh * (dydx[i]   + dyt[i]   + 2.*dym[i]);
    for (i=n; i--;) dotx[i] += hh * (dydx[n+i] + dyt[n+i] + 2.*dym[n+i]);
    return;
} /* end rk4() */

/* free 8*n internal workspace */
void rk4_free ( double h, int n, double *x, double *dotx, double *acc,
                Newtonian_ddotx ddotx )
{
    Free ( yt   );
    Free ( dyt  );
    Free ( dym  );
    Free ( dydx );
    return;
} /* end rk4_free() */


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
static double *x2 = NULL;
static double *x3 = NULL;

/* 4-value Gear predictor-corrector algorithm */
#define gear4_K  1
/* local  truncation error: h^5 */
/* global truncation error: h^4 */

/* Initialize 4-value Gear integrator and allocate 2*n  */
/* internal workspace. No integration is performed yet. */
void gear4_prep ( double h, int n, double *x, double *dotx, double *acc,
                  Newtonian_ddotx ddotx )
{
    register int i;
    double H, *xH, *x_H;
    H = h / GEAR_TEST_H;
    /* integrate backward: error = H^3 */
    MALLOC ( gear4_prep, x_H, 3*n, double );
    xdotx_to_X ( n, x, dotx, x_H );
    symplectic_Ruth83_nofsal ( -H, n, x_H, x_H+n, acc, ddotx );
    /* x(-H),dotx(-H),ddotx(-H) + O(H^3) */
    (*ddotx) ( 0, n, x_H, x_H+2*n );
    /* integrate forward: error = H^3 */
    MALLOC ( gear4_prep, xH, 3*n, double );
    xdotx_to_X ( n, x, dotx, xH );
    symplectic_Ruth83_nofsal ( H, n, xH, xH+n, acc, ddotx );
    /* x(H),dotx(H),ddotx(H) + O(H^3) */
    (*ddotx) ( 0, n, xH, xH+2*n );
    REALLOC ( gear4_prep, x2, n, double );
    (*ddotx) ( 0, n, x, x2 );  /* x(0),dotx(0),ddotx(0) */
    REALLOC ( gear4_prep, x3, n, double );
    for (i=n; i--;)
    { /* This 1st derivative has H^2 truncation error and H^2 int. error */
        x3[i] = ( xH[2*n+i] - x_H[2*n+i] ) / 2. / H
            * pow ( h, 3. ) / factorial(3);
        /* after h^3 -> h^5, which is equal to the Gear4 LTE h^5. */
        x2[i] *= pow ( h, 2. ) / factorial(2);
    }
    free (xH);
    free (x_H);
    return;
} /* end gear4_prep() */


#define gear4_0   ( 1. / 6. )
#define gear4_1   ( 5. / 6. )
#define gear4_3   ( 1. / 3. )
void gear4 ( double h, int n, double *x, double *dotx, double *acc,
             Newtonian_ddotx ddotx )
{
    register int i;
    register double correction;
    for (i=n; i--;)
    {
        /* printf ("%.15e %.15e %.15e %.15e %.15e\n", x[i], dotx[i], */
        /* x2[i]/pow(h,2.), x3[i]/pow(h,3.), x4[i]/pow(h,4.)); */
        x[i]    +=   dotx[i] * h + x2[i] + x3[i];
        dotx[i] += ( 2.*x2[i] + 3.*x3[i] ) / h;
        x2[i]   +=   3.*x3[i];
    }
    /* call k=0: suited for heat current & total energy checksum */
    (*ddotx) ( 0, n, x, acc );
    for (i=n; i--;)
    {
        correction = x2[i] - acc[i] * h * h / 2.;
        x[i]    -= correction * gear4_0;
        dotx[i] -= correction * gear4_1 / h;
        x2[i]   -= correction;
        x3[i]   -= correction * gear4_3;
    }
    return;
} /* end gear4() */

/* free 2*n internal workspace */
void gear4_free ( double h, int n, double *x, double *dotx, double *acc,
                  Newtonian_ddotx ddotx )
{
    Free( x3 );
    Free( x2 );
    return;
} /* end gear4_free() */


/******************************************************/
/* 5-value Gear predictor-corrector algorithm:        */
/* Computer Simulation of Liquids, M.P. Allen and     */
/* D.J. Tildesley, Clarendon Press, (Oxford 1987).    */
/* ftp://ftp.dl.ac.uk/ccp5/ALLEN_TILDESLEY/F.02       */
/* itp.nat.uni-magdeburg.de/~schinner/dip/node22.html */
/******************************************************/

static double *x4 = NULL;

#define gear5_K  1
/* local  truncation error: h^6 */
/* global truncation error: h^5 */

/* Initialize 5-value Gear integrator and allocate 3*n  */
/* internal workspace. No integration is performed yet. */
void gear5_prep ( double h, int n, double *x, double *dotx, double *acc,
                  Newtonian_ddotx ddotx )
{
    register int i;
    double H, *xH, *x2H, *x_H, *x_2H;
    H = h / GEAR_TEST_H;
    /* integrate backward: error = H^5 */
    MALLOC ( gear5_prep, x_H, 3*n, double );
    xdotx_to_X ( n, x, dotx, x_H );
    symplectic_Calvo93_nofsal ( -H, n, x_H, x_H+n, acc, ddotx );
    /* x(-H),dotx(-H),ddotx(-H) + O(H^5) */
    (*ddotx) ( 0, n, x_H, x_H+2*n );
    MALLOC ( gear5_prep, x_2H, 3*n, double );
    for (i=2*n; i--;) x_2H[i] = x_H[i];
    symplectic_Calvo93_nofsal ( -H, n, x_2H, x_2H+n, acc, ddotx );
    /* x(-2H),dotx(-2H),ddotx(-2H) + O(H^5) */
    (*ddotx) ( 0, n, x_2H, x_2H+2*n );
    /* integrate forward: error = H^5 */
    MALLOC ( gear5_prep, xH, 3*n, double );
    xdotx_to_X ( n, x, dotx, xH );
    symplectic_Calvo93_nofsal ( H, n, xH, xH+n, acc, ddotx );
    /* x(H),dotx(H),ddotx(H) + O(H^5) */
    (*ddotx) ( 0, n, xH, xH+2*n );
    MALLOC ( gear5_prep, x2H, 3*n, double );
    for (i=2*n; i--;) x2H[i] = xH[i];
    symplectic_Calvo93_nofsal ( H, n, x2H, x2H+n, acc, ddotx );
    /* x(2H),dotx(2H),ddotx(2H) + O(H^5) */
    (*ddotx) ( 0, n, x2H, x2H+2*n );
    REALLOC ( gear5_prep, x2, n, double );
    (*ddotx) ( 0, n, x, x2 );  /* x(0),dotx(0),ddotx(0) */
    REALLOC ( gear5_prep, x3, n, double );
    REALLOC ( gear5_prep, x4, n, double );
    /* www.orst.edu/instruct/ch490/lessons/lesson11.htm */
    for (i=n; i--;)
    { /* This 1st derivative has H^4 truncation error and H^4 int. error */
        x3[i] = ( 8 * ( xH[2*n+i] - x_H[2*n+i] )
                  -  ( x2H[2*n+i] - x_2H[2*n+i] ) ) / 12. / H
            * pow ( h, 3. ) / factorial(3);
        /* after h^3 -> h^7, which is smaller than Gear LTE h^6. */
        /* This 2nd derivative has H^4 truncation error and H^3 int. error */
        x4[i] = ( 16 * ( (xH[2*n+i] - x2[i]) + (x_H[2*n+i] - x2[i]) )
                  -    ( (x2H[2*n+i] - x2[i]) + (x_2H[2*n+i] - x2[i]) ) )
            / 12 / H / H * pow ( h, 4. ) / factorial(4);
        /* after h^4 -> h^7, which is smaller than Gear LTE h^6. */
        /* !!! */
        /* x4[i] = 0; */
        x2[i] *= pow ( h, 2. ) / factorial(2);
    }
    free (x2H);
    free (xH);
    free (x_2H);
    free (x_H);
    return;
} /* end gear5_prep() */


#define gear5_0   ( 19. / 120. )
#define gear5_1   ( 3.  / 4.   )
#define gear5_3   ( 1.  / 2.   )
#define gear5_4   ( 1.  / 12.  )
void gear5 ( double h, int n, double *x, double *dotx, double *acc,
             Newtonian_ddotx ddotx )
{
    register int i;
    register double correction;
    for (i=n; i--;)
    { /* compare with Notes/Kepler.mws */
        /* printf ("%.15e %.15e %.15e %.15e %.15e\n", x[i], dotx[i], */
        /* x2[i]/pow(h,2.), x3[i]/pow(h,3.), x4[i]/pow(h,4.)); */
        x[i]    +=   dotx[i] * h + x2[i] + x3[i] + x4[i];
        dotx[i] += ( 2.*x2[i] + 3.*x3[i] + 4.*x4[i] ) / h;
        x2[i]   +=   3.*x3[i] + 6.*x4[i];
        x3[i]   +=   4.*x4[i];
    }
    /* call k=0: suited for heat current & total energy checksum */
    (*ddotx) ( 0, n, x, acc );
    for (i=n; i--;)
    {
        correction = x2[i] - acc[i] * h * h / 2.;
        x[i]    -= correction * gear5_0;
        dotx[i] -= correction * gear5_1 / h;
        x2[i]   -= correction;
        x3[i]   -= correction * gear5_3;
        x4[i]   -= correction * gear5_4;
    }
    return;
} /* end gear5() */

/* free 3*n internal workspace */
void gear5_free ( double h, int n, double *x, double *dotx, double *acc,
                  Newtonian_ddotx ddotx )
{
    Free( x4 );
    Free( x3 );
    Free( x2 );
    return;
} /* end gear5_free() */


/* 6-value Gear predictor-corrector algorithm */
static double *x5 = NULL;

#define gear6_K  1
/* local  truncation error: h^7 */
/* global truncation error: h^6 */

/* Initialize 6-value Gear integrator and allocate 4*n  */
/* internal workspace. No integration is performed yet. */
void gear6_prep ( double h, int n, double *x, double *dotx, double *acc,
                  Newtonian_ddotx ddotx )
{
    register int i;
    double H, *xH, *x2H, *x_H, *x_2H;
    H = h / GEAR_TEST_H;
    /* integrate backward: error = H^5 */
    MALLOC ( gear6_prep, x_H, 3*n, double );
    xdotx_to_X ( n, x, dotx, x_H );
    symplectic_Calvo93_nofsal ( -H, n, x_H, x_H+n, acc, ddotx );
    /* x(-H),dotx(-H),ddotx(-H) + O(H^5) */
    (*ddotx) ( 0, n, x_H, x_H+2*n );
    MALLOC ( gear6_prep, x_2H, 3*n, double );
    for (i=2*n; i--;) x_2H[i] = x_H[i];
    symplectic_Calvo93_nofsal ( -H, n, x_2H, x_2H+n, acc, ddotx );
    /* x(-2H),dotx(-2H),ddotx(-2H) + O(H^5) */
    (*ddotx) ( 0, n, x_2H, x_2H+2*n );
    /* integrate forward: error = H^5 */
    MALLOC ( gear6_prep, xH, 3*n, double );
    xdotx_to_X ( n, x, dotx, xH );
    symplectic_Calvo93_nofsal ( H, n, xH, xH+n, acc, ddotx );
    /* x(H),dotx(H),ddotx(H) + O(H^5) */
    (*ddotx) ( 0, n, xH, xH+2*n );
    MALLOC ( gear6_prep, x2H, 3*n, double );
    for (i=2*n; i--;) x2H[i] = xH[i];
    symplectic_Calvo93_nofsal ( H, n, x2H, x2H+n, acc, ddotx );
    /* x(2H),dotx(2H),ddotx(2H) + O(H^5) */
    (*ddotx) ( 0, n, x2H, x2H+2*n );
    REALLOC ( gear6_prep, x2, n, double );
    (*ddotx) ( 0, n, x, x2 );  /* x(0),dotx(0),ddotx(0) */
    REALLOC ( gear6_prep, x3, n, double );
    REALLOC ( gear6_prep, x4, n, double );
    REALLOC ( gear6_prep, x5, n, double );
    /* www.orst.edu/instruct/ch490/lessons/lesson11.htm */
    for (i=n; i--;)
    { /* This 1st derivative has H^4 truncation error and H^4 int. error */
        x3[i] = ( 8 * ( xH[2*n+i] - x_H[2*n+i] )
                  -   ( x2H[2*n+i] - x_2H[2*n+i] ) ) / 12. / H
            * pow ( h, 3. ) / factorial(3);
        /* after h^3 -> h^7, which is equal to Gear6 LTE h^7. */
        /* This 2nd derivative has H^4 truncation error and H^3 int. error */
        x4[i] = ( 16 * ( (xH[2*n+i] - x2[i]) + (x_H[2*n+i] - x2[i]) )
                  -    ( (x2H[2*n+i] - x2[i]) + (x_2H[2*n+i] - x2[i]) ) )
            / 12 / H / H * pow ( h, 4. ) / factorial(4);
        /* after h^4 -> h^7, which is equal to Gear6 LTE h^7. */
        /* This 3rd derivative has H^2 truncation error and H^2 int. error */
        x5[i] = ( ( x2H[2*n+i] - x_2H[2*n+i] ) - 2 *
                  ( xH[2*n+i] - x_H[2*n+i] ) ) / 2 / H / H / H
            * pow ( h, 5. ) / factorial(5);
        /* after h^5 -> h^7, which is equal to Gear6 LTE h^7. */
        x2[i] *= pow ( h, 2. ) / factorial(2);
    }
    free (x2H);
    free (xH);
    free (x_2H);
    free (x_H);
    return;
} /* end gear6_prep() */


#define gear6_02  (   3. /  20. )
#define gear6_12  ( 251. / 360. )
#define gear6_32  (  11. /  18. )
#define gear6_42  (   1. /   6. )
#define gear6_52  (   1. /  60. )
void gear6 ( double h, int n, double *x, double *dotx, double *acc,
             Newtonian_ddotx ddotx )
{
    register int i;
    register double correction;
    for (i=n; i--;)
    {
        x[i]    +=   dotx[i] * h + x2[i] + x3[i] + x4[i] + x5[i];
        dotx[i] += ( 2.*x2[i] + 3.*x3[i] + 4.*x4[i] + 5.*x5[i] ) / h;
        x2[i]   +=   3.*x3[i] + 6.*x4[i] + 10.*x5[i];
        x3[i]   +=   4.*x4[i] + 10.*x5[i];
        x4[i]   +=   5.*x5[i];
    }
    /* call k=0: suited for heat current & total energy checksum */
    (*ddotx) ( 0, n, x, acc );
    for (i=n; i--;)
    {
        correction = x2[i] - acc[i] * h * h / 2.;
        x[i]    -= correction * gear6_02;
        dotx[i] -= correction * gear6_12 / h;
        x2[i]   -= correction;
        x3[i]   -= correction * gear6_32;
        x4[i]   -= correction * gear6_42;
        x5[i]   -= correction * gear6_52;
    }
    return;
} /* end gear6() */

/* free 4*n internal workspace */
void gear6_free ( double h, int n, double *x, double *dotx, double *acc,
                  Newtonian_ddotx ddotx )
{
    Free( x5 );
    Free( x4 );
    Free( x3 );
    Free( x2 );
    return;
} /* end gear6_free() */


/* 7-value Gear predictor-corrector algorithm */
static double *x6 = NULL;

#define gear7_K  1
/* local  truncation error: h^8 */
/* global truncation error: h^7 */

/* Initialize 7-value Gear integrator and allocate 5*n  */
/* internal workspace. No integration is performed yet. */
void gear7_prep ( double h, int n, double *x, double *dotx, double *acc,
                  Newtonian_ddotx ddotx )
{
    register int i;
    double H, *xH, *x2H, *x3H, *x_H, *x_2H, *x_3H;
    H = h / GEAR_TEST_H;
    /* integrate backward: error = H^7 */
    MALLOC ( gear8_prep, x_H, 3*n, double );
    xdotx_to_X ( n, x, dotx, x_H );
    symplectic_Schlier00_6b_nofsal ( -H, n, x_H, x_H+n, acc, ddotx );
    /* x(-H),dotx(-H),ddotx(-H) + O(H^7) */
    (*ddotx) ( 0, n, x_H, x_H+2*n );
    MALLOC ( gear8_prep, x_2H, 3*n, double );
    for (i=2*n; i--;) x_2H[i] = x_H[i];
    symplectic_Schlier00_6b_nofsal ( -H, n, x_2H, x_2H+n, acc, ddotx );
    /* x(-2H),dotx(-2H),ddotx(-2H) + O(H^7) */
    (*ddotx) ( 0, n, x_2H, x_2H+2*n );
    MALLOC ( gear8_prep, x_3H, 3*n, double );
    for (i=2*n; i--;) x_3H[i] = x_2H[i];
    symplectic_Schlier00_6b_nofsal ( -H, n, x_3H, x_3H+n, acc, ddotx );
    /* x(-3H),dotx(-3H),ddotx(-3H) + O(H^7) */
    (*ddotx) ( 0, n, x_3H, x_3H+2*n );
    /* integrate forward: error = H^7 */
    MALLOC ( gear8_prep, xH, 3*n, double );
    xdotx_to_X ( n, x, dotx, xH );
    symplectic_Schlier00_6b_nofsal ( H, n, xH, xH+n, acc, ddotx );
    /* x(H),dotx(H),ddotx(H) + O(H^7) */
    (*ddotx) ( 0, n, xH, xH+2*n );
    MALLOC ( gear8_prep, x2H, 3*n, double );
    for (i=2*n; i--;) x2H[i] = xH[i];
    symplectic_Schlier00_6b_nofsal ( H, n, x2H, x2H+n, acc, ddotx );
    /* x(2H),dotx(2H),ddotx(2H) + O(H^7) */
    (*ddotx) ( 0, n, x2H, x2H+2*n );
    MALLOC ( gear8_prep, x3H, 3*n, double );
    for (i=2*n; i--;) x3H[i] = x2H[i];
    symplectic_Schlier00_6b_nofsal ( H, n, x3H, x3H+n, acc, ddotx );
    /* x(3H),dotx(3H),ddotx(3H) + O(H^7) */
    (*ddotx) ( 0, n, x3H, x3H+2*n );
    REALLOC ( gear7_prep, x2, n, double );
    (*ddotx) ( 0, n, x, x2 );  /* x(0),dotx(0),ddotx(0) */
    REALLOC ( gear7_prep, x3, n, double );
    REALLOC ( gear7_prep, x4, n, double );
    REALLOC ( gear7_prep, x5, n, double );
    REALLOC ( gear7_prep, x6, n, double );
    /* numeric-tools.hypermart.net/docs.htm, or Notes/central_difference.m */
    for (i=n; i--;)
    { /* This 1st derivative has H^6 truncation error and H^6 int. error */
        x3[i] = ( 45 * ( xH[2*n+i] - x_H[2*n+i] ) -
                  9  * ( x2H[2*n+i] - x_2H[2*n+i] ) +
                  ( x3H[2*n+i] - x_3H[2*n+i] ) ) / 60. / H
            * pow ( h, 3. ) / factorial(3);
        /* after h^3 -> h^9, which is smaller than Gear7 LTE h^8. */
        /* This 2nd derivative has H^6 truncation error and H^5 int. error */
        x4[i] = ( 270 * ( (xH[2*n+i] - x2[i]) + (x_H[2*n+i] - x2[i]) ) -
                  27 * ( (x2H[2*n+i] - x2[i]) + (x_2H[2*n+i] - x2[i]) ) +
                  2  * ( (x3H[2*n+i] - x2[i]) + (x_3H[2*n+i] - x2[i]) ) )
            / 180 / H / H * pow ( h, 4. ) / factorial(4);
        /* after h^4 -> h^9, which is smaller than Gear7 LTE h^8. */
        /* This 3rd derivative has H^4 truncation error and H^4 int. error */
        x5[i] = ( -13 * ( xH[2*n+i] - x_H[2*n+i] ) +
                  8  * ( x2H[2*n+i] - x_2H[2*n+i] ) -
                  ( x3H[2*n+i] - x_3H[2*n+i] ) ) / 8. / H / H / H
            * pow ( h, 5. ) / factorial(5);
        /* after h^5 -> h^9, which is smaller than Gear7 LTE h^8. */
        /* This 4th derivative has H^4 truncation error and H^3 int. error */
        x6[i] = ( -39 * ( (xH[2*n+i] - x2[i]) + (x_H[2*n+i] - x2[i]) ) +
                  12 * ( (x2H[2*n+i] - x2[i]) + (x_2H[2*n+i] - x2[i]) ) -
                  ( (x3H[2*n+i] - x2[i]) + (x_3H[2*n+i] - x2[i]) ) )
            / 6 / H / H / H / H * pow ( h, 6. ) / factorial(6);
        /* after h^6 -> h^9, which is smaller than Gear7 LTE h^8. */
        /* !!! */
        /* x6[i] = 0; */
        x2[i] *= pow ( h, 2. ) / factorial(2);
    }
    free (x3H);
    free (x2H);
    free (xH);
    free (x_3H);
    free (x_2H);
    free (x_H);
    return;
} /* end gear7_prep() */


#define gear7_02  ( 863. / 6048. )
#define gear7_12  ( 665. / 1008. )
#define gear7_32  (  25. /   36. )
#define gear7_42  (  35. /  144. )
#define gear7_52  (   1. /   24. )
#define gear7_62  (   1. /  360. )
void gear7 ( double h, int n, double *x, double *dotx, double *acc,
             Newtonian_ddotx ddotx )
{
    register int i;
    register double correction;
    for (i=n; i--;)
    {
        x[i]    +=   dotx[i] * h + x2[i] + x3[i] + x4[i] + x5[i] + x6[i];
        dotx[i] += ( 2.*x2[i] + 3.*x3[i] + 4.*x4[i] + 5.*x5[i] +
                     6.*x6[i] ) / h;
        x2[i]   +=   3.*x3[i] + 6.*x4[i] + 10.*x5[i] + 15.*x6[i];
        x3[i]   +=   4.*x4[i] + 10.*x5[i] + 20.*x6[i];
        x4[i]   +=   5.*x5[i] + 15.*x6[i];
        x5[i]   +=   6.*x6[i];
    }
    /* call k=0: suited for heat current & total energy checksum */
    (*ddotx) ( 0, n, x, acc );
    for (i=n; i--;)
    {
        correction = x2[i] - acc[i] * h * h / 2.;
        x[i]    -= correction * gear7_02;
        dotx[i] -= correction * gear7_12 / h;
        x2[i]   -= correction;
        x3[i]   -= correction * gear7_32;
        x4[i]   -= correction * gear7_42;
        x5[i]   -= correction * gear7_52;
        x6[i]   -= correction * gear7_62;
    }
    return;
} /* end gear7() */

/* free 5*n internal workspace */
void gear7_free ( double h, int n, double *x, double *dotx, double *acc,
                  Newtonian_ddotx ddotx )
{
    Free( x6 );
    Free( x5 );
    Free( x4 );
    Free( x3 );
    Free( x2 );
    return;
} /* end gear7_free() */


/* 8-value Gear predictor-corrector algorithm */
static double *x7 = NULL;

#define gear8_K  1
/* local  truncation error: h^9 */
/* global truncation error: h^8 */

/* Initialize 8-value Gear integrator and allocate 6*n  */
/* internal workspace. No integration is performed yet. */
void gear8_prep ( double h, int n, double *x, double *dotx, double *acc,
                  Newtonian_ddotx ddotx )
{
    register int i;
    double H, *xH, *x2H, *x3H, *x_H, *x_2H, *x_3H;
    H = h / GEAR_TEST_H;
    /* integrate backward: error = H^7 */
    MALLOC ( gear8_prep, x_H, 3*n, double );
    xdotx_to_X ( n, x, dotx, x_H );
    symplectic_Schlier00_6b_nofsal ( -H, n, x_H, x_H+n, acc, ddotx );
    /* x(-H),dotx(-H),ddotx(-H) + O(H^7) */
    (*ddotx) ( 0, n, x_H, x_H+2*n );
    MALLOC ( gear8_prep, x_2H, 3*n, double );
    for (i=2*n; i--;) x_2H[i] = x_H[i];
    symplectic_Schlier00_6b_nofsal ( -H, n, x_2H, x_2H+n, acc, ddotx );
    /* x(-2H),dotx(-2H),ddotx(-2H) + O(H^7) */
    (*ddotx) ( 0, n, x_2H, x_2H+2*n );
    MALLOC ( gear8_prep, x_3H, 3*n, double );
    for (i=2*n; i--;) x_3H[i] = x_2H[i];
    symplectic_Schlier00_6b_nofsal ( -H, n, x_3H, x_3H+n, acc, ddotx );
    /* x(-3H),dotx(-3H),ddotx(-3H) + O(H^7) */
    (*ddotx) ( 0, n, x_3H, x_3H+2*n );
    /* integrate forward: error = H^7 */
    MALLOC ( gear8_prep, xH, 3*n, double );
    xdotx_to_X ( n, x, dotx, xH );
    symplectic_Schlier00_6b_nofsal ( H, n, xH, xH+n, acc, ddotx );
    /* x(H),dotx(H),ddotx(H) + O(H^7) */
    (*ddotx) ( 0, n, xH, xH+2*n );
    MALLOC ( gear8_prep, x2H, 3*n, double );
    for (i=2*n; i--;) x2H[i] = xH[i];
    symplectic_Schlier00_6b_nofsal ( H, n, x2H, x2H+n, acc, ddotx );
    /* x(2H),dotx(2H),ddotx(2H) + O(H^7) */
    (*ddotx) ( 0, n, x2H, x2H+2*n );
    MALLOC ( gear8_prep, x3H, 3*n, double );
    for (i=2*n; i--;) x3H[i] = x2H[i];
    symplectic_Schlier00_6b_nofsal ( H, n, x3H, x3H+n, acc, ddotx );
    /* x(3H),dotx(3H),ddotx(3H) + O(H^7) */
    (*ddotx) ( 0, n, x3H, x3H+2*n );
    REALLOC ( gear8_prep, x2, n, double );
    (*ddotx) ( 0, n, x, x2 );  /* x(0),dotx(0),ddotx(0) */
    REALLOC ( gear8_prep, x3, n, double );
    REALLOC ( gear8_prep, x4, n, double );
    REALLOC ( gear8_prep, x5, n, double );
    REALLOC ( gear8_prep, x6, n, double );
    REALLOC ( gear8_prep, x7, n, double );
    /* numeric-tools.hypermart.net/docs.htm, or Notes/central_difference.m */
    for (i=n; i--;)
    { /* This 1st derivative has H^6 truncation error and H^6 int. error */
        x3[i] = ( 45 * ( xH[2*n+i] - x_H[2*n+i] ) -
                  9  * ( x2H[2*n+i] - x_2H[2*n+i] ) +
                  ( x3H[2*n+i] - x_3H[2*n+i] ) ) / 60. / H
            * pow ( h, 3. ) / factorial(3);
        /* after h^3 -> h^9, which is equal to Gear8 LTE h^9. */
        /* This 2nd derivative has H^6 truncation error and H^5 int. error */
        x4[i] = ( 270 * ( (xH[2*n+i] - x2[i]) + (x_H[2*n+i] - x2[i]) ) -
                  27 * ( (x2H[2*n+i] - x2[i]) + (x_2H[2*n+i] - x2[i]) ) +
                  2  * ( (x3H[2*n+i] - x2[i]) + (x_3H[2*n+i] - x2[i]) ) )
            / 180 / H / H * pow ( h, 4. ) / factorial(4);
        /* after h^4 -> h^9, which is equal to Gear8 LTE h^9. */
        /* This 3rd derivative has H^4 truncation error and H^4 int. error */
        x5[i] = ( -13 * ( xH[2*n+i] - x_H[2*n+i] ) +
                  8  * ( x2H[2*n+i] - x_2H[2*n+i] ) -
                  ( x3H[2*n+i] - x_3H[2*n+i] ) ) / 8. / H / H / H
            * pow ( h, 5. ) / factorial(5);
        /* after h^5 -> h^9, which is equal to Gear8 LTE h^9. */
        /* This 4th derivative has H^4 truncation error and H^3 int. error */
        x6[i] = ( -39 * ( (xH[2*n+i] - x2[i]) + (x_H[2*n+i] - x2[i]) ) +
                  12 * ( (x2H[2*n+i] - x2[i]) + (x_2H[2*n+i] - x2[i]) ) -
                  ( (x3H[2*n+i] - x2[i]) + (x_3H[2*n+i] - x2[i]) ) )
            / 6 / H / H / H / H * pow ( h, 6. ) / factorial(6);
        /* after h^6 -> h^9, which is equal to Gear8 LTE h^9. */
        /* This 5th derivative has H^2 truncation error and H^2 int. error */
        x7[i] = ( 5 * ( xH[2*n+i] - x_H[2*n+i] ) -
                  4 * ( x2H[2*n+i] - x_2H[2*n+i] ) +
                  ( x3H[2*n+i] - x_3H[2*n+i] ) ) / 2. / H / H / H / H / H
            * pow ( h, 7. ) / factorial(7);
        /* after h^7 -> h^9, which is equal to Gear8 LTE h^9. */
        x2[i] *= pow ( h, 2. ) / factorial(2);
    }
    free (x3H);
    free (x2H);
    free (xH);
    free (x_3H);
    free (x_2H);
    free (x_H);
    return;
} /* end gear8_prep() */


#define gear8_02  (  1925. / 14112. )
#define gear8_12  ( 19087. / 30240. )
#define gear8_32  (   137. / 180.   )
#define gear8_42  (     5. / 16.    )
#define gear8_52  (    17. / 240.   )
#define gear8_62  (     1. / 120.   )
#define gear8_72  (     1. / 2520.  )
void gear8 ( double h, int n, double *x, double *dotx, double *acc,
             Newtonian_ddotx ddotx )
{
    register int i;
    register double correction;
    for (i=n; i--;)
    {
        x[i]    +=   dotx[i] * h + x2[i] + x3[i] + x4[i] + x5[i] + x6[i] +
            x7[i];
        dotx[i] += ( 2.*x2[i] + 3.*x3[i] + 4.*x4[i] + 5.*x5[i] +
                     6.*x6[i] + 7.*x7[i] ) / h;
        x2[i]   +=   3.*x3[i] + 6.*x4[i] + 10.*x5[i] + 15.*x6[i] + 21.*x7[i];
        x3[i]   +=   4.*x4[i] + 10.*x5[i] + 20.*x6[i] + 35.*x7[i];
        x4[i]   +=   5.*x5[i] + 15.*x6[i] + 35.*x7[i];
        x5[i]   +=   6.*x6[i] + 21.*x7[i];
        x6[i]   +=   7.*x7[i];
    }
    /* call k=0: suited for heat current & total energy checksum */
    (*ddotx) ( 0, n, x, acc );
    for (i=n; i--;)
    {
        correction = x2[i] - acc[i] * h * h / 2.;
        x[i]    -= correction * gear8_02;
        dotx[i] -= correction * gear8_12 / h;
        x2[i]   -= correction;
        x3[i]   -= correction * gear8_32;
        x4[i]   -= correction * gear8_42;
        x5[i]   -= correction * gear8_52;
        x6[i]   -= correction * gear8_62;
        x7[i]   -= correction * gear8_72;
    }
    return;
} /* end gear8() */

/* free 6*n internal workspace */
void gear8_free ( double h, int n, double *x, double *dotx, double *acc,
                  Newtonian_ddotx ddotx )
{
    Free( x7 );
    Free( x6 );
    Free( x5 );
    Free( x4 );
    Free( x3 );
    Free( x2 );
    return;
} /* end gear8_free() */


Newtonian newtonian [NEWTONIAN_MAX] =
{
    { "Ruth83",        TRUE,   Ruth83_K,         2,
      &symplectic_Ruth83_start,
      &symplectic_Ruth83_fsal,
      &symplectic_Ruth83_stop },
    { "Schlier98_6a",  TRUE,   Schlier98_6a_K,   6,
      &symplectic_Schlier98_6a_start,
      &symplectic_Schlier98_6a_fsal,
      &symplectic_Schlier98_6a_stop },
    { "Tsitouras99",   TRUE,  Tsitouras99_K,    10,
      &symplectic_Tsitouras99_start,
      &symplectic_Tsitouras99_fsal,
      &symplectic_Tsitouras99_stop },

    { "Calvo93",       TRUE,   Calvo93_K,        4,
      &symplectic_Calvo93_start,
      &symplectic_Calvo93_fsal,
      &symplectic_Calvo93_stop },
    { "Schlier00_6b",  TRUE,   Schlier00_6b_K,   6,
      &symplectic_Schlier00_6b_start,
      &symplectic_Schlier00_6b_fsal,
      &symplectic_Schlier00_6b_stop },
    { "Schlier00_8b",  TRUE,  Schlier00_8b_K,    8,
      &symplectic_Schlier00_8b_start,
      &symplectic_Schlier00_8b_fsal,
      &symplectic_Schlier00_8b_stop },
    { "Schlier00_8c",  TRUE,  Schlier00_8c_K,    8,
      &symplectic_Schlier00_8c_start,
      &symplectic_Schlier00_8c_fsal,
      &symplectic_Schlier00_8c_stop },

    { "rk4",          FALSE,   rk4_K,            4,
      &rk4_prep,   &rk4,   &rk4_free },
    { "gear4",        FALSE,   gear4_K,          4,
      &gear4_prep, &gear4, &gear4_free },
    { "gear5",        FALSE,   gear5_K,          5,
      &gear5_prep, &gear5, &gear5_free },
    { "gear6",        FALSE,   gear6_K,          6,
      &gear6_prep, &gear6, &gear6_free },
    { "gear7",        FALSE,   gear7_K,          7,
      &gear7_prep, &gear7, &gear7_free },
    { "gear8",        FALSE,   gear8_K,          8,
      &gear8_prep, &gear8, &gear8_free }
};


/***************************************************************/
/* prep  +  [ firststep, L-1 nextstep, halt ] repeat  +  free  */
/*            |                                                */
/* k=0 property calculation for algorithm 2 and non-symplectic */
/***************************************************************/

/* Allocate and initialize states if necessary, but do NOT integrate. */
void newtonian_integrator_prep
( int id, double g, int n, double *x, double *dotx,
  double *acc, Newtonian_ddotx ddotx )
{
    if ( OUW( id, NEWTONIAN_MAX ) )
        pe ("newtonian_integrator_prep: "
            "there is no integrator of id = %d.\n", id);
    if (! newtonian[id].is_symplectic)
        (*newtonian[id].a)
            (newtonian[id].K*g, n, x, dotx, acc, ddotx);
    return;
} /* end newtonian_integrator_prep() */

/* Get the integrator to run its first full step of a session. */
void newtonian_integrator_firststep
( int id, double g, int n, double *x, double *dotx,
  double *acc, Newtonian_ddotx ddotx )
{
    if ( OUW( id, NEWTONIAN_MAX ) )
        pe ("newtonian_integrator_firststep: "
            "there is no integrator of id = %d.\n", id);
    if (newtonian[id].is_symplectic)
        (*newtonian[id].a)
            (newtonian[id].K*g, n, x, dotx, acc, ddotx);
    else (*newtonian[id].b)
             (newtonian[id].K*g, n, x, dotx, acc, ddotx);
    return;
} /* end newtonian_integrator_firststep() */

/* Run other full steps of a session. */
void newtonian_integrator_nextstep
( int id, double g, int n, double *x, double *dotx,
  double *acc, Newtonian_ddotx ddotx )
{
    if ( OUW( id, NEWTONIAN_MAX ) )
        pe ("newtonian_integrator_nextstep: "
            "there is no integrator of id = %d.\n", id);
    (*newtonian[id].b)
        (newtonian[id].K*g, n, x, dotx, acc, ddotx);
    return;
} /* end newtonian_integrator_nextstep() */

/* Halt the session and do necessary mop-ups. To */
/* start a new session, call _firststep again.   */
void newtonian_integrator_halt
( int id, double g, int n, double *x, double *dotx,
  double *acc, Newtonian_ddotx ddotx )
{
    if ( OUW( id, NEWTONIAN_MAX ) )
        pe ("newtonian_integrator_halt: "
            "there is no integrator of id = %d.\n", id);
    if (newtonian[id].is_symplectic)
        (*newtonian[id].c)
            (newtonian[id].K*g, n, x, dotx, acc, ddotx);
    return;
} /* end newtonian_integrator_halt() */

/* Free workspace if any. No more sessions. To restart, call _prep. */
void newtonian_integrator_free
( int id, double g, int n, double *x, double *dotx,
  double *acc, Newtonian_ddotx ddotx )
{
    if ( OUW( id, NEWTONIAN_MAX ) )
        pe ("newtonian_integrator_free: "
            "there is no integrator of id = %d.\n", id);
    if (! newtonian[id].is_symplectic)
        (*newtonian[id].c)
            (newtonian[id].K*g, n, x, dotx, acc, ddotx);
    return;
} /* end newtonian_integrator_free() */


#ifdef _Kepler_TEST
double energy (int n, double *x, double *dotx)
{
    return ( ( SQUARE(dotx[0]) + SQUARE(dotx[1]) ) / 2 -
             1 / sqrt( SQUARE(x[0]) + SQUARE(x[1]) ) );
} /* end energy() */

void ddotx (int k, int n, double *x, double *acc)
{
    register double r3;
    r3 = pow ( x[0] * x[0] + x[1] * x[1], 1.5 );
    acc[0] = - x[0] / r3;
    acc[1] = - x[1] / r3;
    return;
} /* end ddotx() */

#define n     2
#define e     0.5
#define from  0
#define to   (NEWTONIAN_MAX-1)
int main (int argc, char *argv[])
{
    int id, j, L, periods;
    double xi[2*n] = { 1.-e, 0, 0, sqrt((1.+e)/(1.-e)) };
    double ei, period__g, g, x[3*n], error;
    ei = energy (n, xi, xi+n);
    if (argc != 3)
        pe ("Usage: %s <period/g> <number of periods>\n", argv[0]);
    period__g = atof (argv[1]);
    printf ("period / g ~= %g\n", period__g);
    periods = atoi (argv[2]);
    printf ("number of periods = %d\n", periods);
    printf ("therefore total of ~= %g sub-steps/evaluations.\n",
            periods * period__g);
    printf ("ID  K period/g\tperiod/h\tt_f/h\terror\t\tDE\n");
    for (id=from; id<=to; id++)
    { /* L is the total number of full steps */
        L = rint( periods * period__g / newtonian[id].K );
        g = 2. * PI * periods / L / newtonian[id].K;
        for (j=0; j<2*n; j++) x[j] = xi[j];
        newtonian_integrator_prep      (id, g, n, x, x+n, x+2*n, &ddotx);
        newtonian_integrator_firststep (id, g, n, x, x+n, x+2*n, &ddotx);
        for (j=1; j<L; j++)
            newtonian_integrator_nextstep (id, g, n, x, x+n, x+2*n, &ddotx);
        newtonian_integrator_halt (id, g, n, x, x+n, x+2*n, &ddotx);
        /* for (error=j=0; j<2*n; j++) */
        /* if (error < fabs(x[j]-xi[j])) error = fabs(x[j]-xi[j]); */
        error = sqrt( SQUARE(x[0]-xi[0]) + SQUARE(x[1]-xi[1]) +
                      SQUARE(x[2]-xi[2]) + SQUARE(x[3]-xi[3]) );
        printf ("%2d %2d  %.2f\t%.2f\t\t%d\t%e   \t%e\n", id,
                newtonian[id].K, 2*PI/g, 2*PI/g/newtonian[id].K, L,
                error, energy(n,x,x+n)-ei);
        newtonian_integrator_free (id, g, n, x, x+n, x+2*n, &ddotx);
    }
    return (0);
}
#undef to
#undef from
#undef e
#undef n
#endif /* _Kepler_TEST */
