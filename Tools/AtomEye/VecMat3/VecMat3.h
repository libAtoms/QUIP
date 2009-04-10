/********************************************/
/* libVecMat3: -llapack -lblas -lm          */
/*             -lScalar -lIO                */
/*                                          */
/* Three-dimensional Euclidean space vector */
/* and matrix library.                      */
/*                                          */
/* Nov 11 1999 Ju Li <liju99@mit.edu>       */
/********************************************/

#ifndef _VecMat3_h
#define _VecMat3_h

#include <IO.h>
#include <Scalar.h>

typedef double V3[3];
typedef double M3[3][3];
typedef double (*M3P)[3];

/* V3.c: */
/*****************************/
/* single vector with scalar */
/*****************************/

#define MAX3(a,b,c,tmp) (tmp=MAX(a,b),MAX(tmp,c))
#define MIN3(a,b,c,tmp) (tmp=MIN(a,b),MIN(tmp,c))
#define ABS3(a,b) ((b)[0]=ABS((a)[0]),(b)[1]=ABS((a)[1]),(b)[2]=ABS((a)[2]))
#define FABS3(a,b)((b)[0]=fabs((a)[0]),(b)[1]=fabs((a)[1]),(b)[2]=fabs((a)[2]))

#define ZeroV3 {0.,0.,0.}

/* memory expansion */
#define V3e(a) (a),(a)+1,(a)+2
#define V3E(a) (a)[0],(a)[1],(a)[2]

/* unary vector */

/* a[] := 0; return a[] */
double *V3zero (double a[3]);
#define V3ZERO(a) ((a)[0]=0, (a)[1]=0, (a)[2]=0)
#define V3EYE(a,x) ((x)[0]=(a), (x)[1]=(a), (x)[2]=(a))

/* generate a unit vector: a[i] := 1, a[j!=i] := 0; return a[] */
double *V3unit (double a[3], int i);
#define V3UNIT(a,i) ((a)[0]=0,(a)[1]=0,(a)[2]=0,(a)[i]=1)

#define V3ASSIGN(x0,x1,x2,a) ((a)[0]=(x0),(a)[1]=(x1),(a)[2]=(x2))

/* generate a[] with 3 independent random components on (0.,1.); return a[] */
double *V3frandom (double a[3]);
#define V3Frandom(a) ((a)[0]=Frandom(),(a)[1]=Frandom(),(a)[2]=Frandom())

/* generate a with 3 independent random components on (-0.5,0.5) */
#define V3FRANDOM(a) ((a)[0]=FRANDOM(),(a)[1]=FRANDOM(),(a)[2]=FRANDOM())

/* generate a unit vector of spherically uniform random orientation */
double *V3randomunit (double a[3]);

#define V3SphericalUnit(theta,phi,v) \
  V3ASSIGN(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta),v)

/* b[] := a[]; then return b[] */
double *V3eqv (double a[3], double b[3]);
#define V3EQV(a,b) ((b)[0]=(a)[0], (b)[1]=(a)[1], (b)[2]=(a)[2])

#define V3EQ(a,b) ( ((b)[0]==(a)[0]) && ((b)[1]==(a)[1]) && ((b)[2]==(a)[2]) )
#define V3NE(a,b) ( ((b)[0]!=(a)[0]) || ((b)[1]!=(a)[1]) || ((b)[2]!=(a)[2]) )
#define V3EQZERO(a) ( (0 == (a)[0]) && (0 == (a)[1]) && (0 == (a)[2]) )
#define V3NEZERO(a) ( (0 != (a)[0]) || (0 != (a)[1]) || (0 != (a)[2]) )

#define V3ISTINY(a)  ( V3LENGTH2(a) < TINY*TINY )
#define V3ISSMALL(a) ( V3LENGTH2(a) < SMALL*SMALL )
#define V3ISDELTA(a) ( V3LENGTH2(a) < DELTA*DELTA )

/* If a vector component is tiny, set it to zero. USE WITH CARE! */
#define V3_SET_ZERO_IF_TINY(a) { SET_ZERO_IF_TINY((a)[0]); \
  SET_ZERO_IF_TINY((a)[1]); SET_ZERO_IF_TINY((a)[2]) }

/* b[] := -a[]; then return b[] */
double *V3neg (double a[3], double b[3]);
#define V3NEG(a,b) ((b)[0]=-(a)[0], (b)[1]=-(a)[1], (b)[2]=-(a)[2])
/* a[] := -a[]; then return a[] */
double *V3Neg (double a[3]);
#define V3NeG(a) V3NEG(a,a)

#define V3SWAP(a,b,tmp) ( SWAP((a)[0],(b)[0],tmp), \
  SWAP((a)[1],(b)[1],tmp), SWAP((a)[2],(b)[2],tmp) )

#define V3XIN(v,lower_bound,upper_bound) (XIN((v)[0],lower_bound,upper_bound) \
 && XIN((v)[1],lower_bound,upper_bound) && XIN((v)[2],lower_bound,upper_bound))

#define V3XOU(v,lower_bound,upper_bound) (XOU((v)[0],lower_bound,upper_bound) \
 || XOU((v)[1],lower_bound,upper_bound) || XOU((v)[2],lower_bound,upper_bound))

#define V3NONEED_TRIM(x)  ( (((x)[0]) >= 0.) && (((x)[0]) < 1.) && \
  (((x)[1]) >= 0.) && (((x)[1]) < 1.) && (((x)[2]) >= 0.) && (((x)[2]) < 1.) )
#define V3NEED_TRIM(x)    ( (((x)[0]) <  0.) || (((x)[0]) >= 1.) || \
  (((x)[1]) <  0.) || (((x)[1]) >= 1.) || (((x)[2]) <  0.) || (((x)[2]) >= 1.))
/* change a[] to its own image in [0,1)^3; then return a[] */
double *V3Trim (double a[3]);
#define V3TriM(a) { Trim((a)[0]); Trim((a)[1]); Trim((a)[2]); }
/* faster than V3TRIM() when a[] is around 1 */

/* make b[] an image of a[] in [0,1)^3; then return b[] */
double *V3trim (double a[3], double b[3]);
#define V3TRIM(a,b) (  (b)[0]=TRIM((a)[0]), \
  (b)[1]=TRIM((a)[1]), (b)[2]=TRIM((a)[2]) )
/* faster than V3TriM() when a[] is very large */

#define V3NONEED_IMAGE(x)  ( ((x[0]) >= -0.5) && ((x[0]) < 0.5) && \
  ((x[1]) >= -0.5) && ((x[1]) < 0.5) && ((x[2]) >= -0.5) && ((x[2]) < 0.5) )
#define V3NEED_IMAGE(x)    ( ((x[0]) <  -0.5) || ((x[0]) >= 0.5) || \
  ((x[1]) <  -0.5) || ((x[1]) >= 0.5) || ((x[2]) <  -0.5) || ((x[2]) >= 0.5) )
/* change a[] to its own image in [-0.5,0.5)^3; then return a[] */
double *V3Image (double a[3]);
#define V3ImagE(a) { Image((a)[0]); Image((a)[1]); Image((a)[2]); }
/* faster than V3IMAGE() when a[] is around 1 */

/* make b[] image of a[]'s in [-0.5,0.5)^3; then return b[] */
double *V3image (double a[3], double b[3]);
#define V3IMAGE(a,b) (  (b)[0]=IMAGE((a)[0]), \
  (b)[1]=IMAGE((a)[1]), (b)[2]=IMAGE((a)[2]) )
/* faster than V3ImagE() when a[] is very large */

/* scalar & vector */

/* b[i] := a[i] + x, i=0..2 */
#define V3aDD(a,x,b) ((b)[0]=(a)[0]+(x),(b)[1]=(a)[1]+(x),(b)[2]=(a)[2]+(x))

/* a[i] := a[i] + x, i=0..2 */
#define V3aDd(a,x)   ((a)[0]+=(x),(a)[1]+=(x),(a)[2]+=(x))

/* (a)[0]+(a)[1]+(a)[2] */
#define V3SUM(a)  ((a)[0]+(a)[1]+(a)[2])

/* b[] := multiplier * a[]; then return b[] */
double *V3mul (double multiplier, double a[3], double b[3]);
#define V3MUL(multiplier,a,b) ( (b)[0] = (multiplier)*(a)[0], \
  (b)[1] = (multiplier)*(a)[1], (b)[2] = (multiplier)*(a)[2] )

#define V3muL(multiplier0,multiplier1,multiplier2,a,b)  ( \
  (b)[0] = (multiplier0) * (a)[0], (b)[1] = (multiplier1) * (a)[1], \
  (b)[2] = (multiplier2) * (a)[2] )

/* a[] := multiplier * a[]; then return a[] */
double *V3Mul (double multiplier, double a[3]);
#define V3MuL(multiplier,a) ( (a)[0] *= (multiplier), \
  (a)[1] *= (multiplier),     (a)[2] *= (multiplier) )

#define V3mUL(multiplier0,multiplier1,multiplier2,a) ( \
  (a)[0] *= (multiplier0), (a)[1] *= (multiplier1), (a)[2] *= (multiplier2) )

/* b[] := a[] / divisor; then return b[] */
double *V3div (double a[3], double divisor, double b[3]);
#define V3DIV(a,divisor,b) ( (b)[0] = (a)[0] / (divisor), \
  (b)[1] = (a)[1] / (divisor), (b)[2] = (a)[2] / (divisor) )

#define V3diV(a,divisor0,divisor1,divisor2,b) ( \
  (b)[0] = (a)[0] / (divisor0), (b)[1] = (a)[1] / (divisor1), \
  (b)[2] = (a)[2] / (divisor2) )
/* c[] := a[] ./ b[] */
#define V3ddiV(a,b,c) V3diV(a,b[0],b[1],b[2],c)
/* a[] := a[] ./ b[] */
#define V3ddiv(a,b) ( (a)[0]/=(b)[0], (a)[1]/=(b)[1], (a)[2]/=(b)[2] )

/* a[] := a[] / divisor; then return a[] */
double *V3Div (double a[3], double divisor);
#define V3DiV(a,divisor) \
  ((a)[0]/=(divisor), (a)[1]/=(divisor), (a)[2]/=(divisor))

#define V3dIV(a,divisor0,divisor1,divisor2)  ( \
  (a)[0]/=(divisor0), (a)[1]/=(divisor1), (a)[2]/=(divisor2) )

/* length of a[] */
double V3length (double a[3]);
#define V3LENGTH(a) (sqrt((a)[0]*(a)[0]+(a)[1]*(a)[1]+(a)[2]*(a)[2]))

/* squared length of a[] */
double V3length2 (double a[3]);
#define V3LENGTH2(a) ((a)[0]*(a)[0]+(a)[1]*(a)[1]+(a)[2]*(a)[2])

/* arbitrary exponent-th norm of a[] */
double V3norm (double a[3], double exponent);

/* the largest component of a[] in absolute value */
#define V3INFNORM(a,infnorm) { infnorm = ABS(a[0]); \
  if ( a[1] > infnorm ) infnorm = a[1]; \
  else if ( -a[1] > infnorm ) infnorm = -a[1]; \
  if ( a[2] > infnorm ) infnorm = a[2]; \
  else if ( -a[2] > infnorm ) infnorm = -a[2]; }

/* normalize a[] to unit length; then return a[] */
double *V3normalize (double a[3]);
#define V3NORMALIZE(a,r) (r=V3LENGTH(a),a[0]/=r,a[1]/=r,a[2]/=r)

/* normalize a[] to unit length; then return its original length */
double V3Normalize (double a[3]);

/* b[] := a[] / |a[]|, return b[] */
double *V3direction (double a[3], double b[3]);
#define V3DIRECTION(a,b,r) (r=V3LENGTH(a),V3DIV(a,r,b))

/* b[] := a[] / |a[]|, return |a[]| */
double V3Direction (double a[3], double b[3]);
#define V3DirectioN(a,b,r) (r=V3LENGTH(a),V3DIV(a,r,b),r)

/* assign index array idx[] := 0..2 */
#define SEQ3(idx) ( idx[0]=0, idx[1]=1, idx[2]=2 )

/* return 0..2 permutation idx[] such that length[idx[]] is ascending */
int *V3sort (double length[3], int idx[3]);
#define V3SORT(length,idx,tmpi) { SEQ3(idx); \
  if (length[idx[0]]>length[idx[1]]) SWAP(idx[0],idx[1],tmpi); \
  if (length[idx[1]]>length[idx[2]]) SWAP(idx[1],idx[2],tmpi); \
  if (length[idx[0]]>length[idx[1]]) SWAP(idx[0],idx[1],tmpi); }

/* return 0..2 random permutation index idx[] */
int *I3randompermutation (int idx[3]);
#define I3RANDOMPERMUTATION(idx) ( idx[0] = Fran(0,2), \
  idx[1] = (idx[0]+Fran(1,2))%3, idx[2] = 3 - idx[0] - idx[1] )

/* return 0..2 permutation idx[] such that length[idx[]]  */
/* is non-decreasing. Equal length indices are randomized */
int *V3randomsort (double length[3], int idx[3]);

/* return 0..2 permutation idx[] such that length[idx[]] is     */
/* approximately increasing within epsilon. idx[] is randomized */
int *V3Randomsort (double length[3], int idx[3], double epsilon);

/* Return spherical angles theta=[0,pi] and phi=[0,2*pi) of vector v[]    */
/* in the frame spanned by vectors z[] and x[]. x[] does not have to be   */
/* perpendicular to z[]; its perpendicular component to z[] will be used. */
void V3SphericalAngle (V3 v, V3 z, V3 x, double *theta, double *phi);


/* VV3.c: */

/*******************************/
/* vector & vector with scalar */
/*******************************/


/* vector & vector */

/* c[] := a[] + b[]; then return c[] */
double *V3add (double a[3], double b[3], double c[3]);
#define V3ADD(a,b,c) ( (c)[0] = (a)[0] + (b)[0], \
  (c)[1] = (a)[1] + (b)[1], (c)[2] = (a)[2] + (b)[2] )

/* b[] := a[] + b[]; then return b[] */
double *V3Add (double a[3], double b[3]);
#define V3AdD(a,b) ( (b)[0]+=(a)[0], (b)[1]+=(a)[1], (b)[2]+=(a)[2] )

/* c[] := a[] - b[]; then return c[] */
double *V3sub (double a[3], double b[3], double c[3]);
#define V3SUB(a,b,c) ( (c)[0] = (a)[0] - (b)[0], \
  (c)[1] = (a)[1] - (b)[1], (c)[2] = (a)[2] - (b)[2] )

/* a[] := a[] - b[]; then return a[] */
double *V3Sub (double a[3], double b[3]);
#define V3SuB(a,b) ( (a)[0]-=(b)[0], (a)[1]-=(b)[1], (a)[2]-=(b)[2] )

/* dot product of a[] and b[] */
double V3dot (double a[3], double b[3]);
#define V3DOT(a,b) ((a)[0]*(b)[0]+(a)[1]*(b)[1]+(a)[2]*(b)[2])

/* determine if a[] and b[] are collinear: b[] = c*a[], c=(-inf,inf) */
#define COLLINEAR(a,b,tmp) ( tmp = V3DOT(a,b), \
  tmp *= tmp/V3LENGTH2(a)/V3LENGTH2(b), tmp--, ((tmp>-TINY)&&(tmp<TINY)) )
/* assuming a[] and b[] have non-zero length */

/* cross product of a[] and b[]: c := a x b; then return c[] */
double *V3cross (double a[3], double b[3], double c[3]);
#define V3CROSS(a,b,c) ( (c)[0] = (a)[1]*(b)[2] - (a)[2]*(b)[1], \
  (c)[1] = (a)[2]*(b)[0]-(a)[0]*(b)[2], (c)[2] = (a)[0]*(b)[1]-(a)[1]*(b)[0] )

/* angle between vectors a[],b[] in radian [0,pi] */
double V3angle (double a[3], double b[3]);
#define V3angle_IN_DEGREE(a,b) RADIAN_TO_DEGREE( V3angle(a,b) )
#define V3ANGLE(a,b)           acos( V3DOT(a,b) / V3LENGTH(a) / V3LENGTH(b) )
#define V3ANGLE_IN_DEGREE(a,b) RADIAN_TO_DEGREE( V3ANGLE(a,b) )

/* c[i] := a[i] * b[i], i=0..2; then return c[] */
double *V3prod (double a[3], double b[3], double c[3]);
#define V3PROD(a,b,c) ( (c)[0] = (a)[0]*(b)[0], \
  (c)[1] = (a)[1]*(b)[1], (c)[2] = (a)[2]*(b)[2] )

/* a[i] := a[i] * b[i], i=0..2; then return a[] */
double *V3Prod (double a[3], double b[3]);
#define V3ProD(a,b) ( (a)[0]*=(b)[0], (a)[1]*=(b)[1], (a)[2]*=(b)[2] )


/* vector & vector & scalar */

/* c[] := a[] + multiplier * b[]; then return c[] */
double *V3addmul (double a[3], double multiplier, double b[3], double c[3]);
#define V3ADDMUL(a,multiplier,b,c) ( (c)[0]=(a)[0]+(multiplier)*(b)[0], \
  (c)[1] = (a)[1]+(multiplier)*(b)[1], (c)[2] = (a)[2]+(multiplier)*(b)[2] )

/* c[] := aa * a[] + bb * b[] */
#define V3ADDMULMUL(aa,a,bb,b,c) ( (c)[0] = (aa)*(a)[0]+(bb)*(b)[0], \
  (c)[1] = (aa)*(a)[1]+(bb)*(b)[1], (c)[2] = (aa)*(a)[2]+(bb)*(b)[2] )

/* c[] := a[] + b[] / divisor; then return c[] */
double *V3adddiv (double a[3], double b[3], double divisor, double c[3]);
#define V3ADDDIV(a,b,divisor,c) ( (c)[0] = (a)[0] + (b)[0] / (divisor), \
  (c)[1] = (a)[1] + (b)[1] / (divisor), (c)[2] = (a)[2] + (b)[2] / (divisor) )

/* b[] := b[] + multiplier * a[]; then return b[] */
double *V3ADDmul (double multiplier, double a[3], double b[3]);
#define V3ADDmuL(multiplier,a,b) ( (b)[0] += (multiplier) * (a)[0], \
  (b)[1] += (multiplier) * (a)[1], (b)[2] += (multiplier) * (a)[2] )

/* c[] := c[] + aa * a[] + bb * b[]; then return c[] */
double *V3ADDmulmul (double aa,double a[3],double bb,double b[3], double c[3]);
#define V3ADDmuLmuL(aa,a,bb,b,c) ( (c)[0] += (aa)*(a)[0]+(bb)*(b)[0], \
  (c)[1] += (aa)*(a)[1]+(bb)*(b)[1], (c)[2] += (aa)*(a)[2]+(bb)*(b)[2] )

/* b[] := b[] + a[] / divisor; then return b[] */
double *V3ADDdiv (double a[3], double divisor, double b[3]);
#define V3ADDdiV(a,divisor,b) ( (b)[0] += (a)[0] / (divisor), \
  (b)[1] += (a)[1] / (divisor), (b)[2] += (a)[2] / (divisor) )

/* c[] := a[] - multiplier * b[]; then return c[] */
double *V3submul (double a[3], double multiplier, double b[3], double c[3]);
#define V3SUBMUL(a,multiplier,b,c) ( (c)[0]=(a)[0]-(multiplier)*(b)[0], \
  (c)[1] = (a)[1]-(multiplier)*(b)[1], (c)[2] = (a)[2]-(multiplier)*(b)[2] )

/* c[] := multiplier * a[] - b[]; */
#define V3MULSUB(multiplier,a,b,c) ( (c)[0]=(multiplier)*(a)[0]-(b)[0], \
  (c)[1] = (multiplier)*(a)[1]-(b)[1], (c)[2] = (multiplier)*(a)[2]-(b)[2] )

/* c[] := a[] - b[] / divisor; then return c[] */
double *V3subdiv (double a[3], double b[3], double divisor, double c[3]);
#define V3SUBDIV(a,b,divisor,c) ( (c)[0] = (a)[0] - (b)[0] / (divisor), \
  (c)[1] = (a)[1] - (b)[1] / (divisor), (c)[2] = (a)[2] - (b)[2] / (divisor) )

/* a[] := a[] - multiplier * b[]; then return a[] */
double *V3SUBmul (double a[3], double multiplier, double b[3]);
#define V3SUBmuL(a,multiplier,b) ( (a)[0] -= (multiplier) * (b)[0], \
  (a)[1] -= (multiplier) * (b)[1], (a)[2] -= (multiplier) * (b)[2] )

/* a[] := a[] - b[] / divisor; then return a[] */
double *V3SUBdiv (double a[3], double b[3], double divisor);
#define V3SUBdiV(a,b,divisor) ( (a)[0] -= (b)[0] / (divisor), \
  (a)[1] -= (b)[1] / (divisor), (a)[2] -= (b)[2] / (divisor) )

/* c[] := part of a[] that is perpendicular to b[]; return c[] */
double *V3perpendicular (double a[3], double b[3], double c[3]);
#define V3PERPENDICULAR(a,b,c,tmp) \
  ( tmp=V3DOT(a,b)/V3LENGTH2(b), V3SUBMUL(a, tmp, b, c) )

/* a[] := part of a[] that is perpendicular to b[]; return modified a[] */
double *V3Perpendicular (double a[3], double b[3]);
#define V3PerpendiculaR(a,b,tmp) \
  ( tmp=V3DOT(a,b)/V3LENGTH2(b), V3SUBmuL(a, tmp, b) )


/* vector & vector & vector */

/* d[] := a[] + b[] + c[] */
#define V3ADDADD(a,b,c,d) ( \
  (d)[0] = (a)[0] + (b)[0] + (c)[0], \
  (d)[1] = (a)[1] + (b)[1] + (c)[1], \
  (d)[2] = (a)[2] + (b)[2] + (c)[2] )

/* d[] := a[] + b[] - c[] */
#define V3ADDSUB(a,b,c,d) ( (d)[0] = (a)[0] + (b)[0] - (c)[0], \
  (d)[1] = (a)[1] + (b)[1] - (c)[1], (d)[2] = (a)[2] + (b)[2] - (c)[2] )
/* a[] := a[] + b[] - c[] */
#define V3AdDSUB(a,b,c) V3ADDSUB(a,b,c,a)

/* d[] := a[] - b[] - c[] */
#define V3SUBSUB(a,b,c,d) ( (d)[0] = (a)[0] - (b)[0] - (c)[0], \
  (d)[1] = (a)[1] - (b)[1] - (c)[1], (d)[2] = (a)[2] - (b)[2] - (c)[2] )
/* a[] := a[] - b[] - c[] */
#define V3SuBSUB(a,b,c) V3SUBSUB(a,b,c,a)


/* M3.c: */
#define ZeroM3     {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}
#define IdentityM3 {{1.,0.,0.},{0.,1.,0.},{0.,0.,1.}}
#define GuineaPigM3 \
  {{-0.4326,0.2877,1.1892},{-1.6656,-1.1465,-0.0376},{0.1253,1.1909,0.3273}}
#define GuineaPigSymmetricM3 \
  {{0.3493,-0.7750,0.8397},{-0.7750,4.3664,0.9304},{0.8397,0.9304,0.1186}}
#define GuineaPigSymmetricPositiveDefiniteM3 \
  {{3.2067,2.9391,-0.7604},{2.9391,4.2172,0.6024},{-0.7604,0.6024,1.4385}}
/* O * O' = O' * O = I */
#define GuineaPigOrthonormalM3 \
  {{-0.06937036199895,0.88335767480478,-0.46353745397716}, \
  {0.22707763297048,0.46644129002230,0.85490834103245}, \
  {0.97140285231240,-0.04595368674896,-0.23294797100108}}
/* det|O| = 1 */
#define GuineaPigOrthonormalRightHandedM3 \
  {{0.36172601718872,0.13093569411333,-0.92304394938478}, \
  {-0.29104146616393,0.95646607032430,0.02162224066884}, \
  {0.88569134209494,0.26082273736285,0.38408637858249}}
/* det|O| = -1 */
#define GuineaPigOrthonormalLeftHandedM3 \
  {{0.17625319501775,-0.78998368971498,-0.58724831309280}, \
  {0.89471875807891,-0.12015037695124,0.43016535292863}, \
  {0.41038171887618,0.60124009921394,-0.68563633794414}}

/* memory expansion */
#define M3e(A) A[0][0],A[0][1],A[0][2], \
  A[1][0],A[1][1],A[1][2], A[2][0],A[2][1],A[2][2]

/**********************************************/
/* single matrix & scalar & generating vector */
/**********************************************/

/* matrix */

/* symmetric 3x3 -> 0-5 index conversion */
#define M3OD(i,j) ((5-MIN(i,j))*MIN(i,j)/2+MAX(i,j))

#define ISSYMMETRIC(A) (EQUAL((A)[0][1],(A)[1][0]) && \
  EQUAL((A)[0][2],(A)[2][0]) && EQUAL((A)[2][1],(A)[1][2]))

#define M3ASSIGN(x00,x01,x02,x10,x11,x12,x20,x21,x22,A) ( \
  A[0][0]=(x00),A[0][1]=(x01),A[0][2]=(x02), \
  A[1][0]=(x10),A[1][1]=(x11),A[1][2]=(x12), \
  A[2][0]=(x20),A[2][1]=(x21),A[2][2]=(x22) )

/* A[0][] := x0[], A[1][] := x1[], A[2][] := x2[] */
#define M3Assign(x0,x1,x2,A) \
  M3ASSIGN(x0[0],x0[1],x0[2], x1[0],x1[1],x1[2], x2[0],x2[1],x2[2], A)

/* A[][0] := y0[], A[][1] := y1[], A[][2] := y2[] */
#define M3assign(y0,y1,y2,A) \
  M3ASSIGN(y0[0],y1[0],y2[0], y0[1],y1[1],y2[1], y0[2],y1[2],y2[2], A)

/* memory expansion */
#define V6e(a) (a),(a)+1,(a)+2,(a)+3,(a)+4,(a)+5
#define V6E(a) (a)[0],(a)[1],(a)[2],(a)[3],(a)[4],(a)[5]

/* b[] := a[] */
#define V6EQV(a,b) ( \
  (b)[0]=(a)[0], (b)[1]=(a)[1], (b)[2]=(a)[2],\
  (b)[3]=(a)[3], (b)[4]=(a)[4], (b)[5]=(a)[5] )

/* note that Voigt style indexing is used here */
#define V6_To_M3(a,A) (                         \
        (A)[0][0] = (a)[0],                     \
        (A)[1][1] = (a)[1],                     \
        (A)[2][2] = (a)[2],                     \
        (A)[1][2] = (A)[2][1] = (a)[3],         \
        (A)[2][0] = (A)[0][2] = (a)[4],         \
        (A)[0][1] = (A)[1][0] = (a)[5] )

#define M3_To_V6(A,a) (                         \
        (a)[0] = (A)[0][0],                     \
        (a)[1] = (A)[1][1],                     \
        (a)[2] = (A)[2][2],                     \
        (a)[3] = (A)[1][2],                     \
        (a)[4] = (A)[0][2],                     \
        (a)[5] = (A)[0][1] )

#define V6DiV(a,divisor) (                                              \
        (a)[0]/=(divisor),                                              \
        (a)[1]/=(divisor),                                              \
        (a)[2]/=(divisor),                                              \
        (a)[3]/=(divisor),                                              \
        (a)[4]/=(divisor),                                              \
        (a)[5]/=(divisor) )

/* b[] := multiplier * a[] */
#define V6MUL(multiplier,a,b) ( \
  (b)[0] = (multiplier)*(a)[0], (b)[1] = (multiplier)*(a)[1], \
  (b)[2] = (multiplier)*(a)[2], (b)[3] = (multiplier)*(a)[3], \
  (b)[4] = (multiplier)*(a)[4], (b)[5] = (multiplier)*(a)[5] )

/* a[] := multiplier * a[] */
#define V6MuL(multiplier,a) ( \
  (a)[0] *= (multiplier), (a)[1] *= (multiplier), (a)[2] *= (multiplier) \
  (a)[3] *= (multiplier), (a)[4] *= (multiplier), (a)[5] *= (multiplier) )

/* b[] := a[] / divisor */
#define V6DIV(a,divisor,b) ( \
  (b)[0] = (a)[0] / (divisor), (b)[1] = (a)[1] / (divisor), \
  (b)[2] = (a)[2] / (divisor), (b)[3] = (a)[3] / (divisor), \
  (b)[4] = (a)[4] / (divisor), (b)[5] = (a)[5] / (divisor) )

/* a[] := a[] / divisor */
#define V6DiV(a,divisor) ( \
  (a)[0]/=(divisor), (a)[1]/=(divisor), (a)[2]/=(divisor), \
  (a)[3]/=(divisor), (a)[4]/=(divisor), (a)[5]/=(divisor) )

/* B[][] := multiplier * A[][] */
#define M6MULTIPLY(multiplier,A,B) ( \
  V6MUL(multiplier,A[0],B[0]), V6MUL(multiplier,A[1],B[1]), \
  V6MUL(multiplier,A[2],B[2]), V6MUL(multiplier,A[3],B[3]), \
  V6MUL(multiplier,A[4],B[4]), V6MUL(multiplier,A[5],B[5]) )

/* A[][] := multiplier * A[][] */
#define M6MultiplY(multiplier,A) ( \
  V6MuL(multiplier,A[0]), V6MuL(multiplier,A[1]), V6MuL(multiplier,A[2]), \
  V6MuL(multiplier,A[3]), V6MuL(multiplier,A[4]), V6MuL(multiplier,A[5]) )

/* B[][] := A[][] / divisor */
#define M6DIVIDE(A,divisor,B) ( \
  V6DIV(A[0],divisor,B[0]), V6DIV(A[1],divisor,B[1]), \
  V6DIV(A[2],divisor,B[2]), V6DIV(A[3],divisor,B[3]), \
  V6DIV(A[4],divisor,B[4]), V6DIV(A[5],divisor,B[5]) )

/* A[][] := A[][] / divisor */
#define M6DividE(A,divisor) ( \
  V6DiV(A[0],divisor), V6DiV(A[1],divisor), V6DiV(A[2],divisor), \
  V6DiV(A[3],divisor), V6DiV(A[4],divisor), V6DiV(A[5],divisor) )


/* memory expansion */
#define V9e(a) (a),(a)+1,(a)+2,(a)+3,(a)+4,(a)+5,(a)+6,(a)+7,(a)+8
#define V9E(a) (a)[0],(a)[1],(a)[2],(a)[3],(a)[4],(a)[5],(a)[6],(a)[7],(a)[8]

#define V9_To_M3(a,A) ( \
  (A)[0][0] = (a)[0], \
  (A)[0][1] = (a)[1], \
  (A)[0][2] = (a)[2], \
  (A)[1][0] = (a)[3], \
  (A)[1][1] = (a)[4], \
  (A)[1][2] = (a)[5], \
  (A)[2][0] = (a)[6], \
  (A)[2][1] = (a)[7], \
  (A)[2][2] = (a)[8] )

#define M3_To_V9(A,a) ( \
  (a)[0] = (A)[0][0], \
  (a)[1] = (A)[0][1], \
  (a)[2] = (A)[0][2], \
  (a)[3] = (A)[1][0], \
  (a)[4] = (A)[1][1], \
  (a)[5] = (A)[1][2], \
  (a)[6] = (A)[2][0], \
  (a)[7] = (A)[2][1], \
  (a)[8] = (A)[2][2] )

/* c[] := a[] - b[] */
#define V9SUB(a,b,c) ( \
  (c)[0] = (a)[0] - (b)[0], \
  (c)[1] = (a)[1] - (b)[1], \
  (c)[2] = (a)[2] - (b)[2], \
  (c)[3] = (a)[3] - (b)[3], \
  (c)[4] = (a)[4] - (b)[4], \
  (c)[5] = (a)[5] - (b)[5], \
  (c)[6] = (a)[6] - (b)[6], \
  (c)[7] = (a)[7] - (b)[7], \
  (c)[8] = (a)[8] - (b)[8] )

/* c[] := a[] + b[] */
#define V9ADD(a,b,c) ( \
  (c)[0] = (a)[0] + (b)[0], \
  (c)[1] = (a)[1] + (b)[1], \
  (c)[2] = (a)[2] + (b)[2], \
  (c)[3] = (a)[3] + (b)[3], \
  (c)[4] = (a)[4] + (b)[4], \
  (c)[5] = (a)[5] + (b)[5], \
  (c)[6] = (a)[6] + (b)[6], \
  (c)[7] = (a)[7] + (b)[7], \
  (c)[8] = (a)[8] + (b)[8] )

/* b[] := a[] + b[] */
#define V9AdD(a,b) (    \
        (b)[0]+=(a)[0], \
        (b)[1]+=(a)[1], \
        (b)[2]+=(a)[2], \
        (b)[3]+=(a)[3], \
        (b)[4]+=(a)[4], \
        (b)[5]+=(a)[5], \
        (b)[6]+=(a)[6], \
        (b)[7]+=(a)[7], \
        (b)[8]+=(a)[8] )

/* b[] := multiplier * a[] */
#define V9MUL(multiplier,a,b) (                                         \
        (b)[0] = (multiplier)*(a)[0],                                   \
        (b)[1] = (multiplier)*(a)[1],                                   \
        (b)[2] = (multiplier)*(a)[2],                                   \
        (b)[3] = (multiplier)*(a)[3],                                   \
        (b)[4] = (multiplier)*(a)[4],                                   \
        (b)[5] = (multiplier)*(a)[5],                                   \
        (b)[6] = (multiplier)*(a)[6],                                   \
        (b)[7] = (multiplier)*(a)[7],                                   \
        (b)[8] = (multiplier)*(a)[8] )

/* a[] := multiplier * a[] */
#define V9MuL(multiplier,a) ( \
  (a)[0] *= (multiplier), \
  (a)[1] *= (multiplier), \
  (a)[2] *= (multiplier), \
  (a)[3] *= (multiplier), \
  (a)[4] *= (multiplier), \
  (a)[5] *= (multiplier), \
  (a)[6] *= (multiplier), \
  (a)[7] *= (multiplier), \
  (a)[8] *= (multiplier) )

/* b[] := a[] / divisor */
#define V9DIV(a,divisor,b) ( \
  (b)[0] = (a)[0] / (divisor), (b)[1] = (a)[1] / (divisor), \
  (b)[2] = (a)[2] / (divisor), (b)[3] = (a)[3] / (divisor), \
  (b)[4] = (a)[4] / (divisor), (b)[5] = (a)[5] / (divisor), \
  (b)[6] = (a)[6] / (divisor), (b)[7] = (a)[7] / (divisor), \
  (b)[8] = (a)[8] / (divisor) )

/* a[] := a[] / divisor */
#define V9DiV(a,divisor) ( \
  (a)[0] /= (divisor), \
  (a)[1] /= (divisor), \
  (a)[2] /= (divisor), \
  (a)[3] /= (divisor), \
  (a)[4] /= (divisor), \
  (a)[5] /= (divisor), \
  (a)[6] /= (divisor), \
  (a)[7] /= (divisor), \
  (a)[8] /= (divisor) )

/* a[] := 0 */
#define V9ZERO(a) (\
  (a)[0]=0, (a)[1]=0, (a)[2]=0, \
  (a)[3]=0, (a)[4]=0, (a)[5]=0, \
  (a)[6]=0, (a)[7]=0, (a)[8]=0 )

/* b[] := a[] */
#define V9EQV(a,b) ( \
  (b)[0]=(a)[0], (b)[1]=(a)[1], (b)[2]=(a)[2],\
  (b)[3]=(a)[3], (b)[4]=(a)[4], (b)[5]=(a)[5],\
  (b)[6]=(a)[6], (b)[7]=(a)[7], (b)[8]=(a)[8] )

/* squared length of a[] */
#define V9LENGTH2(a) (\
  (a)[0]*(a)[0]+(a)[1]*(a)[1]+(a)[2]*(a)[2]+\
  (a)[3]*(a)[3]+(a)[4]*(a)[4]+(a)[5]*(a)[5]+\
  (a)[6]*(a)[6]+(a)[7]*(a)[7]+(a)[8]*(a)[8])

/* c[] := a[] + multiplier * b[] */
#define V9ADDMUL(a,multiplier,b,c) (\
  (c)[0] = (a)[0]+(multiplier)*(b)[0], \
  (c)[1] = (a)[1]+(multiplier)*(b)[1], \
  (c)[2] = (a)[2]+(multiplier)*(b)[2], \
  (c)[3] = (a)[3]+(multiplier)*(b)[3], \
  (c)[4] = (a)[4]+(multiplier)*(b)[4], \
  (c)[5] = (a)[5]+(multiplier)*(b)[5], \
  (c)[6] = (a)[6]+(multiplier)*(b)[6], \
  (c)[7] = (a)[7]+(multiplier)*(b)[7], \
  (c)[8] = (a)[8]+(multiplier)*(b)[8] )

/* c[] := a[] - multiplier * b[] */
#define V9SUBMUL(a,multiplier,b,c) (   \
  (c)[0] = (a)[0]-(multiplier)*(b)[0], \
  (c)[1] = (a)[1]-(multiplier)*(b)[1], \
  (c)[2] = (a)[2]-(multiplier)*(b)[2], \
  (c)[3] = (a)[3]-(multiplier)*(b)[3], \
  (c)[4] = (a)[4]-(multiplier)*(b)[4], \
  (c)[5] = (a)[5]-(multiplier)*(b)[5], \
  (c)[6] = (a)[6]-(multiplier)*(b)[6], \
  (c)[7] = (a)[7]-(multiplier)*(b)[7], \
  (c)[8] = (a)[8]-(multiplier)*(b)[8] )

/* b[] := b[] + multiplier * a[] */
#define V9ADDmuL(multiplier,a,b) (                                      \
        (b)[0] += (multiplier) * (a)[0],                                \
        (b)[1] += (multiplier) * (a)[1],                                \
        (b)[2] += (multiplier) * (a)[2],                                \
        (b)[3] += (multiplier) * (a)[3],                                \
        (b)[4] += (multiplier) * (a)[4],                                \
        (b)[5] += (multiplier) * (a)[5],                                \
        (b)[6] += (multiplier) * (a)[6],                                \
        (b)[7] += (multiplier) * (a)[7],                                \
        (b)[8] += (multiplier) * (a)[8] )

/* c[] := a[] + b[] / divisor */
#define V9ADDDIV(a,b,divisor,c) (   \
  (c)[0] = (a)[0]+(b)[0]/(divisor), \
  (c)[1] = (a)[1]+(b)[1]/(divisor), \
  (c)[2] = (a)[2]+(b)[2]/(divisor), \
  (c)[3] = (a)[3]+(b)[3]/(divisor), \
  (c)[4] = (a)[4]+(b)[4]/(divisor), \
  (c)[5] = (a)[5]+(b)[5]/(divisor), \
  (c)[6] = (a)[6]+(b)[6]/(divisor), \
  (c)[7] = (a)[7]+(b)[7]/(divisor), \
  (c)[8] = (a)[8]+(b)[8]/(divisor) )

/* b[] := b[] + a[] / divisor */
#define V9ADDdiV(a,divisor,b) (                                       \
        (b)[0] += (a)[0] / (divisor),                                 \
        (b)[1] += (a)[1] / (divisor),                                 \
        (b)[2] += (a)[2] / (divisor),                                 \
        (b)[3] += (a)[3] / (divisor),                                 \
        (b)[4] += (a)[4] / (divisor),                                 \
        (b)[5] += (a)[5] / (divisor),                                 \
        (b)[6] += (a)[6] / (divisor),                                 \
        (b)[7] += (a)[7] / (divisor),                                 \
        (b)[8] += (a)[8] / (divisor) )

/* c[] := a[] - b[] / divisor */
#define V9SUBDIV(a,b,divisor,c) (   \
  (c)[0] = (a)[0]-(b)[0]/(divisor), \
  (c)[1] = (a)[1]-(b)[1]/(divisor), \
  (c)[2] = (a)[2]-(b)[2]/(divisor), \
  (c)[3] = (a)[3]-(b)[3]/(divisor), \
  (c)[4] = (a)[4]-(b)[4]/(divisor), \
  (c)[5] = (a)[5]-(b)[5]/(divisor), \
  (c)[6] = (a)[6]-(b)[6]/(divisor), \
  (c)[7] = (a)[7]-(b)[7]/(divisor), \
  (c)[8] = (a)[8]-(b)[8]/(divisor) )

/* b[] := b[] - a[] / divisor */
#define V9SUBdiV(a,divisor,b) (                                       \
        (b)[0] -= (a)[0] / (divisor),                                 \
        (b)[1] -= (a)[1] / (divisor),                                 \
        (b)[2] -= (a)[2] / (divisor),                                 \
        (b)[3] -= (a)[3] / (divisor),                                 \
        (b)[4] -= (a)[4] / (divisor),                                 \
        (b)[5] -= (a)[5] / (divisor),                                 \
        (b)[6] -= (a)[6] / (divisor),                                 \
        (b)[7] -= (a)[7] / (divisor),                                 \
        (b)[8] -= (a)[8] / (divisor) )

/* c[] := aa * a[] + bb * b[] */
#define V9ADDMULMUL(aa,a,bb,b,c) (\
  (c)[0] = (aa)*(a)[0]+(bb)*(b)[0], \
  (c)[1] = (aa)*(a)[1]+(bb)*(b)[1], \
  (c)[2] = (aa)*(a)[2]+(bb)*(b)[2], \
  (c)[3] = (aa)*(a)[3]+(bb)*(b)[3], \
  (c)[4] = (aa)*(a)[4]+(bb)*(b)[4], \
  (c)[5] = (aa)*(a)[5]+(bb)*(b)[5], \
  (c)[6] = (aa)*(a)[6]+(bb)*(b)[6], \
  (c)[7] = (aa)*(a)[7]+(bb)*(b)[7], \
  (c)[8] = (aa)*(a)[8]+(bb)*(b)[8] )


/* C[][] := a[]' * b[] */
#define M3ASSIGNV3V3(a,b,C) ( \
  (C)[0][0] = (a)[0] * (b)[0], \
  (C)[0][1] = (a)[0] * (b)[1], \
  (C)[0][2] = (a)[0] * (b)[2], \
  (C)[1][0] = (a)[1] * (b)[0], \
  (C)[1][1] = (a)[1] * (b)[1], \
  (C)[1][2] = (a)[1] * (b)[2], \
  (C)[2][0] = (a)[2] * (b)[0], \
  (C)[2][1] = (a)[2] * (b)[1], \
  (C)[2][2] = (a)[2] * (b)[2] )
#define M3ASSIGNSYMMETRIZEDV3V3(a,b,C) ( \
  (C)[0][0] = (a)[0] * (b)[0], \
  (C)[1][1] = (a)[1] * (b)[1], \
  (C)[2][2] = (a)[2] * (b)[2], \
  (C)[0][1] = (C)[1][0] = ( (a)[0] * (b)[1] + (a)[1] * (b)[0] ) / 2., \
  (C)[0][2] = (C)[2][0] = ( (a)[0] * (b)[2] + (a)[2] * (b)[0] ) / 2., \
  (C)[1][2] = (C)[2][1] = ( (a)[1] * (b)[2] + (a)[2] * (b)[1] ) / 2.  )

/* generate a zero matrix A[][] := 0 */
void M3zero (double A[3][3]);
#define M3ZERO(A) ( A[0][0]=0, A[0][1]=0, A[0][2]=0, A[1][0]=0, \
  A[1][1]=0, A[1][2]=0, A[2][0]=0, A[2][1]=0, A[2][2]=0 )

/* generate an identity matrix A[][] := I[][] */
void M3identity (double A[3][3]);
#define M3IDENTITY(A) ( A[0][0]=1, A[0][1]=0, A[0][2]=0, \
  A[1][0]=0, A[1][1]=1, A[1][2]=0, A[2][0]=0, A[2][1]=0, A[2][2]=1 )

/* generate A[][] with 9 independent random components on (0.,1.) */
void M3frandom (double A[3][3]);
#define M3Frandom(A) ( A[0][0] = Frandom(), A[0][1] = Frandom(), \
  A[0][2] = Frandom(), A[1][0] = Frandom(), A[1][1] = Frandom(), \
  A[1][2] = Frandom(), A[2][0] = Frandom(), A[2][1] = Frandom(), \
  A[2][2] = Frandom())

/* generate A[][] with 9 independent random components on (-0.5,0.5) */
#define M3FRANDOM(A) ( A[0][0] = FRANDOM(), A[0][1] = FRANDOM(), \
  A[0][2] = FRANDOM(), A[1][0] = FRANDOM(), A[1][1] = FRANDOM(), \
  A[1][2] = FRANDOM(), A[2][0] = FRANDOM(), A[2][1] = FRANDOM(), \
  A[2][2] = FRANDOM() )

/* B[][] := A[][] */
void M3eqv (double A[3][3], double B[3][3]);
#define M3EQV(A,B) ( B[0][0] = A[0][0], B[0][1] = A[0][1], \
  B[0][2] = A[0][2], B[1][0] = A[1][0], B[1][1] = A[1][1], \
  B[1][2] = A[1][2], B[2][0] = A[2][0], B[2][1] = A[2][1], \
  B[2][2] = A[2][2] )

#define M3EQ(A,B) ( (B[0][0] == A[0][0]) && (B[0][1] == A[0][1]) && \
  (B[0][2] == A[0][2]) && (B[1][0] == A[1][0]) && (B[1][1] == A[1][1]) && \
  (B[1][2] == A[1][2]) && (B[2][0] == A[2][0]) && (B[2][1] == A[2][1]) && \
  (B[2][2] == A[2][2]) )
#define M3EQZERO(A) ( (0 == A[0][0]) && (0 == A[0][1]) && (0 == A[0][2]) && \
  (0 == A[1][0]) && (0 == A[1][1]) && (0 == A[1][2]) && (0 == A[2][0]) && \
  (0 == A[2][1]) && (0 == A[2][2]) )

#define M3NE(A,B) ( (B[0][0] != A[0][0]) || (B[0][1] != A[0][1]) || \
  (B[0][2] != A[0][2]) || (B[1][0] != A[1][0]) || (B[1][1] != A[1][1]) || \
  (B[1][2] != A[1][2]) || (B[2][0] != A[2][0]) || (B[2][1] != A[2][1]) || \
  (B[2][2] != A[2][2]) )
#define M3NEZERO(A) ( (0 != A[0][0]) || (0 != A[0][1]) || (0 != A[0][2]) || \
  (0 != A[1][0]) || (0 != A[1][1]) || (0 != A[1][2]) || (0 != A[2][0]) || \
  (0 != A[2][1]) || (0 != A[2][2]) )

/* If a matrix element is tiny, set it to zero. USE WITH CARE! */
#define M3_SET_ZERO_IF_TINY(A) { SET_ZERO_IF_TINY(A[0][0]); \
  SET_ZERO_IF_TINY(A[0][1]); SET_ZERO_IF_TINY(A[0][2]); \
  SET_ZERO_IF_TINY(A[1][0]); SET_ZERO_IF_TINY(A[1][1]); \
  SET_ZERO_IF_TINY(A[1][2]); SET_ZERO_IF_TINY(A[2][0]); \
  SET_ZERO_IF_TINY(A[2][1]); SET_ZERO_IF_TINY(A[2][2]) }

#define M3ORTHOGONAL(A) ( ISTINY(V3DOT(A[0],A[1])) && \
  ISTINY(V3DOT(A[0],A[2])) && ISTINY(V3DOT(A[1],A[2])) )
#define M3NONORTHOGONAL(A) ( NOTTINY(V3DOT(A[0],A[1])) || \
  NOTTINY(V3DOT(A[0],A[2])) || NOTTINY(V3DOT(A[1],A[2])) )

/* arbitrary exponent-th norm of A[][] (all components) */
double M3norm (double A[3][3], double exponent);

/* return the largest component of A[][] in absolute value */
double M3infnorm (double A[3][3]);
#define M3INFNORM(A,infnorm) { infnorm = ABS( A[0][0] ); \
  if ( A[0][1] > infnorm ) infnorm = A[0][1]; \
  else if ( -A[0][1] > infnorm ) infnorm = -A[0][1]; \
  if ( A[0][2] > infnorm ) infnorm = A[0][2]; \
  else if ( -A[0][2] > infnorm ) infnorm = -A[0][2]; \
  if ( A[1][0] > infnorm ) infnorm = A[1][0]; \
  else if ( -A[1][0] > infnorm ) infnorm = -A[1][0]; \
  if ( A[1][1] > infnorm ) infnorm = A[1][1]; \
  else if ( -A[1][1] > infnorm ) infnorm = -A[1][1]; \
  if ( A[1][2] > infnorm ) infnorm = A[1][2]; \
  else if ( -A[1][2] > infnorm ) infnorm = -A[1][2]; \
  if ( A[2][0] > infnorm ) infnorm = A[2][0]; \
  else if ( -A[2][0] > infnorm ) infnorm = -A[2][0]; \
  if ( A[2][1] > infnorm ) infnorm = A[2][1]; \
  else if ( -A[2][1] > infnorm ) infnorm = -A[2][1]; \
  if ( A[2][2] > infnorm ) infnorm = A[2][2]; \
  else if ( -A[2][2] > infnorm ) infnorm = -A[2][2]; }

/* return sqrt{ \sum_{i=0..2} \sum_{j=0..2} |A_ij|^2 } */
double M32norm (double A[3][3]);
#define M32NORM2(A) (V3LENGTH2(A[0])+V3LENGTH2(A[1])+V3LENGTH2(A[2]))
#define M32NORM(A) sqrt(V3LENGTH2(A[0])+V3LENGTH2(A[1])+V3LENGTH2(A[2]))

/* B[][] := -A[][] */
void M3neg (double A[3][3], double B[3][3]);
#define M3NEG(A,B) (  B[0][0] = -A[0][0], B[0][1] = -A[0][1], \
  B[0][2] = -A[0][2], B[1][0] = -A[1][0], B[1][1] = -A[1][1], \
  B[1][2] = -A[1][2], B[2][0] = -A[2][0], B[2][1] = -A[2][1], \
  B[2][2] = -A[2][2] )

/* A[][] := -A[][] */
void M3Neg (double A[3][3]);
#define M3NeG(A) M3NEG(A,A)

/* B[][] := A[][]' */
void M3transpose (double A[3][3], double B[3][3]);
/* !! please use M3Transpose() for A[][] := A[][]' !! */
#define M3TRANSPOSE(A,B) ( B[0][0] = A[0][0], B[0][1] = A[1][0], \
  B[0][2] = A[2][0], B[1][0] = A[0][1], B[1][1] = A[1][1], \
  B[1][2] = A[2][1], B[2][0] = A[0][2], B[2][1] = A[1][2], \
  B[2][2] = A[2][2] )
/* !! please use M3Transpose() for A[][] := A[][]' !! */

/* A[][] := A[][]' */
void M3Transpose (double A[3][3]);
#define M3TransposE(A,tmp) ( tmp = A[0][1], A[0][1] = A[1][0], \
  A[1][0] = tmp, tmp = A[0][2], A[0][2] = A[2][0], A[2][0] = tmp, \
  tmp = A[1][2], A[1][2] = A[2][1], A[2][1] = tmp )

/* B[][] := (A[][]+A[][]')/2 */
void M3symmetrize (double A[3][3], double B[3][3]);
/* please use M3Symmetrize() for self-application */
#define M3SYMMETRIZE(A,B) ( B[0][0] = A[0][0], B[1][1] = A[1][1], \
  B[2][2] = A[2][2], B[0][1] = B[1][0] = (A[1][0]+A[0][1])/2., \
  B[0][2]=B[2][0]=(A[2][0]+A[0][2])/2., B[1][2]=B[2][1]=(A[2][1]+A[1][2])/2. )

/* A[][] := (A[][]+A[][]')/2 */
#define M3Symmetrize(A) (A[0][1] = A[1][0] = (A[1][0]+A[0][1])/2., \
  A[0][2]=A[2][0]=(A[2][0]+A[0][2])/2., A[1][2]=A[2][1]=(A[2][1]+A[1][2])/2.)

/* A[][] := A[][]+A[][]' */
#define M3Symmetrize2(A) ( A[0][0]*=2, A[1][1]*=2, A[2][2]*=2, \
  A[0][1]=A[1][0]=A[1][0]+A[0][1], A[0][2]=A[2][0]=A[2][0]+A[0][2], \
  A[1][2]=A[2][1]=A[2][1]+A[1][2] )

/* return the trace of A[][] */
double M3Tr (double A[3][3]);
#define M3TR(A) (A[0][0]+A[1][1]+A[2][2])

/* Tr(AB) */
#define M3TRPROD(A,B) (B[0][0]*A[0][0] + B[0][1]*A[1][0] + \
  B[0][2]*A[2][0] + B[1][0]*A[0][1] + B[1][1]*A[1][1] + \
  B[1][2]*A[2][1] + B[2][0]*A[0][2] + B[2][1]*A[1][2] + B[2][2]*A[2][2] )

/* B[][] := trace(A[][]) / 3 * I[][]; return trace(A) */
double M3trace (double A[3][3], double B[3][3]);
#define M3TRACE(A,B,trace) ( (trace) = A[0][0]+A[1][1]+A[2][2], \
  B[0][0] = (trace)/3., B[0][1] = 0., B[0][2] = 0., \
  B[1][0] = 0., B[1][1] = (trace)/3., B[1][2] = 0., \
  B[2][0] = 0., B[2][1] = 0., B[2][2] = (trace)/3. )

/* A[][] := trace(A[][])/3 * I[][]; return original trace(A) */
double M3Trace (double A[3][3]);
#define M3TracE(A,trace) ( (trace) = A[0][0]+A[1][1]+A[2][2], \
  A[0][0] = (trace)/3., A[0][1] = 0., A[0][2] = 0., \
  A[1][0] = 0., A[1][1] = (trace)/3., A[1][2] = 0., \
  A[2][0] = 0., A[2][1] = 0., A[2][2] = (trace)/3., (trace) )

/* B[][] := A[][] - trace(A[][])/3 * I[][]; return trace(A) */
double M3traceless (double A[3][3], double B[3][3]);
#define M3TRACELESS(A,B,trace) ( (trace) = A[0][0]+A[1][1]+A[2][2], \
  B[0][0] = A[0][0] - (trace)/3., B[0][1] = A[0][1], B[0][2] = A[0][2], \
  B[1][0] = A[1][0], B[1][1] = A[1][1] - (trace)/3., B[1][2] = A[1][2], \
  B[2][0] = A[2][0], B[2][1] = A[2][1], B[2][2] = A[2][2] - (trace)/3. )

/* A[][] := A[][] - trace(A[][])/3 * I[][]; return original trace(A) */
double M3Traceless(double A[3][3]);
#define M3TracelesS(A,trace) ( (trace) = A[0][0]+A[1][1]+A[2][2], \
  A[0][0] -= (trace)/3., A[1][1] -= (trace)/3., A[2][2] -= (trace)/3. )

/* decompose A[][] to b*I[][] + C[][], where b := trace(A)/3; return b */
double M3tracedecompose (double A[3][3], double B[3][3], double C[3][3]);
#define M3TRACEDECOMPOSE(A,B,C,trace) ( (trace) = A[0][0]+A[1][1]+A[2][2], \
  B[0][0] = (trace)/3., B[0][1] = 0., B[0][2] = 0., B[1][0] = 0., \
  B[1][1] = (trace)/3., B[1][2] = 0., B[2][0] = 0., B[2][1] = 0., \
  B[2][2] = (trace)/3., C[0][0] = A[0][0] - (trace)/3., C[0][1] = A[0][1], \
  C[0][2] = A[0][2], C[1][0] = A[1][0], C[1][1] = A[1][1] - (trace)/3., \
  C[1][2] = A[1][2], C[2][0] = A[2][0], C[2][1] = A[2][1], \
  C[2][2] = A[2][2] - (trace)/3. )


/* matrix & scalar */

/* generate an identity matrix A[][] := a x I[][] */
void M3Identity (double a, double A[3][3]);
#define M3IdentitY(a,A) ( A[0][0]=(a), A[0][1]=0., A[0][2]=0., A[1][0] = 0., \
  A[1][1] = (a), A[1][2] = 0., A[2][0] = 0., A[2][1] = 0., A[2][2] = (a) )

/* B[][] := multiplier * A[][] */
void M3multiply (double multiplier, double A[3][3], double B[3][3]);
#define M3MULTIPLY(multiplier,A,B) (B[0][0] = (multiplier) * A[0][0], \
  B[0][1] = (multiplier) * A[0][1], B[0][2] = (multiplier) * A[0][2], \
  B[1][0] = (multiplier) * A[1][0], B[1][1] = (multiplier) * A[1][1], \
  B[1][2] = (multiplier) * A[1][2], B[2][0] = (multiplier) * A[2][0], \
  B[2][1] = (multiplier) * A[2][1], B[2][2] = (multiplier) * A[2][2])

/* A[][] := multiplier * A[][] */
void M3Multiply (double multiplier, double A[3][3]);
#define M3MultiplY(multiplier,A) (A[0][0] *= (multiplier), \
  A[0][1] *= (multiplier), A[0][2] *= (multiplier), A[1][0] *= (multiplier), \
  A[1][1] *= (multiplier), A[1][2] *= (multiplier), A[2][0] *= (multiplier), \
  A[2][1] *= (multiplier), A[2][2] *= (multiplier) )

/* B[][] := A[][] / divisor */
void M3divide (double A[3][3], double divisor, double B[3][3]);
#define M3DIVIDE(A,divisor,B) (  B[0][0] = A[0][0] / (divisor), \
  B[0][1] = A[0][1] / (divisor), B[0][2] = A[0][2] / (divisor), \
  B[1][0] = A[1][0] / (divisor), B[1][1] = A[1][1] / (divisor), \
  B[1][2] = A[1][2] / (divisor), B[2][0] = A[2][0] / (divisor), \
  B[2][1] = A[2][1] / (divisor), B[2][2] = A[2][2] / (divisor) )

/* A[][] := A[][] / divisor */
void M3Divide (double A[3][3], double divisor);
#define M3DividE(A,divisor) ( A[0][0] /= (divisor), A[0][1] /= (divisor), \
  A[0][2] /= (divisor), A[1][0] /= (divisor), A[1][1] /= (divisor), \
  A[1][2] /= (divisor), A[2][0] /= (divisor), A[2][1] /= (divisor), \
  A[2][2] /= (divisor) )

/* B[][] := A[][] + a * I */
void M3adddiag (double A[3][3], double a, double B[3][3]);
#define M3ADDDIAG(A,a,B) ( B[0][0] = A[0][0] + (a), B[0][1] = A[0][1], \
  B[0][2] = A[0][2], B[1][0] = A[1][0], B[1][1] = A[1][1] + (a), \
  B[1][2] = A[1][2], B[2][0] = A[2][0], B[2][1] = A[2][1], \
  B[2][2] = A[2][2] + (a) )

/* B[][] := a * A[][] + x * I */
#define M3MULADDDIAG(a,A,x,B) ( \
  B[0][0] = (a)*A[0][0]+(x), B[0][1] = (a)*A[0][1], B[0][2] = (a)*A[0][2], \
  B[1][0] = (a)*A[1][0], B[1][1] = (a)*A[1][1]+(x), B[1][2] = (a)*A[1][2], \
  B[2][0] = (a)*A[2][0], B[2][1] = (a)*A[2][1], B[2][2] = (a)*A[2][2]+(x) )

/* A[][] := A[][] + a * I */
void M3Adddiag (double A[3][3], double a);
#define M3AdddiaG(A,a) ( A[0][0] += (a), A[1][1] += (a), A[2][2] += (a) )

/* B[][] := A[][] - a * I */
void M3subdiag (double A[3][3], double a, double B[3][3]);
#define M3SUBDIAG(A,a,B) ( B[0][0] = A[0][0] - (a), \
  B[0][1] = A[0][1], B[0][2] = A[0][2], B[1][0] = A[1][0], \
  B[1][1] = A[1][1] - (a), B[1][2] = A[1][2], B[2][0] = A[2][0], \
  B[2][1] = A[2][1], B[2][2] = A[2][2] - (a) )

/* A[][] := A[][] - a * I */
void M3Subdiag (double A[3][3], double a);
#define M3SubdiaG(A,a) ( A[0][0] -= (a), A[1][1] -= (a), A[2][2] -= (a) )


/* matrix & generating vector */

/* generate a diagonal matrix A[i][i] := a_i, others = 0 */
#define M3diagonal(a_0,a_1,a_2,A) ( A[0][0] = (a_0), A[0][1] = 0., \
  A[0][2] = 0., A[1][0] = 0., A[1][1] = (a_1), A[1][2] = 0., \
  A[2][0] = 0., A[2][1] = 0., A[2][2] = (a_2) )
#define M3Diagonal(a,A) M3diagonal(a,a,a,A)

/* generate a diagonal matrix A[i][i] := a[i], others = 0 */
#define M3DIAGONAL(a,A) ( A[0][0] = (a)[0], A[0][1] = 0., \
  A[0][2] = 0., A[1][0] = 0., A[1][1] = (a)[1], A[1][2] = 0., \
  A[2][0] = 0., A[2][1] = 0., A[2][2] = (a)[2] )


/* VM3.c: */

/************************/
/* row space operations */
/************************/

/* rowlength[i] := | A[i][] |; returns rowlength[] */
double *M3rowlengths (double A[3][3], double rowlength[3]);
#define M3ROWLENGTHS(A,rowlength) ( rowlength[0] = V3LENGTH(A[0]), \
  rowlength[1] = V3LENGTH(A[1]),    rowlength[2] = V3LENGTH(A[2]) )

/* returns the maximum Euclidean length of A[0][], A[1][], A[2][] */
double M3maxrowlength (double A[3][3]);

/* c[] := a[] * B[][]; then return c[] */
double *V3mulM3 (double a[3], double B[3][3], double c[3]);
#define V3mM3(a,B,c) ( \
  (c)[0] = (a)[0]*B[0][0] + (a)[1]*B[1][0] + (a)[2]*B[2][0], \
  (c)[1] = (a)[0]*B[0][1] + (a)[1]*B[1][1] + (a)[2]*B[2][1], \
  (c)[2] = (a)[0]*B[0][2] + (a)[1]*B[1][2] + (a)[2]*B[2][2] )

#define V3M3LENGTH2(ds,H,dx) ( V3mM3(ds,H,dx), (dx)[3]=V3LENGTH2(dx) )
#define V3M3LENGTH(ds,H,dx)  ( V3mM3(ds,H,dx), (dx)[3]=V3LENGTH(dx) )

/* x*A*y' */
#define V3ADOT(x,A,y,tmp) (V3mM3(x,A,tmp),V3DOT(tmp,y))

/* a[] := a[] * B[][]; then return a[] */
double *V3MULM3 (double a[3], double B[3][3]);
#define V3MM3(a,B,tmp) (V3mM3(a,B,tmp), V3EQV(tmp,a))

/* a[] := a[] * B[][] * multiplier */
#define V3MM3MUL(a,B,multiplier,tmp) (V3mM3(a,B,tmp), V3MUL(multiplier,tmp,a))

/* a[] := a[] * B[][] / divisor */
#define V3MM3DIV(a,B,divisor,tmp) (V3mM3(a,B,tmp), V3DIV(tmp,divisor,a))

/* d[] := a[] + b[] * C[][]; then return d[] */
double *V3addmulM3 (double a[3], double b[3], double C[3][3], double d[3]);
#define V3addmM3(a,b,C,d) ( \
  (d)[0] = (a)[0] + (b)[0]*C[0][0] + (b)[1]*C[1][0] + (b)[2]*C[2][0], \
  (d)[1] = (a)[1] + (b)[0]*C[0][1] + (b)[1]*C[1][1] + (b)[2]*C[2][1], \
  (d)[2] = (a)[2] + (b)[0]*C[0][2] + (b)[1]*C[1][2] + (b)[2]*C[2][2] )

/* a[] := a[] + b[] * C[][]; then return a[] */
double *V3ADDmulM3 (double a[3], double b[3], double C[3][3]);
#define V3ADDmM3(a,b,C,d) ( \
  (a)[0] += (b)[0]*C[0][0] + (b)[1]*C[1][0] + (b)[2]*C[2][0], \
  (a)[1] += (b)[0]*C[0][1] + (b)[1]*C[1][1] + (b)[2]*C[2][1], \
  (a)[2] += (b)[0]*C[0][2] + (b)[1]*C[1][2] + (b)[2]*C[2][2] )

/* d[] := a[] - b[] * C[][]; then return d[] */
double *V3submulM3 (double a[3], double b[3], double C[3][3], double d[3]);
#define V3submM3(a,b,C,d) ( \
  (d)[0] = (a)[0] - (b)[0]*C[0][0] - (b)[1]*C[1][0] - (b)[2]*C[2][0], \
  (d)[1] = (a)[1] - (b)[0]*C[0][1] - (b)[1]*C[1][1] - (b)[2]*C[2][1], \
  (d)[2] = (a)[2] - (b)[0]*C[0][2] - (b)[1]*C[1][2] - (b)[2]*C[2][2] )

/* d[] := b[] * C[][] - a[]; then return d[] */
double *V3mulsubM3 (double b[3], double C[3][3], double a[3], double d[3]);
#define V3msubM3(b,C,a,d) ( \
  (d)[0] = (b)[0]*C[0][0] + (b)[1]*C[1][0] + (b)[2]*C[2][0] - (a)[0], \
  (d)[1] = (b)[0]*C[0][1] + (b)[1]*C[1][1] + (b)[2]*C[2][1] - (a)[1], \
  (d)[2] = (b)[0]*C[0][2] + (b)[1]*C[1][2] + (b)[2]*C[2][2] - (a)[2] )

/* a[] := a[] - b[] * C[][]; then return a[] */
double *V3SUBmulM3 (double a[3], double b[3], double C[3][3]);
#define V3SUBmM3(a,b,C) ( \
  (a)[0] -= (b)[0]*C[0][0] + (b)[1]*C[1][0] + (b)[2]*C[2][0], \
  (a)[1] -= (b)[0]*C[0][1] + (b)[1]*C[1][1] + (b)[2]*C[2][1], \
  (a)[2] -= (b)[0]*C[0][2] + (b)[1]*C[1][2] + (b)[2]*C[2][2] )


/* MV3.c: */

/***************************/
/* column space operations */
/***************************/

/* columnlength[i] := |A[][i]|; returns columnlength[] */
double *M3columnlengths (double A[3][3], double columnlength[3]);
#define M3COLUMNLENGTHS(A,columnlength) ( \
  (columnlength)[0] = DISTANCE(A[0][0],A[1][0],A[2][0]), \
  (columnlength)[1] = DISTANCE(A[0][1],A[1][1],A[2][1]), \
  (columnlength)[2] = DISTANCE(A[0][2],A[1][2],A[2][2]) )

/* returns the maximum Euclidean length of A[][0], A[][1], A[][2] */
double M3maxcolumnlength (double A[3][3]);

/* column[] := A[][i]; return column[] */
double *M3column (double A[3][3], int i, double column[3]);
#define M3COLUMN(A,i,column) \
  ( (column)[0] = A[0][i], (column)[1] = A[1][i], (column)[2] = A[2][i] )

/* c[] := A[][] * b[]; then return c[] */
double *M3mulV3 (double A[3][3], double b[3], double c[3]);
#define M3mV3(A,b,c) ( \
  (c)[0] = A[0][0]*(b)[0] + A[0][1]*(b)[1] + A[0][2]*(b)[2], \
  (c)[1] = A[1][0]*(b)[0] + A[1][1]*(b)[1] + A[1][2]*(b)[2], \
  (c)[2] = A[2][0]*(b)[0] + A[2][1]*(b)[1] + A[2][2]*(b)[2] )

/* b[] := A[][] * b[]; then return b[] */
double *M3MULV3 (double A[3][3], double b[3]);
#define M3MV3(A,b,tmp) ( M3mV3(A,b,tmp), V3EQV(tmp,b) )

/* d[] := A[][] * b[] + c[]; then return d[] */
double *M3muladdV3 (double A[3][3], double b[3], double c[3], double d[3]);
#define M3maddV3(A,b,c,d) ( \
  (d)[0] = A[0][0]*(b)[0] + A[0][1]*(b)[1] + A[0][2]*(b)[2] + (c)[0], \
  (d)[1] = A[1][0]*(b)[0] + A[1][1]*(b)[1] + A[1][2]*(b)[2] + (c)[1], \
  (d)[2] = A[2][0]*(b)[0] + A[2][1]*(b)[1] + A[2][2]*(b)[2] + (c)[2] )

/* c[] := c[] + A[][] * b[]; then return c[] */
double *M3mulADDV3 (double A[3][3], double b[3], double c[3]);
#define M3mADDV3(A,b,c) ( \
  (c)[0] += A[0][0]*(b)[0] + A[0][1]*(b)[1] + A[0][2]*(b)[2], \
  (c)[1] += A[1][0]*(b)[0] + A[1][1]*(b)[1] + A[1][2]*(b)[2], \
  (c)[2] += A[2][0]*(b)[0] + A[2][1]*(b)[1] + A[2][2]*(b)[2] )

/* d[] := A[][] * b[] - c[]; then return d[] */
double *M3mulsubV3 (double A[3][3], double b[3], double c[3], double d[3]);
#define M3msubV3(A,b,c,d) ( \
  (d)[0] = A[0][0]*(b)[0] + A[0][1]*(b)[1] + A[0][2]*(b)[2] - (c)[0], \
  (d)[1] = A[1][0]*(b)[0] + A[1][1]*(b)[1] + A[1][2]*(b)[2] - (c)[1], \
  (d)[2] = A[2][0]*(b)[0] + A[2][1]*(b)[1] + A[2][2]*(b)[2] - (c)[2] )

/* d[] := c[] - A[][] * b[]; then return d[] */
double *M3submulV3 (double c[3], double A[3][3], double b[3], double d[3]);
#define M3submV3(c,A,b,d) ( \
  (d)[0] = (c)[0] - A[0][0]*(b)[0] - A[0][1]*(b)[1] - A[0][2]*(b)[2], \
  (d)[1] = (c)[1] - A[1][0]*(b)[0] - A[1][1]*(b)[1] - A[1][2]*(b)[2], \
  (d)[2] = (c)[2] - A[2][0]*(b)[0] - A[2][1]*(b)[1] - A[2][2]*(b)[2] )

/* c[] := c[] - A[][] * b[]; then return c[] */
double *M3mulSUBV3 (double A[3][3], double b[3], double c[3]);
#define M3mSUBV3(A,b,c) ( \
  (c)[0] -= A[0][0]*(b)[0] + A[0][1]*(b)[1] + A[0][2]*(b)[2], \
  (c)[1] -= A[1][0]*(b)[0] + A[1][1]*(b)[1] + A[1][2]*(b)[2], \
  (c)[2] -= A[2][0]*(b)[0] + A[2][1]*(b)[1] + A[2][2]*(b)[2] )


/* M3diag.c */

/*********************************************************/
/* Matrix inversion and symmetric matrix diagonalization */
/*********************************************************/

/* B[][] := A[][]^-1; return det(A) */
double M3inv (double A[3][3], double B[3][3]);

/* A[][] := A[][]^-1; return original det(A) */
double M3Inv (double A[3][3]);

#define M3DETERMINANT(A) ( \
  A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1]) + \
  A[0][1]*(A[1][2]*A[2][0]-A[1][0]*A[2][2]) + \
  A[0][2]*(A[1][0]*A[2][1]-A[2][0]*A[1][1]) )
#define M3VOLUME(A) fabs(M3DETERMINANT(A))

/* B[][] := A[][]^-1; determinant := det(A) */
#define M3INV(A,B,determinant) ( \
  B[0][0] = A[1][1]*A[2][2]-A[1][2]*A[2][1], \
  B[1][1] = A[2][2]*A[0][0]-A[2][0]*A[0][2], \
  B[2][2] = A[0][0]*A[1][1]-A[0][1]*A[1][0], \
  B[1][0] = A[1][2]*A[2][0]-A[1][0]*A[2][2], \
  B[2][1] = A[2][0]*A[0][1]-A[2][1]*A[0][0], \
  B[0][2] = A[0][1]*A[1][2]-A[0][2]*A[1][1], \
  B[2][0] = A[1][0]*A[2][1]-A[2][0]*A[1][1], \
  B[0][1] = A[2][1]*A[0][2]-A[0][1]*A[2][2], \
  B[1][2] = A[0][2]*A[1][0]-A[1][2]*A[0][0], \
  (determinant) = A[0][0]*B[0][0] + A[0][1]*B[1][0] + A[0][2]*B[2][0], \
  B[0][0] /= (determinant),B[1][1] /= (determinant),B[2][2] /= (determinant), \
  B[1][0] /= (determinant),B[2][1] /= (determinant),B[0][2] /= (determinant), \
  B[2][0] /= (determinant),B[0][1] /= (determinant),B[1][2] /= (determinant) )

#define M3InV(A,B,volume) { M3INV(A,B,volume); \
  if ((volume) < 0) (volume) = -(volume); }

/* B[][] := A[][]^-T; determinant := det(A) */
#define M3INVTRANSPOSE(A,B,determinant) ( \
  B[0][0] = A[1][1]*A[2][2]-A[1][2]*A[2][1], \
  B[1][1] = A[2][2]*A[0][0]-A[2][0]*A[0][2], \
  B[2][2] = A[0][0]*A[1][1]-A[0][1]*A[1][0], \
  B[0][1] = A[1][2]*A[2][0]-A[1][0]*A[2][2], \
  B[1][2] = A[2][0]*A[0][1]-A[2][1]*A[0][0], \
  B[2][0] = A[0][1]*A[1][2]-A[0][2]*A[1][1], \
  B[0][2] = A[1][0]*A[2][1]-A[2][0]*A[1][1], \
  B[1][0] = A[2][1]*A[0][2]-A[0][1]*A[2][2], \
  B[2][1] = A[0][2]*A[1][0]-A[1][2]*A[0][0], \
  (determinant) = A[0][0]*B[0][0] + A[0][1]*B[0][1] + A[0][2]*B[0][2], \
  B[0][0] /= (determinant),B[1][1] /= (determinant),B[2][2] /= (determinant), \
  B[1][0] /= (determinant),B[2][1] /= (determinant),B[0][2] /= (determinant), \
  B[2][0] /= (determinant),B[0][1] /= (determinant),B[1][2] /= (determinant) )


/***********************************************************************/
/* Find a linear combination of columns A[][0], A[][1], A[][2] so that */
/* A[][]*c[]=0. c[] is normalized and "random" if possible. return c[] */
/***********************************************************************/
double *M3nullvector (double A[3][3], double c[3]);

/******************************************************************/
/* Diagonalize 3x3 real symmetric matrix: A = V^T*Diag(eigval)*V, */
/* eigval[i] will be stored in ascending order, i = 0..2; the     */
/* corresponding eigenvectors V[i][] form a right-handed system.  */
/* Return index i of the eigenvalue of largest absolute value.    */
/******************************************************************/
int M3diag (double A[3][3], double eigval[3], double Q[3][3]);

int FORTRAN_SYMBOL(dsyev)
    (char *jobz, char *uplo, int *n, double *a,
     int *lda, double *w, double *work, int *lwork, 
     int *info);

/* same function as M3diag() but using Lapack */
int M3Diag (double A[3][3], double eigval[3], double V[3][3]);

/* A = M*R, where M is symmetric, R (do not return) is orthogonal */
void M3RightPolarDecompose (M3 A, M3 M);
/* A = M*R, where M is symmetric, R is orthogonal */
void M3RightPolarDECOMPOSE (M3 A, M3 M, M3 R);

/* A = L*M, where L (do not return) is orthogonal, M is symmetric*/
void M3LeftPolarDecompose (M3 A, M3 M);
/* A = L*M, where L is orthogonal, M is symmetric */
void M3LeftPolarDECOMPOSE (M3 A, M3 L, M3 M);

/* MM3.c */

/********************************/
/* matrix & matrix with scalars */
/********************************/

/* matrix & matrix */

/* C[][] := A[][] + B[][] */
void M3add (double A[3][3], double B[3][3], double C[3][3]);
#define M3ADD(A,B,C) (C[0][0] = A[0][0] + B[0][0], \
  C[0][1] = A[0][1] + B[0][1], C[0][2] = A[0][2] + B[0][2], \
  C[1][0] = A[1][0] + B[1][0], C[1][1] = A[1][1] + B[1][1], \
  C[1][2] = A[1][2] + B[1][2], C[2][0] = A[2][0] + B[2][0], \
  C[2][1] = A[2][1] + B[2][1], C[2][2] = A[2][2] + B[2][2])

/* B[][] := A[][] + B[][] */
#define M3AdD(A,B) ( B[0][0] += A[0][0], B[0][1] += A[0][1], \
  B[0][2] += A[0][2], B[1][0] += A[1][0], B[1][1] += A[1][1], \
  B[1][2] += A[1][2], B[2][0] += A[2][0], B[2][1] += A[2][1], \
  B[2][2] += A[2][2] )

/* D[][] := A[][] + b[]' * c[] */
#define M3ADDV3V3(A,b,c,D) ( \
  (D)[0][0] = (A)[0][0] + (b)[0] * (c)[0], \
  (D)[0][1] = (A)[0][1] + (b)[0] * (c)[1], \
  (D)[0][2] = (A)[0][2] + (b)[0] * (c)[2], \
  (D)[1][0] = (A)[1][0] + (b)[1] * (c)[0], \
  (D)[1][1] = (A)[1][1] + (b)[1] * (c)[1], \
  (D)[1][2] = (A)[1][2] + (b)[1] * (c)[2], \
  (D)[2][0] = (A)[2][0] + (b)[2] * (c)[0], \
  (D)[2][1] = (A)[2][1] + (b)[2] * (c)[1], \
  (D)[2][2] = (A)[2][2] + (b)[2] * (c)[2] )

/* C[][] := C[][] + a[]' * b[] */
#define M3AdDV3V3(a,b,C) ( \
  (C)[0][0] += (a)[0] * (b)[0], \
  (C)[0][1] += (a)[0] * (b)[1], \
  (C)[0][2] += (a)[0] * (b)[2], \
  (C)[1][0] += (a)[1] * (b)[0], \
  (C)[1][1] += (a)[1] * (b)[1], \
  (C)[1][2] += (a)[1] * (b)[2], \
  (C)[2][0] += (a)[2] * (b)[0], \
  (C)[2][1] += (a)[2] * (b)[1], \
  (C)[2][2] += (a)[2] * (b)[2] )

/* C[][] := A[][] - B[][] */
void M3sub (double A[3][3], double B[3][3], double C[3][3]);
#define M3SUB(A,B,C) ( C[0][0] = A[0][0] - B[0][0], \
  C[0][1] = A[0][1] - B[0][1], C[0][2] = A[0][2] - B[0][2], \
  C[1][0] = A[1][0] - B[1][0], C[1][1] = A[1][1] - B[1][1], \
  C[1][2] = A[1][2] - B[1][2], C[2][0] = A[2][0] - B[2][0], \
  C[2][1] = A[2][1] - B[2][1], C[2][2] = A[2][2] - B[2][2] )

/* A[][] := A[][] - B[][] */
#define M3SuB(A,B) ( A[0][0] -= B[0][0], A[0][1] -= B[0][1], \
  A[0][2] -= B[0][2], A[1][0] -= B[1][0], A[1][1] -= B[1][1], \
  A[1][2] -= B[1][2], A[2][0] -= B[2][0], A[2][1] -= B[2][1], \
  A[2][2] -= B[2][2] )

/* C[][] := A[][] * B[][] */
void M3mul (double A[3][3], double B[3][3], double C[3][3]);
#define M3MUL(A,B,C) \
( C[0][0] = A[0][0]*B[0][0] + A[0][1]*B[1][0] + A[0][2]*B[2][0], \
  C[0][1] = A[0][0]*B[0][1] + A[0][1]*B[1][1] + A[0][2]*B[2][1], \
  C[0][2] = A[0][0]*B[0][2] + A[0][1]*B[1][2] + A[0][2]*B[2][2], \
  C[1][0] = A[1][0]*B[0][0] + A[1][1]*B[1][0] + A[1][2]*B[2][0], \
  C[1][1] = A[1][0]*B[0][1] + A[1][1]*B[1][1] + A[1][2]*B[2][1], \
  C[1][2] = A[1][0]*B[0][2] + A[1][1]*B[1][2] + A[1][2]*B[2][2], \
  C[2][0] = A[2][0]*B[0][0] + A[2][1]*B[1][0] + A[2][2]*B[2][0], \
  C[2][1] = A[2][0]*B[0][1] + A[2][1]*B[1][1] + A[2][2]*B[2][1], \
  C[2][2] = A[2][0]*B[0][2] + A[2][1]*B[1][2] + A[2][2]*B[2][2] )

/* A[][] := A[][] * B[][] */
void M3Mul (double A[3][3], double B[3][3]);
#define M3MUl(A,B,tmp) ( M3MUL(A,B,tmp), M3EQV(tmp,A) )

/* B[][] := A[][] * B[][] */
void M3muL (double A[3][3], double B[3][3]);
#define M3mUL(A,B,tmp) ( M3MUL(A,B,tmp), M3EQV(tmp,B) )

/* C[][] := A[][] * A'[][] */
#define M3MULT(A,C) \
  (         C[0][0] = A[0][0]*A[0][0] + A[0][1]*A[0][1] + A[0][2]*A[0][2], \
  C[1][0] = C[0][1] = A[0][0]*A[1][0] + A[0][1]*A[1][1] + A[0][2]*A[1][2], \
  C[2][0] = C[0][2] = A[0][0]*A[2][0] + A[0][1]*A[2][1] + A[0][2]*A[2][2], \
            C[1][1] = A[1][0]*A[1][0] + A[1][1]*A[1][1] + A[1][2]*A[1][2], \
  C[2][1] = C[1][2] = A[1][0]*A[2][0] + A[1][1]*A[2][1] + A[1][2]*A[2][2], \
            C[2][2] = A[2][0]*A[2][0] + A[2][1]*A[2][1] + A[2][2]*A[2][2] )

/* C[][] := A'[][] * A[][] */
#define M3TMUL(A,C) \
  (         C[0][0] = A[0][0]*A[0][0] + A[1][0]*A[1][0] + A[2][0]*A[2][0], \
  C[1][0] = C[0][1] = A[0][0]*A[0][1] + A[1][0]*A[1][1] + A[2][0]*A[2][1], \
  C[2][0] = C[0][2] = A[0][0]*A[0][2] + A[1][0]*A[1][2] + A[2][0]*A[2][2], \
            C[1][1] = A[0][1]*A[0][1] + A[1][1]*A[1][1] + A[2][1]*A[2][1], \
  C[2][1] = C[1][2] = A[0][1]*A[0][2] + A[1][1]*A[1][2] + A[2][1]*A[2][2], \
            C[2][2] = A[0][2]*A[0][2] + A[1][2]*A[1][2] + A[2][2]*A[2][2] )

/* D[][] := A[][] * B[][] * C[][] */
void M3mul2 (double A[3][3], double B[3][3], double C[3][3], double D[3][3]);
#define M3MUL2(A,B,C,D,TMP) (M3MUL(A,B,TMP),M3MUL(TMP,C,D))

/* matrix & matrix & scalar */

/* C[] := aa * A[] + bb * B[] */
#define M3ADDMULMUL(aa,A,bb,B,C) (                                      \
        C[0][0] = (aa)*A[0][0]+(bb)*B[0][0],                            \
        C[0][1] = (aa)*A[0][1]+(bb)*B[0][1],                            \
        C[0][2] = (aa)*A[0][2]+(bb)*B[0][2],                            \
        C[1][0] = (aa)*A[1][0]+(bb)*B[1][0],                            \
        C[1][1] = (aa)*A[1][1]+(bb)*B[1][1],                            \
        C[1][2] = (aa)*A[1][2]+(bb)*B[1][2],                            \
        C[2][0] = (aa)*A[2][0]+(bb)*B[2][0],                            \
        C[2][1] = (aa)*A[2][1]+(bb)*B[2][1],                            \
        C[2][2] = (aa)*A[2][2]+(bb)*B[2][2] )


/* strain.c: */

/****************************/
/* Properties of H matrices */
/****************************/

/* calibrated to pure pressure */
#define SymmetricM3HydroInvariant(A) (M3TR(A)/3)
/* calibrated to single shear */
#define SymmetricM3MisesInvariant(A) \
  sqrt( SQUARE(A[0][1]) + SQUARE(A[0][2]) + SQUARE(A[1][2]) + \
  (SQUARE(A[0][0]-A[1][1]) + SQUARE(A[0][0]-A[2][2]) + \
  SQUARE(A[1][1]-A[2][2]))/6. )

/* Crystallographic Notation: http://spot.colorado.edu/~smyth/G30102.html */
typedef struct
{
    double a;
    double b;
    double c;
    double alpha;  /* in degrees */
    double beta;   /* in degrees */
    double gamma;  /* in degrees */
} Crystallographic;

#define CrystallographicAssign(A,B,C,ALPHA,BETA,GAMMA,X) (  (X).a=(A), \
  (X).b=(B), (X).c=(C), (X).alpha=(ALPHA), (X).beta=(BETA), (X).gamma=(GAMMA) )

/* construct crystallographic notation from H[][] */
#define H_to_Crystallographic(H,X) ( \
  (X).a = V3LENGTH(H[0]), (X).b = V3LENGTH(H[1]), (X).c = V3LENGTH(H[2]), \
  (X).alpha = RADIAN_TO_DEGREE( acos(V3DOT(H[1],H[2])/(X).b/(X).c) ), \
  (X).beta  = RADIAN_TO_DEGREE( acos(V3DOT(H[0],H[2])/(X).a/(X).c) ), \
  (X).gamma = RADIAN_TO_DEGREE( acos(V3DOT(H[0],H[1])/(X).a/(X).b) ) )

/* construct H[][] from crystallographic notation */
#define Crystallographic_to_H(X,H) ( \
  H[0][0] = (X).a, H[0][1] = 0, H[0][2] = 0, \
  H[1][0] = cos( DEGREE_TO_RADIAN((X).gamma) ) * (X).b, \
  H[1][1] = sin( DEGREE_TO_RADIAN((X).gamma) ) * (X).b, H[1][2] = 0., \
  H[2][0] = (X).c * cos( DEGREE_TO_RADIAN((X).beta) ), \
  H[2][1] = ( (X).c * (X).b * cos( DEGREE_TO_RADIAN((X).alpha) ) - \
  H[2][0] * H[1][0] ) / H[1][1], \
  H[2][2] = sqrt( SQUARE((X).c) - SQUARE(H[2][0]) - SQUARE(H[2][1]) ) )

/* determine which of the seven vertices of a col-parallelepiped */
/* formed by A[][0], A[][1], A[][2] is the farthest from origin. */
double M3maxcolumnradius (double H[3][3]);

/* determine which of the seven vertices of a row-parallelepiped */
/* formed by A[0][], A[1][], A[2][] is the farthest from origin. */
double M3maxrowradius (double H[3][3]);

/* returns the thickness (>=0) of the parallelepiped formed */
/* by (H[0][], H[1][], H[2][]) in the row i direction.      */
double M3rowthickness (double H[3][3], int i);

/* returns the three thicknesses (>=0) of the parallelepiped */
/* formed by (H[0][], H[1][], H[2][]), in thickness[].       */
double *M3rowthicknesses (double H[3][3], double thickness[3]);

/* Calculate multiplication factors nc[] whereby -nc[0]:nc[0],   */
/* -nc[1]:nc[1], -nc[2]:nc[2] replica of H[0..2][] is guaranteed */
/* to cover sphere of radius R; returns total number of replicas */
int M3rows_to_cover_sphere (double H[3][3], double R, int nc[3]);
#define total_replica(nc) ( (2*(nc)[0]+1) * (2*(nc)[1]+1) * (2*(nc)[2]+1) )

/* return the thickness (>=0) of the parallelepiped formed */
/* by (H[][0], H[][1], H[][2]) in the column i direction.  */
double M3columnthickness (double H[3][3], int i);

/* return the three thicknesses (>=0) of the parallelepiped */
/* formed by (H[][0], H[][1], H[][2]), in thickness[].      */
double *M3columnthicknesses (double H[3][3], double thickness[3]);

/* Calculate multiplication factors nc[] whereby -nc[0]:nc[0],   */
/* -nc[1]:nc[1], -nc[2]:nc[2] replica of H[][0..2] is guaranteed */
/* to cover sphere of radius R; return total number of replicas. */
int M3columns_to_cover_sphere (double H[3][3], double R, int nc[3]);

/* dl^2 = dx_i * (1 + 2 * eta_{ij}) * dx_j, where dx meshes the undeformed */
/* because dx = ds * H0, dx' = ds * H, dl^2 = dx' * (dx')^T => above       */
void Lagrangian_strain (double H0[3][3], double H[3][3], double eta[3][3]);

/* achieve eta without rotation. M := sqrt(1 + 2*eta), H := H0 * M */
void pure_deform (double H0[3][3], double eta[3][3], double H[3][3]);
/* H := H * M */
#define pure_DEFORM(H,eta) pure_deform(H,eta,H)

/* D = H0^-1*H - 1 */
void simple_D (double H0[3][3], double H[3][3], double D[3][3]);

/* H = H0 * (1+D) */
void simple_deform (double H0[3][3], double D[3][3], double H[3][3]);

/* simple strain = (J+J^T)/2-1 = (D+D^T)/2 */
void simple_strain (double H0[3][3], double H[3][3], double eta[3][3]);


/* rotation.c: */

/****************************************/
/* rotation generator & representations */
/****************************************/

/* this is cubic random, NOT spherical random */
#define M3RANDOMROTATION(R,tmp) ( V3FRANDOM(R[0]), \
  V3NORMALIZE(R[0],tmp), V3FRANDOM(R[1]), V3CROSS(R[0],R[1],R[2]), \
  V3NORMALIZE(R[2],tmp), V3CROSS(R[2],R[0],R[1]) )

/* Rotate v[] around axis[] by theta: store and return in w[] */
double *V3axialrotate (double axis[3], double theta, double v[3], double w[3]);

/* Compute rotational matrix R[][] corresponding */
/* to rotation with respect to axis[] by theta.  */
void V3axialrotatematrix (double axis[3], double theta, double R[3][3]);

/* if geodesic ("big-circle") rotation makes a[]->b[], what happens to v[]? */
/* Assume a[] and b[] are already NORMALIZED and NOT collinear; returns v[] */
double *V3geodesic (double a[3], double b[3], double v[3]);

/* Compute the rotational matrix R[][] corresponding to a geodesic */
/* ("big-circle") rotation that makes a[]->b[] as v[] := v[]*R[][] */
/* Assume a[] and b[] are already NORMALIZED and NOT collinear.    */
void M3geodesic (double a[3], double b[3], double R[3][3]);

/* Generate rotational matrix in axis-"i" of angle "theta". */
#define M3axialrmat(i,theta,R) \
  if ((i)==0) \
  M3ASSIGN(1,0,0,0,cos(theta),sin(theta),0,-sin(theta),cos(theta),R); \
  else if ((i)==1) \
  M3ASSIGN(cos(theta),0,sin(theta),0,1,0,-sin(theta),0,cos(theta),R); \
  else \
  M3ASSIGN(cos(theta),sin(theta),0,-sin(theta),cos(theta),0,0,0,1,R)

#endif  /* _VecMat3_h */
