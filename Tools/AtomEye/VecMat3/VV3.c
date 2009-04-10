/********************************************/
/* libVecMat3: -llapck -lblas -lm           */
/*             -lScalar -lIO                */
/*                                          */
/* Three-dimensional Euclidean space vector */
/* and matrix library.                      */
/*                                          */
/* Nov 11 1999 Ju Li <liju99@mit.edu>       */
/********************************************/

#include "VecMat3.h"


/*******************************/
/* vector & vector with scalar */
/*******************************/

/* vector & vector */

/* c[] := a[] + b[]; then return c[] */
double *V3add (double a[3], double b[3], double c[3])
{
    c[0] = a[0] + b[0];
    c[1] = a[1] + b[1];
    c[2] = a[2] + b[2];
    return (c);
} /* end V3add() */


/* b[] := a[] + b[]; then return b[] */
double *V3Add (double a[3], double b[3])
{
    b[0] += a[0];
    b[1] += a[1];
    b[2] += a[2];
    return (b);
} /* end V3Add() */


/* c[] := a[] - b[]; then return c[] */
double *V3sub (double a[3], double b[3], double c[3])
{
    c[0] = a[0] - b[0];
    c[1] = a[1] - b[1];
    c[2] = a[2] - b[2];
    return (c);
} /* end V3sub() */


/* a[] := a[] - b[]; then return a[] */
double *V3Sub (double a[3], double b[3])
{
    a[0] -= b[0];
    a[1] -= b[1];
    a[2] -= b[2];
    return (a);
} /* end V3SUB() */


/* dot product of a[] and b[] */
double V3dot (double a[3], double b[3])
{
    return (a[0]*b[0]+a[1]*b[1]+a[2]*b[2]);
} /* end V3dot() */


/* cross product of a[] and b[]: c := a x b; then return c[] */
double *V3cross (double a[3], double b[3], double c[3])
{
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
    return(c);
} /* end V3cross() */


/* angle between vectors a[],b[] in radian [0,pi] */
double V3angle (double a[3], double b[3])
{
    register double aa, bb, cc;
    aa = V3LENGTH(a);
    bb = V3LENGTH(b);
    if ( (aa == 0) || (bb == 0) ) pe ("V3angle: a[]=0 or b[]=0.\n");
    cc = V3DOT(a,b) / aa / bb;
    if (cc > 1) cc = 1; else if (cc < -1) cc = -1;
    return ( acos( cc ) );
} /* end V3angle() */


/* c[i] := a[i] * b[i], i=0..2; then return c[] */
double *V3prod (double a[3], double b[3], double c[3])
{
    c[0] = a[0] * b[0];
    c[1] = a[1] * b[1];
    c[2] = a[2] * b[2];
    return(c);
} /* end V3prod() */


/* a[i] := a[i] * b[i], i=0..2; then return a[] */
double *V3Prod (double a[3], double b[3])
{
    a[0] *= b[0];
    a[1] *= b[1];
    a[2] *= b[2];
    return(a);
} /* end V3Prod() */


/* vector & vector & scalar */



/* c[] := a[] + multiplier * b[]; then return c[] */
double *V3addmul (double a[3], double multiplier, double b[3], double c[3])
{
    c[0] = a[0] + multiplier * b[0];
    c[1] = a[1] + multiplier * b[1];
    c[2] = a[2] + multiplier * b[2];
    return (c);
} /* end V3addmul() */


/* c[] := a[] + b[] / divisor; then return c[] */
double *V3adddiv (double a[3], double b[3], double divisor, double c[3])
{
    if (divisor == 0.)
    {
        printf ("error: V3adddiv: divisor = %e\n", divisor);
        exit(1);
    }
    c[0] = a[0] + b[0] / divisor;
    c[1] = a[1] + b[1] / divisor;
    c[2] = a[2] + b[2] / divisor;
    return (c);
} /* end V3adddiv() */


/* b[] := b[] + multiplier * a[]; then return b[] */
double *V3ADDmul (double multiplier, double a[3], double b[3])
{
    b[0] += multiplier * a[0];
    b[1] += multiplier * a[1];
    b[2] += multiplier * a[2];
    return (b);
} /* end V3ADDmul() */


/* c[] := c[] + aa * a[] + bb * b[]; then return c[] */
double *V3ADDmulmul (double aa,double a[3], double bb,double b[3], double c[3])
{
    c[0] += aa*a[0] + bb*b[0];
    c[1] += aa*a[1] + bb*b[1];
    c[2] += aa*a[2] + bb*b[2];
    return (c);
} /* end V3ADDmulmul() */


/* b[] := b[] + a[] / divisor; then return b[] */
double *V3ADDdiv (double a[3], double divisor, double b[3])
{
    if (divisor == 0.)
    {
        printf ("error: V3ADDdiv: divisor = %e\n", divisor);
        exit(1);
    }
    b[0] += a[0] / divisor;
    b[1] += a[1] / divisor;
    b[2] += a[2] / divisor;
    return (b);
} /* end V3ADDdiv() */


/* c[] := a[] - multiplier * b[]; then return c[] */
double *V3submul (double a[3], double multiplier, double b[3], double c[3])
{
    c[0] = a[0] - multiplier * b[0];
    c[1] = a[1] - multiplier * b[1];
    c[2] = a[2] - multiplier * b[2];
    return (c);
} /* end V3submul() */


/* c[] := a[] - b[] / divisor; then return c[] */
double *V3subdiv (double a[3], double b[3], double divisor, double c[3])
{
    if (divisor == 0.)
    {
        printf ("error: V3subdiv: divisor = %e\n", divisor);
        exit(1);
    }
    c[0] = a[0] - b[0] / divisor;
    c[1] = a[1] - b[1] / divisor;
    c[2] = a[2] - b[2] / divisor;
    return (c);
} /* end V3subdiv() */


/* a[] := a[] - multiplier * b[]; then return a[] */
double *V3SUBmul (double a[3], double multiplier, double b[3])
{
    a[0] -= multiplier * b[0];
    a[1] -= multiplier * b[1];
    a[2] -= multiplier * b[2];
    return (a);
} /* end V3SUBmul() */


/* a[] := a[] - b[] / divisor; then return a[] */
double *V3SUBdiv (double a[3], double b[3], double divisor)
{
    if (divisor == 0.)
    {
        printf ("error: V3SUBdiv: divisor = %e\n", divisor);
        exit(1);
    }
    a[0] -= b[0] / divisor;
    a[1] -= b[1] / divisor;
    a[2] -= b[2] / divisor;
    return (a);
} /* end V3SUBdiv() */


/* c[] := part of a[] that is perpendicular to b[]; return c[] */
double *V3perpendicular (double a[3], double b[3], double c[3])
{
    double b2 = V3LENGTH2 (b);
    if ( b2 == 0 )
    {
        printf ("error: V3perpendicular: b[] = 0\n");
        exit(1);
    }
    return ( V3submul(a, V3DOT(a,b)/b2, b, c) );
} /* end V3perpendicular() */


/* a[] := part of a[] that is perpendicular to b[]; return modified a[] */
double *V3Perpendicular (double a[3], double b[3])
{
    double b2 = V3LENGTH2 (b);
    if ( b2 == 0 )
    {
        printf ("error: V3Perpendicular: b[] = 0\n");
        exit(1);
    }
    return ( V3SUBmul(a, V3DOT(a,b)/b2, b) );
} /* end V3Perpendicular() */
