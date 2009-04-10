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


/************************/
/* row space operations */
/************************/


/* rowlength[i] := | A[i][] |; returns rowlength[] */
double *M3rowlengths (double A[3][3], double rowlength[3])
{
    rowlength[0] = V3LENGTH (A[0]);
    rowlength[1] = V3LENGTH (A[1]);
    rowlength[2] = V3LENGTH (A[2]);
    return (rowlength);
} /* end M3rowlengths() */


/* returns the maximum Euclidean length of A[0][], A[1][], A[2][] */
double M3maxrowlength (double A[3][3])
{
    double length0, length1, length2, tmp;
    length0 = V3LENGTH (A[0]);
    length1 = V3LENGTH (A[1]);
    length2 = V3LENGTH (A[2]);
    tmp = MAX(length0,length1);
    return (MAX(tmp,length2));
} /* end M3maxrowlength() */


/* c[] := a[] * B[][]; then return c[] */
double *V3mulM3 (double a[3], double B[3][3], double c[3])
{
    c[0] = a[0]*B[0][0] + a[1]*B[1][0] + a[2]*B[2][0];
    c[1] = a[0]*B[0][1] + a[1]*B[1][1] + a[2]*B[2][1];
    c[2] = a[0]*B[0][2] + a[1]*B[1][2] + a[2]*B[2][2];
    return(c);
} /* end V3mulM3() */


/* a[] := a[] * B[][]; then return a[] */
double *V3MULM3 (double a[3], double B[3][3])
{
    double c[3];
    c[0] = a[0];
    c[1] = a[1];
    c[2] = a[2];
    a[0] = c[0]*B[0][0] + c[1]*B[1][0] + c[2]*B[2][0];
    a[1] = c[0]*B[0][1] + c[1]*B[1][1] + c[2]*B[2][1];
    a[2] = c[0]*B[0][2] + c[1]*B[1][2] + c[2]*B[2][2];
    return(a);
} /* end V3MULM3() */


/* d[] := a[] + b[] * C[][]; then return d[] */
double *V3addmulM3 (double a[3], double b[3], double C[3][3], double d[3])
{
    d[0] = a[0] + b[0]*C[0][0] + b[1]*C[1][0] + b[2]*C[2][0];
    d[1] = a[1] + b[0]*C[0][1] + b[1]*C[1][1] + b[2]*C[2][1];
    d[2] = a[2] + b[0]*C[0][2] + b[1]*C[1][2] + b[2]*C[2][2];
    return(d);
} /* end V3addmulM3() */


/* a[] := a[] + b[] * C[][]; then return a[] */
double *V3ADDmulM3 (double a[3], double b[3], double C[3][3])
{
    a[0] += b[0]*C[0][0] + b[1]*C[1][0] + b[2]*C[2][0];
    a[1] += b[0]*C[0][1] + b[1]*C[1][1] + b[2]*C[2][1];
    a[2] += b[0]*C[0][2] + b[1]*C[1][2] + b[2]*C[2][2];
    return(a);
} /* end V3ADDmulM3() */


/* d[] := a[] - b[] * C[][]; then return d[] */
double *V3submulM3 (double a[3], double b[3], double C[3][3], double d[3])
{
    d[0] = a[0] - b[0]*C[0][0] - b[1]*C[1][0] - b[2]*C[2][0];
    d[1] = a[1] - b[0]*C[0][1] - b[1]*C[1][1] - b[2]*C[2][1];
    d[2] = a[2] - b[0]*C[0][2] - b[1]*C[1][2] - b[2]*C[2][2];
    return(d);
} /* end V3submulM3() */


/* d[] := b[] * C[][] - a[]; then return d[] */
double *V3mulsubM3 (double b[3], double C[3][3], double a[3], double d[3])
{
    d[0] = b[0]*C[0][0] + b[1]*C[1][0] + b[2]*C[2][0] - a[0];
    d[1] = b[0]*C[0][1] + b[1]*C[1][1] + b[2]*C[2][1] - a[1];
    d[2] = b[0]*C[0][2] + b[1]*C[1][2] + b[2]*C[2][2] - a[2];
    return(d);
} /* end V3mulsubM3() */


/* a[] := a[] - b[] * C[][]; then return a[] */
double *V3SUBmulM3 (double a[3], double b[3], double C[3][3])
{
    a[0] -= b[0]*C[0][0] + b[1]*C[1][0] + b[2]*C[2][0];
    a[1] -= b[0]*C[0][1] + b[1]*C[1][1] + b[2]*C[2][1];
    a[2] -= b[0]*C[0][2] + b[1]*C[1][2] + b[2]*C[2][2];
    return(a);
} /* end V3SUBmulM3() */
