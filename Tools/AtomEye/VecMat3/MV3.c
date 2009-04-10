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


/***************************/
/* column space operations */
/***************************/


/* columnlength[i] := | A[][i] |; returns columnlength[] */
double *M3columnlengths (double A[3][3], double columnlength[3])
{
    columnlength[0] = DISTANCE (A[0][0],A[1][0],A[2][0]);
    columnlength[1] = DISTANCE (A[0][1],A[1][1],A[2][1]);
    columnlength[2] = DISTANCE (A[0][2],A[1][2],A[2][2]);
    return(columnlength);
} /* end M3columnlengths() */


/* returns the maximum Euclidean length of A[][0], A[][1], A[][2] */
double M3maxcolumnlength (double A[3][3])
{
    double length0, length1, length2, tmp;
    length0 = DISTANCE (A[0][0],A[1][0],A[2][0]);
    length1 = DISTANCE (A[0][1],A[1][1],A[2][1]);
    length2 = DISTANCE (A[0][2],A[1][2],A[2][2]);
    tmp = MAX(length0,length1);
    return (MAX(tmp,length2));
} /* end M3maxcolumnlength() */


/* column[] := A[][i]; return column[] */
double *M3column (double A[3][3], int i, double column[3])
{
    column[0] = A[0][i];
    column[1] = A[1][i];
    column[2] = A[2][i];
    return(column);
} /* end M3column() */


/* c[] := A[][] * b[]; then return c[] */
double *M3mulV3 (double A[3][3], double b[3], double c[3])
{
    c[0] = A[0][0]*b[0] + A[0][1]*b[1] + A[0][2]*b[2];
    c[1] = A[1][0]*b[0] + A[1][1]*b[1] + A[1][2]*b[2];
    c[2] = A[2][0]*b[0] + A[2][1]*b[1] + A[2][2]*b[2];
    return(c);
} /* end M3mulV3() */


/* b[] := A[][] * b[]; then return b[] */
double *M3MULV3 (double A[3][3], double b[3])
{
    double a[3];
    a[0] = b[0];
    a[1] = b[1];
    a[2] = b[2];
    b[0] = A[0][0]*a[0] + A[0][1]*a[1] + A[0][2]*a[2];
    b[1] = A[1][0]*a[0] + A[1][1]*a[1] + A[1][2]*a[2];
    b[2] = A[2][0]*a[0] + A[2][1]*a[1] + A[2][2]*a[2];
    return(b);
} /* end M3MULV3() */


/* d[] := A[][] * b[] + c[]; then return d[] */
double *M3muladdV3 (double A[3][3], double b[3], double c[3], double d[3])
{
    d[0] = A[0][0]*b[0] + A[0][1]*b[1] + A[0][2]*b[2] + c[0];
    d[1] = A[1][0]*b[0] + A[1][1]*b[1] + A[1][2]*b[2] + c[1];
    d[2] = A[2][0]*b[0] + A[2][1]*b[1] + A[2][2]*b[2] + c[2];
    return(d);
} /* end M3muladdV3() */


/* c[] := c[] + A[][] * b[]; then return c[] */
double *M3mulADDV3 (double A[3][3], double b[3], double c[3])
{
    c[0] += A[0][0]*b[0] + A[0][1]*b[1] + A[0][2]*b[2];
    c[1] += A[1][0]*b[0] + A[1][1]*b[1] + A[1][2]*b[2];
    c[2] += A[2][0]*b[0] + A[2][1]*b[1] + A[2][2]*b[2];
    return(c);
} /* end M3mulADDV3() */


/* d[] := A[][] * b[] - c[]; then return d[] */
double *M3mulsubV3 (double A[3][3], double b[3], double c[3], double d[3])
{
    d[0] = A[0][0]*b[0] + A[0][1]*b[1] + A[0][2]*b[2] - c[0];
    d[1] = A[1][0]*b[0] + A[1][1]*b[1] + A[1][2]*b[2] - c[1];
    d[2] = A[2][0]*b[0] + A[2][1]*b[1] + A[2][2]*b[2] - c[2];
    return(d);
} /* end M3mulsubV3() */


/* d[] := c[] - A[][] * b[]; then return d[] */
double *M3submulV3 (double c[3], double A[3][3], double b[3], double d[3])
{
    d[0] = c[0] - A[0][0]*b[0] - A[0][1]*b[1] - A[0][2]*b[2];
    d[1] = c[1] - A[1][0]*b[0] - A[1][1]*b[1] - A[1][2]*b[2];
    d[2] = c[2] - A[2][0]*b[0] - A[2][1]*b[1] - A[2][2]*b[2];
    return(d);
} /* end M3submulV3() */


/* c[] := c[] - A[][] * b[]; then return c[] */
double *M3mulSUBV3 (double A[3][3], double b[3], double c[3])
{
    c[0] -= A[0][0]*b[0] + A[0][1]*b[1] + A[0][2]*b[2];
    c[1] -= A[1][0]*b[0] + A[1][1]*b[1] + A[1][2]*b[2];
    c[2] -= A[2][0]*b[0] + A[2][1]*b[1] + A[2][2]*b[2];
    return(c);
} /* end M3mulSUBV3() */
