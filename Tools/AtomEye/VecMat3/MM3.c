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


/********************************/
/* matrix & matrix with scalars */
/********************************/


/* matrix & matrix */


/* C[][] := A[][] + B[][] */
void M3add (double A[3][3], double B[3][3], double C[3][3])
{
    C[0][0] = A[0][0] + B[0][0];
    C[0][1] = A[0][1] + B[0][1];
    C[0][2] = A[0][2] + B[0][2];
    C[1][0] = A[1][0] + B[1][0];
    C[1][1] = A[1][1] + B[1][1];       
    C[1][2] = A[1][2] + B[1][2];       
    C[2][0] = A[2][0] + B[2][0];           
    C[2][1] = A[2][1] + B[2][1];      
    C[2][2] = A[2][2] + B[2][2];
    return;
} /* end M3add() */


/* C[][] := A[][] - B[][] */
void M3sub (double A[3][3], double B[3][3], double C[3][3])
{
    C[0][0] = A[0][0] - B[0][0];
    C[0][1] = A[0][1] - B[0][1];
    C[0][2] = A[0][2] - B[0][2]; 
    C[1][0] = A[1][0] - B[1][0];         
    C[1][1] = A[1][1] - B[1][1];       
    C[1][2] = A[1][2] - B[1][2];       
    C[2][0] = A[2][0] - B[2][0];           
    C[2][1] = A[2][1] - B[2][1];      
    C[2][2] = A[2][2] - B[2][2];
    return;
} /* end M3sub() */


/* C[][] := A[][] * B[][] */
void M3mul (double A[3][3], double B[3][3], double C[3][3])
{
    C[0][0] = A[0][0]*B[0][0] + A[0][1]*B[1][0] + A[0][2]*B[2][0];
    C[0][1] = A[0][0]*B[0][1] + A[0][1]*B[1][1] + A[0][2]*B[2][1];
    C[0][2] = A[0][0]*B[0][2] + A[0][1]*B[1][2] + A[0][2]*B[2][2];
    C[1][0] = A[1][0]*B[0][0] + A[1][1]*B[1][0] + A[1][2]*B[2][0];
    C[1][1] = A[1][0]*B[0][1] + A[1][1]*B[1][1] + A[1][2]*B[2][1];
    C[1][2] = A[1][0]*B[0][2] + A[1][1]*B[1][2] + A[1][2]*B[2][2];
    C[2][0] = A[2][0]*B[0][0] + A[2][1]*B[1][0] + A[2][2]*B[2][0];
    C[2][1] = A[2][0]*B[0][1] + A[2][1]*B[1][1] + A[2][2]*B[2][1];
    C[2][2] = A[2][0]*B[0][2] + A[2][1]*B[1][2] + A[2][2]*B[2][2];
    return;
} /* end M3mul() */


/* A[][] := A[][] * B[][] */
void M3Mul (double A[3][3], double B[3][3])
{
    M3 tmp;
    M3MUl (A,B,tmp);
    return;
} /* end M3Mul() */


/* B[][] := A[][] * B[][] */
void M3muL (double A[3][3], double B[3][3])
{
    M3 tmp;
    M3mUL (A,B,tmp);
    return;
} /* end M3muL() */


/* D[][] := A[][] * B[][] * C[][] */
void M3mul2 (double A[3][3], double B[3][3], double C[3][3], double D[3][3])
{
    static double E[3][3];
    M3MUL (A, B, E);
    M3MUL (E, C, D);
    return;
} /* end M3mul2() */


/* matrix & matrix & scalar */

