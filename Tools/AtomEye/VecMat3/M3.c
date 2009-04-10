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


/**********************************************/
/* single matrix & scalar & generating vector */
/**********************************************/


/* matrix */


/* generate a zero matrix A[][] := 0 */
void M3zero (double A[3][3])
{
    A[0][0] = 0.;
    A[0][1] = 0.;
    A[0][2] = 0.;
    A[1][0] = 0.;
    A[1][1] = 0.;
    A[1][2] = 0.;
    A[2][0] = 0.;
    A[2][1] = 0.;
    A[2][2] = 0.;
    return;
} /* end M3zero() */


/* generate an identity matrix A[][] := I[][] */
void M3identity (double A[3][3])
{
    A[0][0] = 1.;
    A[0][1] = 0.;
    A[0][2] = 0.;
    A[1][0] = 0.;
    A[1][1] = 1.;
    A[1][2] = 0.;
    A[2][0] = 0.;
    A[2][1] = 0.;
    A[2][2] = 1.;
    return;
} /* end M3identity() */


/* generate A[][] with 9 independent random components on (0.,1.) */
void M3frandom (double A[3][3])
{
    A[0][0] = Frandom();
    A[0][1] = Frandom();
    A[0][2] = Frandom();
    A[1][0] = Frandom();
    A[1][1] = Frandom();
    A[1][2] = Frandom();
    A[2][0] = Frandom();
    A[2][1] = Frandom();
    A[2][2] = Frandom();
    return;
} /* end M3frandom() */


/* B[][] := A[][] */
void M3eqv (double A[3][3], double B[3][3])
{
    B[0][0] = A[0][0];
    B[0][1] = A[0][1];
    B[0][2] = A[0][2];
    B[1][0] = A[1][0];
    B[1][1] = A[1][1];
    B[1][2] = A[1][2];
    B[2][0] = A[2][0];
    B[2][1] = A[2][1];
    B[2][2] = A[2][2];
    return;
} /* end M3eqv() */


/* arbitrary exponent-th norm of A[][] (all components) */
double M3norm (double A[3][3], double exponent)
{
    double tmp;
    if (exponent <= 0.)
        pe ("M3norm: norm exponent = %lf <= 0 is illegal\n", exponent);
    else if ( exponent >= EXPONENT_INFINITY )
    {
	tmp = ABS( A[0][0] );
	if ( A[0][1] > tmp ) tmp = A[0][1];
	else if ( -A[0][1] > tmp ) tmp = -A[0][1];
        if ( A[0][2] > tmp ) tmp = A[0][2];
	else if ( -A[0][2] > tmp ) tmp = -A[0][2];
	if ( A[1][0] > tmp ) tmp = A[1][0];
	else if ( -A[1][0] > tmp ) tmp = -A[1][0];
	if ( A[1][1] > tmp ) tmp = A[1][1];
	else if ( -A[1][1] > tmp ) tmp = -A[1][1];
	if ( A[1][2] > tmp ) tmp = A[1][2];
	else if ( -A[1][2] > tmp ) tmp = -A[1][2];
	if ( A[2][0] > tmp ) tmp = A[2][0];
	else if ( -A[2][0] > tmp ) tmp = -A[2][0];
	if ( A[2][1] > tmp ) tmp = A[2][1];
	else if ( -A[2][1] > tmp ) tmp = -A[2][1];
	if ( A[2][2] > tmp ) tmp = A[2][2];
	else if ( -A[2][2] > tmp ) tmp = -A[2][2];
	return (tmp);
    }
    else
    {
	tmp =
	    pow(ABS(A[0][0]), exponent) +
	    pow(ABS(A[0][1]), exponent) +
	    pow(ABS(A[0][2]), exponent) +
	    pow(ABS(A[1][0]), exponent) +
	    pow(ABS(A[1][1]), exponent) +
	    pow(ABS(A[1][2]), exponent) +
	    pow(ABS(A[2][0]), exponent) +
	    pow(ABS(A[2][1]), exponent) +
	    pow(ABS(A[2][2]), exponent);
	return (pow(tmp, 1./exponent));
    }
    return (0.);
} /* end M3norm() */


/* return the largest component of A[][] in absolute value */
double M3infnorm (double A[3][3])
{
    register double tmp = ABS( A[0][0] );
    if ( A[0][1] > tmp ) tmp = A[0][1];
    else if ( -A[0][1] > tmp ) tmp = -A[0][1];
    if ( A[0][2] > tmp ) tmp = A[0][2];
    else if ( -A[0][2] > tmp ) tmp = -A[0][2];
    if ( A[1][0] > tmp ) tmp = A[1][0];
    else if ( -A[1][0] > tmp ) tmp = -A[1][0];
    if ( A[1][1] > tmp ) tmp = A[1][1];
    else if ( -A[1][1] > tmp ) tmp = -A[1][1];
    if ( A[1][2] > tmp ) tmp = A[1][2];
    else if ( -A[1][2] > tmp ) tmp = -A[1][2];
    if ( A[2][0] > tmp ) tmp = A[2][0];
    else if ( -A[2][0] > tmp ) tmp = -A[2][0];
    if ( A[2][1] > tmp ) tmp = A[2][1];
    else if ( -A[2][1] > tmp ) tmp = -A[2][1];
    if ( A[2][2] > tmp ) tmp = A[2][2];
    else if ( -A[2][2] > tmp ) tmp = -A[2][2];
    return (tmp);
} /* end M3infnorm() */


/* return sqrt{ \sum_{i=0..2} \sum_{j=0..2} |A_ij|^2 } */
double M32norm (double A[3][3])
{
    return (sqrt(V3LENGTH2(A[0])+V3LENGTH2(A[1])+V3LENGTH2(A[2])));
} /* end M32norm() */


/* B[][] := -A[][] */
void M3neg (double A[3][3], double B[3][3])
{
    B[0][0] = -A[0][0];
    B[0][1] = -A[0][1];
    B[0][2] = -A[0][2];
    B[1][0] = -A[1][0];
    B[1][1] = -A[1][1];
    B[1][2] = -A[1][2];
    B[2][0] = -A[2][0];
    B[2][1] = -A[2][1];
    B[2][2] = -A[2][2];
    return;
} /* end M3neg() */


/* A[][] := -A[][] */
void M3Neg (double A[3][3])
{
    A[0][0] = -A[0][0];
    A[0][1] = -A[0][1];
    A[0][2] = -A[0][2];
    A[1][0] = -A[1][0];
    A[1][1] = -A[1][1];
    A[1][2] = -A[1][2];
    A[2][0] = -A[2][0];
    A[2][1] = -A[2][1];
    A[2][2] = -A[2][2];
    return;
} /* end M3Neg() */


/* B[][] := A[][]' */
void M3transpose (double A[3][3], double B[3][3])
{
    B[0][0] = A[0][0];
    B[0][1] = A[1][0];
    B[0][2] = A[2][0];
    B[1][0] = A[0][1];
    B[1][1] = A[1][1];
    B[1][2] = A[2][1];
    B[2][0] = A[0][2];
    B[2][1] = A[1][2];
    B[2][2] = A[2][2];
    return;
} /* end M3transpose() */


/* A[][] := A[][]' */
void M3Transpose (double A[3][3])
{
    double tmp;
    tmp = A[0][1];
    A[0][1] = A[1][0];
    A[1][0] = tmp;
    tmp = A[0][2];
    A[0][2] = A[2][0];
    A[2][0] = tmp;
    tmp = A[1][2];
    A[1][2] = A[2][1];
    A[2][1] = tmp;
    return;
} /* end M3Transpose() */


/* B[][] := (A[][]+A[][]')/2 */
void M3symmetrize (double A[3][3], double B[3][3])
{
    B[0][0] = A[0][0];
    B[1][1] = A[1][1];
    B[2][2] = A[2][2];
    B[0][1] = B[1][0] = (A[1][0]+A[0][1])/2.;
    B[0][2] = B[2][0] = (A[2][0]+A[0][2])/2.;
    B[1][2] = B[2][1] = (A[2][1]+A[1][2])/2.;
    return;
} /* end M3symmetrize() */


/* return the trace of A[][] */
double M3Tr (double A[3][3])
{
    return (A[0][0]+A[1][1]+A[2][2]);
} /* end M3Tr() */


/* B[][] := trace(A[][]) / 3 * I[][]; return trace(A) */
double M3trace (double A[3][3], double B[3][3])
{
    double trace = A[0][0]+A[1][1]+A[2][2];
    B[0][0] = trace/3.;
    B[0][1] = 0.;
    B[0][2] = 0.;
    B[1][0] = 0.;
    B[1][1] = trace/3.;
    B[1][2] = 0.;
    B[2][0] = 0.;
    B[2][1] = 0.;
    B[2][2] = trace/3.;
    return (trace);
} /* return M3trace() */


/* A[][] := trace(A[][])/3 * I[][]; return original trace(A) */
double M3Trace (double A[3][3])
{
    double trace = A[0][0]+A[1][1]+A[2][2];
    A[0][0] = trace/3.;
    A[0][1] = 0.;
    A[0][2] = 0.;
    A[1][0] = 0.;
    A[1][1] = trace/3.;
    A[1][2] = 0.;
    A[2][0] = 0.;
    A[2][1] = 0.;
    A[2][2] = trace/3.;
    return (trace);
} /* return M3Trace() */


/* B[][] := A[][] - trace(A[][])/3 * I[][]; return trace(A) */
double M3traceless (double A[3][3], double B[3][3])
{
    double trace = A[0][0]+A[1][1]+A[2][2];
    B[0][0] = A[0][0] - trace/3.;
    B[0][1] = A[0][1];
    B[0][2] = A[0][2];
    B[1][0] = A[1][0];
    B[1][1] = A[1][1] - trace/3.;
    B[1][2] = A[1][2];
    B[2][0] = A[2][0];
    B[2][1] = A[2][1];
    B[2][2] = A[2][2] - trace/3.;
    return (trace);
} /* return M3traceless() */


/* A[][] := A[][] - trace(A[][])/3 * I[][]; return original trace(A) */
double M3Traceless (double A[3][3])
{
    double trace = A[0][0]+A[1][1]+A[2][2];
    A[0][0] -= trace/3.;
    A[1][1] -= trace/3.;
    A[2][2] -= trace/3.;
    return (trace);
} /* return M3Traceless() */


/* decompose A[][] to b*I[][] + C[][], where b := trace(A)/3; return b */
double M3tracedecompose (double A[3][3], double B[3][3], double C[3][3])
{
    double trace = A[0][0] + A[1][1] + A[2][2];
    B[0][0] = trace/3.;
    B[0][1] = 0.;
    B[0][2] = 0.;
    B[1][0] = 0.;
    B[1][1] = trace/3.;
    B[1][2] = 0.;
    B[2][0] = 0.;
    B[2][1] = 0.;
    B[2][2] = trace/3.;
    C[0][0] = A[0][0] - trace/3.;
    C[0][1] = A[0][1];
    C[0][2] = A[0][2];
    C[1][0] = A[1][0];
    C[1][1] = A[1][1] - trace/3.;
    C[1][2] = A[1][2];
    C[2][0] = A[2][0];
    C[2][1] = A[2][1];
    C[2][2] = A[2][2] - trace/3.;
    return (trace);
} /* end M3tracedecompose() */



/* matrix & scalar */



/* generate an identity matrix A[][] := a x I[][] */
void M3Identity (double a, double A[3][3])
{
    A[0][0] = a;
    A[0][1] = 0.;
    A[0][2] = 0.;
    A[1][0] = 0.;
    A[1][1] = a;
    A[1][2] = 0.;
    A[2][0] = 0.;
    A[2][1] = 0.;
    A[2][2] = a;
    return;
} /* end M3Identity() */


/* B[][] := multiplier * A[][] */
void M3multiply (double multiplier, double A[3][3], double B[3][3])
{
    B[0][0] = multiplier * A[0][0];
    B[0][1] = multiplier * A[0][1];
    B[0][2] = multiplier * A[0][2];
    B[1][0] = multiplier * A[1][0];
    B[1][1] = multiplier * A[1][1];
    B[1][2] = multiplier * A[1][2];
    B[2][0] = multiplier * A[2][0];
    B[2][1] = multiplier * A[2][1];
    B[2][2] = multiplier * A[2][2];
    return;
} /* end M3multiply() */


/* A[][] := multiplier * A[][] */
void M3Multiply (double multiplier, double A[3][3])
{
    A[0][0] *= multiplier;
    A[0][1] *= multiplier;
    A[0][2] *= multiplier;
    A[1][0] *= multiplier;
    A[1][1] *= multiplier;
    A[1][2] *= multiplier;
    A[2][0] *= multiplier;
    A[2][1] *= multiplier;
    A[2][2] *= multiplier;
    return;
} /* end M3Multiply() */


/* B[][] := A[][] / divisor */
void M3divide (double A[3][3], double divisor, double B[3][3])
{
    if (divisor == 0.) pe ("M3divide: divisor = %e\n", divisor);
    B[0][0] = A[0][0] / divisor;
    B[0][1] = A[0][1] / divisor;
    B[0][2] = A[0][2] / divisor;
    B[1][0] = A[1][0] / divisor;
    B[1][1] = A[1][1] / divisor;
    B[1][2] = A[1][2] / divisor;
    B[2][0] = A[2][0] / divisor;
    B[2][1] = A[2][1] / divisor;
    B[2][2] = A[2][2] / divisor;
    return;
} /* end M3divide() */


/* A[][] := A[][] / divisor */
void M3Divide (double A[3][3], double divisor)
{
    if (divisor == 0.) pe ("M3Divide: divisor = %e\n", divisor);
    A[0][0] /= divisor;
    A[0][1] /= divisor;
    A[0][2] /= divisor;
    A[1][0] /= divisor;
    A[1][1] /= divisor;
    A[1][2] /= divisor;
    A[2][0] /= divisor;
    A[2][1] /= divisor;
    A[2][2] /= divisor;
    return;
} /* end M3Divide() */


/* B[][] := A[][] + a * I */
void M3adddiag (double A[3][3], double a, double B[3][3])
{
    B[0][0] = A[0][0] + a;
    B[0][1] = A[0][1];
    B[0][2] = A[0][2];
    B[1][0] = A[1][0];
    B[1][1] = A[1][1] + a;
    B[1][2] = A[1][2];
    B[2][0] = A[2][0];
    B[2][1] = A[2][1];
    B[2][2] = A[2][2] + a;
    return;
} /* end M3adddiag() */


/* A[][] := A[][] + a * I */
void M3Adddiag (double A[3][3], double a)
{
    A[0][0] += a;
    A[1][1] += a;
    A[2][2] += a;
    return;
} /* end M3Adddiag() */


/* B[][] := A[][] - a * I */
void M3subdiag (double A[3][3], double a, double B[3][3])
{
    B[0][0] = A[0][0] - a;
    B[0][1] = A[0][1];
    B[0][2] = A[0][2];
    B[1][0] = A[1][0];
    B[1][1] = A[1][1] - a;
    B[1][2] = A[1][2];
    B[2][0] = A[2][0];
    B[2][1] = A[2][1];
    B[2][2] = A[2][2] - a;
    return;
} /* end M3subdiag() */


/* A[][] := A[][] - a * I */
void M3Subdiag (double A[3][3], double a)
{
    A[0][0] -= a;
    A[1][1] -= a;
    A[2][2] -= a;
    return;
} /* end M3Subdiag() */



/* matrix & generating vector */

