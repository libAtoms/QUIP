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


/*****************************/
/* single vector with scalar */
/*****************************/

/* unary vector */

/* a[] := 0; return a[] */
double *V3zero (double a[3])
{
    a[0] = 0;
    a[1] = 0;
    a[2] = 0;
    return(a);
} /* end V3zero() */


/* generate a unit vector: a[i] := 1, a[j!=i] := 0; return a[] */
double *V3unit (double a[3], int i)
{
    a[0] = 0;
    a[1] = 0;
    a[2] = 0;
    a[i] = 1;
    return(a);
} /* end V3unit() */


/* generate a[] with 3 independent random components on (0.,1.); return a[] */
double *V3frandom (double a[3])
{
    a[0] = Frandom();
    a[1] = Frandom();
    a[2] = Frandom();
    return(a);
} /* end V3frandom() */


/* generate a unit vector of spherically uniform random orientation */
double *V3randomunit (double a[3])
{
    double length;
    while(1)
    {
        a[0] = FRANDOM();
        a[1] = FRANDOM();
        a[2] = FRANDOM();
        length = V3LENGTH(a);
        if ( (length > SMALL) && (length < 0.5 - SMALL) )
        { /* ensure numerical accuracy */
            a[0] /= length;
            a[1] /= length;
            a[2] /= length;
            return(a);
        }
    }
} /* end V3randomunit() */


/* b[] := a[]; then return b[] */
double *V3eqv (double a[3], double b[3])
{
    b[0] = a[0];
    b[1] = a[1];
    b[2] = a[2];
    return(b);
} /* end V3eqv() */


/* b[] := -a[]; then return b[] */
double *V3neg (double a[3], double b[3])
{
    b[0] = -a[0];
    b[1] = -a[1];
    b[2] = -a[2];
    return(b);
} /* end V3neg() */


/* a[] := -a[]; then return a[] */
double *V3Neg (double a[3])
{
    a[0] = -a[0];
    a[1] = -a[1];
    a[2] = -a[2];
    return(a);
} /* end V3Neg() */


/* make b[] an image of a[] in [0,1)^3; then return b[] */
double *V3trim (double a[3], double b[3])
{
    b[0] = TRIM(a[0]);
    b[1] = TRIM(a[1]);
    b[2] = TRIM(a[2]);
    return(b);
} /* end V3trim() */


/* change a[] to its own image in [0,1)^3; then return a[] */
double *V3Trim (double a[3])
{
    a[0] = TRIM(a[0]);
    a[1] = TRIM(a[1]);
    a[2] = TRIM(a[2]);
    return(a);
} /* end V3Trim() */


/* make b[] image of a[]'s in [-0.5,0.5)^3; then return b[] */
double *V3image (double a[3], double b[3])
{
    b[0] = IMAGE(a[0]);
    b[1] = IMAGE(a[1]);
    b[2] = IMAGE(a[2]);
    return(b);
} /* end V3image() */


/* change a[] to its own image in [-0.5,0.5)^3; then return a[] */
double *V3Image (double a[3])
{
    a[0] = IMAGE(a[0]);
    a[1] = IMAGE(a[1]);
    a[2] = IMAGE(a[2]);
    return(a);
} /* end V3Image() */


/* scalar & vector */


/* b[] := multiplier * a[]; then return b[] */
double *V3mul (double multiplier, double a[3], double b[3])
{
    b[0] = multiplier * a[0];
    b[1] = multiplier * a[1];
    b[2] = multiplier * a[2];
    return(b);
} /* end V3mul() */


/* a[] := multiplier * a[]; then return a[] */
double *V3Mul (double multiplier, double a[3])
{
    a[0] *= multiplier;
    a[1] *= multiplier;
    a[2] *= multiplier;
    return(a);
} /* end V3Mul() */


/* b[] := a[] / divisor; then return b[] */
double *V3div (double a[3], double divisor, double b[3])
{
    if (divisor == 0.)
        pe ("V3div: divisor = %e\n", divisor);
    b[0] = a[0] / divisor;
    b[1] = a[1] / divisor;
    b[2] = a[2] / divisor;
    return(b);
} /* end V3div() */


/* a[] := a[] / divisor; then return a[] */
double *V3Div (double a[3], double divisor)
{
    if (divisor == 0.)
        pe ("V3DIV: divisor = %e\n", divisor);
    a[0] /= divisor;
    a[1] /= divisor;
    a[2] /= divisor;
    return(a);
} /* end V3DIV() */


/* length of a[] */
double V3length (double a[3])
{
    return (sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]));
} /* end V3length() */


/* squared length of a[] */
double V3length2 (double a[3])
{
    return (a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
} /* end V3length2() */


/* arbitrary exponent-th norm of a[] */
double V3norm (double a[3], double exponent)
{
    double tmp;
    if (exponent <= 0.)
        pe ("V3norm: norm exponent = %lf <= 0 is illegal\n", exponent);
    else if (exponent >= EXPONENT_INFINITY)
    {
        tmp = MAX(fabs(a[0]),fabs(a[1]));
        return (MAX(tmp,fabs(a[2])));
    }
    else
    {
        tmp = pow(fabs(a[0]), exponent) +
            pow(fabs(a[1]), exponent) +
            pow(fabs(a[2]), exponent); 
        return (pow(tmp, 1./exponent));
    }
    return (0.);
} /* end V3norm() */


/* normalize a[] to unit length; then return a[] */
double *V3normalize (double a[3])
{
    double length = sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
    if (length == 0.)
    {
        pr ("warning: V3normalize: length = 0,\n"
            "will not normalize input vector\n");
        return (a);
    }
    a[0] /= length;
    a[1] /= length;
    a[2] /= length;
    return (a);
} /* end V3normalize() */


/* normalize a[] to unit length; then return its original length */
double V3Normalize (double a[3])
{
    double length = sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
    if (length == 0.)
    {
        pr ("warning: V3NORMALIZE: length = 0,"
            "will not normalize input vector\n");
        return (length);
    }
    a[0] /= length;
    a[1] /= length;
    a[2] /= length;
    return (length);
} /* end V3Normalize() */


/* b[] := a[] / |a[]|, return b[] */
double *V3direction (double a[3], double b[3])
{
    double length = sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
    if (length == 0.)
    {
        pr ("warning: V3direction: length = 0,\n"
            "will not normalize input vector\n");
        return (b);
    }
    b[0] = a[0] / length;
    b[1] = a[1] / length;
    b[2] = a[2] / length;
    return (b);
} /* end V3direction() */


/* b[] := a[] / |a[]|, return |a[]| */
double V3Direction (double a[3], double b[3])
{
    double length = sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
    if (length == 0.)
    {
        pr ("warning: V3Direction: length = 0,\n"
            "will not normalize input vector\n");
        return (0);
    }
    b[0] = a[0] / length;
    b[1] = a[1] / length;
    b[2] = a[2] / length;
    return(length);
} /* end V3Direction() */


/* return 0..2 permutation idx[] such that length[idx[]] is ascending */
int *V3sort (double length[3], int idx[3])
{
    register int tmp;
    SEQ3(idx);
    if (length[idx[0]]>length[idx[1]]) SWAP(idx[0],idx[1],tmp);
    if (length[idx[1]]>length[idx[2]]) SWAP(idx[1],idx[2],tmp);
    if (length[idx[0]]>length[idx[1]]) SWAP(idx[0],idx[1],tmp);
    return (idx);
} /* end V3sort() */


/* return 0..2 random permutation index idx[] */
int *I3randompermutation (int idx[3])
{
    idx[0] = Fran(0,2);
    idx[1] = (idx[0]+Fran(1,2)) % 3;
    idx[2] = 3 - idx[0] - idx[1];
    return(idx);
} /* end I3randompermutation() */


/* return 0..2 permutation idx[] such that length[idx[]]  */
/* is non-decreasing. Equal length indices are randomized */
int *V3randomsort (double length[3], int idx[3])
{
    int tmp;
    I3RANDOMPERMUTATION (idx);
    if (length[idx[0]]>length[idx[1]]) SWAP(idx[0],idx[1],tmp);
    if (length[idx[1]]>length[idx[2]]) SWAP(idx[1],idx[2],tmp);
    if (length[idx[0]]>length[idx[1]]) SWAP(idx[0],idx[1],tmp);
    return (idx);
} /* end V3randomsort() */


/* return 0..2 permutation idx[] such that length[idx[]] is     */
/* approximately increasing within epsilon. idx[] is randomized */
int *V3Randomsort (double length[3], int idx[3], double epsilon)
{
    int tmp;
    I3RANDOMPERMUTATION (idx);
    if (length[idx[0]]>length[idx[1]]) SWAP(idx[0],idx[1],tmp);
    if (length[idx[1]]>length[idx[2]]) SWAP(idx[1],idx[2],tmp);
    if (length[idx[0]]>length[idx[1]]) SWAP(idx[0],idx[1],tmp);
    BEPOSITIVE(epsilon);
    if ( SMALLSEPARATION(length[idx[0]],length[idx[1]],epsilon)
         && HALF_DECISION() )  SWAP(idx[0],idx[1],tmp);
    if ( SMALLSEPARATION(length[idx[1]],length[idx[2]],epsilon)
         && HALF_DECISION() )  SWAP(idx[1],idx[2],tmp);
    if ( SMALLSEPARATION(length[idx[0]],length[idx[1]],epsilon)
         && HALF_DECISION() )  SWAP(idx[0],idx[1],tmp);
    return (idx);
} /* end V3Randomsort() */


/* Return spherical angles theta=[0,pi] and phi=[0,2*pi) of vector v[]    */
/* in the frame spanned by vectors z[] and x[]. x[] does not have to be   */
/* perpendicular to z[]; its perpendicular component to z[] will be used. */
void V3SphericalAngle (V3 v, V3 z, V3 x, double *theta, double *phi)
{
    V3 vv, zz, xx, yy;
    register double cc;
    if ( V3EQZERO(v) )
    {
        *theta = 0;
        *phi   = 0;
        return;
    }
    V3DIRECTION (v, vv, cc);
    if (V3EQZERO(z)) pe ("V3SphericalAngle: z is zero vector.\n");
    V3DIRECTION (z, zz, cc);
    cc = V3DOT(vv, zz);
    *theta = acos(cc);
    V3SUBmuL (vv, cc, zz);
    if (V3EQZERO(vv))
    {
        *phi = 0;
        return;
    }
    else V3NORMALIZE (vv,cc);
    if (V3EQZERO(x)) pe ("V3SphericalAngle: x is zero vector.\n");
    if (COLLINEAR(z,x,cc)) pe ("V3SphericalAngle: x and z are collinear.\n");
    V3CROSS (zz,x,yy);
    V3NORMALIZE (yy,cc);
    V3CROSS (yy,zz,xx);
    *phi = acos(V3DOT(vv,xx));
    if ( V3DOT(vv,yy) < 0 ) *phi = 2*PI-*phi;
    return;
} /* end V3SphericalAngle() */

#ifdef _V3SphericalAngle_TEST
#define TESTS 100
int main (int argc, char *argv[])
{
    int i;
    double theta=2.46018, phi=5.01675;
    V3 v, x={1,0,0}, z={0,0,1};
    for (i=TESTS; i--;)
    {
        theta = Frandom() * PI;
        phi = Frandom() * 2 * PI;
        printf ("%g %g ", theta, phi);
        V3SphericalUnit (theta, phi, v);
        V3SphericalAngle (v,z,x,&theta,&phi);
        printf ("%g %g\n", theta, phi);
    }
    return (0);
}
#undef TESTS
#endif /* _V3SphericalAngle_TEST */
