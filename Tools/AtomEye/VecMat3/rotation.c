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


/****************************************/
/* rotation generator & representations */
/****************************************/


/* Rotate v[] around axis[] by theta: store and return in w[] */
double *V3axialrotate (double axis[3], double theta, double v[3], double w[3])
{
    double axislength, cosine, sine, vertical[3];
    if ( V3EQZERO(axis) ) pe ("V3axialrotate: axis=0\n");
    axislength = V3LENGTH(axis);
    V3CROSS (axis, v, vertical);
    V3DiV (vertical, axislength);
    cosine = cos(theta);
    sine   = sin(theta);
    V3ADDMULMUL (cosine,v, sine,vertical, w);
    return (w);
} /* end V3axialrotate() */


/* Compute rotational matrix R[][] corresponding */
/* to rotation with respect to axis[] by theta.  */
void V3axialrotatematrix (double axis[3], double theta, double R[3][3])
{
    M3 E = IdentityM3;
    if ( V3EQZERO(axis) ) pe ("V3axialrotatematrix: axis=0\n");    
    V3axialrotate (axis, theta, &E[0][0], &R[0][0]);
    V3axialrotate (axis, theta, &E[1][0], &R[1][0]);
    V3axialrotate (axis, theta, &E[2][0], &R[2][0]);
    return;
} /* end V3axialrotatematrix() */


/* If geodesic ("big-circle") rotation makes a[]->b[], what happens to v[]? */
/* Assume a[] and b[] are already NORMALIZED and NOT collinear; returns v[] */
double *V3geodesic (double a[3], double b[3], double v[3])
{
    double c[3], ca[3], cb[3], vc, va, vca;
    V3CROSS(a,b,c);
    V3NORMALIZE(c,cb[0]);
    V3CROSS(c,a,ca);
    V3CROSS(c,b,cb);
    vc  = V3DOT(v,c);
    va  = V3DOT(v,a);
    vca = V3DOT(v,ca);
    /* recombine in the new coordinate frame */
    v[0] = vc*c[0] + va*b[0] + vca*cb[0];
    v[1] = vc*c[1] + va*b[1] + vca*cb[1];
    v[2] = vc*c[2] + va*b[2] + vca*cb[2];
    return(v);
} /* end V3geodesic() */


/* Compute the rotational matrix R[][] corresponding to a geodesic */
/* ("big-circle") rotation that makes a[]->b[] as v[] := v[]*R[][] */
/* Assume a[] and b[] are already NORMALIZED and NOT collinear.    */
void M3geodesic (double a[3], double b[3], double R[3][3])
{
    double c[3], ca[3], cb[3];
    /* M3identity (1.,R); */
    V3CROSS(a,b,c);
    V3NORMALIZE(c,cb[0]);
    V3CROSS(c,a,ca);
    V3CROSS(c,b,cb);
    /* V3geodesic(a,b,&R[0][0]); */
    R[0][0] = c[0]*c[0] + a[0]*b[0] + ca[0]*cb[0];
    R[0][1] = c[0]*c[1] + a[0]*b[1] + ca[0]*cb[1];
    R[0][2] = c[0]*c[2] + a[0]*b[2] + ca[0]*cb[2];
    /* V3geodesic(a,b,&R[1][0]); */
    R[1][0] = c[1]*c[0] + a[1]*b[0] + ca[1]*cb[0];
    R[1][1] = c[1]*c[1] + a[1]*b[1] + ca[1]*cb[1];
    R[1][2] = c[1]*c[2] + a[1]*b[2] + ca[1]*cb[2];
    /* V3geodesic(a,b,&R[2][0]); */
    R[2][0] = c[2]*c[0] + a[2]*b[0] + ca[2]*cb[0];
    R[2][1] = c[2]*c[1] + a[2]*b[1] + ca[2]*cb[1];
    R[2][2] = c[2]*c[2] + a[2]*b[2] + ca[2]*cb[2];
    return;
} /* end M3geodesic() */
