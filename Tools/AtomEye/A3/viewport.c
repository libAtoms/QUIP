/*******************************************/
/* Threaded Atomistic Configuration Viewer */
/*                                         */
/* Fri Mar 31 2000 Ju Li <liju99@mit.edu>  */
/*******************************************/

#include "A.h"

/* Change the distance between viewpoint and anchor: if  */
/* delta >=0, d -> d/(1+delta), else d -> d * (1-delta). */
bool advance (int iw, double delta)
{
    double tmp[3];
    if (delta >= 0) delta = 1/(1+delta);
    else delta = (1-delta);
    if (n[iw].anchor >= 0)
    {
        V3SUB (AX_3D[iw].x, B->BALL[n[iw].anchor].x, tmp);
        V3ADDMUL (B->BALL[n[iw].anchor].x, delta, tmp, AX_3D[iw].x);
    }
    else
    {
        V3SUB (AX_3D[iw].x, n[iw].hook, tmp);
        V3ADDMUL (n[iw].hook, delta, tmp, AX_3D[iw].x);
    }
    return (TRUE);
} /* end advance() */


bool pointer_advance (int iw, int to_x, int to_y)
{
    double delta, tmp[3];
    if (n[iw].ly == to_y) return (FALSE);
    delta = exp(2*(n[iw].ly-to_y)/n[iw].mgs_radius);
    if (n[iw].anchor >= 0)
    {
        V3SUB (AX_3D[iw].x, B->BALL[n[iw].anchor].x, tmp);
        V3ADDMUL (B->BALL[n[iw].anchor].x, delta, tmp, AX_3D[iw].x);
    }
    else
    {
        V3SUB (AX_3D[iw].x, n[iw].hook, tmp);
        V3ADDMUL (n[iw].hook, delta, tmp, AX_3D[iw].x);
    }
    n[iw].lx = to_x;
    n[iw].ly = to_y;
    return (TRUE);
} /* end pointer_advance() */


/* Translate the viewport along axis-i by amount d*lengthscale */
bool translate (int iw, int i, double d)
{
    d *= lengthscale;
    V3ADDmuL (-d, AX_3D[iw].V[i], AX_3D[iw].x);
    return (TRUE);
} /* end translate() */


/* use pointer device to translate the viewport */
bool pointer_translate (int iw, int to_x, int to_y)
{
    double z, tmp[3];
    if (n[iw].anchor >= 0) V3SUB (B->BALL[n[iw].anchor].x, AX_3D[iw].x, tmp);
    else V3SUB (n[iw].hook, AX_3D[iw].x, tmp);
    z = V3DOT (AX_3D[iw].V[2], tmp);
    tmp[0] = (n[iw].lx - to_x) / AX_3D[iw].k * z;
    V3ADDmuL (tmp[0], AX_3D[iw].V[0], AX_3D[iw].x);
    tmp[1] = (n[iw].ly - to_y) / AX_3D[iw].k * z;
    V3ADDmuL (tmp[1], AX_3D[iw].V[1], AX_3D[iw].x);
    n[iw].lx = to_x;
    n[iw].ly = to_y;
    return (TRUE);
} /* end pointer_translate() */


/* Rotate viewport as V' = R*V but letting the anchor stay at the same place */
void rotate (int iw, double R[3][3])
{
    double V0[3][3], s[3], tmp[3];
    if (n[iw].anchor >= 0)
    {
        V3SUB (B->BALL[n[iw].anchor].x, AX_3D[iw].x, tmp);
        M3mV3 (AX_3D[iw].V, tmp, s);
    }
    else
    {
        V3SUB (n[iw].hook, AX_3D[iw].x, tmp);
        M3mV3 (AX_3D[iw].V, tmp, s);
    }
    M3EQV (AX_3D[iw].V, V0);
    V3mM3 (R[0], V0, AX_3D[iw].V[0]);
    V3mM3 (R[1], V0, AX_3D[iw].V[1]);
    V3NORMALIZE (AX_3D[iw].V[1], tmp[0]);
    /* ensure the orthogonality numerically */
    V3CROSS (AX_3D[iw].V[0], AX_3D[iw].V[1], AX_3D[iw].V[2]);
    V3NORMALIZE (AX_3D[iw].V[2], tmp[0]);
    V3CROSS (AX_3D[iw].V[1], AX_3D[iw].V[2], AX_3D[iw].V[0]);
    if (n[iw].anchor >= 0)
    {
        V3mM3 (s, AX_3D[iw].V, tmp);
        V3SUB (B->BALL[n[iw].anchor].x, tmp, AX_3D[iw].x);
    }
    else
    {
        V3mM3 (s, AX_3D[iw].V, tmp);
        V3SUB (n[iw].hook, tmp, AX_3D[iw].x);
    }
    return;
} /* end rotate() */


/* magic sphere: map a plane onto an upper half-sphere surface */
void mgs (double x, double y, double r, double a[3])
{
    double a2;
    if ( (a2=x*x+y*y) < r*r )
    { 
        a[0] = x / r;
        a[1] = y / r;
        a[2] = -sqrt(1.-a[0]*a[0]-a[1]*a[1]);
    }
    else
    {
        a[0] = x / sqrt(a2);
        a[1] = y / sqrt(a2);
        a[2] = 0.;
    }
    return;
} /* end mgs() */


#define MIN_RADIAN (2*PI/720)
/* use pointer device to rotate (via magic sphere) */
bool pointer_rotate (int iw, int to_x, int to_y)
{
    double a[3], b[3], c[3], R[3][3];
    mgs ( n[iw].lx - AX_size[iw].width/2.,
          n[iw].ly - AX_size[iw].height/2., n[iw].mgs_radius, a );
    mgs ( to_x - AX_size[iw].width/2.,
          to_y - AX_size[iw].height/2., n[iw].mgs_radius, b );
    V3CROSS (a, b, c);
    if (V3LENGTH2(c) < MIN_RADIAN*MIN_RADIAN) return (FALSE);
    M3geodesic (b, a, R);
    rotate (iw, R);
    n[iw].lx = to_x;
    n[iw].ly = to_y;
    return (TRUE);
} /* end pointer_rotate() */


/* rotate around axis "i" by "theta" radians */
bool axis_rotate (int iw, int i, double theta)
{
    double R[3][3];
    int j=(i+1)%3, k=(i+2)%3;
    R[i][i]=1; R[i][j]=R[i][k]=0;
    R[j][i]=0; R[j][j]=cos(theta);  R[j][k]=-sin(theta);
    R[k][i]=0; R[k][j]=sin(theta); R[k][k]=cos(theta);
    rotate (iw, R);
    return (TRUE);
} /* end axis_rotate() */
