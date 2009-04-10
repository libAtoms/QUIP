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

/****************************/
/* Properties of H matrices */
/****************************/

#ifdef _Crystallographic_to_H_TEST
int main (int argc, char *argv[])
{
    double H[3][3];
    Crystallographic X;
    CrystallographicAssign(1.2,1.3,1.4,60,70,20,X);
    Crystallographic_to_H(X,H);
    S3PR("H = %M\n", H);
    H_to_Crystallographic(H,X);
    printf ("%f %f %f %f %f %f\n", X.a, X.b, X.c, X.alpha, X.beta, X.gamma);
    return (0);
}
#endif /* _Crystallographic_to_H_TEST */


/* determine which of the seven vertices of a col-parallelepiped */
/* formed by A[][0], A[][1], A[][2] is the farthest from origin. */
double M3maxcolumnradius (double H[3][3])
{
    int i;
    double dx[3], r2, r2max;
    double ds[7][3] = {{1,0,0},{0,1,0},{0,0,1},{1,1,0},
                       {0,1,1},{1,0,1},{1,1,1}};
    for (r2max=i=0; i<7; i++)
    {
        M3mV3 (H, ds[i], dx);
        r2 = V3DOT (dx,dx);
        if (r2 > r2max) r2max = r2;
    }
    return (sqrt(r2max));
} /* end M3maxcolumnradius() */


/* determine which of the seven vertices of a row-parallelepiped */
/* formed by A[0][], A[1][], A[2][] is the farthest from origin. */
double M3maxrowradius (double H[3][3])
{
    int i;
    double dx[3], r2, r2max;
    double ds[7][3] = {{1,0,0},{0,1,0},{0,0,1},{1,1,0},
                       {0,1,1},{1,0,1},{1,1,1}};
    for (r2max=i=0; i<7; i++)
    {
        V3mM3 (ds[i], H, dx);
        r2 = V3DOT (dx,dx);
        if (r2 > r2max) r2max = r2;
    }
    return (sqrt(r2max));
} /* end M3maxrowradius() */


/* returns the thickness (>=0) of the parallelepiped formed */
/* by (H[0][], H[1][], H[2][]) in the row i direction.      */
double M3rowthickness (double H[3][3], int i)
{
    static double normal[3];
    if (i==0)
        V3Normalize(V3cross(&H[1][0],&H[2][0],normal));
    else if (i==1)
        V3Normalize(V3cross(&H[2][0],&H[0][0],normal));
    else if (i==2)
        V3Normalize(V3cross(&H[0][0],&H[1][0],normal));
    else
    {
        printf ("error: M3rowthickness: i = %d illegal\n", i);
        exit(1);
    }
    return (fabs(V3DOT(normal,&H[i][0])));
} /* end M3rowthickness() */


/* returns the three thicknesses (>=0) of the parallelepiped */
/* formed by (H[0][], H[1][], H[2][]), in thickness[].       */
double *M3rowthicknesses (double H[3][3], double thickness[3])
{
    static double normal[3];
    V3Normalize(V3cross(&H[1][0],&H[2][0],normal));
    thickness[0] = fabs(V3DOT(normal,&H[0][0]));
    V3Normalize(V3cross(&H[2][0],&H[0][0],normal));
    thickness[1] = fabs(V3DOT(normal,&H[1][0]));
    V3Normalize(V3cross(&H[0][0],&H[1][0],normal));
    thickness[2] = fabs(V3DOT(normal,&H[2][0]));
    return (thickness);
} /* end M3rowthicknesses() */


/* Calculate multiplication factors nc[] whereby -nc[0]:nc[0],   */
/* -nc[1]:nc[1], -nc[2]:nc[2] replica of H[0..2][] is guaranteed */
/* to cover sphere of radius R; returns total number of replicas */
int M3rows_to_cover_sphere (double H[3][3], double R, int nc[3])
{
    int i;
    static double thickness[3];
    M3rowthicknesses (H, thickness);
    for (i=0; i<3; i++)
        if (thickness[i] == 0.)
        {
            printf ("error: M3rows_to_cover_sphere: "
                    "thickness[%d] = 0.\n", i);
            exit(1);
        }
        else nc[i] =(int)ceil(ABS(R)/thickness[i]);
    return ((2*nc[0]+1)*(2*nc[1]+1)*(2*nc[2]+1));
} /* end M3rows_to_cover_sphere() */


/* return the thickness (>=0) of the parallelepiped formed */
/* by (H[][0], H[][1], H[][2]) in the column i direction.  */
double M3columnthickness (double H[3][3], int i)
{
    static double HT[3][3];
    M3TRANSPOSE (H, HT);
    return (M3rowthickness(HT,i));
} /* end M3columnthickness() */


/* return the three thicknesses (>=0) of the parallelepiped */
/* formed by (H[][0], H[][1], H[][2]), in thickness[].      */
double *M3columnthicknesses (double H[3][3], double thickness[3])
{
    static double HT[3][3];
    M3TRANSPOSE (H, HT);
    return (M3rowthicknesses(HT,thickness));
} /* end M3columnthicknesses() */


/* Calculate multiplication factors nc[] whereby -nc[0]:nc[0],   */
/* -nc[1]:nc[1], -nc[2]:nc[2] replica of H[][0..2] is guaranteed */
/* to cover sphere of radius R; return total number of replicas. */
int M3columns_to_cover_sphere (double H[3][3], double R, int nc[3])
{
    static double HT[3][3];
    M3TRANSPOSE (H, HT);
    return (M3rows_to_cover_sphere(HT,R,nc));
} /* end M3columns_to_cover_sphere() */


/* dl^2 = dx_i * (1 + 2 * eta_{ij}) * dx_j, where dx meshes the undeformed */
/* because dx = ds * H0, dx' = ds * H, dl^2 = dx' * (dx')^T => above       */
void Lagrangian_strain (double H0[3][3], double H[3][3], double eta[3][3])
{
    double determinant;
    M3 H0I,J;
    M3INV (H0,H0I,determinant);
    if (determinant == 0)
        pe ("Lagrangian_strain: H0 is singular.\n");
    M3MUL (H0I,H,J);
    M3MULT (J,eta);
    M3SubdiaG (eta, 1);
    M3DividE (eta, 2);
    return;
} /* end Lagrangian_strain() */


/* achieve eta without rotation. M := sqrt(1 + 2*eta), H := H0 * M */
void pure_deform (double H0[3][3], double eta[3][3], double H[3][3])
{
    int i;
    static double eigval[3],L[3][3]={{0}},Q[3][3],QT[3][3],M[3][3],tmp[3][3];
    M3Diag (eta, eigval, Q);
    for (i=0; i<3; i++)
    {
        eigval[i] = 1. + 2 * eigval[i];
        if (eigval[i] < 0)
        {
            printf ("error: pure_deform: Lagrangian strain infeasible\n");
            S3pr ("\neta = %M\n ", eta[0]);
            exit(1);
        }
        L[i][i] = sqrt(eigval[i]);
    }
    M3TRANSPOSE (Q, QT);
    M3MUL2 (QT, L, Q, M, tmp);
    if (H[0] == H0[0])
    { /* H and H0 are the same matrix */
        M3MUL (H0, M, tmp);
        M3EQV (tmp, H);
    }
    else M3MUL (H0, M, H);
    return;
} /* end pure_deform() */


#ifdef _pure_deform_TEST
int main (int argc, char *argv[])
{
    double H0[3][3], H[3][3], eta[3][3];
    M3FRANDOM(H0);
    S3PR("H0 = %M\n", H0);
    M3FRANDOM(eta);
    M3Symmetrize(eta);
    M3MultiplY(0.2,eta);
    S3PR("eta = %M\n", eta);
    pure_deform(H0,eta,H);
    S3PR("H = %M\n", H);
    Lagrangian_strain(H0,H,eta);
    S3PR("eta = %M\n", eta);
    return (0);
}
#endif /* _pure_deform_TEST */


/* D = H0^-1*H - 1 */
void simple_D (double H0[3][3], double H[3][3], double D[3][3])
{
    double determinant;
    M3 H0I;
    M3INV (H0,H0I,determinant);
    if (determinant == 0)
        pe ("simple_D: H0 is singular.\n");
    M3MUL (H0I,H,D);
    M3SubdiaG (D, 1);
    return;
} /* end simple_D() */


/* H = H0 * (1+D) */
void simple_deform (double H0[3][3], double D[3][3], double H[3][3])
{
    M3 J;
    M3ADDDIAG(D,1,J);
    M3MUL (H0,J,H);
    return;
} /* end simple_deform() */


/* simple strain = (J+J^T)/2-1 = (D+D^T)/2 */
void simple_strain (double H0[3][3], double H[3][3], double eta[3][3])
{
    double determinant;
    M3 H0I;
    M3INV (H0,H0I,determinant);
    if (determinant == 0)
        pe ("simple_strain: H0 is singular.\n");
    M3MUL (H0I,H,eta);
    M3Symmetrize(eta);
    M3SubdiaG (eta, 1);
    return;
} /* end simple_strain() */


#ifdef _simpledeform_TEST
int main (int argc, char *argv[])
{
    M3 H0 = GuineaPigM3, D, H, eta;
    M3FRANDOM(D);
    M3MultiplY(0.01,D);
    S3PR("D = %M\n ", D);
    simple_deform (H0, D, H);
    simple_D (H0, H, D);
    S3PR("D = %M\n ", D);
    simple_strain (H0, H, eta);
    S3PR("simple_strain = %M\n ", eta);
    Lagrangian_strain (H0, H, eta);
    S3PR("Lagrangian_strain = %M\n", eta);
    return (0);
}
#endif /* _simpledeform_TEST */
