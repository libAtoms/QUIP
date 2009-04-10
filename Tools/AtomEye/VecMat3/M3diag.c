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


/*********************************************************/
/* Matrix inversion and symmetric matrix diagonalization */
/*********************************************************/


/* B[][] := A[][]^-1; return det(A) */
double M3inv (double A[3][3], double B[3][3])
{
    double determinant;
    B[0][0] = A[1][1]*A[2][2]-A[1][2]*A[2][1];
    B[1][1] = A[2][2]*A[0][0]-A[2][0]*A[0][2];
    B[2][2] = A[0][0]*A[1][1]-A[0][1]*A[1][0];
    B[1][0] = A[1][2]*A[2][0]-A[1][0]*A[2][2];
    B[2][1] = A[2][0]*A[0][1]-A[2][1]*A[0][0];
    B[0][2] = A[0][1]*A[1][2]-A[0][2]*A[1][1];
    B[2][0] = A[1][0]*A[2][1]-A[2][0]*A[1][1];
    B[0][1] = A[2][1]*A[0][2]-A[0][1]*A[2][2];
    B[1][2] = A[0][2]*A[1][0]-A[1][2]*A[0][0];
    determinant = A[0][0]*B[0][0] + A[0][1]*B[1][0] + A[0][2]*B[2][0];
    if (determinant == 0.)
    {
        printf ("error: M3inv: determinant = %e\n"
                "matrix is singular\n", determinant);
        exit(1);
    }
    B[0][0] /= determinant;
    B[1][1] /= determinant;
    B[2][2] /= determinant;
    B[1][0] /= determinant;
    B[2][1] /= determinant;
    B[0][2] /= determinant;
    B[2][0] /= determinant;
    B[0][1] /= determinant;
    B[1][2] /= determinant;
    return (determinant);
} /* end M3inv() */


#ifdef _M3inv_TEST
int main (int argc, char *argv[])
{
    double A[3][3], B[3][3], C[3][3];
    M3FRANDOM(A);
    S3PR("\nA = %M\n ", A);
    printf ("determinant = %f\n", M3inv(A,B));
    S3PR("\nB = %M\n ", B);
    M3mul(A,B,C);
    S3PR("A * B = %M\n ", C);
    return(0);
}
#endif  /* _M3inv_TEST */


/* A[][] := A[][]^-1; return original det(A) */
double M3Inv (double A[3][3])
{
    double determinant, D11, D22, D33, D12, D23, D31, D13, D21, D32;
    D11 = A[1][1]*A[2][2]-A[1][2]*A[2][1];
    D22 = A[2][2]*A[0][0]-A[2][0]*A[0][2];
    D33 = A[0][0]*A[1][1]-A[0][1]*A[1][0];
    D12 = A[1][2]*A[2][0]-A[1][0]*A[2][2];
    D23 = A[2][0]*A[0][1]-A[2][1]*A[0][0];
    D31 = A[0][1]*A[1][2]-A[0][2]*A[1][1];
    D13 = A[1][0]*A[2][1]-A[2][0]*A[1][1];
    D21 = A[2][1]*A[0][2]-A[0][1]*A[2][2];
    D32 = A[0][2]*A[1][0]-A[1][2]*A[0][0];
    determinant = A[0][0]*D11+A[0][1]*D12+A[0][2]*D13;
    A[0][0] = D11/determinant;
    A[1][1] = D22/determinant;
    A[2][2] = D33/determinant;
    A[0][1] = D21/determinant;
    A[1][2] = D32/determinant;
    A[2][0] = D13/determinant;
    A[1][0] = D12/determinant;
    A[2][1] = D23/determinant;
    A[0][2] = D31/determinant;
    return (determinant);
} /* end M3Inv() */


#ifdef _M3Inv_TEST
int main (int argc, char *argv[])
{
    double A[3][3], B[3][3], C[3][3];
    M3FRANDOM(A);
    S3PR("\nA = %M\n ", A);
    printf ("determinant = %f\n", M3Inv(A));
    S3PR("\nAfter M3Inv, A = %M\n ", A);
    return(0);
}
#endif  /* _M3Inv_TEST */


/***********************************************************************/
/* Find a linear combination of columns A[][0], A[][1], A[][2] so that */
/* A[][]*c[]=0. c[] is normalized and "random" if possible. return c[] */
/***********************************************************************/
#define M3diag_TINY (100*EPS)
#define M3diag_COLLINEAR(a,b) \
(fabs(V3DOT(a,b)/V3LENGTH(a)/V3LENGTH(b)-1) < M3diag_TINY)
double *M3nullvector (double A[3][3], double c[3])
{
    double scale;
    static int idx[3], tmp;
    static double B[3][3], rowlength[3];
    V3Frandom(c);
    M3INFNORM(A,scale);  /* largest component of A[][] in absolute value */
    if ( scale == 0. )
    { /* perfectly zero matrix */
        V3NORMALIZE(c,scale);
        return (c);
    }
    M3DIVIDE (A, scale, B);  /* set up the zoom: safe to do so */
    M3ROWLENGTHS (B, rowlength);
    /* V3Randomsort (rowlength, idx, TINY); */
    /* V3randomsort (rowlength, idx); */
    V3SORT (rowlength, idx, tmp);
    /* the longest row vector must be non-zero */
    if ( rowlength[idx[1]] < M3diag_TINY )  /* indeterministic */
        return(V3normalize(V3Perpendicular (c, B[idx[2]])));
    /* the second longest row vector is now non-zero */
    if (! M3diag_COLLINEAR(B[idx[2]], B[idx[1]]))  /* deterministic */
        return(V3normalize(V3cross(B[idx[2]],B[idx[1]],c)));
    /* now that the first two rows are indeed collinear */
    if ( rowlength[idx[0]] < M3diag_TINY )    /* indeterministic */
        return(V3normalize(V3Perpendicular (c, B[idx[2]])));
    /* the third row vector is now non-zero */
    if (! M3diag_COLLINEAR(B[idx[2]], B[idx[0]]))  /* deterministic */
        return(V3normalize(V3cross(B[idx[2]],B[idx[0]],c)));
    /* all three rows are collinear now */
    return(V3normalize(V3Perpendicular (c, B[idx[2]])));
} /* end M3nullvector() */


#define M3diag_ACC  (1000*EPS)

/******************************************************************/
/* Diagonalize 3x3 real symmetric matrix: A = V^T*Diag(eigval)*V, */
/* eigval[i] will be stored in ascending order, i = 0..2; the     */
/* corresponding eigenvectors V[i][] form a right-handed system.  */
/* Return index i of the eigenvalue of largest absolute value.    */
/******************************************************************/
int M3diag (double A[3][3], double eigval[3], double V[3][3])
{
    int i, j, idx[3];
    static double B[3][3], VV[3][3], C[3][3], VT[3][3], E[3][3], e[3];
    double scale, a, p, q, Q, R, theta, tmp;
    /* find the largest component of A[][] in absolute value */
    scale = M3infnorm(A);
    if ( scale == 0. )
    { /* this is a perfectly zero matrix */
        M3IDENTITY (V);
        V3ZERO (eigval);
        return (0);
    }
    /* set up the zoom: safe to do so */
    M3DIVIDE (A, scale, B);
    /* we really go after numerical accuracy */
    if ( ISSYMMETRIC(B) ) M3Symmetrize(B);
    else
    {
        printf ("error: M3diag: matrix is not symmetric\n");
        Mprintf ("\n%3M|| %15.8le %15.8le %15.8le |\n ", A[0]);
        exit (1);
    }
    /* eigenvalue equation: x^3 + a*x^2 + b*x + c = 0: */
    a = - B[0][0] - B[1][1] - B[2][2];
    /* b = B[0][0] * B[1][1] - B[0][1] * B[0][1]+ */
    /* B[0][0] * B[2][2] - B[0][2] * B[0][2]+ */
    /* B[1][1] * B[2][2] - B[1][2] * B[1][2]; */
    /* c = B[0][0] * B[1][2] * B[1][2] + B[1][1] * B[0][2] * B[0][2] + */
    /* B[2][2] * B[0][1] * B[0][1] - B[0][0] * B[1][1] * B[2][2] - */
    /* 2 * B[0][1] * B[0][2] * B[1][2]; */
    /* cubic equation standard form: x^3 + p * x = q */
    p = ( B[0][0] * (B[1][1] - B[0][0]) +
          B[2][2] * (B[0][0] - B[2][2]) +
          B[1][1] * (B[2][2] - B[1][1]) ) / 3
        - B[0][1] * B[0][1] - B[0][2] * B[0][2] - B[1][2] * B[1][2];
    q = ( -3 * B[0][0] * B[0][0] * B[1][1] +
          9 * B[0][0] * B[0][1] * B[0][1] -
          3 * B[0][0] * B[0][0] * B[2][2] +
          9 * B[0][0] * B[0][2] * B[0][2] -
          3 * B[0][0] * B[1][1] * B[1][1] +
          9 * B[1][1] * B[0][1] * B[0][1] -
          3 * B[1][1] * B[1][1] * B[2][2] -
          18 * B[0][0] * B[1][2] * B[1][2] -
          18 * B[1][1] * B[0][2] * B[0][2] -
          18 * B[2][2] * B[0][1] * B[0][1] +
          12 * B[0][0] * B[1][1] * B[2][2] +
          54 * B[0][1] * B[0][2] * B[1][2] +
          2 * B[0][0] * B[0][0] * B[0][0] +
          2 * B[1][1] * B[1][1] * B[1][1] +
          2 * B[2][2] * B[2][2] * B[2][2] +
          9 * B[1][1] * B[1][2] * B[1][2] -
          3 * B[0][0] * B[2][2] * B[2][2] +
          9 * B[2][2] * B[0][2] * B[0][2] -
          3 * B[1][1] * B[2][2] * B[2][2] +
          9 * B[2][2] * B[1][2] * B[1][2] ) / 27;
    /* Vieta's substitution: x = w - Q / w */
    Q = p / 3;
    R = q / 2;
    /* w's solution: w^3  = R \pm sqrt{R * R + Q * Q *Q} */
    /* D = R * R  +  Q * Q * Q; */
    /* discriminant */
    /**********************************************************/
    /* D > 0: one root is real and two are complex conjugates */
    /* D = 0: all roots are real and at least two are equal   */
    /* D < 0, all roots are real and unequal                  */
    /* cf. mathworld.wolfram.com/CubicEquation.html           */
    /**********************************************************/
    eigval[0] = eigval[1] = eigval[2] = -a/3;
    if (Q >= 0)
    { /* this can't happen in reality; must be numerical artifact */
        M3IDENTITY (VV);
        goto exit;
    } /* here on Q < 0 */
    tmp = R / sqrt(-Q) / (-Q);
    if (tmp >= 1.)  tmp = 1.;   /* numerical artifact */
    if (tmp <= -1.) tmp = -1.;  /* numerical artifact */
    theta = acos(tmp);
    eigval[0] += 2 * sqrt(-Q) * cos( (theta)      / 3 );
    eigval[1] += 2 * sqrt(-Q) * cos( (theta+2*PI) / 3 );
    eigval[2] += 2 * sqrt(-Q) * cos( (theta+4*PI) / 3 );
    ABS3 (eigval,e);
    V3SORT (e,idx,i);
    /* all above is just a preconditioner, next I do relaxation */
    for (i=0; i<100; i++)
    {
        j = 1;
        do
        { /* calculate the eigenvectors */
            if (j%3 == 0) I3RANDOMPERMUTATION(idx);
            M3SUBDIAG (B, eigval[idx[2]], C);
            M3nullvector (C, VV[2]);
            M3SUBDIAG (B, eigval[idx[1]], C);
            M3nullvector (C, VV[1]);
            V3CROSS (VV[1], VV[2], VV[0]);
            if ( (++j) == 100 )
            {
                printf ("error: M3diag: program error 1\n");
                Mprintf ("\n%3M|| %15.8le %15.8le %15.8le |\n ", A[0]);
                exit (1);
            }
        } while (V3length(VV[0]) < SMALL);
        V3NORMALIZE (VV[0], VT[0][0]);
        V3CROSS (VV[2], VV[0], VV[1]);
        V3NORMALIZE (VV[1], VT[0][0]);
        /* test our accuracy */
        M3TRANSPOSE (VV, VT);
        M3MUL2(VV, B, VT, E, C);
        eigval[0] = E[0][0];
        eigval[1] = E[1][1];
        eigval[2] = E[2][2];
        if (DISTANCE(E[0][1],E[0][2],E[1][2]) < M3diag_ACC) goto exit;
        I3RANDOMPERMUTATION(idx);
    }
    printf ("error: M3diag: program error 2\n");
    Mprintf ("\n%3M|| %15.8le %15.8le %15.8le |\n ", A[0]);
    Mprintf ("%3M|| %15.8le %15.8le %15.8le |\n ", E[0]);
    exit (1);
  exit:  /* zoom back */
    V3MUL (scale,eigval,e);
    /* make eigval[0] <= eigval[1] <= eigval[2] */
    V3SORT (e, idx, i);
    eigval[0] = e[idx[0]];
    V3EQV (VV[idx[0]], V[0]); 
    eigval[1] = e[idx[1]];
    V3EQV (VV[idx[1]], V[1]);
    eigval[2] = e[idx[2]];
    V3EQV (VV[idx[2]], V[2]);
    ABS3 (eigval,e);
    V3SORT (e,idx,i);
    return (idx[2]);
} /* end M3diag() */

#ifdef _M3diag_TEST
int main (int argc, char *argv[])
{
    double A[3][3],eigval[3],Q[3][3],QT[3][3],L[3][3],QI[3][3];
    M3FRANDOM (A);
    M3Symmetrize (A);
    M3IDENTITY (A);
    A[0][1] += TINY;
    A[1][0] += TINY;
    A[2][2] -= TINY/2;
    /* A[2][2] = 0.00001; */
    /* A[1][1] = 1.1; */
    /* A[0][0] = 0.3; */
    M3diag (A,eigval,Q);
    S3PR ("\nA = %M\n ", A);
    V3pr ("eigval[] = %M\n ", eigval);
    printf ("%25.18le %25.18le %25.18le\n",eigval[0], eigval[1], eigval[2]); 
    S3PR ("\nQ = %M\n ", Q);
    M3TRANSPOSE (Q,QT);
    M3MUL2 (Q,A,QT,L,QI);
    S3PR ("Q * A * Q^T = %M\n ", L);
    M3MUL (Q,QT,L);
    S3PR("Q * Q^T = %M\n ", L);
    printf ("det|Q| = %f\n\n", M3DETERMINANT(Q));
    return(0);
}
#endif  /* _M3diag_TEST */


/* same function as M3diag() but using Lapack */
#define M3DIAG_LWORK 216
int M3Diag (double A[3][3], double eigval[3], double V[3][3])
{
    int info, N=3, lwork=M3DIAG_LWORK;
    static int idx[3];
    static double work[M3DIAG_LWORK];
    if ( ! ISSYMMETRIC(A) )
    {
        printf ("error: M3Diag: matrix is not symmetric\n");
        Mprintf ("\n%3M|| %15.8le %15.8le %15.8le |\n ", A[0]);
        exit (1);
    }
    M3EQV (A, V);
    FORTRAN_SYMBOL(dsyev)
        ("V", "U", &N, V[0], &N, eigval, work, &lwork, &info);
    if ( info > 0 )
    {
        printf ("error: M3Diag: matrix failed to converge.\n");
        Mprintf ("\n%3M|| %15.8le %15.8le %15.8le |\n ", A[0]);
        exit (1);
    }
    ABS3 (eigval,work);
    V3SORT (work,idx,info);
    return (idx[2]);
} /* end M3Diag() */

#ifdef _M3Diag_TEST
int main (int argc, char *argv[])
{
    double A[3][3],eigval[3],Q[3][3],QT[3][3],L[3][3],QI[3][3];
    M3FRANDOM (A);
    M3Symmetrize (A);
    /* M3IDENTITY (A); */
    /* A[2][2] = 0.00001; */
    /* A[1][1] = 1.1; */
    /* A[0][0] = 0.3; */
    M3Diag (A,eigval,Q);
    S3PR ("\nA = %M\n ", A);
    V3pr ("eigval[] = %M\n ", eigval);
    S3PR ("Q = %M\n ", Q);
    M3TRANSPOSE (Q,QT);
    M3MUL2 (Q,A,QT,L,QI);
    S3PR ("Q * A * Q^T = %M\n ", L);
    M3MUL (Q,QT,L);
    S3PR("Q * Q^T = %M\n ", L);
    printf ("det|Q| = %f\n\n", M3inv(Q,QI));
    return(0);
}
#endif  /* _M3Diag_TEST */


#ifdef _M3diagContest_TEST
#define M3diagContest_TRIALS  100000
#include <Timer.h>
int main(int argc, char *argv[])
{
    int i;
    double A[3][3],eigval[3],Q[3][3],QT[3][3],L[3][3],TMP[3][3];
    start_chronometer();
    for (i=0; i<M3diagContest_TRIALS; i++)
    { /* control group sugar pill */
        M3FRANDOM (A);
        M3Symmetrize (A);
        M3EQV(A,Q);
        M3TRANSPOSE (Q,QT);
        M3MUL2 (Q,A,QT,L,TMP);
    }
    stop_chronometer();
    printf ("empty run: %d 3x3 matrix manipulations took %s\n",
            M3diagContest_TRIALS, earthtime(stopped_usertime()));
    start_chronometer();
    for (i=0; i<M3diagContest_TRIALS; i++)
    {
        M3FRANDOM (A);
        M3Symmetrize (A);
        M3diag (A,eigval,Q);
        M3TRANSPOSE (Q,QT);
        M3MUL2 (Q,A,QT,L,TMP);
        if (DISTANCE(L[0][1],L[0][2],L[1][2]) > M3diag_ACC)
        { /* highest level quality control */
            Mprintf ("\n%3M|| %15.8le %15.8le %15.8le |\n ", A[0]);
            Mprintf ("%3M|| %15.8le %15.8le %15.8le |\n ", L[0]);
            exit(1);
        }
    }
    stop_chronometer();
    printf ("Ju Li's M3diag: %d 3x3 diagonalizations took %s\n",
            M3diagContest_TRIALS, earthtime(stopped_usertime()));
    start_chronometer();
    for (i=0; i<M3diagContest_TRIALS; i++)
    {
        M3FRANDOM (A);
        M3Symmetrize (A);
        M3Diag (A,eigval,Q);
        M3TRANSPOSE (Q,QT);
        M3MUL2 (Q,A,QT,L,TMP);
        if (DISTANCE(L[0][1],L[0][2],L[1][2]) > M3diag_ACC)
        { /* highest level quality control */
            Mprintf ("\n%3M|| %15.8le %15.8le %15.8le |\n ", A[0]);
            Mprintf ("%3M|| %15.8le %15.8le %15.8le |\n ", L[0]);
            exit(1);
        }
    }
    stop_chronometer();
    printf ("LAPACK M3Diag: %d 3x3 diagonalizations took %s\n",
            M3diagContest_TRIALS, earthtime(stopped_usertime()));
    return(0);
}
#endif  /* _M3diagContest_TEST */


/* A = M*R, where M is symmetric, R (do not return) is orthogonal */
void M3RightPolarDecompose (M3 A, M3 M)
{
    V3 eigval;
    M3 A2,V,E,VT,TMP;
    if (M3EQZERO(A))
    {
        M3ZERO(M);
        return;
    }
    M3MULT(A,A2);
    M3Diag(A2,eigval,V);
    M3diagonal(sqrt(eigval[0]),sqrt(eigval[1]),sqrt(eigval[2]),E);
    M3TRANSPOSE(V,VT);
    M3MUL2 (VT,E,V,M,TMP);
    return;
} /* end M3RightPolarDecompose() */


/* A = M*R, where M is symmetric, R is orthogonal */
void M3RightPolarDECOMPOSE (M3 A, M3 M, M3 R)
{
    M3 MI;
    M3RightPolarDecompose (A, M);
    if (M3DETERMINANT(M) == 0.)
        pe("M3RightPolarDECOMPOSE: A is singular.\n");
    M3inv (M,MI);
    M3MUL (MI,A,R);
    return;
} /* end M3RightPolarDECOMPOSE() */

#ifdef _M3RightPolarDECOMPOSE_TEST
int main (int argc, char *argv[])
{
    M3 A = GuineaPigM3, M,R,RT,RRT, RESIDUE;
    M3RightPolarDECOMPOSE (A, M, R);
    M3MUL (M, R, RESIDUE);
    M3SuB (RESIDUE, A);
    S3PR("residue = %M", RESIDUE);
    S3PR("R = %M", R);
    M3TRANSPOSE (R,RT);
    M3MUL (R,RT,RRT);
    S3PR("R * R^T = %M", RRT);
    return (0);
}
#endif /* _M3RightPolarDECOMPOSE_TEST */


/* A = L*M, where L (do not return) is orthogonal, M is symmetric */
void M3LeftPolarDecompose (M3 A, M3 M)
{
    V3 eigval;
    M3 A2,V,E,VT,TMP;
    if (M3EQZERO(A))
    {
        M3ZERO(M);
        return;
    }
    M3TMUL(A,A2);
    M3Diag(A2,eigval,V);
    M3diagonal(sqrt(eigval[0]),sqrt(eigval[1]),sqrt(eigval[2]),E);
    M3TRANSPOSE(V,VT);
    M3MUL2 (VT,E,V,M,TMP);
    return;
} /* end M3LeftPolarDecompose() */


/* A = L*M, where L is orthogonal, M is symmetric */
void M3LeftPolarDECOMPOSE (M3 A, M3 L, M3 M)
{
    M3 MI;
    M3LeftPolarDecompose (A, M);
    if (M3DETERMINANT(M) == 0.)
        pe("M3LeftPolarDECOMPOSE: A is singular.\n");
    M3inv (M,MI);
    M3MUL (A,MI,L);
    return;
} /* end M3LeftPolarDECOMPOSE() */

#ifdef _M3LeftPolarDECOMPOSE_TEST
int main (int argc, char *argv[])
{
    M3 A = GuineaPigM3, L,M,LT,LLT, RESIDUE;
    M3LeftPolarDECOMPOSE (A, L, M);
    M3MUL (L, M, RESIDUE);
    M3SuB (RESIDUE, A);
    S3PR("residue = %M", RESIDUE);
    S3PR("L = %M", L);
    M3TRANSPOSE (L,LT);
    M3MUL (L,LT,LLT);
    S3PR("L * L^T = %M", LLT);
    return (0);
}
#endif /* _M3LeftPolarDECOMPOSE_TEST */
