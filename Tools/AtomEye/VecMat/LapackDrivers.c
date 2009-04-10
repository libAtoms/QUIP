/*********************************************/
/* libVecMat: -llapck -lblas -lm             */
/*            -lScalar -lIO                  */
/*                                           */
/* General Vector & Matrix Operation Library */
/*                                           */
/* Nov 19 1999 Ju Li <liju99@mit.edu>        */
/*********************************************/

#include "VecMat.h"


/***************************************************************/
/* BLAS and LAPACK library drivers with simplified interface   */
/* and temporary memory hierarchy (static..dynamic) management */
/* invisible to the end-user.                                  */
/***************************************************************/

/* B := inv(A), both A and B are real matrices of rank x rank */
void Minv (int rank, double A[], double B[])
{
    int info, *ipiv;
    double *a, *work;
    a = Mmem(rank*rank);
    work = Vmem(rank);
    ipiv = VImem(rank);
    memcpy (a, A, (long)rank*rank*sizeof(double));
    FORTRAN_SYMBOL(dgetrf) (&rank, &rank, a, &rank, ipiv, &info);
    if ( info != 0 )
    {
	printf ("error: Minv: matrix of rank %d is singular.\n",
		rank);
	mump(rank,A);
	exit(1);
    }
    FORTRAN_SYMBOL(dgetri) (&rank, a, &rank, ipiv, work, &rank, &info);
    memcpy (B, a, (long)rank*rank*sizeof(double));
    VFREE (ipiv);
    Vfree (work);
    Mfree (a);
    return;
} /* end Minv() */

#ifdef _Minv_TEST
#define DIM 3
int main ()
{
    double A[DIM][DIM], B[DIM][DIM];
    Vfrandom(DIM*DIM, A[0]);
    Minv (DIM, A[0], B[0]);
    S3pr("\nA = %M\n ", A[0]);
    S3pr("B = %M\n ", B[0]);
    return(0);
} /* end main() */
#endif


/*******************************************************************/
/* Diagonalize a complex Hermitian matrix described by ap[], order */
/* its (real) eigenvalues in ascending order, and compute all      */
/* (complex) eigenvectors associated with them.                    */
/*                                                                 */
/* ap=double[rank*(rank+1)]: upper triangle of A packed columnwise */
/* ap[j*(j+1)+2*i] = Re(A_{ij}), ap[j*(j+1)+2*i+1] = Img(A_{ij})   */
/* i, j = 0 .. rank-1, i <= j.                                     */
/*                                                                 */
/* eigenvalues = double[rank], eigenvectors = double[2*rank*rank]. */
/*******************************************************************/

void MDiag (int rank, double ap[], double eigenvalues[], double eigenvectors[])
{
    int info;
    double *work, *rwork;
    work  = Vmem(4*rank);
    rwork = Vmem(3*rank);
    info = 0;
    FORTRAN_SYMBOL(zhpev) ("V", "U", &rank, ap, eigenvalues, eigenvectors,
                           &rank, work, rwork, &info);
    if ( info > 0 )
    {
	printf ("error: MDiag: rank %d matrix failed to converge.\n",
		rank);
	exit(1);
    }
    Vfree (rwork);
    Vfree (work);
} /* end MDiag() */


