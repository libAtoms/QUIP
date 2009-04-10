/***************************************/
/* libMin: Minimization Library        */
/*         -lm -lScalar -lVecMat       */
/*                                     */
/* Dec. 6, 1999 Ju Li <liju99@mit.edu> */
/***************************************/

#include "Min.h"

#ifdef _mina_TEST
#define DIM 20
double L[DIM], U[DIM], xguess[DIM], xmin[DIM];
double fe (double x[])
{
    static int count = 0;
    int i;
    double error;
    for (error=i=0; i<DIM; i++) error += SQUARE(x[i]-i);
    ++count;
    if (count%500 == 1)
	printf ("Count = %d  error = %e\n", count, error);
    return (1.1+error);
}

int main()
{
    int i, n=DIM, ndiv=20;
    double del=0.5, emin;
    for (i=0; i<DIM; i++)
    {
	L[i] = i - Frandom();
	U[i] = i + Frandom();
	xguess[i] = (L[i] + U[i])/2;
    }
    emin = FORTRAN_SYMBOL(minaf77) (fe, &n, &ndiv, &del, L, U, xguess, xmin);
    printf ("\nMinimization complete (emin = %f):\n\n", emin);
    for (i=0; i<DIM; i++) printf ("%d %f\n", i, xmin[i]);
    return(1);
}
#endif


#ifdef _simann_TEST
#define DIM 20
double L[DIM], U[DIM], xguess[DIM], xmin[DIM];
double fe (double x[])
{
    static int count = 0;
    int i;
    double error;
    for (error=i=0; i<DIM; i++) error += SQUARE(x[i]-i);
    ++count;
    if (count%5000 == 1)
	printf ("Count = %d  error = %e\n", count, error);
    return (1.1+error);
}

int main()
{
    int i, n=DIM, ISEED=20;
    double T, RT, emin;
    for (i=0; i<DIM; i++)
    {
	L[i] = i - Frandom();
	U[i] = i + Frandom();
	xguess[i] = (L[i] + U[i])/2;
    }
    T = 1.0;
    RT = 0.99;
    ISEED = 1;
    emin = FORTRAN_SYMBOL(simannf77)
        (fe, &n, &T, &RT, &ISEED,L, U, xguess, xmin);
    printf ("\nMinimization complete (emin = %f):\n\n", emin);
    for (i=0; i<DIM; i++) printf ("%d %f\n", i, xmin[i]);
    return(1);
}
#endif
