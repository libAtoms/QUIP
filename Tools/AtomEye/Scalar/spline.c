/*******************************************/
/* libScalar:                              */
/*                                         */
/* Scalar Arithmetic Operations and Macros */
/*                                         */
/* Nov 11 1999 Ju Li <liju99@mit.edu>      */
/*******************************************/

#include "Scalar.h"


/* Allocate and initialize Spline structure based on */
/* lookup table x[0..N-1],y[0..N-1] and spline type. */
Spline *spline_initialize (int N, double *x, double *y, int spline_type)
{
    register int i;
    Spline *sp;
    MALLOC (spline_initialize, sp, 1, Spline);
    switch (spline_type)
    {
        case SPLINE_CUBIC_NATURAL:
            sp->spline = gsl_spline_alloc (gsl_interp_cspline, N);
            break;
        case SPLINE_CUBIC_PERIODIC:
            sp->spline = gsl_spline_alloc (gsl_interp_cspline_periodic, N);
            break;
        case SPLINE_AKIMA_NATURAL:
            sp->spline = gsl_spline_alloc (gsl_interp_akima, N);
            break;
        case SPLINE_AKIMA_PERIODIC:
            sp->spline = gsl_spline_alloc (gsl_interp_akima_periodic, N);
            break;
        default:
            pe("spline_initialize: "
               "Spline type %d does not exist.\n", spline_type);
    }
    sp->acc = gsl_interp_accel_alloc();
    gsl_spline_init (sp->spline, x, y, N);
    sp->xmin =  DOUBLE_PRECISION_INFINITY;
    sp->xmax = -DOUBLE_PRECISION_INFINITY;
    sp->ymin =  DOUBLE_PRECISION_INFINITY;
    sp->ymax = -DOUBLE_PRECISION_INFINITY;
    for (i=N; i--;)
    {
        if (x[i] < sp->xmin) sp->xmin=x[i];
        if (x[i] > sp->xmax) sp->xmax=x[i];
        if (y[i] < sp->ymin) sp->ymin=y[i];
        if (y[i] > sp->ymax) sp->ymax=y[i];
    }
    return (sp);
} /* end spline_initialize() */


/* Evaluate spline value. */
double spline_value (Spline *sp, double x)
{
    return(gsl_spline_eval (sp->spline, x, sp->acc));
} /* end spline_value() */


/* Evaluate first-order derivative. */
double spline_deriv (Spline *sp, double x)
{
    return(gsl_spline_eval_deriv (sp->spline, x, sp->acc));
} /* end spline_deriv() */


/* Evaluate second-order derivative. */
double spline_deriv2 (Spline *sp, double x)
{
    return(gsl_spline_eval_deriv2 (sp->spline, x, sp->acc));
} /* end spline_deriv2() */


/* Free all data structure. */
void spline_finalize (Spline *sp)
{
    gsl_spline_free (sp->spline);
    gsl_interp_accel_free (sp->acc);
    Free (sp);
    return;
} /* end spline_finalize() */


#ifdef _spline_TEST
#define N 10
int main (int argc, char *argv[])
{
    int i;
    double xi, x[N], y[N];
    Spline *sp;
    printf ("#m=0,S=2\n");
    for (i = 0; i < N; i++)
    {
        x[i] = i + 0.5 * sin (i);
        y[i] = i + cos (i * i);
        printf ("%g %g\n", x[i], y[i]);
    }
    printf ("#m=1,S=0\n");
    sp = spline_initialize (N, x, y, SPLINE_CUBIC_NATURAL);
    for (xi = x[0]; xi < x[9]; xi += 0.01)
        printf ("%g %g\n", xi, spline_value(sp,xi));
    spline_finalize (sp);
    return 0;
}
#undef N
/*
  ./spline_TEST > interp.dat
  graph -T ps < interp.dat > interp.ps
  g interp.ps
*/
#endif /* _spline_TEST */
