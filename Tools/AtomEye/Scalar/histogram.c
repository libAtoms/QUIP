/*******************************************/
/* libScalar:                              */
/*                                         */
/* Scalar Arithmetic Operations and Macros */
/*                                         */
/* Nov 11 1999 Ju Li <liju99@mit.edu>      */
/*******************************************/

#include "Scalar.h"


/* Given finite xmin < xmax and bins >= 1, reallocate */
/* data structure in h and reset all counters.        */
void histogram_reinitialize
(double xmin, double xmax, int bins, Histogram *h)
{
    register int i;
    if (xmin >= xmax)
        pe( "histogram_reinitialize: xmin=%g >= xmax=%g.\n", xmin, xmax);
    if (bins < 1)
        pe( "histogram_reinitialize: bins=%d < 1.\n", bins);
    h->xmin  = xmin;
    h->xmax  = xmax;
    h->bins = bins;
    h->xbin  = (xmax - xmin) / bins;
    h->totalcount = 0;
    h->domaincount = 0;
    REALLOC_ZERO( histogram_reinitialize, h->bincount,        bins, double );
    REALLOC     ( histogram_reinitialize, h->x_sample,        bins, double );
    for (i=0; i<h->bins; i++)
        h->x_sample[i]   = h->xmin + (i + 0.5) * h->xbin;
    REALLOC_ZERO( histogram_reinitialize, h->density_profile, bins, double );
    return;
} /* end histogram_reinitialize() */


/* Update histogram by data array x[0..n-1]; return domaincount */
double histogram_intercept (int n, double *x, Histogram *h)
{
    register int i, j;
    register double d;
    for (i=n; i--;)
    {
        d = (x[i] - h->xmin) / h->xbin;
        if ( (d < 0) || (d > h->bins) ) continue;
        else if (d == 0)
        {
            h->bincount[0] += 0.5;
            h->domaincount += 0.5;
            continue;
        }
        else if (d == h->bins)
        {
            h->bincount[h->bins-1] += 0.5;
            h->domaincount += 0.5;
            continue;
        }
        j = (int) d;
        if (d == j)
        {
            h->bincount[j]   += 0.5;
            h->bincount[j-1] += 0.5;
            h->domaincount ++;
            continue;
        }
        h->bincount[j]++;
        h->domaincount++;
    }
    h->totalcount += n;
    return( h->domaincount );
} /* end histogram_intercept() */


/* Return the interception rate = domaincount / totalcount */
double histogram_interception_rate (Histogram *h)
{
    if ( h->totalcount <= 0 ) return(0);
    return( h->domaincount / h->totalcount );
} /* end histogram_interception_rate() */


/* Update the density profile; return the interception rate */
double histogram_update_density_profile (Histogram *h)
{
    register int i;
    if ( h->totalcount <= 0 ) return(0);
    for (i=0; i<h->bins; i++)
        h->density_profile[i] = h->bincount[i] / h->totalcount / h->xbin;
    return( h->domaincount / h->totalcount );
} /* end histogram_update_density_profile() */


/* Save histogram in matlab script "fname" as  */
/* symbol "token" and notify I/O stream "info" */
double histogram_save_as_matlab
(Histogram *h, FILE *info, char *token, char *fname)
{
    register int i;
    FILE *fp;
    char *tok, *fn, *p;
    if ( h->totalcount <= 0 )
    {
        Fprintf(info, "histogram_save_as_matlab: "
                "there is no data in \"%s\".\n", token);
        return(0);
    }
    tok = IOClone(token);
    for (p=tok; *p!=EOS; p++) if (*p==' ') *p='_';
    histogram_update_density_profile (h);
    fn = IOClone(fname);
    for (p=fn; *p!=EOS; p++) if (*p==' ') *p='_';
    fp = wOpen(fn);
    fprintf (fp, "%% totalcount = %d, domaincount = %g, "
             "interception rate = %.3g%%\n",
             h->totalcount, h->domaincount,
             100. * h->domaincount / h->totalcount);
    fprintf (fp, "%s = [\n", tok);
    for (i=0; i<h->bins; i++)
        fprintf (fp, "%e %e\n", h->x_sample[i], h->density_profile[i]);
    fprintf (fp, "];\n\n");
    fprintf (fp, "clf;\n");
    fprintf (fp, "plot(%s(:,1),%s(:,2));\n", tok, tok);
    fprintf (fp, "v = axis;\n");
    fprintf (fp, "axis([%g %g 0 v(4)]);\n", h->xmin, h->xmax);
    fprintf (fp, "xlabel('%s'); ylabel('probability density of %s');\n\n",
             token, token);
    fclose( fp );
    Fprintf(info, "histogram of \"%s\" is saved on \"%s\".\n",
            tok, fn);
    free (tok);
    free (fn);
    return( h->domaincount / h->totalcount );
} /* end histogram_save_as_matlab() */


/* Clear all counters as if it is just the beginning */
void histogram_reset_all_counters (Histogram *h)
{
    bzero (VOIDP(h->bincount),        h->bins*sizeof(double));
    bzero (VOIDP(h->density_profile), h->bins*sizeof(double));
    h->totalcount = 0;
    h->domaincount = 0;
    return;
} /* end histogram_clear_all_counters() */


/* Free the histogram data structure and reset all counters */
void histogram_free (Histogram *h)
{
    histogram_reset_all_counters (h);
    h->xmin  = 0;
    h->xmax  = 0;
    h->bins = 0;
    Free ( h->bincount );
    Free ( h->x_sample );
    Free ( h->density_profile );
} /* end histogram_Free() */


#ifdef _histogram_TEST
#define N  1024
#define M (1024 * 8)
#define XMIN -0.5
#define XMAX  1.5
#define BINS  100
int main (int argc, char *argv[])
{
    int i,j;
    double x[N];
    Histogram h[1] = {{0}};
    histogram_reinitialize (XMIN, XMAX, BINS, h);
    TimeRandomize();
    for (i=0; i<M; i++)
    {
        for (j=0; j<N; j++) x[j] = FRANDNORM();
        histogram_intercept (N, x, h);
    }
    printf ("interception rate = %g (%g).\n",
            histogram_save_as_MATLAB(h,stdout,"a"),
            (erf(XMAX/SQRT2)-erf(XMIN/SQRT2))/2);
    histogram_free (h);
    return (0);
}
/* a_h; hold on; plot( a(:,1), exp(-a(:,1).^2/2) / sqrt(2*pi), 'r'); */
#endif /* _histogram_TEST */
