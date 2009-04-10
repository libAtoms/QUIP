/***********************************************/
/* Spectral analyses and correlation functions */
/*                                             */
/* Jun. 16, 2000 Ju Li <liju99@mit.edu>.       */
/***********************************************/


/********************************************************/
/* http://fftw.org/doc/fftw_1.html#SEC1                 */
/* "FFTW works most efficiently for arrays whose size   */
/* can be factored into small primes (2, 3, 5, and 7)." */
/********************************************************/

#define FFTW_MAX_SMALL_PRIMES  4
static const int FFTW_small_primes [FFTW_MAX_SMALL_PRIMES]
= {2,3,5,7};

/* test if n can be expressed as 2^i * 3^j * 5^k * 7^l ... */
static bool FFTW_small_prime_factorizable (int n)
{
    register int i;
    for (i=0; i<FFTW_MAX_SMALL_PRIMES; i++)
        while (n % FFTW_small_primes[i] == 0)
            n /= FFTW_small_primes[i];
    if (n==1) return (TRUE);
    else return (FALSE);
} /* end FFTW_small_prime_factorizable() */


/* find m<=n but which is small_prime_factorizable */
static int FFTW_fast_n (int n)
{
    if (n > 0)
        while (!FFTW_small_prime_factorizable(n)) n--;
    return (n);
} /* end FFTW_fast_n() */


#ifdef _rfftw_TEST
#define N 10
int main (int argc, char *argv[])
{
    int i;
    fftw_real in[N], out[N];
    rfftw_plan p;
    for (i=0; i<N; i++)
        printf ("%f\n", in[i] = FRANDOM());
    p = rfftw_create_plan (N, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE);
    rfftw_one (p, in, out);
    rfftw_destroy_plan (p);
    p = rfftw_create_plan (N, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE);
    rfftw_one (p, out, in);
    rfftw_destroy_plan (p);
    printf ("\n");
    for (i=0; i<N; i++) printf ("%f\n", in[i]/N);
    return (0);
}
#undef N
#endif /* _rfftw_TEST */


/* write x to binary current in appropriate floating precision */
void spectral_fwrite (double x, FILE *fp)
{
    fftw_real cc = x;
    fwrite (&cc, sizeof(fftw_real), 1, fp);
    return;
} /* end spectral_fwrite() */


/* save "function" array to "fn" with delta being the first element */
void spectral_fsave
(fftw_real delta, int cutoff, fftw_real *function, char *fn, FILE *info)
{
    register int i;
    FILE *fp;
    if (Finvalid(fn)) return;
    fp = wopen(fn);
    fprintf (fp, "%e\n", delta);
    for (i=0; i<cutoff; i++) fprintf (fp, "%e\n", function[i]);
    fclose (fp);
    Fprintf (info, "Saving %d values on \"%s\".\n", cutoff, fn);
    return;
}  /* end spectral_fsave() */


/* summing the correlation for a rough estim. of the transport coefficient */
fftw_real spectral_sumcorr (fftw_real delta, int N, int cutoff,
                            fftw_real *function, int method,
                            fftw_real parameter, FILE *info)
{
    register int i, j;
    fftw_real sum;
    if (method == SPECTRAL_CORRELATION_SUM_METHOD_FRACTION)
    {
        if (parameter == 0) parameter = SPECTRAL_CORRELATION_SUM_DEF_PARAM;
        else if (parameter < 1) parameter = 1;
        j = ceil ( N / parameter );
        if ( j < cutoff )
            Fprintf (info, "Summing corr. to %d / %g = %d steps (%e s),\n",
                     N, parameter, j, j*delta);
        else
        {
            Fprintf (info, "We'd like to sum corr. to %d / %g = "
                     "%d steps,\n", N, parameter, j);
            Fprintf (info, "but the correlation accumulated is not "
                     "that long\n");
            j = cutoff;
            Fprintf (info, "and we have to stop at %d steps (%e s),\n",
                     j, j*delta);
        }
    }
    else if (method == SPECTRAL_CORRELATION_SUM_METHOD_DIP)
    {
        for ( j=0; j<cutoff; j++ ) if (function[j] < 0) break;
        if ( j < cutoff )
            Fprintf (info, "Summing corr. to first dip at %d steps (%e s),\n",
                     j, j*delta);
        else Fprintf (info, "Summing corr. (sees no dip!) to end at %d steps "
                      "(%e s),\n", j, j*delta);
    }
    else pe ("spectral_sumcorr: method %d does not exist.\n", method);
    for (sum=i=0; i<j; i++) sum += function[i];
    Fprintf (info, "we get transport coefficient of %e.\n", sum);
    return (sum);
}  /* end spectral_sumcorr() */


/* Self-correlation analysis of a single current (with builtin accumulation) */
fftw_real spectral_selfcorr (char *fn, Spectrum *r, FILE *info)
{
    register int i, j, k;
    FILE *fp;
    fftw_real *in, *out, cc;
    rfftw_plan p;

    fp = ropen(fn);
    Fprintf (info, "Opening current file \"%s\",\n", fn);
    fread (&cc, sizeof(fftw_real), 1, fp);
    Fprintf (info, "found delta = %e s,\n", cc);
    if ((r->delta != 0) && (r->delta != cc))
        pe ("spectral_selfcorr: new delta does not agree with old\n"
            "%e, cannot accumulate.", r->delta);
    r->delta = cc;

    j = fpsize(fp) / sizeof(fftw_real) - 1;
    i = FFTW_fast_n (j);
    Fprintf (info, "duration = %d (%d lost) steps = %e s,\n",
             i, j-i, i*r->delta);
    if ( (r->N != 0) && (r->N != i) )
    {
        j = i;
        Fprintf (info, "** warning: spectral_selfcorr: new duration\n"
                 "%d does not agree with old duration %d,\n"
                 "truncate to %d.\n", j, r->N, i=MIN(i,r->N));
    }
    r->N = i;

    Fprintf (info, "Loading \"%s\"...\n", fn);
    MALLOC ( spectral_selfcorr, in,  r->N, fftw_real );
    MALLOC ( spectral_selfcorr, out, r->N, fftw_real );
    fread (in, sizeof(fftw_real), r->N, fp );

    Fprintf (info, "raw data loaded, FFT this current...\n");
    p = rfftw_create_plan (r->N, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE);
    rfftw_one (p, in, out);
    rfftw_destroy_plan (p);

    Fprintf (info, "FFT complete, getting the power spectrum...\n");
    out[0] = 0; /* DC component set to 0 */
    j = (r->N + 1) / 2; /* N/2 rounded up */
    for (k=1; k<j; k++)
    {
        out[k] = SQUARE(out[k]/r->N) + SQUARE(out[r->N-k]/r->N);
        out[r->N-k] = 0;
    }
    if (r->N % 2 == 0)  /* N is even */
        out[r->N/2] = SQUARE(out[r->N/2]/r->N);  /* Nyquist freq. */

    if (r->power_spectrum_cutoff == 0)
        r->power_spectrum_cutoff = SPECTRAL_POWER_SPECTRUM_DEF_CUTOFF;
    r->power_spectrum_cutoff = MIN(j, r->power_spectrum_cutoff);
    RECALLOC (spectral_selfcorr, r->power_spectrum,
              r->power_spectrum_cutoff, fftw_real);
    for (k=0; k<r->power_spectrum_cutoff; k++)
    {
        in[k] = out[k] * r->N / 2;
        r->power_spectrum[k] += in[k];
    }
    spectral_fsave (r->delta, r->power_spectrum_cutoff, in,
                    r->fn_power_spectrum, info);

    Fprintf (info, "Inverse FFT the power spectrum...\n");
    p = rfftw_create_plan (r->N, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE);
    rfftw_one (p, out, in);
    rfftw_destroy_plan (p);
    Fprintf (info, "complete, accumulate to counters.\n");
    
    r->transport_coefficient =
        spectral_sumcorr (r->delta, r->N, r->N, in, r->correlation_sum_method,
                          r->correlation_sum_parameter, info);
    if (r->correlation_cutoff_parameter == 0)
        r->correlation_cutoff_parameter =
            SPECTRAL_CORRELATION_CUTOFF_DEF_PARAM;
    if ( r->correlation_cutoff_method ==
         SPECTRAL_CORRELATION_CUTOFF_METHOD_FRACTION )
        i = ceil(r->N / r->correlation_cutoff_parameter);
    else if ( r->correlation_cutoff_method ==
              SPECTRAL_CORRELATION_CUTOFF_METHOD_DIP )
    {
        for (j=0; j<r->N; j++) if (in[j] < 0) break;
        i = MIN( r->N, ceil(r->correlation_cutoff_parameter * j) );
    }
    else pe ("spectral_selfcorr: cutoff method %d does not exist.\n",
             r->correlation_cutoff_method);
    spectral_fsave (r->delta, i, in, r->fn_correlation, info);

    if (r->correlation_cutoff == 0) r->correlation_cutoff = i;
    else r->correlation_cutoff = MIN(i, r->correlation_cutoff);

    RECALLOC ( spectral_selfcorr, r->correlation,
               r->correlation_cutoff, fftw_real );
    for (k=0; k<r->correlation_cutoff; k++)
        r->correlation[k] += in[k];

    free (in);
    free (out);
    return (r->transport_coefficient);
} /* end spectral_selfcorr() */


#ifdef _spectral_selfcorr_TEST
#define L       100
#define N      (L*10000)
#define FNAME  "/tmp/j.out"
int main (int argc, char *argv[])
{
    register int i;
    register double a, factor, correct;
    Spectrum r[1] = {0};
    FILE *fp;
    factor = exp(-1./L);
    correct = sqrt(1./12./(1-SQUARE(factor))/(1-factor));
    fp = wopen(FNAME);
    spectral_fwrite (1., fp);
    for (a=i=0; i<N; i++)
    {
        a = a * factor + FRANDOM() / correct;
        spectral_fwrite (a, fp);
    }
    fclose(fp);
    /* for this unimodal (simple exponential, no oscillation) correlation */
    r->correlation_sum_method = SPECTRAL_CORRELATION_SUM_METHOD_DIP;
    /* function, we override the default FRACTION method and use the DIP. */
    r->fn_power_spectrum = "/tmp/freq.out";
    r->fn_correlation    = "/tmp/corr.out";
    printf ("result = %e\n", spectral_selfcorr(FNAME, r, stdout));
    return (0);
}
#undef FNAME
#undef N
#undef L
#endif /* _spectral_selfcorr_TEST */


#ifdef _selfcorr2_TEST
#define FNAME  "/asm/home/tmp/jx.out"
int main (int argc, char *argv[])
{
    register int i;
    Spectrum r[1] = {0};
    r->correlation_sum_method = SPECTRAL_CORRELATION_SUM_METHOD_DIP;
    r->fn_power_spectrum = "/tmp/freq.out";
    r->fn_correlation    = "/tmp/corr.out";
    printf ("result = %e\n", spectral_selfcorr(FNAME, r, stdout));
    return (0);
}
#endif /* _selfcorr2_TEST */


/* Check out "power_spectrum" and "correlation" after "navg" averages */
double spectral_checkout (int navg, Spectrum *r, FILE *info)
{
    register int i;
    if ( (navg <= 0) || (r->delta <= 0) || (r->N <= 0) )
        pe ("spectral_checkout: no average seems possible.\n");
    if (r->power_spectrum_cutoff > 0)
    {
        Fprintf (info, "Averaging the power spectrum...\n");
        for (i=0; i<r->power_spectrum_cutoff; i++)
            r->power_spectrum[i] /= navg;
        spectral_fsave (r->delta, r->power_spectrum_cutoff,
                        r->power_spectrum, r->fn_power_spectrum, info);
    }
    if (r->correlation_cutoff > 0)
    {
        Fprintf (info, "Averaging the correlation...\n");
        for (i=0; i<r->correlation_cutoff; i++)
            r->correlation[i] /= navg;
        spectral_fsave (r->delta, r->correlation_cutoff,
                        r->correlation, r->fn_correlation, info);
        r->transport_coefficient =
            spectral_sumcorr (r->delta, r->N, r->correlation_cutoff,
                              r->correlation, r->correlation_sum_method,
                              r->correlation_sum_parameter, info);
    }
    else r->transport_coefficient = 0;
    spectral_free (r);
    return (r->transport_coefficient);
} /* end spectral_checkout() */


#define DIM 1
char *spectral_fnflux1[SPECTRAL_TYPE_MAX][DIM] = {{"j.out"}, {"s.out"}};
FILE *spectral_fpflux1[SPECTRAL_TYPE_MAX][DIM] = {0};
char *spectral_fnfreq1[SPECTRAL_TYPE_MAX] = {"freq.out", "sfreq.out"};
char *spectral_fncorr1[SPECTRAL_TYPE_MAX] = {"corr.out", "scorr.out"};

/* open binary current handles and fill in delta_in_S as first element */
void spectral_openflux1 (double delta_IN_S, int type)
{
    int i;
    for (i=0; i<DIM; i++)
    {
        spectral_fpflux1[type][i] = wopen ( spectral_fnflux1[type][i] );
        spectral_fwrite ( delta_IN_S, spectral_fpflux1[type][i] );
    }
    return;
} /* end spectral_openflux1() */


void spectral_writeflux1 (double *flux, double factor_to_SI, int type)
{
    int i;
    for (i=0; i<DIM; i++)
        if (spectral_fpflux1[type][i])
            spectral_fwrite ( factor_to_SI * flux[i],
                              spectral_fpflux1[type][i] );
    return;
} /* end spectral_writeflux1() */


void spectral_flushflux1 (int type)
{
    int i;
    for (i=0; i<DIM; i++)
        if (spectral_fpflux1[type][i])
            fflush (spectral_fpflux1[type][i]);
    return;
} /* end spectral_flushflux1() */


/* close binary current handles */
void spectral_closeflux1 (int type)
{
    int i;
    for (i=0; i<DIM; i++)
        if (spectral_fpflux1[type][i])
        {
            fclose (spectral_fpflux1[type][i]);
            spectral_fpflux1[type][i] = NULL;
        }
    return;
} /* end spectral_closeflux1() */


/* 1D self-correlation analysis */
double smile1 (int type, FILE *info)
{
    int i;
    Spectrum r[1] = {0};
    if ( OUW(type, SPECTRAL_TYPE_MAX) )
        pe ("smile1: invalid flux type %d.\n", type);
    for (i=0; i<DIM; i++)
        spectral_selfcorr (spectral_fnflux1[type][i], r, info);
    /* r->correlation_sum_method = SPECTRAL_CORRELATION_SUM_METHOD_DIP; */
    r->fn_power_spectrum = spectral_fnfreq1[type];
    r->fn_correlation    = spectral_fncorr1[type];
    spectral_checkout (DIM, r, info);
    return (r->transport_coefficient);
} /* end smile1() */


#ifdef _smile1_TEST
#define L      10
#define N     (L*10000)
#define TYPE   SPECTRAL_TYPE_HEAT
int main (int argc, char *argv[])
{
    register int i,j;
    double flux[DIM], factor, correct;
    factor = exp(-1./L);
    correct = sqrt(1./12./(1-SQUARE(factor))/(1-factor));
    spectral_openflux1 (1., TYPE);
    flux[0] = 0;
    for (i=0; i<N; i++)
    {
        for (j=0; j<DIM; j++)
            flux[j] = flux[j] * factor + FRANDOM() / correct;
        spectral_writeflux1 (flux, 1., TYPE);
    }
    spectral_closeflux1 (TYPE);
    smile1 (TYPE, stdout);
    return(0);
} /* end main() */
#undef TYPE
#undef N
#undef L
#endif
#undef DIM


char *spectral_fnflux[SPECTRAL_TYPE_MAX][DIMENSION] =
{{"jx.out", "jy.out", "jz.out"}, {"syz.out", "sxz.out", "sxy.out"}};
FILE *spectral_fpflux[SPECTRAL_TYPE_MAX][DIMENSION] = {0};
char *spectral_fnfreq[SPECTRAL_TYPE_MAX] = {"freq.out", "sfreq.out"};
char *spectral_fncorr[SPECTRAL_TYPE_MAX] = {"corr.out", "scorr.out"};

/* open binary current handles and fill in delta_in_S as first element */
void spectral_openflux (double delta_IN_S, int type)
{
    int i;
    for (i=0; i<DIMENSION; i++)
    {
        spectral_fpflux[type][i] = wopen ( spectral_fnflux[type][i] );
        spectral_fwrite ( delta_IN_S, spectral_fpflux[type][i] );
    }
    return;
} /* end spectral_openflux() */


void spectral_writeflux (double *flux, double factor_to_SI, int type)
{
    int i;
    for (i=0; i<DIMENSION; i++)
        if (spectral_fpflux[type][i])
            spectral_fwrite ( factor_to_SI * flux[i],
                              spectral_fpflux[type][i] );
    return;
} /* end spectral_writeflux() */


void spectral_flushflux (int type)
{
    int i;
    for (i=0; i<DIMENSION; i++)
        if (spectral_fpflux[type][i])
            fflush (spectral_fpflux[type][i]);
    return;
} /* end spectral_flushflux() */


/* close binary current handles */
void spectral_closeflux (int type)
{
    int i;
    for (i=0; i<DIMENSION; i++)
        if (spectral_fpflux[type][i])
        {
            fclose (spectral_fpflux[type][i]);
            spectral_fpflux[type][i] = NULL;
        }
    return;
} /* end spectral_closeflux() */


/* directional averaged self-correlation analysis */
double smile (int type, FILE *info)
{
    int i;
    Spectrum r[1] = {0};
    if ( OUW(type, SPECTRAL_TYPE_MAX) )
        pe ("smile: invalid flux type %d.\n", type);
    for (i=0; i<DIMENSION; i++)
        spectral_selfcorr (spectral_fnflux[type][i], r, info);
    /* r->correlation_sum_method = SPECTRAL_CORRELATION_SUM_METHOD_DIP; */
    r->fn_power_spectrum = spectral_fnfreq[type];
    r->fn_correlation    = spectral_fncorr[type];
    spectral_checkout (DIMENSION, r, info);
    return (r->transport_coefficient);
} /* end smile() */


#ifdef _smile_TEST
#define L      10
#define N     (L*10000)
#define TYPE   SPECTRAL_TYPE_HEAT
int main (int argc, char *argv[])
{
    register int i,j;
    double flux[DIMENSION], factor, correct;
    factor = exp(-1./L);
    correct = sqrt(1./12./(1-SQUARE(factor))/(1-factor));
    spectral_openflux (1., TYPE);
    V3ZERO (flux);
    for (i=0; i<N; i++)
    {
        for (j=0; j<DIMENSION; j++)
            flux[j] = flux[j] * factor + FRANDOM() / correct;
        spectral_writeflux (flux, 1., TYPE);
    }
    spectral_closeflux (TYPE);
    smile (TYPE, stdout);
    return(0);
} /* end main() */
#undef TYPE
#undef N
#undef L
#endif


/* get self-correlation at a certain integer step by direct multiplication */
double poke1 (char *fname, int correlation_step, FILE *info)
{
    FILE *fp;
    fftw_real delta,*memory,*memory_end,*mp;
    int i,total_step;
    double datum,selfcorr_sum,correlation_sum;

    /* time-invariant */
    if (correlation_step<0) correlation_step=-correlation_step;
    fp = ropen(fname);
    Fprintf (info, "Opening current file \"%s\",\n", fname);
    fread (&delta, sizeof(fftw_real), 1, fp);
    Fprintf (info, "found delta = %e s.\n", delta);

    total_step = fpsize(fp) / sizeof(fftw_real) - 1;
    if (total_step < correlation_step+1)
        pe ("Sorry, \"%s\" contains only %d data,\n"
            "cannot do %d-step correlation.\n",
            fname, total_step, correlation_step);

    if (correlation_step > 0)
    {
        MALLOC(main,memory,correlation_step,fftw_real);
        memory_end = memory + correlation_step;
        mp = memory;
        fread(memory,sizeof(fftw_real),correlation_step,fp);
        selfcorr_sum = correlation_sum = 0;
        for (i=0; i<total_step-correlation_step; i++)
        {
            datum = *mp;
            selfcorr_sum += datum*datum;
            fread(mp,sizeof(fftw_real),1,fp);
            correlation_sum += datum*(*mp);
            mp++;
            if (mp == memory_end) mp = memory;
        }
        free(memory);
    }
    else
    {
        selfcorr_sum = 0;
        for (i=0; i<total_step-correlation_step; i++)
        {
            fread(&datum,sizeof(fftw_real),1,fp);
            selfcorr_sum += datum*datum;
        }
        correlation_sum = selfcorr_sum;
    }

    Fprintf (info, "Summing corr. for %d steps (%e s),\n",
             i, i*delta);
    Fprintf (info, "<x^2> = %e,\n", selfcorr_sum/i);
    Fprintf (info, "<x_0x_%d> = <x(0)x(%e s)> = %e,\n",
             correlation_step, correlation_step*delta,
             correlation_sum/i);
    Fprintf (info, "<x(0)x(%e s)> / <x^2> = %e.\n",
             correlation_step*delta,correlation_sum/selfcorr_sum);
    return(correlation_sum/i);
} /* end poke1() */
