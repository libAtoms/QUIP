/*******************************************/
/* libDSpectral: -ldrfftw -ldfftw -lIO -lm */
/*                                         */
/* Double-precision spectral analyses and  */
/* correlation functions.                  */
/*                                         */
/* Jun. 16, 2000 Ju Li <liju99@mit.edu>.   */
/*******************************************/

#ifndef _DSpectral_h
#define _DSpectral_h

#include <math.h>
#include <Atoms.h>
#include <Scalar.h>
#include <IO.h>
#include <drfftw.h>

/*********************************************************************/
/* Suppose the concerned transport coefficient is K = C \int_0^\inf  */
/* dt <x(0)x(t)>, but we have just x(0), x(d), x(2d) .. x((N-1)d)    */
/* [everything is in SI unit]. Let us define y_i = \sqrt(Cd)x(id),   */
/* then clearly K = \sum_{j=0}^\inf <y_0y_j>. Define z_k =           */
/* \sum_{j=0}^{N-1} y_j e_kj, the unnormalized FFT, then y_j = 1/N   */
/* \sum_{k=-N/2}^{N/2} z_k e^kj, and it is easy to see that <y_0y_j> */
/* = 1/N^2 \sum_{k=-N/2}^{N/2} |z_k|^2 e^kj. Furthermore, there is   */
/* |z_k|^2 = \sum_{j,j'=0}^{N-1} y_j y_j' e_k(j-j') \approx = 2NK    */
/* if k is small compared to the correlation but large to N.         */
/* The first element of a current file stores d in seconds, followed */
/* by y_0, y_1, .. y_{N-1} in \sqrt(SI of the K). The value of d is  */
/* irrelevant to the K result & is used only for labeling purposes.  */
/*********************************************************************/

/* Spectral.c: */

/* write x to binary current in appropriate floating precision */
void spectral_fwrite (double x, FILE *fp);

/* save "function" array to "fn" with delta being the first element */
void spectral_fsave
(fftw_real delta, int cutoff, fftw_real *function, char *fn, FILE *info);

/* sum correlation function up to a certain fraction of the total length */
#define SPECTRAL_CORRELATION_SUM_METHOD_FRACTION  0
#define SPECTRAL_CORRELATION_SUM_DEF_PARAM        40.
/* sum correlation function up to the point of first dip */
#define SPECTRAL_CORRELATION_SUM_METHOD_DIP       1

/* summing the correlation for a rough estim. of the transport coefficient */
fftw_real spectral_sumcorr (fftw_real delta, int N, int cutoff,
                            fftw_real *function, int method,
                            fftw_real parameter, FILE *info);

/* how many channels of the power spectrum is saved */
#define SPECTRAL_POWER_SPECTRUM_DEF_CUTOFF                2048
/* save correlation function up to a certain fraction of the total length */
#define SPECTRAL_CORRELATION_CUTOFF_METHOD_FRACTION       0
#define SPECTRAL_CORRELATION_CUTOFF_DEF_PARAM             10.
/* save correlation function up to a certain multiple of the first dip */
#define SPECTRAL_CORRELATION_CUTOFF_METHOD_DIP            1

typedef struct
{
    int N;  /* total number of time steps */
    fftw_real delta;  /* time-step label in seconds */
    int power_spectrum_cutoff;  /* if 0, then use the default */
    fftw_real *power_spectrum;  /* |z_k|^2 / 2N */
    char *fn_power_spectrum;  /* if not NULL then save to file */
    int correlation_sum_method;  /* just a rough estimate */
    double correlation_sum_parameter;  /* if 0., then use the default */
    /* one should really use Matlab to view the correlation function */
    fftw_real transport_coefficient;
    int correlation_cutoff_method;
    double correlation_cutoff_parameter; /* if 0., then use the default */
    int correlation_cutoff;
    fftw_real *correlation;
    char *fn_correlation;  /* if not NULL then save to file */
} Spectrum;

/* Self-correlation analysis of a single current (with builtin accumulation) */
fftw_real spectral_selfcorr (char *fn, Spectrum *r, FILE *info);

#define spectral_clear(r)  bzero(r, sizeof(Spectrum))
#define spectral_free(r) { Free((r)->power_spectrum); Free((r)->correlation); }

/* Check out "power_spectrum" and "correlation" after "navg" averages */
double spectral_checkout (int navg, Spectrum *r, FILE *info);

#define SPECTRAL_TYPE_HEAT     0
#define SPECTRAL_TYPE_STRESS   1
#define SPECTRAL_TYPE_MAX      2

extern char *spectral_fnflux1[SPECTRAL_TYPE_MAX][1];
extern FILE *spectral_fpflux1[SPECTRAL_TYPE_MAX][1];
extern char *spectral_fnfreq1[SPECTRAL_TYPE_MAX];
extern char *spectral_fncorr1[SPECTRAL_TYPE_MAX];

/* open binary current handles and fill in delta_in_S as first element */
void spectral_openflux1 (double delta_IN_S, int type);

void spectral_writeflux1 (double *flux, double factor_to_SI, int type);

void spectral_flushflux1 (int type);

#define spectral_flushwriteflux1(flux,factor_to_SI,type) ( \
  spectral_writeflux(flux,factor_to_SI,type), spectral_flushflux(type) )

/* close binary current handles */
void spectral_closeflux1 (int type);

/* 1D self-correlation analysis */
double smile1 (int which, FILE *output);

extern char *spectral_fnflux[SPECTRAL_TYPE_MAX][3];
extern FILE *spectral_fpflux[SPECTRAL_TYPE_MAX][3];
extern char *spectral_fnfreq[SPECTRAL_TYPE_MAX];
extern char *spectral_fncorr[SPECTRAL_TYPE_MAX];

/* open binary current handles and fill in delta_in_S as first element */
void spectral_openflux (double delta_IN_S, int type);

void spectral_writeflux (double *flux, double factor_to_SI, int type);

void spectral_flushflux (int type);

#define spectral_flushwriteflux(flux,factor_to_SI,type) ( \
  spectral_writeflux(flux,factor_to_SI,type), spectral_flushflux(type) )

/* close binary current handles */
void spectral_closeflux (int type);

/* directional averaged self-correlation analysis */
double smile (int which, FILE *output);

/* get self-correlation at a certain integer step by direct multiplication */
double poke1 (char *fname, int correlation_step, FILE *info);

#endif  /* _DSpectral_h */
