/*******************************************/
/* libPhonon: -llapack -lblas -lm          */
/*            -lVecMat -lScalar -lIO       */
/*                                         */
/* Force constant, phonon diagonalization, */
/* and thermodynamic functions.            */
/*                                         */
/* Jul 9 2000 Ju Li <liju99@mit.edu>.      */
/*******************************************/

#ifndef _Phonon_h
#define _Phonon_h

#include <Atoms.h>

/* Phonon.c: */
/*******************************************************/
/* Input: a configuration with proper Chemtab, Tp, and */
/* a Neighborlist with no atoms sitting right on edge, */
/* and which is nc0 x nc1 x nc2 replica of a unit cell */
/* whose atom indices is sequence [0, np/nc0/nc1/nc2). */
/*                                                     */
/* Output: force constant & specification in SI units. */
/*******************************************************/
#define PHONON_DEF_TEST_VOLUMETRIC_STRAIN        (1e-5)
#define PHONON_DEF_TEST_DISPLACEMENT_IN_A        (1e-5)
#define PHONON_DEF_MIN_FORCE_CONSTANT_IN_N__M    (1e-3)
#define PHONON_DEF_FORCE_CONSTANT_MIN_COMPONENT  (1e-5)
#define PHONON_DEF_FCONST_FN_BASENAME            "fconst"
typedef struct
{ /* if zero then default will be used */
    double test_volumetric_strain;       /* to evaluate bulk modulus */
    double test_displacement_in_A;       /* to evaluate force constants */
    double min_force_constant_in_N__M;   /* (9-component) magnitude cutoff */
    double force_constant_min_component; /* component relative ratio cutoff */
    char *fconst_fn_basename;  /* disk load/save interface */
    /* EVERYTHING BELOW ARE IN SI UNITS */
    M3 H;          /* H of unit cell [m] */
    int npa;       /* number of atoms in unit cell */
    double *mass;  /* mass of atoms [kg] */
    char *symbol;  /* atom symbols */
    double B;      /* bulk modulus [Pa] */
    int *cn;       /* how many nonzero interactions for 0 <= i<=j < npa */
    int total;     /* total sum of cn[] */
    /* separation measured in unit cell s, and force constant [N/m] */
    double *fc;
} Phonon_force_constant;

void phonon_force_constant_calc
(Alib_Declare_Config, Chemtab *ct, Tp *tp, Neighborlist *N,
 int nc0, int nc1, int nc2, StaticPotential potential,
 Phonon_force_constant *P, FILE *info);
#define PHONON_force_constant_calc(Config_Alib_to_Alib,ct,tp,N,nc,p,P,info) \
  phonon_force_constant_calc(Config_Alib_to_Alib,ct,tp,N,nc,nc,nc,p,P,info)
#define PHONON_force_constant_CALC(Config_Alib_to_Alib,ct,tp,N,nc,p,P) \
  phonon_force_constant_calc(Config_Alib_to_Alib,ct,tp,N,nc,nc,nc,p,P,stdout)

/* free force constant resident memory */
void phonon_force_constant_free (Phonon_force_constant *P);

/* save force constants to P->fconst_fn_basename */
void phonon_force_constant_save (Phonon_force_constant *P, FILE *info);
#define phonon_force_constant_SAVE(P) phonon_force_constant_save(P,stdout)

/* Load force constants from P->fconst_fn_basename */
void phonon_force_constant_load (Phonon_force_constant *P, FILE *info);
#define phonon_force_constant_LOAD(P) phonon_force_constant_load(P,stdout)

typedef struct
{
    bool gamma_only;
    bool calculate_eigenvector;
    /* workspace */
    double *ap;
    double *work;
    double *rwork;
    /* output: the eigenvalues are stored in ascending order */
    double *w2;
    double *e;
} Phonon_diagonalization;

/* reallocate memory for phonon diagonalization */
void phonon_diagonalization_realloc
(Phonon_force_constant *P, Phonon_diagonalization *D, FILE *info);
#define phonon_diagonalization_REALLOC(P,D) \
  phonon_diagonalization_realloc(P,D,stdout)

/* free all resident space */
void phonon_diagonalization_free (Phonon_diagonalization *D, FILE *info);
#define phonon_diagonalization_FREE(D) phonon_diagonalization_free(D,stdout)

/* Diagonalize Cartesian k=2*PI*ks*(H^-T), or k*x^T=2*PI*ks*s^T */
void phonon_dk (Phonon_force_constant *P, Phonon_diagonalization *D, V3 ks);

/* Diagonalize k=0 using faster routine and less memory */
void phonon_dgamma (Phonon_force_constant *P, Phonon_diagonalization *D);


/* Debye.c: */

/* D(x) = 3 * x^3 * \int_0^{1/x}dy y^4 * e^y / (e^y-1)^2. See Notes. */
double Debye_function (double x);

/*************************************************************/
/* Given Treal and corresponding dTmd__dTreal (heat capacity */
/* normalized to 1),  deduce the Debye temperature.          */
/*************************************************************/
double Debye_temperature (double Treal, double dTmd__dTreal);

/**************************************************************/
/* Given Debye temperature and its associated DOS which leads */
/* to certain quantum average energy at T, find Tmd whereby   */
/* a classical system would have the same average energy.     */
/**************************************************************/
double Debye_Tmd (double TD, double T);

/* TD <= 0 is the signal for no temperature rescaling */
double Debye_Treal_to_Tmd (double TD, double Treal, double *dTmd__dTreal);

/* TD <= 0 is the signal for no temperature rescaling */
double Debye_Tmd_to_Treal (double TD, double Tmd, double *dTmd__dTreal);

#endif  /* _Phonon_h */
