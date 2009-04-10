/*******************************************/
/* libPhonon: -llapack -lblas -lm          */
/*            -lVecMat -lScalar -lIO       */
/*                                         */
/* Force constant, phonon diagonalization, */
/* and thermodynamic functions.            */
/*                                         */
/* Jul 9 2000 Ju Li <liju99@mit.edu>.      */
/*******************************************/

#include "Phonon.h"

/*******************************************************/
/* Input: a configuration with proper Chemtab, Tp, and */
/* a Neighborlist with no atoms sitting right on edge, */
/* and which is nc0 x nc1 x nc2 replica of a unit cell */
/* whose atom indices is sequence [0, np/nc0/nc1/nc2). */
/*                                                     */
/* Output: force constant & specification in SI units. */
/*******************************************************/

void phonon_force_constant_calc
(Alib_Declare_Config, Chemtab *ct, Tp *tp, Neighborlist *N,
 int nc0, int nc1, int nc2, StaticPotential potential,
 Phonon_force_constant *P, FILE *info)
{
    int i, j, k, l, m, n, w, *cn, *idx;
    double cc, delta, *f, *force, *fc;
    V3 old_s;
    M3 stress, HH, HI;

    /* fill in defaults */
    if (P->test_volumetric_strain == 0)
        P->test_volumetric_strain = PHONON_DEF_TEST_VOLUMETRIC_STRAIN;
    if (P->test_displacement_in_A == 0)
        P->test_displacement_in_A = PHONON_DEF_TEST_DISPLACEMENT_IN_A;
    if (P->min_force_constant_in_N__M == 0)
        P->min_force_constant_in_N__M = PHONON_DEF_MIN_FORCE_CONSTANT_IN_N__M;
    if (P->force_constant_min_component == 0)
        P->force_constant_min_component =
            PHONON_DEF_FORCE_CONSTANT_MIN_COMPONENT;

    i = nc0 * nc1 * nc2;
    if ((*np) % i != 0)
        pe ("phonon_force_constant_calc: you really screwed up:\n"
            "np = %d is NOT an integer multiple of %d x %d x %d.\n",
            (*np), nc0, nc1, nc2);
    /* unit cell properties in SI */
    V3MUL ( ulength_IN_M / nc0, H[0], P->H[0] );
    V3MUL ( ulength_IN_M / nc1, H[1], P->H[1] );
    V3MUL ( ulength_IN_M / nc2, H[2], P->H[2] );
    P->npa = (*np) / i;

    MALLOC (phonon_force_constant_calc, f, DIMENSION*(*np), double );
    (*potential) ((*np), H, (*s), tp, N, f, stress);
    VLENGTH (DIMENSION*(*np), f, i, cc);
    cc *= uforce_IN_N / sqrt(DOUBLE(*np));
    if ( cc > P->min_force_constant_in_N__M *
         P->test_displacement_in_A * A_IN_M )
        Fprintf (info, "phonon_force_constant_calc: warning: residual\n"
                 "|force| = %e N/m is big, you may have a problem.\n", cc);

    MALLOC ( phonon_force_constant_calc, force, DIMENSION2*(*np), double );
    MALLOC ( phonon_force_constant_calc, idx, (*np), int );
    i = P->npa * (P->npa + 1) / 2;
    REALLOC ( phonon_force_constant_calc, P->cn, i, int );

    /* calculate the unrelaxed bulk modulus */
    delta = P->test_volumetric_strain;
    M3MULTIPLY ( 1-delta/6, H, HH );
    (*potential) ((*np), HH, (*s), tp, N, f, stress);
    cc = M3TR(stress) / 3;
    M3MULTIPLY ( 1+delta/6, H, HH );
    (*potential) ((*np), HH, (*s), tp, N, f, stress);
    P->B = (cc-M3TR(stress)/3) / delta * ustress_IN_PA;

    /* calculate the force constants */
    delta = P->test_displacement_in_A / ulength_IN_A;
    M3INV (H, HI, cc);
    for (cn=P->cn,fc=P->fc,j=0; j<P->npa; j++)
    {
        V3EQV( &(*s)[DIMENSION*j], old_s );
        for (n=0; n<DIMENSION; n++)
        {
            V3SUBMUL (old_s, delta, HI[n], &(*s)[DIMENSION*j]);
            (*potential) ( (*np), H, (*s), tp, N, f, stress );
            V3ADDMUL (old_s, delta, HI[n], &(*s)[DIMENSION*j]);
            (*potential)
                ( (*np), H, (*s), tp, N, force+n*DIMENSION*(*np), stress );
            VSuB (DIMENSION*(*np), force+n*DIMENSION*(*np), f, k);
        }
        V3EQV ( old_s, &(*s)[DIMENSION*j] );
        cc = - 0.5 * uforce_IN_N / ulength_IN_M / delta;
        VMuL ( DIMENSION2*(*np), cc, force, k );
        for (l=k=0; k<(*np); k++)
        {
            f[k] = sqrt ( V3LENGTH2(force+DIMENSION*k) +
                          V3LENGTH2(force+DIMENSION*(k+(*np))) +
                          V3LENGTH2(force+DIMENSION*(k+2*(*np))) );
            if ((k % P->npa <= j) && (f[k] > P->min_force_constant_in_N__M))
            { /* the sorter gives ascending order, so... */
                f[k] = -f[k];
                idx[l++] = k;
            }
        }
        qsort_numerical_recipes ( l, f, idx, USE_OLD_IDX );
        k = fc - P->fc;
        REALLOC ( phonon_force_constant_calc, P->fc,
                  k + (DIMENSION+DIMENSION2) * l, double );
        fc = P->fc + k;
        for (i=0; i<=j; i++,cn++)
        {
            *cn = 0;
            for (m=0; m<l; m++)
            {
                k = idx[m];
                if (k % P->npa != i) continue;
                V3SUB ( &(*s)[DIMENSION*j], &(*s)[DIMENSION*k], fc );
                V3ImagE ( fc );
                V3mUL ( nc0, nc1, nc2, fc );
                fc += DIMENSION;
                for (n=0; n<DIMENSION; n++)
                    for (w=0; w<DIMENSION; w++,fc++)
                    {
                        *fc = force [n*DIMENSION*(*np) + DIMENSION*k + w];
                        if ( fabs(*fc / f[k]) <
                             P->force_constant_min_component ) *fc = 0;
                    }
                (*cn)++;
            }
        }
    }
    P->total = (fc - P->fc) / (DIMENSION+DIMENSION2);
    VRECLONEMul (P->npa, umass_IN_KG, (*mass), P->mass, i);
    RECLONE ((*symbol), P->npa * SYMBOL_SIZE, char, P->symbol);
    Free (idx);
    Free (force);
    Free (f);
    return;
} /* end phonon_force_constant_calc() */


/* free resident memory */
void phonon_force_constant_free (Phonon_force_constant *P)
{
    Free (P->symbol);
    Free (P->mass);
    Free (P->cn);
    Free (P->fc);
    return;
} /* end phonon_force_constant_save() */


/* save force constants to P->fconst_fn_basename */
void phonon_force_constant_save (Phonon_force_constant *P, FILE *info)
{
    int i, j, k, w, *cn;
    double *fc;
    TermString fname;
    FILE *fp;

    if (P->fconst_fn_basename == NULL)
        P->fconst_fn_basename = PHONON_DEF_FCONST_FN_BASENAME;
    sprintf (fname, "%s.spec", P->fconst_fn_basename);
    fp = wopen(fname);
    for (i=0; i<DIMENSION; i++)
    {
        for (j=0; j<DIMENSION; j++) fprintf(fp, "%g ", P->H[i][j]);
        fprintf(fp, "\n");
    }
    fprintf(fp, "%d\n", P->npa);
    for (i=0; i<P->npa; i++)
        fprintf(fp, "%2s %e\n", Symbol(P->symbol,i),P->mass[i]);
    fprintf (fp, "%e\n", P->B);
    for (cn=P->cn; cn<P->cn+P->npa*(P->npa+1)/2; cn++)
        fprintf (fp, "%d\n", *cn);
    fclose (fp);
    Fprintf (info, "Force constant specification is saved in "
             "\"%s\".\n", fname);
    sprintf (fname, "%s.out",  P->fconst_fn_basename);
    fp = wopen(fname);
    for (cn=P->cn,fc=P->fc,i=0; i<P->npa; i++)
        for (j=i; j<P->npa; j++,cn++)
        {
            for (k=0; k<*cn; k++)
            {
                for (w=0; w<DIMENSION+DIMENSION2; w++,fc++)
                    fprintf (fp, "%g ", *fc);
                fprintf (fp, "\n");
            }
            Fprintf (info, "%d(%s)-%d(%s): %d\n", i,Symbol(P->symbol,i),
                     j,Symbol(P->symbol,j), *cn);
        }
    fclose(fp);
    Fprintf (info, "Force constant values are saved in "
             "\"%s\".\n", fname);
    return;
} /* end phonon_force_constant_save() */


/* Load force constants from P->fconst_fn_basename */
void phonon_force_constant_load (Phonon_force_constant *P, FILE *info)
{
    int i, j, k, w, *cn;
    double *fc;
    char as[SYMBOL_SIZE+1]={' '};
    TermString fname;
    FILE *fp;

    if (P->fconst_fn_basename == NULL)
        P->fconst_fn_basename = PHONON_DEF_FCONST_FN_BASENAME;
    sprintf (fname, "%s.spec", P->fconst_fn_basename);
    fp = ropen(fname);
    for (i=0; i<DIMENSION; i++)
    {
        for (j=0; j<DIMENSION; j++)
            fscanf(fp, "%lf ", &P->H[i][j]);
        fscanf(fp, "\n");
    }
    fscanf(fp, "%d\n", &P->npa);
    REALLOC ( phonon_force_constant_load, P->mass,   P->npa, double );
    REALLOC ( phonon_force_constant_load, P->symbol, P->npa * SYMBOL_SIZE,
              char );
    for (i=0; i<P->npa; i++)
    {
        fscanf(fp, "%s %lf\n", as+1, &P->mass[i]);
        for (j=2; as[j]!=EOS; j++);
        for (k=j-SYMBOL_CHAR; k<=j; k++)
            Symbol(P->symbol,i)[k+SYMBOL_CHAR-j] = as[k];
    }
    fscanf (fp, "%lf\n", &P->B);
    i = P->npa * (P->npa + 1) / 2;
    REALLOC ( phonon_force_constant_load, P->cn, i, int );
    for (P->total=0,cn=P->cn; cn<P->cn+i; cn++)
    {
        fscanf (fp, "%d\n", cn);
        P->total += *cn;
    }
    fclose (fp);
    Fprintf (info, "Force constant specification loaded from "
             "\"%s\".\n", fname);
    Fprintf (info, "Bulk modulus = %e Pa.\n", P->B);
    Fprintf (info, "%d atoms in unit cell:\n", P->npa);
    REALLOC ( phonon_force_constant_load, P->fc,
              (DIMENSION+DIMENSION2) * P->total, double );
    sprintf (fname, "%s.out",  P->fconst_fn_basename);
    fp = ropen(fname);
    for (cn=P->cn,fc=P->fc,i=0; i<P->npa; i++)
        for (j=i; j<P->npa; j++,cn++)
        {
            for (k=0; k<*cn; k++)
            {
                for (w=0; w<DIMENSION+DIMENSION2; w++,fc++)
                    fscanf (fp, "%lf ", fc);
                fscanf (fp, "\n");
            }
            Fprintf (info, "%d(%s)-%d(%s): %d\n", i,Symbol(P->symbol,i),
                     j,Symbol(P->symbol,j), *cn);
        }
    fclose(fp);
    Fprintf (info, "Force constant values loaded from \"%s\".\n", fname);
    return;
} /* end phonon_force_constant_load() */


/* reallocate memory for phonon diagonalization */
void phonon_diagonalization_realloc
(Phonon_force_constant *P, Phonon_diagonalization *D, FILE *info)
{
    int i=0, j;

    j = 3 * P->npa * (3 * P->npa + 1);
    if (D->gamma_only) j /= 2;
    REALLOC ( phonon_diagonalization_realloc, D->ap, j, double ); i+=j;

    j = 4 * 3 * P->npa;
    if (D->gamma_only) j = 0;
    REALLOC ( phonon_diagonalization_realloc, D->work, j, double ); i+=j;

    j = 3 * 3 * P->npa;
    REALLOC ( phonon_diagonalization_realloc, D->rwork, j, double ); i+=j;

    j = 3 * P->npa;
    REALLOC ( phonon_diagonalization_realloc, D->w2, j, double ); i+=j;

    if (D->calculate_eigenvector)
    {
        j = 2 * 3 * P->npa * 3 * P->npa;
        if (D->gamma_only) j /= 2;
        REALLOC ( phonon_diagonalization_realloc, D->e, j, double ); i+=j;
    }
    Fprintf (info, "Phonon diagonalization now occupies %d bytes memory.\n",
             i * sizeof(double));
    return;
} /* end phonon_diagonalization_realloc() */


/* free all resident space */
void phonon_diagonalization_free (Phonon_diagonalization *D, FILE *info)
{
    Free (D->ap);
    Free (D->work);
    Free (D->rwork);
    Free (D->w2);
    Free (D->e);
    Fprintf (info, "Phonon diagonalization memory freed.\n");
    return;
} /* end phonon_diagonalization_free() */


/* Diagonalize Cartesian k=2*PI*ks*(H^-T), or k*x^T=2*PI*ks*s^T */
void phonon_dk (Phonon_force_constant *P, Phonon_diagonalization *D, V3 ks)
{
    int i, j, k, *cn;
    register double *fc, *dp, cc, dd;
    double ee, ff;
    static double pie [2*DIMENSION2];
    for (cn=P->cn,fc=P->fc,j=0; j<P->npa; j++)
    {
        for (i=0; i<j; i++,cn++)
        {
            pie[0] = 0;
            pie[1] = 0;
            pie[2] = 0;
            pie[3] = 0;
            pie[4] = 0;
            pie[5] = 0;
            pie[6] = 0;
            pie[7] = 0;
            pie[8] = 0;
            pie[9] = 0;
            pie[10] = 0;
            pie[11] = 0;
            pie[12] = 0;
            pie[13] = 0;
            pie[14] = 0;
            pie[15] = 0;
            pie[16] = 0;
            pie[17] = 0;
            ff = sqrt(P->mass[i]*P->mass[j]);
            for (k=*cn; k--;)
            {
                ee = 2 * PI * V3DOT( ks, fc );
                fc += DIMENSION;
                cc = cos(ee) / ff;
                dd = sin(ee) / ff;
                pie[0] += cc * (*fc);
                pie[1] += dd * (*(fc++));
                pie[2] += cc * (*fc);
                pie[3] += dd * (*(fc++));
                pie[4] += cc * (*fc);
                pie[5] += dd * (*(fc++));
                pie[6] += cc * (*fc);
                pie[7] += dd * (*(fc++));
                pie[8] += cc * (*fc);
                pie[9] += dd * (*(fc++));
                pie[10] += cc * (*fc);
                pie[11] += dd * (*(fc++));
                pie[12] += cc * (*fc);
                pie[13] += dd * (*(fc++));
                pie[14] += cc * (*fc);
                pie[15] += dd * (*(fc++));
                pie[16] += cc * (*fc);
                pie[17] += dd * (*(fc++));
            }
            dp = D->ap + 3*j * (3*j + 1) + 6*i;
            dp[0] = pie[0];
            dp[1] = pie[1];
            dp[2] = pie[2];
            dp[3] = pie[3];
            dp[4] = pie[4];
            dp[5] = pie[5];
            dp += 6 * j + 2;
            dp[0] = pie[6];
            dp[1] = pie[7];
            dp[2] = pie[8];
            dp[3] = pie[9];
            dp[4] = pie[10];
            dp[5] = pie[11];
            dp += 6 * j + 4;
            dp[0] = pie[12];
            dp[1] = pie[13];
            dp[2] = pie[14];
            dp[3] = pie[15];
            dp[4] = pie[16];
            dp[5] = pie[17];
        }
        pie[0] = 0;
        pie[1] = 0;
        pie[6] = 0;
        pie[7] = 0;
        pie[8] = 0;
        pie[9] = 0;
        pie[12] = 0;
        pie[13] = 0;
        pie[14] = 0;
        pie[15] = 0;
        pie[16] = 0;
        pie[17] = 0;
        for (k=*cn; k--;)
        {
            ee = 2 * PI * V3DOT( ks, fc );
            fc += DIMENSION;
            cc = cos(ee) / P->mass[j];
            dd = sin(ee) / P->mass[j];
            pie[0] += cc * (*fc);
            pie[1] += dd * (*fc);
            fc += 3;
            pie[6] += cc * (*fc);
            pie[7] += dd * (*(fc++));
            pie[8] += cc * (*fc);
            pie[9] += dd * (*fc);
            fc += 2;
            pie[12] += cc * (*fc);
            pie[13] += dd * (*(fc++));
            pie[14] += cc * (*fc);
            pie[15] += dd * (*(fc++));
            pie[16] += cc * (*fc);
            pie[17] += dd * (*(fc++));
        }
        dp = D->ap + 3*j * (3*j + 1) + 6*i;
        dp[0] = pie[0];
        dp[1] = pie[1];
        dp += 6 * j + 2;
        dp[0] = pie[6];
        dp[1] = pie[7];
        dp[2] = pie[8];
        dp[3] = pie[9];
        dp += 6 * j + 4;
        dp[0] = pie[12];
        dp[1] = pie[13];
        dp[2] = pie[14];
        dp[3] = pie[15];
        dp[4] = pie[16];
        dp[5] = pie[17];
        cn++;
    }
    i = 0;
    j = 3 * P->npa;
    FORTRAN_SYMBOL(zhpev)
        (D->calculate_eigenvector ? "V" : "N",  "U",
         &j, D->ap, D->w2, D->e, &j, D->work, D->rwork, &i);
    if (i < 0) pe("phonon_dk: %d-th argument had an illegal value.\n", -i);
    else if (i > 0) pe ("phonon_dk: LAPACK failed to converge (%d).\n", i);
    return;
} /* end phonon_dk() */


#ifdef _dk
int main (int argc, char *argv[])
{
    int i;
    Phonon_force_constant  P[1]={0};
    Phonon_diagonalization D[1]={0};
    V3 ks;
    if (argc == 5)
    {
        P->fconst_fn_basename = argv[1];
        ks[0] = atof(argv[2]);
        ks[1] = atof(argv[3]);
        ks[2] = atof(argv[4]);
    }
    else if (argc == 4)
    {
        ks[0] = atof(argv[1]);
        ks[1] = atof(argv[2]);
        ks[2] = atof(argv[3]);
    }
    else
    {
        printf ("Usage: %s [force const file basename] ks0 ks1 ks2\n",
                argv[0]);
        return (1);
    }
    phonon_force_constant_LOAD (P);
    phonon_diagonalization_REALLOC (P,D);
    phonon_dk (P, D, ks);
    for (i=0; i<3*P->npa; i++)
    {
        printf ("%e ", D->w2[i]);
        if (D->w2[i] < 0) printf ("!!\n");
        else printf (" -> %.8g THz\n", sqrt(D->w2[i])/2/PI*HZ_IN_THZ);
    }
    return (0);
}
#endif /* _dk */


#ifdef _cv
#define CV_REPORT_FREQ  50000
int main (int argc, char *argv[])
{
    int i, w_sample;
    double a, b, c, T, CV;
    V3 ks, k;
    M3 HI;
    Phonon_force_constant  P[1]={0};
    Phonon_diagonalization D[1]={0};

    if (argc == 3)
    {
        P->fconst_fn_basename = argv[1];
        T = atof(argv[2]);
    }
    else if (argc == 2) T = atof(argv[1]);
    else
    {
        printf ("Usage: %s [force const file basename] T\n", argv[0]);
        return (1);
    }
    phonon_force_constant_LOAD (P);
    phonon_diagonalization_REALLOC (P,D);
    for (CV=w_sample=0; ;)
    {
        V3FRANDOM (ks);
        phonon_dk (P, D, ks);
        for (i=0; i<3*P->npa; i++)
            if (D->w2[i] < 0)
            {
                printf ("** Found w2 = %e [THz^2] for\n", D->w2[i]);
                V3pr ("ks = %M, or", ks);
                M3INV (P->H, HI, a);
                /* k = 2*PI*ks*(H^-T) */
                M3mV3 (HI, ks, k);
                /* V3MuL (2*PI, k); */
                V3MuL (A_IN_M, k);
                V3pr("Cartesian k = %M [2*PI/Angstrom].", k);
            }
            else
            {
                a = HBAR_IN_J_S__R * sqrt(D->w2[i]);
                b = a / BOLZ_IN_J__K / T;
                c = exp(b);
                CV += SQUARE(b) * c / SQUARE(c-1);
                w_sample++;
                if (w_sample % CV_REPORT_FREQ == 0)
                {
                    c = CV / w_sample;
                    printf ("Averaging %d eigen-frequencies,\n", w_sample);
                    printf ("cv = %g [k_B/oscillator] = %g [J/mol/K]\n",
                            c, 3 * AVO * c * BOLZ_IN_J__K);
                    printf ("-> Debye temperature = %g [K].\n",
                            Debye_temperature(T,c));
                }
            }
    }
    return (0);
}
#undef CV_REPORT_FREQ
#endif /* _cv */


/* Diagonalize k=0 using faster routine and less memory */
void phonon_dgamma (Phonon_force_constant *P, Phonon_diagonalization *D)
{
    int i, j, k, *cn;
    register double *fc, *dp, ff;
    static double pie [DIMENSION2];
    for (cn=P->cn,fc=P->fc,j=0; j<P->npa; j++)
    {
        for (i=0; i<j; i++,cn++)
        {
            pie[0] = 0;
            pie[1] = 0;
            pie[2] = 0;
            pie[3] = 0;
            pie[4] = 0;
            pie[5] = 0;
            pie[6] = 0;
            pie[7] = 0;
            pie[8] = 0;
            ff = sqrt(P->mass[i]*P->mass[j]);
            for (k=*cn; k--;)
            {
                fc += DIMENSION;
                pie[0] += (*(fc++)) / ff;
                pie[1] += (*(fc++)) / ff;
                pie[2] += (*(fc++)) / ff;
                pie[3] += (*(fc++)) / ff;
                pie[4] += (*(fc++)) / ff;
                pie[5] += (*(fc++)) / ff;
                pie[6] += (*(fc++)) / ff;
                pie[7] += (*(fc++)) / ff;
                pie[8] += (*(fc++)) / ff;
            }
            dp = D->ap + 3*j * (3*j + 1) / 2 + 3*i;
            dp[0] = pie[0];
            dp[1] = pie[1];
            dp[2] = pie[2];
            dp += 3 * j + 1;
            dp[0] = pie[3];
            dp[1] = pie[4];
            dp[2] = pie[5];
            dp += 3 * j + 2;
            dp[0] = pie[6];
            dp[1] = pie[7];
            dp[2] = pie[8];
        }
        pie[0] = 0;
        pie[3] = 0;
        pie[4] = 0;
        pie[6] = 0;
        pie[7] = 0;
        pie[8] = 0;
        ff = P->mass[j];
        for (k=*cn; k--;)
        {
            fc += DIMENSION;
            pie[0] += (*fc) / ff;
            fc += 3;
            pie[3] += (*(fc++)) / ff;
            pie[4] += (*fc) / ff;
            fc += 2;
            pie[6] += (*(fc++)) / ff;
            pie[7] += (*(fc++)) / ff;
            pie[8] += (*(fc++)) / ff;
        }
        dp = D->ap + 3*j * (3*j + 1) / 2 + 3*i;
        dp[0] = pie[0];
        dp += 3 * j + 1;
        dp[0] = pie[3];
        dp[1] = pie[4];
        dp += 3 * j + 2;
        dp[0] = pie[6];
        dp[1] = pie[7];
        dp[2] = pie[8];
        cn++;
    }
    i = 0;
    j = 3 * P->npa;
    FORTRAN_SYMBOL(dspev)
        (D->calculate_eigenvector ? "V" : "N",  "U",
         &j, D->ap, D->w2, D->e, &j, D->rwork, &i);
    if (i < 0) pe("phonon_dgamma: %d-th argument had an illegal value.\n", -i);
    else if (i > 0) pe ("phonon_dgamma: LAPACK failed to converge (%d).\n", i);
    return;
} /* end phonon_dgamma() */


#ifdef _dgamma
int main (int argc, char *argv[])
{
    int i;
    Phonon_force_constant  P[1]={0};
    Phonon_diagonalization D[1]={0};
    if (argc == 2) P->fconst_fn_basename = argv[1];
    else if (argc != 1)
    {
        printf ("Usage: %s [force const file basename]\n", argv[0]);
        return (1);
    }
    phonon_force_constant_LOAD (P);
    phonon_diagonalization_REALLOC (P,D);
    phonon_dgamma (P, D);
    for (i=0; i<3*P->npa; i++)
    {
        printf ("%e ", D->w2[i]);
        if (D->w2[i] < 0) printf (" !\n");
        else printf (" -> %.8g THz\n", sqrt(D->w2[i])/2/PI*HZ_IN_THZ);
    }
    return (0);
}
#endif /* _dgamma */
