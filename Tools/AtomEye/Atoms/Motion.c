/*************************************************/
/* Atoms: -llapack -lblas -lm                    */
/*        -lVecMat3 -lVecMat -lScalar -lIO       */
/*                                               */
/* Physical constants, macros, and configuration */
/* operators for atomistic simulations.          */
/*                                               */
/* Dec. 12, 1999  Ju Li <liju99@mit.edu>         */
/*************************************************/

#include "Atoms.h"

/***************************************/
/* Mobility, Velocity & Kinetic Energy */
/***************************************/

/* Average drift, total kinetic energy, system T & printouts */
ConfigMotion *Config_analyze_motion (Alib_Declare_Config, FILE *out)
{
    register int i;
    static ConfigMotion m[1];
    V3 ds1, v;
    m->nmobile = 0;
    m->mobilemass = 0.;
    V3ZERO( m->s1_avg );
    for (i=(*np); i--;)
        if (MOBILE(i))
        {
            m->nmobile ++;
            m->mobilemass += (*mass)[i];
            V3ADDmuL ( (*mass)[i], &((*s1)[DIMENSION*i]), m->s1_avg );
        }
    if (m->nmobile > 0) V3DiV (m->s1_avg, m->mobilemass);
    m->kine = 0;
    for (i=(*np); i--;)
        if (MOBILE(i))
        {
            V3SUB ( &((*s1)[DIMENSION*i]), m->s1_avg, ds1 );
            V3mM3 ( ds1, H, v );
            m->kine += (*mass)[i] * V3LENGTH2(v);
        }
    if (m->nmobile > 1)
        m->T_in_K = m->kine / DIMENSION / (m->nmobile-1) *
            uenergy_IN_J / BOLZ_IN_J__K;
    else m->T_in_K = 0;
    m->kine /= 2;
    if ( out && (m->nmobile > 0) )
    {
        fprintf (out, "%d mobile atoms,\n", m->nmobile);
        V3mM3 ( m->s1_avg, H, v );
        V3fpr (out, "average velocity = %M,", v);
        fprintf (out, "total kinetic energy = %.16g = %.16g [eV]\n",
                 m->kine, m->kine * uenergy_IN_EV);
        fprintf (out, "avg. atomic kinetic energy = %.16g = %.16g [eV]\n"
                 "= %.16g [kJ/mol]\n", m->kine / m->nmobile,
                 m->kine * uenergy_IN_EV / m->nmobile,
                 m->kine / m->nmobile * uenergy_IN_KJ__MOL);
        fprintf (out, "MD temperature = %.16g = %.16g [K].\n",
                 BOLZ_IN_J__K * m->T_in_K / uenergy_IN_J, m->T_in_K);
    }
    return (m);
} /* end Config_analyze_motion() */


/* Set average mobile s1 drift to zero */
ConfigMotion *Config_motion_zero_drift (ConfigMotion *m, Alib_Declare_Config)
{
    register int i;
    if (m == NULL)
        m = Config_analyze_motion (Config_Alib_to_Alib, NULL);
    for (i=(*np); i--;)
        if (MOBILE(i))
            V3SuB ( &((*s1)[DIMENSION*i]), m->s1_avg );
    V3ZERO ( m->s1_avg );
    return (m);
} /* end Config_motion_zero_drift() */


/* Multiply the current kinetic energy and temperature by "factor" */
ConfigMotion *Config_motion_scale_T
(double factor, ConfigMotion *m, Alib_Declare_Config)
{
    register int i;
    m = Config_motion_zero_drift (m, Config_Alib_to_Alib);
    m->kine *= factor;
    m->T_in_K *= factor;
    factor = sqrt(factor);
    for (i=(*np); i--;)
        if (MOBILE(i))
            V3MuL ( factor, &((*s1)[DIMENSION*i]) );
    return (m);
} /* end Config_motion_scale_T() */


/* Assign random thermal velocities such that temperature=T_in_K */
ConfigMotion *Config_motion_initialize (double T_in_K, Alib_Declare_Config)
{
    register int i;
    ConfigMotion *m;
    V3 v;
    M3 HI;
    M3inv (H, HI);
    for (i=(*np); i--;)
        if (MOBILE(i))
        { /* assign Gaussian with <v^2> ~ 1/m */
            Vfrandnorm (DIMENSION, 0., 1./(*mass)[i], v);
            V3mM3 ( v, HI, &((*s1)[DIMENSION*i]) );
        }
    m = Config_analyze_motion (Config_Alib_to_Alib, NULL);
    m = Config_motion_scale_T (T_in_K / m->T_in_K, m, Config_Alib_to_Alib);
    return (m);
} /* end Config_motion_initialize() */


#ifdef _iniT
int main (int argc, char *argv[])
{
    Aapp_Define_Config;
    double T_in_K;
    if (argc != 4)
    {
        printf ("Purpose: assign kinetic energy to configuration.\n");
        printf ("Usage: %s in_file T_in_K out_file\n", argv[0]);
        return (1);
    }
    printf ("Loading \"%s\"...\n\n", argv[1]);
    Config_Load (argv[1], stdout, Config_Aapp_to_Alib);
    T_in_K = atof(argv[2]);
    Config_motion_initialize (T_in_K, Config_Aapp_to_Alib);
    Config_save (Config_Aapp_to_Alib, TRUE, argv[3]);
    printf ("Tmd = %g [K] -> \"%s\".\n", T_in_K, argv[3]);
    return (0);
}
#endif /* _iniT */


/* "m" & the kinetic contribution to unsymmetrized/unnormalized stress */
void Config_motion_stress (Alib_Declare_Config, ConfigMotion *m, M3 stress)
{
    register int i;
    V3 ds1, v;
    m->nmobile = 0;
    m->mobilemass = 0.;
    V3ZERO( m->s1_avg );
    for (i=(*np); i--;)
        if (MOBILE(i))
        {
            m->nmobile ++;
            m->mobilemass += (*mass)[i];
            V3ADDmuL ( (*mass)[i], &((*s1)[DIMENSION*i]), m->s1_avg );
        }
    if (m->nmobile > 0) V3DiV (m->s1_avg, m->mobilemass);
    m->kine = 0;
    SYMMAT_ZERO(stress);
    for (i=(*np); i--;)
        if (MOBILE(i))
        {
            V3SUB ( &((*s1)[DIMENSION*i]), m->s1_avg, ds1 );
            V3mM3 ( ds1, H, v );
            m->kine += (*mass)[i] * V3LENGTH2(v);
            SYMMAT_accumulate (v, (*mass)[i], v, stress);
        }
    if (m->nmobile > 1)
        m->T_in_K = m->kine / DIMENSION / (m->nmobile-1) *
            uenergy_IN_J / BOLZ_IN_J__K;
    else m->T_in_K = 0;
    m->kine /= 2;
    return;
} /* end Config_motion_stress() */


/* Compute "m", "v", and kinetic contribution to "stress" & "ep" */
void Config_motion_stress_v_ep
(Alib_Declare_Config, ConfigMotion *m, M3 stress, double *v, double *ep)
{
    register int i;
    V3 ds1;
    m->nmobile = 0;
    m->mobilemass = 0.;
    V3ZERO( m->s1_avg );
    for (i=(*np); i--;)
        if (MOBILE(i))
        {
            m->nmobile ++;
            m->mobilemass += (*mass)[i];
            V3ADDmuL ( (*mass)[i], &((*s1)[DIMENSION*i]), m->s1_avg );
        }
    if (m->nmobile > 0) V3DiV (m->s1_avg, m->mobilemass);
    m->kine = 0;
    SYMMAT_ZERO (stress);
    for (i=(*np); i--;)
        if (MOBILE(i))
        {
            V3SUB ( &((*s1)[DIMENSION*i]), m->s1_avg, ds1 );
            V3mM3 ( ds1, H, &(v[DIMENSION*i]) );
            ep[i] = (*mass)[i] * V3LENGTH2(&(v[DIMENSION*i])) / 2;
            m->kine += ep[i];
            SYMMAT_accumulate ( &(v[DIMENSION*i]), (*mass)[i],
                                &(v[DIMENSION*i]), stress );
        }
        else
        {
            V3ZERO ( &(v[DIMENSION*i]) );
            ep[i] = 0;
        }
    if (m->nmobile > 1)
        m->T_in_K = 2 * m->kine / DIMENSION / (m->nmobile-1) *
            uenergy_IN_J / BOLZ_IN_J__K;
    else m->T_in_K = 0;
    return;
} /* end Config_motion_stress_v_ep() */


/* TRUE if all atoms in the configuration are mobile */
bool Config_all_mobile ( Alib_Declare_Config )
{
    register int i;
    for (i=(*np); i--;)
        if (IMMOBILE(i)) return (FALSE);
    return (TRUE);
} /* end Config_all_mobile() */


#ifdef _label_immobile
int main (int argc, char *argv[])
{
    register int i, k;
    Aapp_Define_Config;
    if (argc != 4)
    {
        printf ("Purpose: label all immobile atoms by certain "
                "chemical symbol.\n");
        printf ("Usage: %s in_file symbol out_file\n", argv[0]);
        return (1);
    }
    printf ("Loading \"%s\"...\n\n", argv[1]);
    Config_Load (argv[1], stdout, Config_Aapp_to_Alib);
    k = 0;
    for (i=np; i--;)
        if (IMMOBILE_ATOM(i))
        {
            safe_symbol (argv[2], SYM(i));
            k++;
        }
    Config_save (Config_Aapp_to_Alib, TRUE, argv[3]);
    printf ("%d immobile atoms labelled \"%s\" -> \"%s\".\n",
            k, argv[2], argv[3]);
    return (0);
}
#endif /* _label_immobile */
