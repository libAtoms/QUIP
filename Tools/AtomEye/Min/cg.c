/***************************************/
/* libMin: Minimization Library        */
/*         -lm -lScalar -lVecMat       */
/*                                     */
/* Dec. 6, 1999 Ju Li <liju99@mit.edu> */
/***************************************/

#include "Min.h"

/*********************************/
/* Conjugate Gradient Minimizers */
/*********************************/

CG_Option CG_option =
{
    {
        MACOPT_END_IF_SMALL_GRADIENT, /* if TRUE use small-step criterion */
        TRUE, /* if TRUE the program runs a bit slower but more stable */
        20,    /* in maclinmin: Probably not worth fiddling with */
        0.01,  /* typical step length */
        2.0,   /* factors for growing and shrinking the interval */
        1.25,  /* don't fiddle unless you really mean it */
        0.5    /* don't fiddle unless you really mean it */
    },
    HUGE_INT,  /* do not stop at maximum number of evaluations */
    EPS,       /* gradient_square_sum_tolerance */
    EPS,       /* step_square_sum_tolerance */
    0.01       /* rough_potential_decrease_expectation */
};

CG_State CG_state;


static void macopt_restart (int start) 
{
    register int j;
    register double *g, *h, *xi;
    g = CG_state.macopt.g;
    h = CG_state.macopt.h;
    xi = CG_state.macopt.xi;
    if (start == 0)
        CG_state.macopt.lastx = CG_option.macopt.lastx_default;
    for (j=CG_state.macopt.N; j--;)
    {
        if (CG_state.macopt.restart != 2) g[j] = -xi[j];
        xi[j] = h[j] = g[j];
    }
    CG_state.macopt.restart = 0;
    return;
} /* end macopt_restart() */


static double macprod (double *p , double *gy , double y)
{
    register double *pt = CG_state.macopt.pt;
    register double *xi = CG_state.macopt.xi;
    register int i;
    register double s = 0;
    for (i=CG_state.macopt.N; i--;)
        pt[i] = p[i] + y * xi[i];
    CG_state.potential_min = CG_state.potential (pt, gy);
    CG_state.number_of_evaluations++;
    for (i=CG_state.macopt.N; i--;)
        s += gy[i] * xi[i];
    return (s);
} /* end macprod() */


static double maclinmin (double *p)
{
    double x, y;
    double s, t, m;
    register int its = 0, i;
    double step, tmpd; 
    register double *gx, *gy;
    gx = CG_state.macopt.gx;
    gy = CG_state.macopt.gy;
    x = CG_state.macopt.lastx / CG_state.macopt.gtyp;
    s = macprod(p, gx, x);
    if ( s < 0 )
    {  /* we need to go further */
        do
        {
            y = x * CG_option.macopt.linmin_g1;
            t = macprod (p, gy, y);
            if ( t >= 0.0 ) break;
            x = y;
            s = t;
            CG_state.macopt.gunused = gx;
            gx = gy;
            gy = CG_state.macopt.gunused; 
            its++;
        }
        while (its < CG_option.macopt.linmin_maxits);
    }
    else if ( s > 0 )
    { /* need to step back inside interval */
        do
        {
            y = x * CG_option.macopt.linmin_g3;
            t = macprod (p, gy, y); 
            if ( t <= 0.0 ) break;
            x = y;
            s = t;
            CG_state.macopt.gunused = gx;
            gx = gy;
            gy = CG_state.macopt.gunused;
            its ++;
        } while (its < CG_option.macopt.linmin_maxits);
    }
    else
    { /* hole in one s = 0.0 */
        t = 1.0;
        y = x;
    }
    if (its > CG_option.macopt.linmin_maxits)
    {
        fprintf (stderr, "Warning! maclinmin overran" );
        fprintf (stderr, "- inner product at 0 = %9.4g\n",
                 tmpd = macprod (p, gy, 0.0) ); 
        if ( (tmpd > 0) && (!CG_option.macopt.rich) )
        {
            fprintf (stderr, "setting rich to 1\n" );
            CG_option.macopt.rich = TRUE;
        }
        if ( tmpd > 0 ) CG_state.macopt.restart = 1; 
    }
    if ( s < 0.0 ) s = - s;
    if ( t < 0.0 ) t = - t;
    m = ( s + t );
    s /= m;
    t /= m;
    m =  s * y + t * x;
    for (step=0,i=CG_state.macopt.N; i--;)
    {
        tmpd = m * CG_state.macopt.xi[i];
        p[i] += tmpd;
        step += tmpd * tmpd;
        CG_state.macopt.xi[i] = s * gy[i] + t * gx[i];
    }
    CG_state.macopt.lastx = m * CG_option.macopt.linmin_g2 *
        CG_state.macopt.gtyp;
    return (step); 
} /* end maclinmin() */


/* http://131.111.48.24/mackay/c/macopt.html; mackay@mrao.cam.ac.uk */
double macopt (double (*potential) (double *x, double *gradient),
               int N, double *p)
{
    register int j;
    register double *g, *h, *xi;
    double gg, gam, dgg;
    double step, tmpd;

    CG_state.macopt.N = N;
    CG_state.potential = potential;
    CG_state.macopt.lastx = CG_option.macopt.lastx_default;
    CG_state.error = FALSE;

    MALLOC( macopt, CG_state.macopt.g, 6*N, double );
    g = CG_state.macopt.g;
    h = CG_state.macopt.h = g + N;
    xi = CG_state.macopt.xi = h + N;
    CG_state.macopt.pt = xi + N;
    CG_state.macopt.gx = CG_state.macopt.pt + N;
    CG_state.macopt.gy = CG_state.macopt.gx + N;

    CG_state.potential_min = potential(p, xi);
    CG_state.number_of_evaluations = 1;
    macopt_restart(1);

    while (CG_state.number_of_evaluations <
           CG_option.max_potential_evaluations)
    {
        VLENGTH2 (N, g, j, gg);
        CG_state.macopt.gtyp = sqrt ( gg / (double)(N) );
        if ( (CG_option.macopt.convergence_condition ==
              MACOPT_END_IF_SMALL_GRADIENT) &&
             (gg < CG_option.gradient_square_sum_tolerance))
            goto exit;
        step = maclinmin(p);
        if (CG_state.macopt.restart == 0)
            if ((CG_option.macopt.convergence_condition ==
                 MACOPT_END_IF_SMALL_STEP) &&
                (step < CG_option.step_square_sum_tolerance))
                goto exit;
        if (CG_option.macopt.rich || CG_state.macopt.restart)
        {
            CG_state.potential_min = potential(p, xi);
            CG_state.number_of_evaluations++;
        }
        if ( CG_state.macopt.restart )
        {
            fprintf(stderr, "Restarting macopt\n");
            macopt_restart(0);
        }
        else
        {
            dgg = 0.0;
            for (j=N; j--;) dgg += ( xi[j] + g[j] ) * xi[j];
            gam = dgg / gg;
            for ( tmpd = 0.0, j=N; j--; )
            {
                g[j] = -xi[j];
                xi[j] = h[j] = g[j] + gam * h[j];
                tmpd -= xi[j] * g[j];
            }
            if ( tmpd > 0.0 )
                fprintf (stderr,"new line search has inner prod %9.4g\n",
                         tmpd ); 
            if ( tmpd > 0.0 )
            {
                if (!CG_option.macopt.rich)
                {
                    fprintf (stderr, "Setting rich to TRUE; " ); 
                    CG_option.macopt.rich = TRUE; 
                }
                CG_state.macopt.restart = 2;
                fprintf(stderr,"Restarting macopt (2)\n" );
                macopt_restart (0);
            }
        }
    }
    CG_state.error = TRUE;
    CG_state.error_mesg =
        "macopt: iteration is terminated because potential evaluation\n"
        "limit CG_option.max_potential_evaluations is reached.\n";
    fprintf(stderr,"Reached iteration limit in macopt; continuing.\n");
  exit:
    free (CG_state.macopt.g);
    return (CG_state.potential_min);
} /* end macopt() */


static void zxcgr_dummy (int *N, double *x, double *value, double *gradient)
{
    *value = CG_state.potential(x,gradient);
    CG_state.number_of_evaluations++;
    return;
} /* end zxcgr_dummy() */

/* Driver of IMSL conjugate gradient minimizer zxcgr_() */
double zxcgr (double (*potential) (double *x, double *gradient),
              int N, double *x)
{
    CG_state.potential = potential;
    MALLOC (zxcgr, CG_state.gradient, 7*N, double);
    CG_state.number_of_evaluations = 0;
    FORTRAN_SYMBOL(zxcgrf77)
        (zxcgr_dummy,
         &N,
         &CG_option.gradient_square_sum_tolerance,
         &CG_option.max_potential_evaluations,
         &CG_option.rough_potential_decrease_expectation,
         x,
         CG_state.gradient,
         &CG_state.potential_min,
         CG_state.gradient+N,
         &CG_state.error);
    free (CG_state.gradient);
    switch (CG_state.error)
    {
        case 0:  /* success */
            break;
        case 129:
            CG_state.error_mesg =
                "zxcgr: LINE SEARCH OF AN INTEGRATION WAS ABANDONED.\n"
                "THIS ERROR MAY BE CAUSED BY AN ERROR IN THE GRADIENT.\n";
            break;
        case 130:
            CG_state.error_mesg =
                "zxcgr: CALCULATION CANNOT CONTINUE BECAUSE THE\n"
                "SEARCH DIRECTION IS UPHILL.\n";
            break;
        case 131:
            CG_state.error_mesg =
                "zxcgr: ITERATION WAS TERMINATED BECAUSE POTENTIAL EVALUATION"
                "\nLIMIT CG_option.max_potential_evaluations WAS EXCEEDED.\n";
            break;
        case 132:
            CG_state.error_mesg =
                "zxcgr: CALCULATION WAS TERMINATED BECAUSE TWO\n"
                "CONSECUTIVE ITERATIONS FAILED TO REDUCE F.\n";
            break;
        default:
            break;
    }
    return(CG_state.potential_min);
} /* end zxcgr() */


#ifdef _cg_contest_TEST
#define N 20000
double xanswer[N];
double test_potential (double *x, double *gradient)
{
    int i;
    double value;
    for (value=i=0; i<N; i++)
    {
        value += SQUARE(x[i] - xanswer[i]) / (i+1);
        gradient[i] = 2 * (x[i] - xanswer[i]) / (i+1);
    }
    return (value);
} /* end test_potential() */

int main (int argc, char *argv[])
{
    register int i,j;
    double xstart[N], xmin[N], value, gradient[N], g2, xdiff2;
    struct List
    {
        char *name;
        double (*cg) (double (*potential) (double *x, double *gradient),
                      int n, double *x);
    } funs[] =
      {{"macopt", &macopt},
       {"IMSL zxcgr", &zxcgr},};
    Vfrandom (N, xanswer);
    Vfrandom (N, xstart);
    for (i=0; i<sizeof(funs)/sizeof(struct List); i++)
    {
        VEQV(N, xstart, xmin);
        printf ("Testing %s...", funs[i].name);
        value = funs[i].cg (test_potential, N, xmin);
        if (CG_state.error) pe(CG_state.error_mesg);
        printf (" %d evaluations\n", CG_state.number_of_evaluations);
        if (value!=test_potential(xmin,gradient))
            printf("warning: value incongruence of %e vs %e\n",
                   value, test_potential(xmin,gradient));
        VLENGTH2 (N, gradient, j, g2);
        VDIFF2 (N, xmin, xanswer, j, xdiff2);
        printf ("|g|^2=%e, |Dx|^2=%e, e=%e\n", g2, xdiff2, value);
    }
    return (0);
}
#endif /* _cg_contest_TEST */
