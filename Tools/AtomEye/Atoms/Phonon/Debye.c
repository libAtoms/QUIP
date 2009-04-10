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

/**************************************************/
/* Debye DOS and average energy model. See Notes. */
/**************************************************/

static double Debye_function_integrand (double y)
{
    register double a, b;
    if (y <= 0.) return (0.);
    a = exp(y);
    b = a - 1;
    return ( QUAD(y) * a / SQUARE(b) );
} /* end Debye_function_integrand() */


/* D(x) = 3 * x^3 * \int_0^{1/x}dy y^4 * e^y / (e^y-1)^2 */
double Debye_function (double x)
{
    if (x <= 0.) return(0.);
    return( 3 * CUBE(x) * qromb ( &Debye_function_integrand, 0., 1./x ) );
} /* end Debye_function() */

#ifdef _Debye_function_TEST
int main (int argc, char *argv[])
{
    printf ("%e\n", Debye_function(1.) );
    return (0);
}
#endif /* _Debye_function_TEST */


static double Debye_function_derivative (double x)
{
    if (x <= 0.) return(0.);
    return( 9 * SQUARE(x) * qromb ( &Debye_function_integrand, 0., 1./x ) -
            3 * x * Debye_function_integrand(1./x) );
} /* end Debye_function_derivative() */


static double Debye_inverse_function (double D)
{
    double x, a, b, newx;
    if (D <= 0.) return(0.);
    else if (D >= 1.) return(SINGLE_PRECISION_INFINITY);
    for (x=1.; ; x=newx)
    {
        a = Debye_function (x);
        if (fabs(a-D) < SMALL) break;
        b = Debye_function_derivative (x);
        newx = x - (a - D) / b;
        if (newx <= 0) newx = x / 2;
    }
    return(x);
} /* end Debye_inverse_function() */


/*************************************************************/
/* Given Treal and corresponding dTmd__dTreal (heat capacity */
/* normalized to 1),  deduce the Debye temperature.          */
/*************************************************************/
double Debye_temperature (double Treal, double dTmd__dTreal)
{
    double Treal__TD;
    Treal__TD = Debye_inverse_function (dTmd__dTreal);
    if (Treal__TD == 0) return(SINGLE_PRECISION_INFINITY);
    return( Treal / Treal__TD );
} /* end Debye_function() */


static double Debye_Tmd_integrand (double y)
{
    if (y <= 0.) return (0.);
    return ( CUBE(y) * (0.5 + 1./(exp(y)-1)) );
} /* end Debye_Tmd_integrand() */

/**************************************************************/
/* Given Debye temperature and its associated DOS which leads */
/* to certain quantum average energy at T, find Tmd whereby   */
/* a classical system would have the same average energy.     */
/**************************************************************/
double Debye_Tmd (double TD, double T)
{
    register double TD__T;
    if (T<=0) return (3./8.*TD);
    TD__T = TD / T;
    return( 3 * TD / QUAD(TD__T) * qromb(&Debye_Tmd_integrand,0.,TD__T) );
} /* end Debye_Tmd() */


/* TD <= 0 is the signal for no temperature rescaling */
double Debye_Treal_to_Tmd (double TD, double Treal, double *dTmd__dTreal)
{
    double Tmd;
    if (TD <= 0)
    {  /* signal for no temperature rescaling */
        Tmd = Treal;
        *dTmd__dTreal = 1.0;
        return (Tmd);
    }
    Tmd = Debye_Tmd (TD, Treal);
    *dTmd__dTreal = Debye_function (Treal / TD);
    return (Tmd);
} /* end Debye_Treal_to_Tmd() */


/* TD <= 0 is the signal for no temperature rescaling */
double Debye_Tmd_to_Treal (double TD, double Tmd, double *dTmd__dTreal)
{
    double Treal, a;
    if (TD <= 0)
    { /* signal for no temperature rescaling */
        Treal = Tmd;
        *dTmd__dTreal = 1.0;
        return (Treal);
    }
    if (Tmd <= 3./8.*TD)
    {
        Treal = 0;
        *dTmd__dTreal = 0;
        return (Treal);
    }
    /* Newton's method: we like super-linear convergence */
    for (Treal=Tmd; ;Treal-=(a-Tmd)/(*dTmd__dTreal))
        if (fabs((a=Debye_Treal_to_Tmd(TD,Treal,dTmd__dTreal))-Tmd) < SMALL)
            break;
    return (Treal);
} /* end Debye_Tmd_to_Treal() */


#ifdef _Debye_TEST
int main (int argc, char *argv[])
{
    double TD=1000, Treal=77, Tmd, dTmd__dTreal[1];
    Tmd = Debye_Treal_to_Tmd (TD, Treal, dTmd__dTreal);
    printf ("TD=%g Treal=%g Tmd=%g dTmd__dTreal=%g\n",
            TD, Treal, Tmd, dTmd__dTreal[0]);
    Treal = Debye_Tmd_to_Treal (TD, Tmd, dTmd__dTreal);
    printf ("TD=%g Treal=%g Tmd=%g dTmd__dTreal=%g\n",
            TD, Treal, Tmd, dTmd__dTreal[0]);
    TD = Debye_temperature (Treal, dTmd__dTreal[0]);
    printf ("TD=%g Treal=%g Tmd=%g dTmd__dTreal=%g\n",
            TD, Treal, Tmd, dTmd__dTreal[0]);
    return (0);
}
#endif /* _Debye_TEST */
