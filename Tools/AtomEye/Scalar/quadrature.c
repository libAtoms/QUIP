/*******************************************/
/* libScalar:                              */
/*                                         */
/* Scalar Arithmetic Operations and Macros */
/*                                         */
/* Nov 11 1999 Ju Li <liju99@mit.edu>      */
/*******************************************/

#include "Scalar.h"


/* Integration of one-dimensional function: */

#define FUNC(x) ((*func)(x))
static double trapzd (double (*func)(double), double a, double b, int n)
{
    double x,tnm,sum,del;
    static double s;
    int it,j;
    
    if (n == 1)
    {
        return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
    }
    else
    {
        for (it=1,j=1;j<n-1;j++) it <<= 1;
        tnm=it;
        del=(b-a)/tnm;
        x=a+0.5*del;
        for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
        s=0.5*(s+(b-a)*sum/tnm);
        return s;
    }
} /* end trapzd() */
#undef FUNC


/* Gaussian quadrature with W(x)=1, N=10, NR pp.148 */
double qgaus (double (*func)(double), double a, double b)
{
    int j;
    double xr,xm,dx,s;
    static double x[]={0.0,0.1488743389,0.4333953941,
                      0.6794095682,0.8650633666,0.9739065285};
    static double w[]={0.0,0.2955242247,0.2692667193,
                      0.2190863625,0.1494513491,0.0666713443};
    xm=0.5*(b+a);
    xr=0.5*(b-a);
    s=0;
    for (j=1;j<=5;j++)
    {
        dx=xr*x[j];
        s += w[j]*((*func)(xm+dx)+(*func)(xm-dx));
    }
    return s *= xr;
} /* end qgaus() */


/* Romberg Quadrature: */

double qromb_tolerance = QROMB_DEFAULT_TOLERANCE;

static void qromb_polint
(double xa[], double ya[], int n, double x, double *y, double *dy)
{
    int i,m,ns=1;
    double den,dif,dift,ho,hp,w;
    double c[QROMB_K+1],d[QROMB_K+1];
    dif=fabs(x-xa[1]);
    for (i=1;i<=n;i++)
    {
        if ( (dift=fabs(x-xa[i])) < dif)
        {
            ns=i;
            dif=dift;
        }
        c[i]=ya[i];
        d[i]=ya[i];
    }
    *y=ya[ns--];
    for (m=1;m<n;m++)
    {
        for (i=1;i<=n-m;i++)
        {
            ho=xa[i]-x;
            hp=xa[i+m]-x;
            w=c[i+1]-d[i];
            if ( (den=ho-hp) == 0.0)
                pe("qromb_polint: error in routine qromb_polint.\n");
            den=w/den;
            d[i]=hp*den;
            c[i]=ho*den;
        }
        *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
    }
    return;
} /* end qromb_polint() */


double qromb (double (*func)(double), double a, double b)
{
    double ss, dss;
    double s[QROMB_JMAX+1], h[QROMB_JMAX+1+1];
    int j;
    h[1] = 1.0;
    for (j=1; j<=QROMB_JMAX; j++)
    {
        s[j] = trapzd (func, a, b, j);
        if (j >= QROMB_K)
        {
            qromb_polint
                (&h[j-QROMB_K], &s[j-QROMB_K], QROMB_K, 0.0, &ss, &dss);
            if (fabs(dss) <= qromb_tolerance*fabs(ss)) return(ss);
        }
        /* s[j+1] = s[j]; */
        h[j+1] = 0.25 * h[j];
    }
    fprintf (stderr, "qromb: Too many steps in routine qromb.\n");
    /* pe("qromb: Too many steps in routine qromb.\n"); */
    return (0);
} /* end qromb() */


#ifdef _quadrature_contest_TEST
int main (int argc, char *argv[])
{
    register int i,j;
    double value;
    struct List
    {
        char *name;
        double (*quadrature) (double (*func)(double), double a, double b);
    } funs[] =
      {{"qgaus", &qgaus},
       {"qromb", &qromb},};
    for (i=0; i<sizeof(funs)/sizeof(struct List); i++)
    {
        printf ("Testing %s...", funs[i].name);
        value = funs[i].quadrature (sin, 0, PI/2);
        printf (" value=%.15e\n", value);
    }
    return (0);
}
#endif /* _quadrature_contest_TEST */
