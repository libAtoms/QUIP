/*********************************************/
/* libVecMat:                                */
/*                                           */
/* General Vector & Matrix Operation Library */
/*                                           */
/* Nov 19 1999 Ju Li <liju99@mit.edu>        */
/*********************************************/

#include "VecMat.h"


/* double type vector operations */

/* x[0:n-1] := Frandom() [uniform on (0,1)]; then return x[] */
double *Vfrandom (int n, double x[])
{
    int i;
    for (i=n; i--;) x[i] = Frandom();
    return (x);
} /* end Vfrandom() */


/* x[0:n-1] := FRANDOM() [uniform on (-0.5,0.5)]; then return x[] */
double *VFRANDOM (int n, double x[])
{
    int i;
    for (i=n; i--;) x[i] = FRANDOM();
    return (x);
} /* end VFRANDOM() */


/* x[0:n-1] := Frandnorm(E,sigma2); then return x[] */
double *Vfrandnorm (int n, double E, double sigma2, double x[])
{
    int i;
    double s1, s2, sigma;
    sigma = sqrt(sigma2);
    for (i=0; i<n; i+=2)
    {
	s1 = sqrt(-2.*log(Frandom()));
	s2 = 2*acos(-1.)*Frandom();
	x[i] = s1*sigma*cos(s2) + E;
	if (i+1 < n) x[i+1] = s1*sigma*sin(s2) + E;
    }
    return (x);
} /* end Vfrandnorm() */


/* x[0:n-1] := Frandnorm(0,1); then return x[] */
double *VFRANDNORM (int n, double x[])
{
    int i;
    double s1, s2;
    for (i=0; i<n; i+=2)
    {
	s1 = sqrt(-2.*log(Frandom()));
	s2 = 2*acos(-1.)*Frandom();
	x[i] = s1*cos(s2);
	if (i+1 < n) x[i+1] = s1*sin(s2);
    }
    return (x);
} /* end VFRANDNORM() */


/* return the Euclidean norm squared := |x|^2 */
double Vlength2 (int n, double x[])
{
    int i;
    register double sum;
    sum = 0;
    for (i=n; i--;) sum += SQUARE(x[i]);
    return (sum);
} /* end Vlength2() */


/* return the Euclidean distance squared := |a-b|^2 */
double Vdistance2 (int n, double a[], double b[])
{
    int i;
    register double d, sum;
    sum = 0;
    for (i=n; i--;)
    {
        d = a[i] - b[i];
        sum += SQUARE(d);
    }
    return (sum);
} /* end Vdistance2() */


/* normalize x[] to unit length; then return x[] */
double *Vnormalize (int n, double x[])
{
    int i;
    register double sum;
    sum = 0;
    for (i=n; i--;) sum += SQUARE(x[i]);
    sum = sqrt(sum);
    if (sum > 0) 
        for (i=n; i--;)
            x[i] /= sum;
    return (x);
} /* end Vnormalize() */


/* b[] := -a[]; then return b[] */
double *Vneg (int n, double *a, double *b)
{
    register int i;
    for (i=n; i--;) b[i] = -a[i];
    return(b);
} /* end Vneg() */


/* a[] := -a[]; then return a[] */
double *VNeg (int n, double *a)
{
    register int i;
    for (i=n; i--;) a[i] = -a[i];
    return(a);
} /* end VNeg() */


/* make b[] an image of a[] in [0,1)^3; then return b[] */
double *Vtrim (int n, double a[], double b[])
{
    int i;
    for (i=n; i--;)
	b[i] = TRIM(a[i]);
    return(b);
} /* end Vtrim() */


/* change a[] to its own image in [0,1)^3; then return a[] */
double *VTRIM (int n, double a[])
{
    int i;
    for (i=n; i--;)
	a[i] = TRIM(a[i]);
    return(a);
} /* end VTRIM() */


/* make b[] image of a[]'s in [-0.5,0.5)^3; then return b[] */
double *Vimage (int n, double a[], double b[])
{
    int i;
    for (i=n; i--;)
	b[i] = IMAGE(a[i]);
    return(b);
} /* end Vimage() */


/* change a[] to its own image in [-0.5,0.5)^3; then return a[] */
double *VIMAGE (int n, double a[])
{
    int i;
    for (i=n; i--;)
	a[i] = IMAGE(a[i]);
    return(a);
} /* end VIMAGE() */


/* returns max(a[]) */
double VMAX (int n, double a[])
{
    int i;
    double vmax;
    for (vmax=a[0], i=1; i<n; i++)
	if (a[i]>vmax) vmax = a[i];
    return(vmax);
} /* end VMAX() */


/* returns index imax, such that a[imax] >= a[] */
int Vmax (int n, double a[])
{
    int i, imax;
    double vmax;
    for (vmax=a[0], imax=0, i=1; i<n; i++)
	if (a[i] > vmax)
	{
	    imax = i;
	    vmax = a[i];
	}
    return(imax);
} /* end Vmax() */


/* *vmax := max(a[]); and returns the index imax, s.t. a[imax] = max(a[]) */
int VMax (int n, double a[], double *vmax)
{
    int i, imax;
    for (*vmax=a[0], imax=0, i=1; i<n; i++)
	if (a[i] > *vmax)
	{
	    imax = i;
	    *vmax = a[i];
	}
    return(imax);
} /* end VMax() */


/* returns min(a[]) */
/* double VMIN (int n, double a[]) */
/* { */
/* int i; */
/* double vmin; */
/* for (vmin=a[0], i=1; i<n; i++) */
/* if (a[i]<vmin) vmin = a[i]; */
/* return(vmin); */
/* } */
/* end VMIN() */


/* returns index imin, such that a[imin] <= a[] */
int Vmin (int n, double a[])
{
    int i, imin;
    double vmin;
    for (vmin=a[0], imin=0, i=1; i<n; i++)
	if (a[i] < vmin)
	{
	    imin = i;
	    vmin = a[i];
	}
    return(imin);
} /* end Vmin() */


/* *vmin := min(a[]); and returns the index imin, s.t. a[imin] = min(a[]) */
int VMin (int n, double a[], double *vmin)
{
    int i, imin;
    for (*vmin=a[0], imin=0, i=1; i<n; i++)
	if (a[i] < *vmin)
	{
	    imin = i;
	    *vmin = a[i];
	}
    return(imin);
} /* end VMin() */


/* summing all elements in a[] */
double Vsum (int n, double a[])
{
    register int i;
    register double vsum=0;
    for (i=n; i--;) vsum += a[i];
    return(vsum);
} /* end Vsum() */


/* b[] := multiplier * a[]; then return b[] */
double *Vmul (int n, double multiplier, double a[], double b[])
{
    int i;
    for (i=n; i--;) b[i] = multiplier * a[i];
    return(b);
} /* end Vmul() */


/* a[] := multiplier * a[]; then return a[] */
double *VMUL (int n, double multiplier, double a[])
{
    int i;
    for (i=n; i--;) a[i] *= multiplier;
    return(a);
} /* end VMUL() */


/* b[] := a[] / divisor; then return b[] */
double *Vdiv (int n, double a[], double divisor, double b[])
{
    int i;
    if (divisor == 0.)
    {
	printf ("error: Vdiv: divisor = %e\n", divisor);
	exit(1);
    }
    for (i=n; i--;) b[i] *= a[i] / divisor;
    return(b);
} /* end Vdiv() */


/* a[] := a[] / divisor; then return a[] */
double *VDIV (int n, double a[], double divisor)
{
    int i;
    if (divisor == 0.)
    {
	printf ("error: VDIV: divisor = %e\n", divisor);
	exit(1);
    }
    for (i=n; i--;) a[i] /= divisor;
    return(a);
} /* end VDIV() */


/* c[] := a[] + b[]; then return c[] */
double *Vadd (int n, double a[], double b[], double c[])
{
    int i;
    for (i=n; i--;) c[i] = a[i] + b[i];
    return (c);
} /* end Vadd() */


/* b[] := a[] + b[]; then return b[] */
double *VADD (int n, double a[], double b[])
{
    int i;
    for (i=n; i--;) b[i] += a[i];
    return (b);
} /* end VADD() */


/* c[] := a[] - b[]; then return c[] */
double *Vsub (int n, double a[], double b[], double c[])
{
    int i;
    for (i=n; i--;) c[i] = a[i] - b[i];
    return (c);
} /* end Vsub() */


/* a[] := a[] - b[]; then return a[] */
double *VSUB (int n, double a[], double b[])
{
    int i;
    for (i=n; i--;) a[i] -= b[i];
    return (a);
} /* end VSUB() */


/* dot product of a[] and b[] */
double Vdot (int n, double a[], double b[])
{
    int i;
    double product = 0;
    for (i=n; i--;) product += a[i]*b[i];
    return(product);
} /* end Vdot() */


/* VECTOR & (SCALAR x VECTOR) OPERATIONS */

/* c[] := a[] + multiplier * b[]; then return c[] */
double *Vaddmul (int n, double a[], double multiplier, double b[], double c[])
{
    int i;
    for (i=n; i--;) 
	c[i] = a[i] + multiplier * b[i];
    return (c);
} /* end Vaddmul() */


/* c[] := aa * a[] + bb * b[]; then return c[] */
double *Vaddmulmul
(int n, double aa, double a[], double bb, double b[], double c[])
{
    register int i;
    for (i=n; i--;) c[i] = aa * a[i] + bb * b[i];
    return (c);
} /* end Vaddmulmul() */


/* c[] := a[] + b[] / divisor; then return c[] */
double *Vadddiv (int n, double a[], double b[], double divisor, double c[])
{
    int i;
    if (divisor == 0.)
    {
	printf ("error: Vadddiv: divisor = %e\n", divisor);
	exit(1);
    }
    for (i=n; i--;) c[i] = a[i] + b[i] / divisor;
    return (c);
} /* end Vadddiv() */


/* b[] := b[] + multiplier * a[]; then return b[] */
double *VADDmul (int n, double multiplier, double a[], double b[])
{
    int i;
    for (i=n; i--;) b[i] += multiplier * a[i];
    return (b);
} /* end VADDmul() */


/* b[] := b[] + a[] / divisor; then return b[] */
double *VADDdiv (int n, double a[], double divisor, double b[])
{
    int i;
    if (divisor == 0.)
    {
	printf ("error: VADDdiv: divisor = %e\n", divisor);
	exit(1);
    }
    for (i=n; i--;) b[i] += a[i] / divisor;
    return (b);
} /* end VADDdiv() */


/* c[] := a[] - multiplier * b[]; then return c[] */
double *Vsubmul (int n, double a[], double multiplier, double b[], double c[])
{
    int i;
    for (i=n; i--;) c[i] = a[i] - multiplier * b[i];
    return (c);
} /* end Vsubmul() */


/* c[] := a[] - b[] / divisor; then return c[] */
double *Vsubdiv (int n, double a[], double b[], double divisor, double c[])
{
    int i;
    if (divisor == 0.)
    {
	printf ("error: Vsubdiv: divisor = %e\n", divisor);
	exit(1);
    }
    for (i=n; i--;) c[i] = a[i] - b[i] / divisor;
    return (c);
} /* end Vsubdiv() */


/* a[] := a[] - multiplier * b[]; then return a[] */
double *VSUBmul (int n, double a[], double multiplier, double b[])
{
    int i;
    for (i=n; i--;) a[i] -= multiplier * b[i];
    return (a);
} /* end VSUBmul() */


/* a[] := a[] - b[] / divisor; then return a[] */
double *VSUBdiv (int n, double a[], double b[], double divisor)
{
    int i;
    if (divisor == 0.)
    {
	printf ("error: VSUBdiv: divisor = %e\n", divisor);
	exit(1);
    }
    for (i=n; i--;) a[i] -= b[i] / divisor;
    return (a);
} /* end VSUBdiv() */

