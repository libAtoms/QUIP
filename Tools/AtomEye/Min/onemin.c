/***************************************/
/* libMin: Minimization Library        */
/*         -lm -lScalar -lVecMat       */
/*                                     */
/* Dec. 6, 1999 Ju Li <liju99@mit.edu> */
/***************************************/

#include "Min.h"

/************************************************************************/
/* Minimization of single-variable f(x) on interval [a,b]: do parabolic */
/* interpolation when applicable; returns min; xmin is stored is *xmin  */
/* which is guaranteed to be within tol from the real minimum on [a,b]. */
/************************************************************************/
double onemin (double (*f)(double), double xa, double xb, double tol,
               double *xmin)
{
    double a,b,c,d,e,eps,xm,p,q,r,tol1,t2,u,v,w,fu,fv,fw, fx,x,tol3;
    if (xa > xb) SWAP(xa,xb,a); /* added by Ju Li */
    c=0.5*(3.0-sqrt(5.0)); /* squared inverse of the Golden Ratio */
    /* eps: approximately square root of relative machine precision */
    eps = SQUARE(EPS); /* added by Ju Li to enhance the precision */
    tol1=eps+1.0;
    eps=sqrt(eps);
    a=xa;
    b=xb;
    v=a+c*(b-a);
    w=v;
    x=v;
    e=0.0;
    fx=f(x);
    fv=fx;
    fw=fx;
    tol3=tol/3.0;
label_20: /* main loop starts here */
    xm=0.5*(a+b);
    tol1=eps*fabs(x)+tol3;
    t2=2.0*tol1;
    if (fabs(x-xm)<=(t2-0.5*(b-a))) goto label_190; /* stopping criterion */
    p=0.0;
    q=0.0;
    r=0.0;
    if (fabs(e)<=tol1) goto label_50;
    r=(x-w)*(fx-fv); /* fit parabola */
    q=(x-v)*(fx-fw);
    p=(x-v)*q-(x-w)*r;
    q=2.0*(q-r);
    if (q<=0.0) goto label_30;
    p=-p;
    goto label_40;
label_30:
    q=-q;
label_40:
    r=e;
    e=d;
label_50:
    if ((fabs(p)>=fabs(0.5*q*r))||(p<=q*(a-x))
	||(p>=q*(b-x))) goto label_60;
    d=p/q; /* a parabolic-interpolation step */
    u=x+d; /* f must not be evaluated too close to xa or xb */
    if (((u-a)>=t2)&&((b-u)>=t2)) goto label_90;
    d=tol1;
    if (x>=xm) d=-d;
    goto label_90;
label_60: /* a golden-section step */
    if (x>=xm) goto label_70;
    e=b-x;
    goto label_80;
label_70:
    e=a-x;
label_80:
    d=c*e;
label_90: /* f must not be evaluated too close to x */
    if (fabs(d)<tol1) goto label_100;
    u=x+d;
    goto label_120;
label_100:
    if (d<=0.0) goto label_110;
    u=x+tol1;
    goto label_120;
label_110:
    u=x-tol1;
label_120:
    fu=f(u); /* update  a, b, v, w, and x */
    if (fx>fu) goto label_140;
    if (u>=x) goto label_130;
    a=u;
    goto label_140;
label_130:
    b=u;
label_140:
    if (fu>fx) goto label_170;
    if (u>=x) goto label_150;
    b=x;
    goto label_160;
label_150:
    a=x;
label_160:
    v=w;
    fv=fw;
    w=x;
    fw=fx;
    x=u;
    fx=fu;
    goto label_20;
label_170:
    if ((fu>fw)&&(w!=x)) goto label_180;
    v=w;
    fv=fw;
    w=u;
    fw=fu;
    goto label_20;
label_180:
    if ((fu>fv)&&(v!=x)&&(v!=w)) goto label_20;
    v=u;
    fv=fu;
    goto label_20; /* end of the main loop */
label_190:
    *xmin = x;
    return (fx);
} /* end onemin() */


#ifdef _onemin_TEST
double func (double x)
{
    /* return (fabs(x-1)); */
    return (SQUARE(x)+CUBE(x)/4+QUAD(x)/17);
}

void main()
{
    double xmin, ymin;
    /* ymin = onemin (func, 11.5, 9., EPS, &xmin); */
    ymin = onemin (func, Frandom(), -Frandom(), EPS, &xmin);
    printf ("xmin = %16.13f  fmin = %16.13f\n", xmin, ymin);
    return;
}
#endif
