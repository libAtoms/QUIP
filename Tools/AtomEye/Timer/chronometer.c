/**************************************/
/* libTimer: -lm                      */
/*                                    */
/* Time and Schedule Manager          */
/*                                    */
/* Nov 19 1999 Ju Li <liju99@mit.edu> */
/**************************************/

#include "Timer.h"

/* Accurate usertime benchmark of a single task */

struct rusage chronometer_start, chronometer_stop;

#ifdef _mem_TEST
#define MEMMAX 1000000
#define MEMRPT 100
int main (int argc, char *argv[])
{
    int a[MEMMAX], b[MEMMAX];
    register int i,j;
    start_chronometer();
    for (j=0; j<MEMRPT; j++)
        for (i=0; i<MEMMAX; i++) b[i] = a[i];
    stop_chronometer();
    printf ("%f\n", stopped_usertime());
    start_chronometer();
    for (j=0; j<MEMRPT; j++)
        for (i=MEMMAX; i--;) b[i] = a[i];
    stop_chronometer();
    printf ("%f\n", stopped_usertime());
    return (0);
} /* end main() */
#endif  /* _mem_TEST */


#ifdef _chronometer_TEST
#define LOOP 2000000L
#define OP   50
int main(int argc, char *argv[])
{
    register int i;
    register double a=0., b=0., c=0., d=0., e=0., f=0.;
    start_chronometer();
    for (i=0; i<LOOP; i++)
    {
	a = b + c; /* 1 */
	d = c + f; /* 2 */
	e = b + f; /* 3 */
	a = f + b; /* 4 */
	c = b + c; /* 5 */
	e = b + e; /* 6 */
	a = b + e; /* 7 */
	d = f + c; /* 8 */
	b = e + c; /* 9 */
	c = a + d; /* 10 */
	a = b + c; /* 1 */
	d = c + f; /* 2 */
	e = b + f; /* 3 */
	a = f + b; /* 4 */
	c = b + c; /* 5 */
	e = b + e; /* 6 */
	a = b + e; /* 7 */
	d = f + c; /* 8 */
	b = e + c; /* 9 */
	c = a + d; /* 10 */
	a = b + c; /* 1 */
	d = c + f; /* 2 */
	e = b + f; /* 3 */
	a = f + b; /* 4 */
	c = b + c; /* 5 */
	e = b + e; /* 6 */
	a = b + e; /* 7 */
	d = f + c; /* 8 */
	b = e + c; /* 9 */
	c = a + d; /* 10 */
	a = b + c; /* 1 */
	d = c + f; /* 2 */
	e = b + f; /* 3 */
	a = f + b; /* 4 */
	c = b + c; /* 5 */
	e = b + e; /* 6 */
	a = b + e; /* 7 */
	d = f + c; /* 8 */
	b = e + c; /* 9 */
	c = a + d; /* 10 */
	a = b + c; /* 1 */
	d = c + f; /* 2 */
	e = b + f; /* 3 */
	a = f + b; /* 4 */
	c = b + c; /* 5 */
	e = b + e; /* 6 */
	a = b + e; /* 7 */
	d = f + c; /* 8 */
	b = e + c; /* 9 */
	c = a + d; /* 10 */
    }
    stop_chronometer();
    printf("%.1f %.1f %.1f %.1f %.1f %.1f\n", a, b, c, d, e, f);
    printf("In a %d Mflop test, user time spent = %f seconds\n",
	   (int)rint(LOOP*1E-6*OP), stopped_usertime());
    printf("this machine's peak floprate = %d Mflop/s\n",
	   (int)rint(LOOP/stopped_usertime()*OP*1E-6));
    return(0);
}
#endif

/* cc -O3 -D_chronometer_TEST chronometer.c -lm; time a.out */
