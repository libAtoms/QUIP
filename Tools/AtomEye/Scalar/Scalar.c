/*******************************************/
/* libScalar:                              */
/*                                         */
/* Scalar Arithmetic Operations and Macros */
/*                                         */
/* Nov 11 1999 Ju Li <liju99@mit.edu>      */
/*******************************************/

#include "Scalar.h"


/* greatest common divisor of a and b */
long gcd (long a, long b)
{ /* Euclid's algorithm, perhaps the first algorithm */
    long tmp;
    a = ABS(a);
    b = ABS(b);
    ENSURE(b,a,tmp);
    for (;;)
    {
	if (b == 0) return a;
	else if (b == 1) return b;
	else
	{
	    tmp = b;
	    b = a % b;
	    a = tmp;
	}
    }
} /* end gcd() */


/* find k,l such that k * x + l * y = gcd(x,y) */
void diophantine2 (int x, int y, int *k, int *l)
{
    register int a, b, u, v, p, q;
    a = abs(x);
    b = abs(y);
    (*k) = 0;
    (*l) = 1;
    u = 1;
    v = 0;
    while (1)
    {
        if (a % b == 0)
        {
            if (x < 0) (*k) = -(*k);
            if (y < 0) (*l) = -(*l);
            return;
        }
        q = a / b;
        p = b;
        b = a - q * b;
        a = p;
        p = (*k);
        (*k) = u - q * (*k);
        u = p;
        p = (*l);
        (*l) = v - q * (*l);
        v = p;
    }
    return;
} /* end diophantine2() */


#ifdef _diophantine2_TEST
int main (int argc, char *argv[])
{
    int x, y, k, l;
    printf ("x y: ");
    scanf ("%d %d", &x, &y);
    diophantine2 (x, y, &k, &l);
    printf ("(%d) x (%d) + (%d) * (%d) = (%d) = (%d)\n",
            k, x, l, y, k*x+l*y, gcd(x,y));
    return (0);
}
#endif /* _diophantine2_TEST */


/******************************************************************/
/* set Bitmap "b", starting from the bit "offset", every "period" */
/* bits to the lowest "period" bits of "value", for "n" periods.  */
/******************************************************************/
Bmap *Bperiod (Bmap *b, int offset, int n, int period, int value)
{
    /* Use type long int as chunk carrier, which has LONG_IN_BITS bits */
    /* therefore the period repeats in LCM(period,LONG_IN_BITS) bits,  */
    /* and is commensurate with array of long ints. Since period <=    */
    /* number of bits in value (int), LCM <= INT_IN_BITS*LONG_IN_BITS  */
    /* which means the long int array needs most INT_IN_BITS elements. */
    static long pattern[INT_IN_BITS];
    long *start;
    int i, j, k, m, h, longperiod, num_longperiods;
    if (period > INT_IN_BITS)
    {
	printf ("error: Bperiod: period = %d > %d bits of int,\n"
		"use BPERIOD() instead.\n", period, (int)INT_IN_BITS);
	exit(1);
    }
    i = offset;  /* bit index in b */
    j = 0;       /* bit index in value */
#define val(j) BIT(value,(j)%period)
    while ( (!ISFACTOR(LONG_IN_BITS,i)) && (i<offset+n*period) )
    { /* align to longs starting from b */
	BASSIGN (b,i,val(j));
	i++;
	j++;
    }
    if (i==offset+n*period) return(b);
    longperiod = LCM (period, LONG_IN_BITS);  /* in bits */
    for (k=0; k<longperiod; k++) BASSIGN(pattern,k,val(j+k));
    num_longperiods = ((long)n*period - j) / longperiod;
    start = ((long *)b) + i / LONG_IN_BITS;
    for (h=k=0; k<num_longperiods; k++)
	for (m=0; m<longperiod/LONG_IN_BITS; m++,h++)
	    start[h] = pattern[m];
    i += num_longperiods * longperiod;
    while ( i < offset + n*period )
    {
	BASSIGN (b,i,val(j));
	i++;
	j++;
    }
    return(b);
} /* end Bperiod() */

#ifdef _Bperiod_TEST
#define size_in_bits    256
#define value           6
#define value_period    3
#define period_repeated 20
#define offset          7
void main()
{
    Bmap *b = VBmem(size_in_bits);
    Bdump(b, 0, size_in_bits);
    Bperiod(b,offset,period_repeated,value_period,value);
    Bdump(b, 0, size_in_bits);
}
#endif  /* _Bperiod_TEST */


/* relative error of x compared to target */
double Ferr(double x, double target)
{
    if (x == target) return (0.);
    if (target == 0.) return (LARGE);
    return (x/target-1);
} /* end Ferr() */


#ifdef _Frandexp_TEST
void main()
{
    printf ("%lf\n", Frandom());
    printf ("Frandexp(1.) = %lf\n", Frandexp(1.));
}
#endif  /* _Frandexp_TEST */


/* calculates |_ a/b _| with mathematically correct floor */
int FloorDiv (int a, int b)
{
    int floor;
    if (b>0)
    {
        if (a>0)
        {
            return a /b;
        }
        else
        {
            floor = -((-a)/b);
            if ((-a)%b != 0) floor--;
        }
        return floor;
    }
    else
    {
        if (a>0)
        {
            floor = -(a/(-b));
            if (a%(-b) != 0) floor--;
            return floor;
        }
        else
        {
            return (-a)/(-b);
        }
    }
}  /* end FloorDiv() */


/* calculates |^ a/b ^| with mathamatically correct ceiling */
int CeilDiv (int a,int b)
{
    if (b>0)
        return FloorDiv(a-1,b)+1;
    else
        return FloorDiv(-a-1,-b)+1;
} /* end CeilDiv() */


/* return (n!) = GAMMA(n+1) in double precision */
double factorial (int n)
{
    register int i;
    register double ans;
    for (ans=1, i=2; i<=n; i++) ans *= i;
    return (ans);
} /* end factorial() */


/******************************************************************/
/* Take average on structures whose first element is the counter  */
/* to avoid alignment trouble please define the counter as double */
/******************************************************************/

/* set zero all memory */
void avg_clear (int n, double *s)
{
    bzero ((void *)(s), (n)*sizeof(double));
    return;
} /* end avg_clear() */

/* accumulate instances */
int avg_add (int n, double *s, ...)
{
    int i, counter;
    va_list ap;
    counter = ++(*(s++));
    va_start (ap, s);
    for (i=0; i<n-1; i++)
        (*(s++)) += va_arg (ap, double);
    va_end (ap);
    return(counter);
} /* end avg_add() */

/* do the average */
int avg_done (int n, double *s)
{
    int i, counter;
    counter = (*(s++));
    if (counter==0) return(counter);
    for (i=0; i<n-1; i++)
        (*(s++)) /= counter;
    return(counter);
} /* end avg_done() */

/* undo the effects of avg_done(): please do not use this too often */
int avg_recover (int n, double *s)
{
    int i, counter;
    counter = (*(s++));
    for (i=0; i<n-1; i++)
        (*(s++)) *= counter;
    return(counter);
} /* end avg_recover() */

#ifdef _avg_TEST
typedef struct
{
    double counter;  /* first element must be the counter */
    double a[2];
} Avg;

int main (int argc, char *argv[])
{
    Avg b;
    avg_clear ( Avg_what(Avg,b) );
    avg_add ( Avg_what(Avg,b), 1., 2. );
    avg_add ( Avg_what(Avg,b), 3., 4. );
    avg_clear ( Avg_what(Avg,b) );
    avg_add ( Avg_what(Avg,b), 11., 0. );
    avg_clear ( Avg_what(Avg,b) );
    avg_done ( Avg_what(Avg,b) );
    printf ("%f %f %f\n",  b.counter, b.a[0], b.a[1]);
    return (0);
}
#endif /* _avg_TEST */


/* Generate a string of n-1 random base64 characters ended by EOS. */
char *RandomBase64String (int n, char *a)
{
    register int i;
    for (i=0; i<n-1; i++) a[i] = Base64Alphabet[fran(0,64)];
    a[i] = EOS;
    return(a);
} /* end RandomBase64String() */

#ifdef _RandomBase64String_TEST
#define SIZE 32
int main (int argc, char *argv[])
{
    char str[SIZE];
    printf ("%s\n", RandomBase64String(SIZE,str));
    return (0);
}
#undef SIZE
#endif /* _RandomBase64String_TEST */


/* (-1)%3=-1, but positive_remainder(-1,3)=2 */
int positive_remainder(int a, int b)
{
    register int c = a % b;
    if (c < 0) c += b;
    return (c);
} /* end positive_remainder() */
