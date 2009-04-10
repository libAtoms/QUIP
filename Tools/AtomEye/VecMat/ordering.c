/*********************************************/
/* libVecMat:                                */
/*                                           */
/* General Vector & Matrix Operation Library */
/*                                           */
/* Nov 19 1999 Ju Li <liju99@mit.edu>        */
/*********************************************/

#include "VecMat.h"

/********************************************/
/* elementary order, swapping & permutation */
/********************************************/

/* assign index array idx[] := a..b (1..10, 10..1); return idx[] */
int *sequentially_index (int idx[], int a, int b)
{
    register int i=0, c=a;
    if (a <= b)
    {
      loop1:
        REPEAT256( idx[i++]=(c++); if (c>b) goto exit; );
        goto loop1;
    }
    else
    {
      loop2:
        REPEAT256( idx[i++]=(c--); if (c<b) goto exit; );
        goto loop2;
    }
  exit:
    return (idx);
} /* end sequentially_index() */


/* assign index array idx[] := a:b:c; return idx[] */
int *sequentially_Index (int idx[], int a, int b, int c)
{
    register int i=0, d=a;
    if (a <= c)
    {
      loop1:
        REPEAT256( idx[i++]=d; d+=b; if (d>c) goto exit; );
        goto loop1;
    }
    else
    {
      loop2:
        REPEAT256( idx[i++]=d; d-=b; if (d<c) goto exit; );
        goto loop2;
    }
  exit:
    return (idx);
} /* end sequentially_Index() */


/* assign index array idx[] := 0..N-1 */
void Sequentially_index (int N, int idx[])
{
    register int i=0;
    if (i>=N) goto exit; idx[i] = i; i++;
  loop:
    REPEAT256(if (i==N) goto exit; idx[i] = i; i++;);
    goto loop;
  exit:
    return;
} /* end Sequentially_index() */


/* assign index array idx[] := a, a+b, a+2*b, .., a+(N-1)*b; a=0..b-1 */
void Sequentially_Index (int N, int a, int b, int idx[])
{
    register int i=0, j=a;
    if (a>=b) pe("Sequentially_Index: a=%d >= b=%d.\n", a,b);
    if (i>=N) goto exit; idx[i++] = j; j+=b;
  loop:
    REPEAT256(if (i==N) goto exit; idx[i++] = j; j+=b;);
    goto loop;
  exit:
    return;
} /* end Sequentially_Index() */


/* memset() implementation in glibc, at least, is */
/* not in raw assembly: try looking for memset.c. */


/* swap contents of *a and *b, each of size_in_bytes bytes */
void swapb (char *a, char *b, size_t size_in_bytes)
{
    register char c, *d=a+size_in_bytes;
    for (; a<d; a++,b++) SWAP(*a, *b, c);
    return;
} /* end swapb() */


/* swap contents of *a and *b, each of size_in_ints integers */
void swapi (int *a, int *b, size_t size_in_ints)
{
    register int c, *d=a+size_in_ints;
    for (; a<d; a++,b++) SWAP(*a, *b, c);
    return;
} /* end swapi() */


/* swap contents of *a and *b, each of size_in_doubles doubles */
void swapd (double *a, double *b, size_t size_in_doubles)
{
    register double c, *d=a+size_in_doubles;
    for (; a<d; a++,b++) SWAP(*a, *b, c);
    return;
} /* end swapd() */


/* randomly permute array a[min..max] (e.g. a[0..n-1]),  */
/* each element, like a[max], is of size_in_bytes bytes. */
char *r_permuteb (char *a, int min, int max, size_t size_in_bytes)
{
    int i, j;
    ENSURE(min,max,j);
    for (i=min; i<max; i++)
    { /* randomly select a member in i..max */
	j = Fran(i,max);
	swapb (a+i*size_in_bytes, a+j*size_in_bytes, size_in_bytes);
    }
    return (a);
} /* end r_permuteb() */


/* randomly permute array a[min..max] (e.g. a[0..n-1]), */
/* each element, like a[max], is of size_in_ints ints.  */
int *r_permutei (int *a, int min, int max, size_t size_in_ints)
{
    int i, j;
    ENSURE(min,max,j);
    for (i=min; i<max; i++)
    { /* randomly select a member in i..max */
	j = Fran(i, max);
	swapi (a+i*size_in_ints, a+j*size_in_ints, size_in_ints);
    }
    return (a);
} /* end r_permutei() */


/* randomly permute array a[min..max] (e.g. a[0..n-1]), */
/* each element, like a[max], is of size_in_doubles doubles. */
double *r_permuted
(double *a, int min, int max, size_t size_in_doubles)
{
    int i, j;
    ENSURE(min,max,j);
    for (i=min; i<max; i++)
    { /* randomly select a member in i..max */
	j = Fran(i, max);
	swapd (a+i*size_in_doubles, a+j*size_in_doubles, size_in_doubles);
    }
    return (a);
} /* end r_permuted() */


/* randomly permute index array a[min..max] (e.g. a[0..n-1]) */
int *r_permute (int *a, int min, int max)
{
    int i, j, tmp;
    ENSURE(min,max,j);
    for (i=min; i<max; i++)
    { /* randomly select a member in i..max */
	j = Fran(i, max);
	SWAP(a[i], a[j], tmp);
    }
    return (a);
} /* end r_permute() */


/* determine if values in idx[] is a permutation of min..max */
int is_permutation (int idx[], int min, int max)
{
    int i, n;
    static BMAP(a,is_permutation_STATIC_BITS);
    Bmap *b;
    ENSURE(min,max,i);
    n = max-min+1;
    b = Bmem(n,a,VBmem);
    BZERO(b,n);
    for (i=0; i<n; i++)
	if ((idx[i]>=min)&&(idx[i]<=max)) BSET(b,idx[i]-min);
	else
	{
	    FREE(b,a,VFREE);
	    return (FALSE);
	}
    for (i=0; i<n; i++)
	if (!BVAL(b,i))
	{
	    FREE(b,a,VFREE);
	    return (FALSE);
	}
    FREE(b,a,VFREE);
    return(TRUE);
} /* end is_permutation() */

#ifdef _is_permutation_TEST
#define min  15
#define max  370
int main (int argc, char *argv[])
{
    int idx [max-min+1];
    AX_sequentially_index (idx, min, max);
    AX_r_permute(idx,0,max-min);
    if (! AX_is_permutation(idx,min,max))
	printf ("one of your subroutines is wrong.\n");
    return(0);
}
#endif  /* _is_permutation_TEST */


/* b[i] := a[idx[i]], i=0..N-1, each a[i] is object of size_in_chars chars */
void rearrange (int N, size_t size_in_chars, char a[], int idx[], char b[])
{
    register int i, j;
    for (i=0; i<N; i++)
	for (j=0; j<size_in_chars; j++)
	    b[i*size_in_chars+j] = a[idx[i]*size_in_chars+j];
    return;
} /* end rearrange() */


/* a[] := a[idx[]], each a[i] is an object of size_in_chars chars */
void Rearrange (int N, size_t size_in_chars, char a[], int idx[])
{
    register int i, j;
    char *b = MCmem(size_in_chars);
    for (i=0; i<N; i++)
	for (j=0; j<size_in_chars; j++)
	    b[i*size_in_chars+j] = a[idx[i]*size_in_chars+j];
    memcpy ((void *)a, (void *)b, (size_t)N*size_in_chars*sizeof(char));
    MFREE(b);
    return;
} /* end Rearrange() */


/* b[i] := a[idx[i]], i=0..N-1, each a[i] is object of size_in_ints ints */
void rearrangei (int N, size_t size_in_ints, int a[], int idx[], int b[])
{
    register int i, j;
    for (i=0; i<N; i++)
	for (j=0; j<size_in_ints; j++)
	    b[i*size_in_ints+j] = a[idx[i]*size_in_ints+j];
    return;
} /* end rearrangei() */


/* a[] := a[idx[]], each a[i] is an object of size_in_ints ints */
void Rearrangei (int N, size_t size_in_ints, int a[], int idx[])
{
    register int i, j;
    int *b = MImem(size_in_ints);
    for (i=0; i<N; i++)
	for (j=0; j<size_in_ints; j++)
	    b[i*size_in_ints+j] = a[idx[i]*size_in_ints+j];
    memcpy ((void *)a, (void *)b, (size_t)N*size_in_ints*sizeof(int));
    MFREE(b);
    return;
} /* end Rearrangei() */


/* b[i] := a[idx[i]], i=0..N-1, each a[i] is object of size_in_dbles doubles */
void rearranged
(int N, size_t size_in_dbles, double a[], int idx[], double b[])
{
    register int i, j;
    for (i=0; i<N; i++)
	for (j=0; j<size_in_dbles; j++)
	    b[i*size_in_dbles+j] = a[idx[i]*size_in_dbles+j];
    return;
} /* end rearranged() */


/* a[] := a[idx[]], each a[i] is an object of size_in_dbles doubles */
void Rearranged (int N, size_t size_in_dbles, double a[], int idx[])
{
    register int i, j;
    double *b = Mmem(size_in_dbles);
    for (i=0; i<N; i++)
	for (j=0; j<size_in_dbles; j++)
	    b[i*size_in_dbles+j] = a[idx[i]*size_in_dbles+j];
    for (i=0; i<N*size_in_dbles; i++) a[i] = b[i];
    memcpy ((void *)a, (void *)b, (size_t)N*size_in_dbles*sizeof(double));
    MFREE(b);
    return;
} /* end Rearranged() */


/* returns max(a[]) */
int IMAX (int n, int a[])
{
    int i, max;
    for (max=a[0], i=1; i<n; i++)
        if (a[i]>max) max = a[i];
    return(max);
} /* end IMAX() */


/* returns index imax, such that a[imax] >= a[] */
int Imax (int n, int a[])
{
    int i, max, imax;
    for (max=a[0], imax=0, i=1; i<n; i++)
        if (a[i] > max)
        {
            imax = i;
            max = a[i];
        }
    return(imax);
} /* end Imax() */


/* returns min(a[]) */
int IMIN (int n, int a[])
{
    int i, min;
    for (min=a[0], i=1; i<n; i++)
        if (a[i] < min) min = a[i];
    return(min);
} /* end IMIN() */


/* returns index imin, such that a[imin] <= a[] */
int Imin (int n, int a[])
{
    int i, imin, min;
    for (min=a[0], imin=0, i=1; i<n; i++)
        if (a[i] < min)
        {
            imin = i;
            min = a[i];
        }
    return(imin);
} /* end Imin() */


/* summing all elements in a[] */
int Isum (int n, int a[])
{
    register int i, isum=0;
    for (i=n; i--;) isum += a[i];
    return(isum);
} /* end Isum() */


/* see sorters.c: qsort_glibc */
#define MAXTHRESH       2
#define PUSH(low,high)  ((*(top++))=(low),(*(top++))=(high))
#define	POP(low,high)   (high=(*(--top)),low=(*(--top)))
/* modified from /usr/local/src/glibc-2.1.2/stdlib/Iqsort.c */
void Iqsort_glibc (int N,  int x[], int idx[], int option)
{
    register int tmp, pivot;
    register int *base_ptr=idx, *left_ptr, *right_ptr;
    static int **top, *lo, *hi, *mid, *tmp_ptr, *end_ptr,
	*run_ptr, *thresh, *trav;
    static int *stack[16*sizeof(unsigned long)];
    
    if (option==USE_NEW_IDX) Sequentially_index(N,idx);
    
    if (N > MAXTHRESH)
    {
        lo = idx;
        hi = &lo[N-1];
        top = stack + 2;
        while (stack < top)
        {
            mid = lo + ((hi - lo) >> 1);
            if (x[*mid] < x[*lo]) SWAP(*mid, *lo, tmp);
            if (x[*hi] < x[*mid]) SWAP(*mid, *hi, tmp);
	    else goto jump_over;
            if (x[*mid] < x[*lo]) SWAP(*mid, *lo, tmp);
        jump_over:;
            pivot = *mid;
            left_ptr  = lo + 1;
            right_ptr = hi - 1;
            do
            {
                while (x[pivot] < x[*right_ptr]) right_ptr--;
                while (x[*left_ptr] < x[pivot])  left_ptr++;
                if (left_ptr <= right_ptr)
                {
		    tmp = *left_ptr;
		    *(left_ptr++) = *right_ptr;
		    *(right_ptr--) = tmp;
                }
            } while (left_ptr <= right_ptr);
            if ((right_ptr - lo) <= MAXTHRESH)
            {
                if ((hi - left_ptr) <= MAXTHRESH)
		    POP(lo, hi);
                else lo = left_ptr;
            }
            else if ( (hi - left_ptr) <= MAXTHRESH)
		hi = right_ptr;
            else if ((right_ptr - lo) > (hi - left_ptr))
            {
                PUSH(lo, right_ptr);
                lo = left_ptr;
            }
            else
            {
                PUSH(left_ptr, hi);
                hi = right_ptr;
            }
        }
    }
    end_ptr = &idx[N-1];
    tmp_ptr = base_ptr;
    thresh = MIN(end_ptr, base_ptr+MAXTHRESH);
    for (run_ptr = tmp_ptr+1; run_ptr <= thresh; run_ptr++)
	if (x[*run_ptr] < x[*tmp_ptr]) tmp_ptr = run_ptr;
    if (tmp_ptr != base_ptr) SWAP(*tmp_ptr, *base_ptr, tmp);
    run_ptr = base_ptr+1;
    while (++run_ptr <= end_ptr)
    {
	tmp_ptr = run_ptr-1;
	while ( x[*run_ptr] < x[*tmp_ptr] ) tmp_ptr--;
	tmp_ptr++;
	if (tmp_ptr != run_ptr)
	{
	    trav = run_ptr+1;
	    while (--trav >= run_ptr)
	    {
		tmp = *trav;
		for (hi = lo = trav; --lo >= tmp_ptr; hi = lo)
		    *hi = *lo;
		*hi = tmp;
	    }
	}
    }
}
#undef POP
#undef PUSH
#undef MAXTHRESH
/* end Iqsort_glibc() */


/* see sorters.c: qsort_numerical_recipes */
#define QUICKSORT_M  7
void Iqsort_numerical_recipes (int N,  int x[], int idx[], int option)
{
    register int i, j, tmp;
    int idxt, ir=N, k, l=1, jstack=0;
    static int istack[128];
    register int a;
    if (option==USE_NEW_IDX) Sequentially_index(N,idx);
    idx--;
    for (;;)
    {
	if (ir-l < QUICKSORT_M)
	{
	    for (j=l+1; j<=ir; j++)
	    {
		idxt = idx[j];
		a = x[idxt];
		for (i=j-1; i>=l; i--)
		{
		    if (x[idx[i]] <= a) break;
		    idx[i+1]=idx[i];
		}
		idx[i+1]=idxt;
	    }
	    if (jstack == 0) break;
	    ir = istack[jstack--];
	    l = istack[jstack--];
	}
	else
	{
	    k = (l+ir) >> 1;
	    SWAP(idx[k],idx[l+1],tmp);
	    if (x[idx[l]] > x[idx[ir]]) SWAP(idx[l],idx[ir],tmp);
	    if (x[idx[l+1]] > x[idx[ir]]) SWAP(idx[l+1],idx[ir],tmp);
	    if (x[idx[l]] > x[idx[l+1]]) SWAP(idx[l],idx[l+1],tmp);
	    i = l+1;
	    j = ir;
	    idxt = idx[l+1];
	    a = x[idxt];
	    for (;;)
	    {
		do i++; while (x[idx[i]] < a);
		do j--; while (x[idx[j]] > a);
		if (j < i) break;
		SWAP (idx[i],idx[j],tmp);
	    }
	    idx[l+1]=idx[j];
	    idx[j]=idxt;
	    jstack += 2;
	    if (ir-i+1 >= j-l)
	    {
		istack[jstack] = ir;
		istack[jstack-1] = i;
		ir = j-1;
	    }
	    else
	    {
		istack[jstack] = j-1;
		istack[jstack-1] = l;
		l = i;
	    }
	}
    }
}
#undef QUICKSORT_M
/* end Iqsort_numerical_recipes() */


/* determine whether x[idx[i]], i=0..N-1, is non-decreasing */
int Is_nondecreasing (int N, int x[], int idx[])
{
    register int i;
    for (i=0; i<N-1; i++)
        if ( x[idx[i]] > x[idx[i+1]] )
            return (FALSE);
    return (TRUE);
} /* end Is_nondecreasing() */


/***************************************************/
/* Fastest int array index sort subroutine contest */
/***************************************************/
#ifdef _IsortersContest_TEST
#include <Timer.h>
#define N 1000001
int main()
{
    struct List
    {
	char *name;
	void (*fun) (int M,  int x[], int idx[], int option);
    } funs[] = {
        {"glibc Iqsort", &Iqsort_glibc},
        {"numerical recipes Iqsort", &Iqsort_numerical_recipes},
    };
    static int x[N];
    static int idx[N];
    int i;
    printf ("%ld random numbers are generated...\n", (long)N);
    for (i=0; i<N; i++) x[i] = Fran(0, 10*N);
    for (i=0; i<sizeof(funs)/sizeof(struct List); i++)
    {
	printf ("Testing %s..", funs[i].name);
	start_chronometer();
	funs[i].fun (N, x, idx, USE_NEW_IDX);
	stop_chronometer();
	if ( (!is_permutation(idx,0,N-1)) || (!Is_nondecreasing(N,x,idx)) )
	    printf (" it gave wrong order\n");
	else printf (" %f seconds\n", stopped_usertime());
    }
    start_chronometer();
    is_permutation(idx,0,N-1);
    Is_nondecreasing(N,x,idx);
    stop_chronometer();
    printf ("correctness check.. %f seconds\n", stopped_usertime());
    return (0);
}
#endif  /* _IsortersContest_TEST */


#ifdef _IqSort_TEST
#define N 4
int main (int argc, char *argv[])
{
    int i, idx[N];
    int s[3*N] = {0, 2,  0,
                  0, 1,  0,
                  0, 3,  0,
                  0, -1, 0};
    IqSort_numerical_recipes(N,s,idx,1,3);
    for (i=0; i<N; i++) printf ("%d\n", s[idx[i]]);
    return (0);
}
#endif /* _IqSort_TEST */


/* Find out if there are identical elements in an integer array */
bool Iall_different (int N,  int x[])
{
    int i, *idx;
    MALLOC (Iall_different, idx, N, int);
    Iqsort_numerical_recipes (N, x, idx, USE_NEW_IDX);
    for (i=0; i<N-1; i++)
        if ( x[idx[i]] == x[idx[i+1]] )
        {
            Free (idx);
            return (FALSE);
        }
    Free (idx);
    return (TRUE);
} /* end Iall_different() */


#ifdef _Iall_different_TEST
#define N 4
int main (int argc, char *argv[])
{
    int x[] = {0, 2, -10};
    printf ("It %s all different!\n",
            Iall_different(sizeof(x)/sizeof(int), x) ? "is" : "is not");
    return (0);
}
#endif /* _Iall_different_TEST */

