/*********************************************/
/* libVecMat: -llapck -lblas -lm             */
/*            -lScalar -lIO                  */
/*                                           */
/* General Vector & Matrix Operation Library */
/*                                           */
/* Nov 19 1999 Ju Li <liju99@mit.edu>        */
/*********************************************/

#include "VecMat.h"
#ifdef _sortersContest_TEST
#include <Timer.h>  /* internal profiling */
#endif

/**********************************************/
/* fast sorting routines for array of doubles */
/**********************************************/

/* determine whether x[idx[i]], i=0..N-1, is non-decreasing */
int is_nondecreasing (int N, double x[], int idx[])
{
    register int i;
    for (i=0; i<N-1; i++)
	if ( x[idx[i]] > x[idx[i+1]] )
	    return (FALSE);
    return (TRUE);
} /* end is_nondecreasing() */


/* rearrange array of doubles: y[i] := x[idx[i]], i=0..N-1 */
void Vrearrange (int N, double x[], int idx[], double y[])
{
    register int i;
    for (i=0; i<N; i++)  y[i] = x[idx[i]];
    return;
} /* end Vrearrange() */


/* rearrange array of doubles: x[] := x[idx[]] */
void VRearrange (int N, double x[], int idx[])
{
    register int i;
    double *y = Mmem(N);
    for (i=0; i<N; i++) y[i] = x[idx[i]];
    memcpy ((void *)x, (void *)y, (size_t)N*sizeof(double));
    MFREE (y);
    return;
} /* end VRearrange() */


/* qsort_glibc() */
#define MAXTHRESH       2
#define PUSH(low,high)  ((*(top++))=(low),(*(top++))=(high))
#define	POP(low,high)   (high=(*(--top)),low=(*(--top)))
/* modified from /usr/local/src/glibc-2.1.2/stdlib/qsort.c */
void qsort_glibc (int N,  double x[], int idx[], int option)
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
/* end qsort_glibc() */


/* qsort_numerical_recipes() */
#define QUICKSORT_M  7
void qsort_numerical_recipes (int N,  double x[], int idx[], int option)
{
    register int i, j, tmp;
    int idxt, ir=N, k, l=1, jstack=0;
    static int istack[128];
    register double a;
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
/* end qsort_numerical_recipes() */


/* msort_Lee() */
static int *msort_Lee_out;
static double *msort_Lee_x;
/* http://www.ontek.com/mikey/flogger.tar.uu */
#define ensure(i,j,tmp) if (msort_Lee_x[i]>msort_Lee_x[j]) SWAP(i,j,tmp);
static void msort_Lee_mergeIt (int N, int idx[])
{
    register int split, tmp;
    register int *a, *b, *c, *d;
    if (N <= 4)
    {
	if (N == 4)
	{
	    ensure(idx[0],idx[1],tmp);
	    ensure(idx[1],idx[2],tmp);
	    ensure(idx[2],idx[3],tmp);
	    ensure(idx[0],idx[1],tmp);
	    ensure(idx[1],idx[2],tmp);
	    ensure(idx[0],idx[1],tmp);
	    goto exit;
	}
	if (N == 3)
	{
	    ensure(idx[0],idx[1],tmp);
	    ensure(idx[1],idx[2],tmp);
	    ensure(idx[0],idx[1],tmp);
	    goto exit;
	}
	if (N == 2) ensure(idx[0],idx[1],tmp);
	goto exit;
    }
    split = (N + 1) >> 1;
    msort_Lee_mergeIt (split, idx);
    msort_Lee_mergeIt (N-split, idx+split);
    if (msort_Lee_x[idx[split]] >= msort_Lee_x[idx[split-1]]) goto exit;
    a = idx;
    b = idx + split;
    c = msort_Lee_out;
    while (a < idx + split)
	if ( (b >= idx + N) || (msort_Lee_x[*a] <= msort_Lee_x[*b]) )
	    *(c++) = *(a++);
	else
	    *(c++) = *(b++);
    d = msort_Lee_out;
loop:
    REPEAT256( if (d==c) goto exit; *(idx++) = *(d++); );
    goto loop;
exit:
    return;
} /* end msort_Lee_mergeIt() */

void msort_Lee (int N,  double x[], int idx[], int option)
{
    if (option==USE_NEW_IDX) Sequentially_index (N, idx);
    msort_Lee_out = MImem(N);
    msort_Lee_x = x;
    msort_Lee_mergeIt (N, idx);
    MFREE(msort_Lee_out);
    return;
} /* end msort_Lee() */


/* very fast index sort of four numbers in fastsort_x */
static double *fastsort_x;
static void fastsort4 (int idx[])
{
    register int itmp;
    register double x0, x1, x2, x3, xtmp;
    x0 = fastsort_x[idx[0]];
    x1 = fastsort_x[idx[1]];
    x2 = fastsort_x[idx[2]];
    x3 = fastsort_x[idx[3]];
    if (x0 > x1) { SWAP(x1,x0,xtmp); SWAP(idx[1],idx[0],itmp); }
    if (x2 > x3) { SWAP(x2,x3,xtmp); SWAP(idx[2],idx[3],itmp); }
    if (x1 > x2) { SWAP(x1,x2,xtmp); SWAP(idx[1],idx[2],itmp); }
    else goto exit;
    if (x0 > x1) { x1 = x0; SWAP(idx[0],idx[1],itmp); }
    if (x2 > x3) { x2 = x3; SWAP(idx[2],idx[3],itmp); }
    if (x1 > x2) SWAP(idx[1],idx[2],itmp);
  exit:
    return;
}  /* end fastsort4() */


/* very fast index sort of eight numbers in fastsort_x */
static void fastsort8 (int idx[])
{
    register int *a, *b;
    register int itmp;
    static int swap[7];
    register double x0, x1, x2, x3, x4, x5, xtmp;
    
    x0 = fastsort_x[idx[0]];
    x1 = fastsort_x[idx[1]];
    x2 = fastsort_x[idx[2]];
    x3 = fastsort_x[idx[3]];
    /* #define SWAP(x,y,tmp) ((tmp)=(x),(x)=(y),(y)=(tmp)) */
    if (x0 > x1) { SWAP(x1,x0,xtmp); SWAP(idx[1],idx[0],itmp); }
    if (x2 > x3) { SWAP(x2,x3,xtmp); SWAP(idx[2],idx[3],swap[0]); }
    if (x1 > x2) { SWAP(x1,x2,xtmp); SWAP(idx[1],idx[2],swap[1]); }
    else goto exit1;
    /* keep the bottom of this stack in x0 */
    if (x0 > x1) { SWAP(x0,x1,xtmp); SWAP(idx[0],idx[1],swap[2]); }
    /* keep the top of this stack in x3 */
    if (x2 > x3) { SWAP(x2,x3,xtmp); SWAP(idx[2],idx[3],swap[3]); }
    if (x1 > x2)  SWAP(idx[1],idx[2],swap[4]);
  exit1:
    
    x1 = fastsort_x[idx[4]];
    x2 = fastsort_x[idx[5]];
    x4 = fastsort_x[idx[6]];
    x5 = fastsort_x[idx[7]];
    if (x1 > x2) { SWAP(x2,x1,xtmp); SWAP(idx[5],idx[4],swap[5]); }
    if (x4 > x5) { SWAP(x4,x5,xtmp); SWAP(idx[6],idx[7],swap[0]); }
    if (x2 > x4) { SWAP(x2,x4,xtmp); SWAP(idx[5],idx[6],swap[1]); }
    else goto exit2;
    /* keep the bottom of this stack in x1 */
    if (x1 > x2) { SWAP(x1,x2,xtmp); SWAP(idx[4],idx[5],swap[5]); }
    /* keep the top of this stack in x5 */
    if (x4 > x5) { SWAP(x4,x5,xtmp); SWAP(idx[6],idx[7],itmp); }
    if (x2 > x4)  SWAP(idx[5],idx[6],swap[6]);
  exit2:

    /* stack-stack boundary test */
    if (x1 >= x3) goto exit;
    if (x0 >= x5)
    {
        SWAP(idx[0],idx[4],itmp);
        SWAP(idx[1],idx[5],swap[0]);
        SWAP(idx[2],idx[6],swap[1]);
        SWAP(idx[3],idx[7],swap[2]);
        goto exit;
    }
    /* assert there is going to be mixing */
    
    a = idx;    /* A stack pointer */
    b = idx+4;  /* B stack pointer */
    
    if (x0 <= x1) swap[0] = *(a++);  else swap[0] = *(b++);
    /* Because mixing must happen, there cannot be spent */
    /* A or B stack until the fifth in the swap space    */
    if (fastsort_x[*a] <= fastsort_x[*b]) swap[1] = *(a++);
    else swap[1] = *(b++);
    if (fastsort_x[*a] <= fastsort_x[*b]) swap[2] = *(a++);
    else swap[2] = *(b++);
    if (fastsort_x[*a] <= fastsort_x[*b]) swap[3] = *(a++);
    else swap[3] = *(b++);
    if (fastsort_x[*a] <= fastsort_x[*b]) swap[4] = *(a++);
    else swap[4] = *(b++);

    /* now one of A or B stack could be spent */
    if (a == idx+4)
    {
        idx[0] = swap[0];
        idx[1] = swap[1];
        idx[2] = swap[2];
        idx[3] = swap[3];
        idx[4] = swap[4];
        goto exit;
    }
    else if (b == idx+8)
    {
        idx[0] = swap[0];
        idx[5] = idx[1];
        idx[6] = idx[2];
        idx[7] = idx[3];
        idx[1] = swap[1];
        idx[2] = swap[2];
        idx[3] = swap[3];
        idx[4] = swap[4];
        goto exit;
    }

    if (fastsort_x[*a] <= fastsort_x[*b]) swap[5] = *(a++);
    else swap[5] = *(b++);
    
    if (a == idx+4)
    {
        idx[0] = swap[0];
        idx[1] = swap[1];
        idx[2] = swap[2];
        idx[3] = swap[3];
        idx[4] = swap[4];
        idx[5] = swap[5];
        goto exit;
    }
    else if (b == idx+8)
    {
        idx[0] = swap[0];
        idx[1] = swap[1];
        idx[6] = idx[2];
        idx[7] = idx[3];
        idx[2] = swap[2];
        idx[3] = swap[3];
        idx[4] = swap[4];
        idx[5] = swap[5];
        goto exit;
    }

    if (fastsort_x[*a] <= fastsort_x[*b]) swap[6] = *(a++);
    else swap[6] = *(b++);
    
    if (a == idx+4)
    {
        idx[0] = swap[0];
        idx[1] = swap[1];
        idx[2] = swap[2];
        idx[3] = swap[3];
        idx[4] = swap[4];
        idx[5] = swap[5];
        idx[6] = swap[6];
        goto exit;
    }
    else if (b == idx+8)
    {
        idx[0] = swap[0];
        idx[1] = swap[1];
        idx[2] = swap[2];
        idx[7] = idx[3];
        idx[3] = swap[3];
        idx[4] = swap[4];
        idx[5] = swap[5];
        idx[6] = swap[6];
    }
  exit:
    return;
} /* end fastsort8() */


static void fastsort16 (int idx[])
{
    register int *a, *b, *c;
    static int swap[16];
    fastsort8 (idx);
    fastsort8 (idx+8);
    if (fastsort_x[*(idx+8)] >= fastsort_x[*(idx+7)])
        goto exit;
    c = swap;
    a = idx;
    b = idx+8;
    REPEAT8( if (fastsort_x[*a] <= fastsort_x[*b]) *(c++) = *(a++); \
             else *(c++) = *(b++); );
    while (a < idx+8)
        if ( (b == idx+16) || (fastsort_x[*a] <= fastsort_x[*b]) )
            *(c++) = *(a++);
        else *(c++) = *(b++);
    a = idx;
    b = swap;
    REPEAT16 ( if (b==c) goto exit; *(a++) = *(b++); );
  exit:
    return;
}  /* end fastsort16() */


static void fastsort32 (int idx[])
{
    register int *a, *b, *c;
    static int swap[32];
    fastsort16 (idx);
    fastsort16 (idx+16);
    if (fastsort_x[*(idx+16)] >= fastsort_x[*(idx+15)])
        goto exit;
    c = swap;
    a = idx;
    b = idx+16;
    REPEAT16( if (fastsort_x[*a] <= fastsort_x[*b]) *(c++) = *(a++); \
             else *(c++) = *(b++); );
    while (a < idx+16)
        if ( (b == idx+32) || (fastsort_x[*a] <= fastsort_x[*b]) )
            *(c++) = *(a++);
        else *(c++) = *(b++);
    a = idx;
    b = swap;
    REPEAT32 ( if (b==c) goto exit; *(a++) = *(b++); );
  exit:
    return;
}  /* end fastsort32() */


static void fastsort64 (int idx[])
{
    register int *a, *b, *c;
    static int swap[64];
    fastsort32 (idx);
    fastsort32 (idx+32);
    if (fastsort_x[*(idx+32)] >= fastsort_x[*(idx+31)])
        goto exit;
    c = swap;
    a = idx;
    b = idx+32;
    REPEAT32( if (fastsort_x[*a] <= fastsort_x[*b]) *(c++) = *(a++); \
             else *(c++) = *(b++); );
    while (a < idx+32)
        if ( (b == idx+64) || (fastsort_x[*a] <= fastsort_x[*b]) )
            *(c++) = *(a++);
        else *(c++) = *(b++);
    a = idx;
    b = swap;
    REPEAT64 ( if (b==c) goto exit; *(a++) = *(b++); );
  exit:
    return;
}  /* end fastsort64() */


static void fastsort128 (int idx[])
{
    register int *a, *b, *c;
    static int swap[128];
    fastsort64 (idx);
    fastsort64 (idx+64);
    if (fastsort_x[*(idx+64)] >= fastsort_x[*(idx+63)])
        goto exit;
    c = swap;
    a = idx;
    b = idx+64;
    REPEAT64( if (fastsort_x[*a] <= fastsort_x[*b]) *(c++) = *(a++); \
              else *(c++) = *(b++); );
    while (a < idx+64)
        if ( (b == idx+128) || (fastsort_x[*a] <= fastsort_x[*b]) )
            *(c++) = *(a++);
        else *(c++) = *(b++);
    a = idx;
    b = swap;
    REPEAT128 ( if (b==c) goto exit; *(a++) = *(b++); );
  exit:
    return;
}  /* end fastsort128() */


static void fastsort256 (int idx[])
{
    register int *a, *b, *c;
    static int swap[256];
    fastsort128 (idx);
    fastsort128 (idx+128);
    if (fastsort_x[*(idx+128)] >= fastsort_x[*(idx+127)])
        goto exit;
    c = swap;
    a = idx;
    b = idx+128;
    while (a < idx+128)
        if ( (b == idx+256) || (fastsort_x[*a] <= fastsort_x[*b]) )
            *(c++) = *(a++);
        else *(c++) = *(b++);
    a = idx;
    b = swap;
    REPEAT256 ( if (b==c) goto exit; *(a++) = *(b++); );
  exit:
    return;
}  /* end fastsort256() */


/* smallest cell size: use fast routine for that */
#define BASIC 256
/* non-recursive merge-sort by Ju Li */
void msort_Ju (int N,  double x[], int idx[], int option)
{
    register int i, cell;
    static int *swap;
    register int *a, *b, *C, *c;
    
    /* buffer space for merge sort */
    swap = MImem(N);
    if (option==USE_NEW_IDX) Sequentially_index (N, idx);
    
    fastsort_x = x;
    for ( a = idx; a < idx+N/BASIC*BASIC; a+=BASIC )
        fastsort256 (a);
    
    /* mop up the rest with a reasonable sorter */
    if (a != idx+N)
        qsort_numerical_recipes( idx+N-a, x, a, USE_OLD_IDX);
    
    cell = BASIC;
    while ( cell < N )
    {
        a = idx;         /* A stack pointer */
        b = idx + cell;  /* B stack pointer */
        C = idx + MIN((cell << 1), N);  /* B stack end */
        /* i goes over number of groups: if cells is odd, all groups */
        /* are complete; if cells is even, the last group may be     */
        /* incomplete. In both cases we only need to rearrange each group */
        for (i=0; i < (((N-1)/cell+1)>>1); i++)
        { /* there is still competition */
            if ( x[idx[i*(cell<<1)+cell]] >= x[idx[i*(cell<<1)+cell-1]] )
                goto exit;
            c = swap;      /* storage stack pointer */
            while (a < idx+i*(cell<<1)+cell)
                if ( (b == C) || (x[*a] <= x[*b]) ) *(c++) = *(a++);
                else *(c++) = *(b++);
            /* a = idx+i*(cell<<1); */
            /* b = swap; */
            /* loop: */
            /* REPEAT16 ( if (b==c) goto exit; *(a++) = *(b++); ); */
            /* goto loop; */
            memcpy (idx+i*(cell<<1), swap, (c-swap)*sizeof(int));
          exit:
            a = C;
            b = a + cell;
            C = MIN( C+(cell<<1), idx+N);
        }
        cell <<= 1;
    }
    MFREE (swap);
    return;
} /* end msort_Ju() */
#undef BASIC


/******************************************************/
/* Fastest double array index sort subroutine contest */
/******************************************************/
#ifdef _sortersContest_TEST
#include <Timer.h>
#define N 1000001
void main()
{
    struct List
    {
	char *name;
	void (*fun) (int M,  double x[], int idx[], int option);
    } funs[] = {
        {"glibc qsort", &qsort_glibc},
        {"numerical recipes qsort", &qsort_numerical_recipes},
        {"recursive mergesort by Michael E. Lee", &msort_Lee},
	{"non-recursive mergesort by Ju Li", &msort_Ju},
    };
    static double x[N];
    static int idx[N];
    int i, j;
    printf ("%ld random numbers are generated...\n", (long)N);
    /* for (i=0; i<N; i++) x[i] = (double)FRANDOM()*10000.; */
    for (i=0; i<N; i++) x[i] = (double) i + 2 * FRANDOM();
    /* for (i=0; i<N; i++) x[i] = (double) i; */
    /* for (i=0; i<N; i++) x[i] = -i; */
    for (i=0; i<sizeof(funs)/sizeof(struct List); i++)
    {
	printf ("Testing %s..", funs[i].name);
	start_chronometer();
	funs[i].fun (N, x, idx, USE_NEW_IDX);
	stop_chronometer();
	if ( (!is_permutation(idx,0,N-1)) || (!is_nondecreasing(N,x,idx)) )
	    printf (" it gave wrong order\n");
	else printf (" %lf seconds\n", stopped_usertime());
    }
    start_chronometer();
    is_permutation(idx,0,N-1);
    is_nondecreasing(N,x,idx);
    stop_chronometer();
    printf ("correctness check.. %lf seconds\n", stopped_usertime());
}
#endif  /* _sortersContest_TEST */

/*******************************************************************/
/* Quicksort references from Akira Wada:                           */
/* http://www.mars.dti.ne.jp/~a-wada/qsortlib.html                 */
/*                                                                 */
/* 1. R.Sedgewick, Algorithms in C++                               */
/* 2. D.E.Knuth, The Art of Computer Programming, Vol.3            */
/* 3. ftp://ring.aist.go.jp/pub/TeX/CTAN/indexing/                 */
/*          makeindex/src/qsort.c                                  */
/* 4. ftp://ftp.sra.co.jp/pub/cmd/postgres/6.2.1/postgresql-6.2.1/ */
/*          src/backend/lib/qsort.c                                */
/* 5. ftp://ftp.jp.freebsd.org/pub/FreeBSD/FreeBSD-current/src/    */
/*         /sys/libkern/qsort.c                                    */
/* 6. ftp://ftp.nizkor.vex.net/local/src/darcy/common/qsort.c      */
/* 7. http://www.tokuyama.ac.jp/home/~kawamura/ qs6.c & qs7.c      */
/* 8. ftp://ftp.groggs.group.cam.ac.uk/pub/software/qsort.c-1.12   */
/* 9. ftp://ftp.linux.org.uk/pub/linux/libc/libc-5.4.44.tar.gz     */
/*         bc/stdlib/_quicksort.c  /libc/stdlib/qsor               */
/*******************************************************************/


#ifdef _qSort_TEST
#define N 4
int main (int argc, char *argv[])
{
    int i, idx[N];
    double s[3*N] = {0, 2,  0,
                     0, 1,  0,
                     0, 3,  0,
                     0, -1, 0};
    qSort_numerical_recipes(N,s,idx,1,3);
    for (i=0; i<N; i++) printf ("%g\n", s[idx[i]]);
    return (0);
}
#endif /* _qSort_TEST */
