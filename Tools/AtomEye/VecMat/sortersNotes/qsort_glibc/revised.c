/* Modified from Glibc stdlib source qsort.c */

#include "match.h"

#define MAXTHRESH 2
#define PUSH(low,high) (((top->lo=(low)),(top->hi=(high)),++top))
#define	POP(low,high) ((--top,(low = top->lo),(high = top->hi)))
struct IdxStack
{
    IdxType *lo, *hi;
};
void glibc (QUICKSORT_INPUT_ARGS)
{
    register IdxType tmp, pivot;
    register IdxType *base_ptr = idx, *left_ptr, *right_ptr;
    static IdxType *lo, *hi, *mid, *tmp_ptr, *end_ptr,
	*run_ptr, *thresh, *trav;
    static struct IdxStack stack[8*sizeof(unsigned long int)];
    static struct IdxStack *top;
    
    if (N == 0) return;
    for (tmp=0; tmp<N; tmp++) idx[tmp] = tmp;
    if (N > MAXTHRESH)
    {
        lo = idx;
        hi = &lo[N-1];
        top = stack + 1;
        while (stack < top)
        {
            mid = lo + ((hi - lo) >> 1);
            if(x[*mid] < x[*lo]) SWAP(*mid, *lo, tmp);
            if(x[*hi] < x[*mid]) SWAP(*mid, *hi, tmp);
	    else goto jump_over;
            if(x[*mid] < x[*lo]) SWAP(*mid, *lo, tmp);
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
            }
            while (left_ptr <= right_ptr);
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
    thresh = MIN(end_ptr,base_ptr+MAXTHRESH);
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
#undef MIN
#undef SWAP
#undef POP
#undef PUSH
#undef MAXTHRESH
