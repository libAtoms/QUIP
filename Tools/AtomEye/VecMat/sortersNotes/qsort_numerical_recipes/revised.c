/* Numerical Recipes */
#include "match.h"
#define QUICKSORT_M 7
#define QUICKSORT_S 50
void nr (QUICKSORT_INPUT_ARGS)
{
    register IdxType i, j, tmp;
    IdxType idxt, ir=N, k, l=1, jstack=0;
    static IdxType istack[QUICKSORT_S];
    register DataType a;
    x--;
    idx--;
    for (i=1; i<=N; i++) idx[i] = i;
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
	    if (jstack > QUICKSORT_S)
	    {
		printf ("error: quicksort: QUICKSORT_S too small");
                exit(1);
	    }
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
    for (i=1; i<=N; i++) idx[i]--;
}
#undef QUICKSORT_S
#undef QUICKSORT_M
