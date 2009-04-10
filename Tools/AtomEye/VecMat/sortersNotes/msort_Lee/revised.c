#include "match.h"

#define MAXN (1024L*1024L)
#define MAYBE_THRESHOLD 8

static IdxType out[MAXN];

static void mergeIt (QUICKSORT_INPUT_ARGS)
{
    int split;
    int a, b;
    IdxType *c;
    int in_order;
    split = (N + 1) >> 1;
    if (split > 1) mergeIt (split, x, idx);
    if (N - split > 1) mergeIt (N-split, x, idx+split);
    if (N>MAYBE_THRESHOLD && x[idx[split]]>x[idx[split-1]]) return;
    a = 0; 
    b = split;
    c = out;
    while (a < split)
    {
	if (b >= N) in_order = 1;
	else in_order = (x[idx[a]] <= x[idx[b]]);
	if (in_order)
	{
	    *c = *(idx + a);
	    c ++;
	    a ++;
	}
	else
	{
	    *c = *(idx + b);
	    c ++;
	    b ++;
	}
    }
    for (a = 0; a < c - out; a ++) idx[a] = out[a];
}

void merge (QUICKSORT_INPUT_ARGS)
{
    IdxType tmp;
    if (N > MAXN)
    {
	printf ("error: mergesort: N = %ld exceed MAXN = %ld\n",
		(long)N, (long)MAXN);
	exit(1);
    }
    for (tmp=0; tmp<N; tmp++) idx[tmp] = tmp;
    mergeIt (N, x, idx);
}
