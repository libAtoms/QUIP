/*********************************************/
/* libVecMat: -llapck -lblas -lm             */
/*            -lScalar -lIO                  */
/*                                           */
/* General Vector & Matrix Operation Library */
/*                                           */
/* Nov 19 1999 Ju Li <liju99@mit.edu>        */
/*********************************************/

#include "VecMat.h"

/******************************************************************/
/* 1D memory tectonics to represent 2D structure of neighbor list */
/******************************************************************/


/* clear all records in an uncompressed list */
void Iclearlist (int *idx, int n)
{
    register int i;
    for (i=0; i<n; i++) idx[2*i+1] = idx[2*i];
    return;
} /* end Iclearlist() */


/* return a one-dimensional list[] evenly indexed by (*idx)[] */
/* where list[(*idx)[2*i]<=k<(*idx)[2*i+1]] belong to owner i */
/* total n owners, each "maximally" keeps likely_size entries */
int *Icreatelist(int **idx, int n, int likely_size)
{
    int i, *list;
    /* cache-friendly layout */
    MALLOC( Icreatelist, *idx, n*2+1+n*likely_size, int );
    list = *idx + n*2 + 1;
    /* evenly spread out with 0 initial width */
    for (i=0; i<n; i++)
    {
	(*idx)[2*i]   = i * likely_size;
	(*idx)[2*i+1] = (*idx)[2*i];
    }
    /* the stopper element */
    (*idx)[2*n] = n * likely_size; 
    return(list);
} /* end Icreatelist() */


/* Same as Icreatelist() except using realloc() instead of malloc() */
int *Irecreatelist(int **idx, int n, int likely_size)
{
    int i, *list;
    /* cache-friendly layout */
    REALLOC( Irecreatelist, *idx, n*2+1+n*likely_size, int );
    list = *idx + n*2 + 1;
    /* evenly spread out with 0 initial width */
    for (i=0; i<n; i++)
    {
        (*idx)[2*i]   = i * likely_size;
	(*idx)[2*i+1] = (*idx)[2*i];
    }
    /* the stopper element */
    (*idx)[2*n] = n * likely_size; 
    return(list);
} /* end Irecreatelist() */


/* Kill memory holes, free extra space, and change */
/* [2*i]-[2*i+1] notation to [i]-[i+1] notation.   */
int *Icompresslist (int **idx, int n)
{
    register int i, j, *list=(*idx)+2*n+1;
    for (i=1; i<n; i++)
    {
        (*idx)[i] = (*idx)[i-1] + (*idx)[2*i-1] - (*idx)[2*i-2];
        for (j=(*idx)[2*i]; j<(*idx)[2*i+1]; j++)
            list[(*idx)[i]+j-(*idx)[2*i]] = list[j];
    }
    (*idx)[n] = (*idx)[n-1] + (*idx)[2*n-1] - (*idx)[2*n-2];
    memmove( (*idx)+n+1, list, (*idx)[n]*sizeof(int) );
    REALLOC ( Icompresslist, *idx, n+1+(*idx)[n], int );
    return ((*idx)+n+1);
} /* end Icompresslist() */


/* Insert memory holes of size "pad_each" to each entry and  */
/* change [i]-[i+1] notation back to [2*i]-[2*i+1] notation. */
int *Idecompresslist (int **idx, int n, int pad_each)
{
    register int i, j, *new_idx, *new_list, *list=(*idx)+n+1;
    MALLOC ( Idecompresslist, new_idx, 2*n+1+(*idx)[n]+n*pad_each, int );
    new_list = new_idx + 2*n + 1;
    new_idx[0] = 0;
    for (i=0; i<n; i++)
    {
        for (j=(*idx)[i]; j<(*idx)[i+1]; j++)
            new_list[new_idx[2*i]+j-(*idx)[i]] = list[j];
        new_idx[2*i+1] = new_idx[2*i]+j-(*idx)[i];
        new_idx[2*i+2] = new_idx[2*i+1] + pad_each;
    }
    free (*idx);
    *idx = new_idx;
    return (new_list);
} /* end Idecompresslist() */


/* free memory allocated by Icreatelist()/Irecreatelist() and set NULL */
void Ifreelist(int **idx, int **list)
{
    Free (*idx);
    *list = NULL;
    return;
} /* end Ifreelist() */


/***************************************************************/
/* Above are for combined idx/list allocation models which can */
/* be regarded as a subclass. Next are more general operators. */
/***************************************************************/


/* print allocation and occupation statistics */
void Ilist_statistics
(char *name, int idx[], int n, bool compressed, FILE *out)
{
    register int i, j;
    int occupied, max, min;
    double mean, fluc, avg_alloc;
    if (!out) return;
    fprintf (out, "%s %s list: %d entries, %ld bytes allocated,\n",
             compressed?"Compressed":"Uncompressed",
             name, n, compressed?LONG(n+1+idx[n])*sizeof(int):
             LONG(2*n+1+idx[2*n])*sizeof(int));
    if (compressed)
    {
        for (i=fluc=max=0,min=HUGE_INT; i<n; i++)
        {
            j = idx[i+1] - idx[i];
            fluc += SQUARE(j);
            if (j > max) max = j;
            if (j < min) min = j;
        }
        avg_alloc = mean = DOUBLE(idx[n]) / n;
        fluc = sqrt(fluc/n - SQUARE(mean));
        fprintf(out, "max=%d, min=%d, avg=%.2f, std.dev.=%.2f (%.1f%%).\n",
                max, min, mean, fluc, 100*fluc/((avg_alloc==0)?1:avg_alloc) );
    }
    else
    {
        for (i=occupied=fluc=max=0,min=HUGE_INT; i<n; i++)
        {
            j = idx[2*i+1] - idx[2*i];
            occupied += j;
            fluc += SQUARE(j);
            if (j > max) max = j;
            if (j < min) min = j;
        }
        mean = DOUBLE(occupied) / n;
        avg_alloc = DOUBLE(idx[2*n]) / n;
        fluc = sqrt(fluc/n - SQUARE(mean));
        fprintf( out, "max=%d, min=%d, avg=%.2f (%.1f%%), "
                 "std.dev.=%.2f (%.1f%%).\n", max, min, mean,
                 100*mean/avg_alloc, fluc,
                 100*fluc/((avg_alloc==0)?1:avg_alloc) );
    }
    return;
} /* end Ilist_statistics() */


/* Without changing idx[2*min] and idx[2*max], solve the memory */
/* conflict at idx[2*i+1], idx[2*i+2]. You must be sure there   */
/* is conflict idx[2*i+1]==idx[2*i+2] before calling Imakespace */
/* Imakespace (faster than IMAKESPACE) does NOT preserve order. */
void Imakespace (int idx[], int list[], int i, int min, int max)
{
    int j, k, a, left_shift; 
    /* look for left free space */
    for (j=i; (j>min)&&(idx[2*j]==idx[2*j-1]); j--);
    /* look for right free space */
    for (k=i+1; (k<max)&&(idx[2*k+1]==idx[2*k+2]); k++);
    if ((j > min) && (k < max))
	left_shift = /* to prevent systematic squeeze */
	     (idx[2*i+1]-idx[2*j] <  idx[2*k+1]-idx[2*i+2]) ||
	    ((idx[2*i+1]-idx[2*j] == idx[2*k+1]-idx[2*i+2]) &&
	     (Frandom() < 0.5));
    else if ((j > min) && (k >= max)) left_shift = 1;
    else if ((j <= min) && (k < max)) left_shift = 0;
    else pe ("Imakespace: min=%d max=%d jammed between "
             "i=%d and i+1.\n", min, max, i);
    /* lazy: shift only by 1 unit */
    if (left_shift)
    {
	for (a=j;   a<=i;     a++) list[idx[2*a]-1] = list[idx[2*a+1]-1];
	for (a=2*j; a<=2*i+1; a++) idx[a]--;
    }
    else
    {
	for (a=k;     a>=i+1;   a--) list[idx[2*a+1]] = list[idx[2*a]];
	for (a=2*k+1; a>=2*i+2; a--) idx[a]++;
    }
    return;
} /* end Imakespace() */


/* Safely append value to the end of owner i's (unordered) list; */
/* original entry order of this and other owners' may be altered */
void Iappend (int idx[], int list[], int value, int i, int min, int max)
{
    int conflict_count = 0;
    if (idx[2*i+1] < idx[2*i+2])
    {
	list[idx[2*i+1]++] = value;
	return;
    }
    else if (idx[2*i+1] == idx[2*i+2])
    {
	conflict_count++;
	Imakespace (idx, list, i, min, max);
	list[idx[2*i+1]++] = value;
	return;
    }
    pe ("Iappend: detects overflow at i=%d and i+1,\n"
        "min=%d, max=%d, conflict_count=%d\n",
        i, min, max, conflict_count);
    return;
} /* end Iappend() */


/* Without changing idx[2*min] and idx[2*max], solve the memory */
/* conflict at idx[2*i+1], idx[2*i+2]. You must be sure there   */
/* is conflict idx[2*i+1]==idx[2*i+2] before calling IMAKESPACE */
/* IMAKESPACE (slower than Imakespace) preserves list[] order   */
void IMAKESPACE (int idx[], int list[], int i, int min, int max)
{
    int j, k, a, left_shift; 
    /* look for left free space */
    for (j=i; (j>min)&&(idx[2*j]==idx[2*j-1]); j--);
    /* look for right free space */
    for (k=i+1; (k<max)&&(idx[2*k+1]==idx[2*k+2]); k++);
    if ((j > min) && (k < max))
	left_shift = /* to prevent systematic squeeze */
	     (idx[2*i+1]-idx[2*j] <  idx[2*k+1]-idx[2*i+2]) ||
	    ((idx[2*i+1]-idx[2*j] == idx[2*k+1]-idx[2*i+2]) &&
	     (Frandom() < 0.5));
    else if ((j > min) && (k >= max)) left_shift = 1;
    else if ((j <= min) && (k < max)) left_shift = 0;
    else pe ("IMAKESPACE: min=%d max=%d jammed between "
             "i=%d and i+1.\n", min, max, i);
    /* lazy: shift only by 1 unit */
    if (left_shift)
    {
	for (a=idx[2*j]; a<idx[2*i+1]; a++) list[a-1] = list[a];
	for (a=2*j; a<=2*i+1; a++) idx[a]--;
    }
    else
    {
	for (a=idx[2*k+1]; a>idx[2*i+2]; a--) list[a] = list[a-1];
	for (a=2*k+1; a>=2*i+2; a--) idx[a]++;
    }
    return;
} /* end IMAKESPACE() */


/* safely append value to the end of owner i's (ordered) list */
void IAPPEND (int idx[], int list[], int value, int i, int min, int max)
{
    int conflict_count = 0;
    if (idx[2*i+1] < idx[2*i+2])
    {
	list[idx[2*i+1]++] = value;
	return;
    }
    else if (idx[2*i+1] == idx[2*i+2])
    {
	conflict_count++;
	IMAKESPACE (idx, list, i, min, max);
	list[idx[2*i+1]++] = value;
	return;
    }
    pe ("IAPPEND: detects overflow at i=%d and i+1,\n"
        "min=%d, max=%d, conflict_count=%d\n",
        i, min, max, conflict_count);
    return;
} /* end IAPPEND() */
