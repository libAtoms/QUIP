/*********************************************/
/* libVecMat: -llapck -lblas -lm             */
/*            -lScalar -lIO                  */
/*                                           */
/* General Vector & Matrix Operation Library */
/*                                           */
/* Nov 19 1999 Ju Li <liju99@mit.edu>        */
/*********************************************/

#include "VecMat.h"


/*******************************************/
/* Matrix scale memory, a.k.a. mega-memory */
/*                                         */
/* scenarios expected:                     */
/*                                         */
/* a) Hermitian matrix diagonalization.    */
/* b) Complex matrix inversion.            */
/* c) Huge array merge sort indexing.      */
/*******************************************/

/* fast matrix (static.. dynamic) memory manager */

static double Pool[MMEM_SEGS][MMEM_UNIT];
static bool this_seg_allocated[MMEM_SEGS]={FALSE};
static double *client[MMEM_MAX_CLIENTS]={NULL};
static int client_segs[MMEM_MAX_CLIENTS];

double *Mmem (size_t size_in_dbles)
{
    register int i, j, k, slot;
    for (slot=0; slot<MMEM_MAX_CLIENTS; slot++)
	if (client[slot] == NULL) break;
    if ( slot == MMEM_MAX_CLIENTS )
    {
	printf ("error: Mmem: MMEM_MAX_CLIENTS = %d "
		"reached, no open slots\n", MMEM_MAX_CLIENTS);
	exit(1);
    }
    for (i=0; i<MMEM_SEGS; i=j+1)
	for (j=i; (j<MMEM_SEGS)&&(!this_seg_allocated[j]); j++)
	    if ( j-i+1 >= size_in_dbles / MMEM_UNIT )
	    {
		client[slot] = &Pool[i][0];
		client_segs[slot] = j-i+1;
		/* unroll the innermost k loop: */
		/* for (k=i; k<=j; k++) this_seg_allocated[k] = TRUE; */
		this_seg_allocated[i] = TRUE; if (i==j) goto exit; k=i+1;
	    loop:
		this_seg_allocated[k] = TRUE; if (k==j) goto exit; k++;
		this_seg_allocated[k] = TRUE; if (k==j) goto exit; k++;
		this_seg_allocated[k] = TRUE; if (k==j) goto exit; k++;
		this_seg_allocated[k] = TRUE; if (k==j) goto exit; k++;
		this_seg_allocated[k] = TRUE; if (k==j) goto exit; k++;
		this_seg_allocated[k] = TRUE; if (k==j) goto exit; k++;
		this_seg_allocated[k] = TRUE; if (k==j) goto exit; k++;
		this_seg_allocated[k] = TRUE; if (k==j) goto exit; k++;
		goto loop;
	    }
    client[slot] = (double *)malloc(size_in_dbles*sizeof(double));
    client_segs[slot] = 0;
exit:
    return (client[slot]);
} /* end Mmem() */


void Mfree (double *ptr)
{
    register int i;
    register bool *p, *q;
    if (ptr == NULL) return;
    for (i=0; i<MMEM_MAX_CLIENTS; i++)
	if (ptr == client[i])
	{
	    if ( client_segs[i] )
	    {
		/* unroll the innermost k loop: */
		/* for (k=0; k<client_segs[i]; k++) */
		/* this_seg_allocated[(ptr-&Pool[0][0])/MMEM_UNIT+k] = 0; */
		p = this_seg_allocated + (ptr-&Pool[0][0]) / MMEM_UNIT;
		q = p + client_segs[i] - 1;
	    loop:
		*p = FALSE; if (p==q) goto exit; p++;
		*p = FALSE; if (p==q) goto exit; p++;
		*p = FALSE; if (p==q) goto exit; p++;
		*p = FALSE; if (p==q) goto exit; p++;
		*p = FALSE; if (p==q) goto exit; p++;
		*p = FALSE; if (p==q) goto exit; p++;
		*p = FALSE; if (p==q) goto exit; p++;
		*p = FALSE; if (p==q) goto exit; p++;
		goto loop;
	    }
	    else free ((void *)client[i]);
	exit:
	    client[i] = NULL;
	    return;
	}
} /* end Mfree() */


void Mfreeall()
{
    register int i;
    for (i=0; i<MMEM_SEGS; i++)
	this_seg_allocated[i] = FALSE;
    for (i=0; i<MMEM_MAX_CLIENTS; i++)
	if (client[i] != NULL)
	{
	    if (!client_segs[i]) free ((void *)client[i]);
	    client[i] = NULL;
	}
    return;
} /* end Mfreeall() */


/* clone an array of size_in_dbles using Mmem() memory */
/* remember to free it later with Mfree()              */
double *Mclone (double *ptr, size_t size_in_dbles)
{
    double *c = Mmem(size_in_dbles);
    memcpy((void *)c, (void *)ptr, size_in_dbles*sizeof(double));
    return (c);
} /* end Mclone() */


#ifdef _Mmem_TEST
#define TRIALS  100000
#define DBLES   MMEM_UNIT
#define TIMES   MMEM_SEGS
#include <Timer.h>
void main()
{
    int i, j;
    double *a[TIMES];

    printf ("\nMmem/Mfree configuration:\n");
    printf ("static unit: MMEM_UNIT = %d [dbles] (%ld bytes)\n",
	    MMEM_UNIT, MMEM_UNIT*sizeof(double));
    printf ("number of static units: MMEM_SEGS = %d\n", MMEM_SEGS);
    printf ("maximal static & dynamic clients: MMEM_MAX_CLIENTS = %d\n\n",
	    MMEM_MAX_CLIENTS);
    
    start_chronometer();
    for (i=0; i<TRIALS; i++)
    {
	for (j=0; j<TIMES; j++) a[j] = NULL;
	for (j=0; j<TIMES; j++) Mfree(a[j]);
    }
    stop_chronometer();
    printf ("empty loop: %d times %d [%d dbles] took %s\n",
	    TRIALS, TIMES, DBLES, earthtime(stopped_usertime()));
    
    start_chronometer();
    for (i=0; i<TRIALS; i++)
    {
	for (j=0; j<TIMES; j++) a[j] = Mmem(DBLES);
	for (j=0; j<TIMES; j++) Mfree(a[j]);
    }
    stop_chronometer();
    printf ("Mmem/Mfree: %d times %d [%d dbles] took %s\n",
	    TRIALS, TIMES, DBLES, earthtime(stopped_usertime()));
    
    start_chronometer();
    for (i=0; i<TRIALS; i++)
    {
	for (j=0; j<TIMES; j++) a[j] = Mmem(DBLES);
	Mfreeall();
    }
    stop_chronometer();
    printf ("Mmem/Mfreeall: %d times %d [%d dbles] took %s\n",
	    TRIALS, TIMES, DBLES, earthtime(stopped_usertime()));
    
    start_chronometer();
    for (i=0; i<TRIALS; i++)
    {
	for (j=0; j<TIMES; j++) a[j] = (double *)malloc(DBLES*sizeof(double));
	for (j=0; j<TIMES; j++) free(a[j]);
    }
    stop_chronometer();
    printf ("malloc/free: %d times %d [%d dbles] took %s\n\n",
	    TRIALS, TIMES, DBLES, earthtime(stopped_usertime()));
}
#endif
