/*********************************************/
/* libVecMat: -llapck -lblas -lm             */
/*            -lScalar -lIO                  */
/*                                           */
/* General Vector & Matrix Operation Library */
/*                                           */
/* Nov 19 1999 Ju Li <liju99@mit.edu>        */
/*********************************************/

#include "VecMat.h"


/* double type vector memory allocation and manager */

/* return pointer to an array of n doubles */
double *Valloc (int n)
{
    void *ptr;
    ptr = malloc(n*sizeof(double));
    if (ptr == NULL)
    {
        printf ("error: Valloc: attempt to malloc %d "
                "doubles, or %ld bytes, failed", n,
                LONG(n)*sizeof(double));
        exit(1);
    }
    return((double *)ptr);
} /* end Valloc() */


/* return pointer to an array of n doubles which are all cleared to 0 */
double *VALLOC (int n)
{
    void *ptr;
    ptr = calloc(n, sizeof(double));
    if (ptr == NULL)
    {
        printf ("error: VALLOC: attempt to calloc %d "
                "doubles, or %ld bytes, failed", n,
                LONG(n)*sizeof(double));
        exit(1);
    }
    return((double *)ptr);
} /* end VALLOC() */


/* return pointer to an array of n doubles which are all set to "value" */
double *VAlloc (int n, double value)
{
    double *ptr;
    int i;
    ptr = (double *)malloc(n*sizeof(double));
    if (ptr == NULL)
    {
        printf ("error: VAlloc: attempt to malloc %d doubles\n"
                "(set to %e), or %ld bytes, failed", n,
                value, LONG(n)*sizeof(double));
        exit(1);
    }
    for (i=0; i<n; i++) ptr[i] = value;
    return(ptr);
} /* end VAlloc() */


/* return pointer to an array of n doubles (reallocated) */
double *Vrealloc (double *a, int n)
{
    void *ptr;
    ptr = realloc( (void *)a, n*sizeof(double) );
    if (ptr == NULL)
    {
        printf ("error: Vrealloc: attempt to realloc %d\n"
                "doubles, or %ld bytes, failed", n, LONG(n)*sizeof(double));
        exit(1);
    }
    return((double *)ptr);
} /* end Vrealloc() */


/****************************************************************/
/* (eigen)Vector scale memory, or, sub-mega memory.             */
/*                                                              */
/* scenarios expected:                                          */
/*                                                              */
/* a) single eigenvector from Hermitian matrix diagonalization. */
/* b) small matrices, say, k-space dynamical matrix of quartz.  */
/* c) medium sized bitmap for ordinary tagging purposes.        */
/****************************************************************/

/* fast vector (static.. dynamic) memory manager */

static double Pool[VMEM_SEGS][VMEM_UNIT];
static bool this_seg_allocated[VMEM_SEGS]={FALSE};
static double *client[VMEM_MAX_CLIENTS]={NULL};
static int client_segs[VMEM_MAX_CLIENTS];

double *Vmem (size_t size_in_dbles)
{
    register int i, j, k, slot;
    for (slot=0; slot<VMEM_MAX_CLIENTS; slot++)
        if (client[slot] == NULL) break;
    if ( slot == VMEM_MAX_CLIENTS )
    {
        printf ("error: Vmem: VMEM_MAX_CLIENTS = %d "
                "reached, no open slots\n", VMEM_MAX_CLIENTS);
        exit(1);
    }
    for (i=0; i<VMEM_SEGS; i=j+1)
        for (j=i; (j<VMEM_SEGS)&&(!this_seg_allocated[j]); j++)
            if ( j-i+1 >= size_in_dbles / VMEM_UNIT )
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
} /* end Vmem() */


void Vfree (double *ptr)
{
    register int i;
    register bool *p, *q;
    if (ptr == NULL) return;
    for (i=0; i<VMEM_MAX_CLIENTS; i++)
        if (ptr == client[i])
        {
            if ( client_segs[i] )
            {
                /* unroll the innermost k loop: */
                /* for (k=0; k<client_segs[i]; k++) */
                /* this_seg_allocated[(ptr-&Pool[0][0])/VMEM_UNIT+k] = 0; */
                p = this_seg_allocated + (ptr-&Pool[0][0]) / VMEM_UNIT;
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
} /* end Vfree() */


void Vfreeall()
{
    register int i;
    for (i=0; i<VMEM_SEGS; i++)
        this_seg_allocated[i] = FALSE;
    for (i=0; i<VMEM_MAX_CLIENTS; i++)
        if (client[i] != NULL)
        {
            if (!client_segs[i]) free ((void *)client[i]);
            client[i] = NULL;
        }
    return;
} /* end Vfreeall() */


/* clone an array of size_in_dbles using Vmem() memory */
/* remember to free it later with Vfree()              */
double *Vclone (double *ptr, size_t size_in_dbles)
{
    double *c = Vmem(size_in_dbles);
    memcpy((void *)c, (void *)ptr, size_in_dbles*sizeof(double));
    return (c);
} /* end Vclone() */


#ifdef _Vmem_TEST
#define TRIALS  100000
#define DBLES   VMEM_UNIT
#define TIMES   (VMEM_SEGS/2)
#include <Timer.h>
void main()
{
    int i, j;
    double *a[TIMES];

    printf ("\nVmem/Vfree configuration:\n");
    printf ("static unit: VMEM_UNIT = %d [dbles] (%ld bytes)\n",
            VMEM_UNIT, VMEM_UNIT*sizeof(double));
    printf ("number of static units: VMEM_SEGS = %d\n", VMEM_SEGS);
    printf ("maximal static & dynamic clients: VMEM_MAX_CLIENTS = %d\n\n",
            VMEM_MAX_CLIENTS);
    
    start_chronometer();
    for (i=0; i<TRIALS; i++)
    {
        for (j=0; j<TIMES; j++) a[j] = NULL;
        for (j=0; j<TIMES; j++) Vfree(a[j]);
    }
    stop_chronometer();
    printf ("empty loop: %d times %d [%d dbles] took %s\n",
            TRIALS, TIMES, DBLES, earthtime(stopped_usertime()));
    
    start_chronometer();
    for (i=0; i<TRIALS; i++)
    {
        for (j=0; j<TIMES; j++) a[j] = Vmem(DBLES);
        for (j=0; j<TIMES; j++) Vfree(a[j]);
    }
    stop_chronometer();
    printf ("Vmem/Vfree: %d times %d [%d dbles] took %s\n",
            TRIALS, TIMES, DBLES, earthtime(stopped_usertime()));
    
    start_chronometer();
    for (i=0; i<TRIALS; i++)
    {
        for (j=0; j<TIMES; j++) a[j] = Vmem(DBLES);
        Vfreeall();
    }
    stop_chronometer();
    printf ("Vmem/Vfreeall: %d times %d [%d dbles] took %s\n",
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
