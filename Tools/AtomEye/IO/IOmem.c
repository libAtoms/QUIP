/***************************************/
/* libIO:                              */
/*                                     */
/* Strings, parsers & files            */
/*                                     */
/* Dec.12, 1999 Ju Li <liju99@mit.edu> */
/***************************************/

#include "IO.h"


/* return pointer to an array of n chars */
char *IOalloc (int n)
{
    void *ptr;
    ptr = malloc(n*sizeof(char));
    if (ptr == NULL)
	pe("IOalloc: attempt to malloc %d "
           "chars, or %ld bytes, failed", n,
           (long)(n*sizeof(char)) );
    return((char *)ptr);
} /* end IOalloc() */


/* return pointer to an array of n chars which are all cleared to 0 */
char *IOALLOC (int n)
{
    void *ptr;
    ptr = calloc(n, sizeof(char));
    if (ptr == NULL)
	pe("IOALLOC: attempt to calloc %d "
           "chars, or %ld bytes, failed", n,
           (long)(n*sizeof(char)) );
    return((char *)ptr);
} /* end IOALLOC() */


/* clone a string using malloc() memory. Remember to free it with free()! */
char *IOClone (char *str)
{
    char *c, *d;
    for (c=str; *c!=EOS; c++);
    d = IOalloc(c-str+1);
    for (c=str; *c!=EOS; c++) d[c-str] = *c;
    d[c-str] = EOS;
    return (d);
} /* end IOClone() */


/* fast IO (static.. dynamic) memory manager for string manipulation */

static char Pool[IOMEM_SEGS][IOMEM_UNIT];
static bool this_seg_allocated[IOMEM_SEGS]={FALSE};
static char *client[IOMEM_MAX_CLIENTS]={NULL};
static int client_segs[IOMEM_MAX_CLIENTS];

char *IOmem (size_t size_in_chars)
{
    register int i, j, k, slot;
    for (slot=0; slot<IOMEM_MAX_CLIENTS; slot++)
	if (client[slot] == NULL) break;
    if ( slot == IOMEM_MAX_CLIENTS )
	pe("IOmem: IOMEM_MAX_CLIENTS = %d "
           "reached, no open slots\n", IOMEM_MAX_CLIENTS);
    for (i=0; i<IOMEM_SEGS; i=j+1)
	for (j=i; (j<IOMEM_SEGS)&&(!this_seg_allocated[j]); j++)
	    if ( j-i+1 >= size_in_chars / IOMEM_UNIT )
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
    client[slot] = (char *)malloc(size_in_chars*sizeof(char));
    client_segs[slot] = 0;
exit:
    return (client[slot]);
} /* end IOmem() */


void IOfree (char *ptr)
{
    register int i;
    register bool *p, *q;
    if (ptr == NULL) return;
    for (i=0; i<IOMEM_MAX_CLIENTS; i++)
	if (ptr == client[i])
	{
	    if ( client_segs[i] )
	    {
		/* unroll the innermost k loop: */
		/* for (k=0; k<client_segs[i]; k++) */
		/* this_seg_allocated[(ptr-&Pool[0][0])/IOMEM_UNIT+k] = 0; */
		p = this_seg_allocated + (ptr-&Pool[0][0]) / IOMEM_UNIT;
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
} /* end IOfree() */


void IOfreeall()
{
    register int i;
    for (i=0; i<IOMEM_SEGS; i++)
	this_seg_allocated[i] = FALSE;
    for (i=0; i<IOMEM_MAX_CLIENTS; i++)
	if (client[i] != NULL)
	{
	    if (!client_segs[i]) free ((void *)client[i]);
	    client[i] = NULL;
	}
    return;
} /* end IOfreeall() */


/* Clone a string using IOmem() memory. Remember to free it with IOfree()! */
char *IOclone (char *str)
{
    char *c, *d;
    for (c=str; *c!=EOS; c++);
    d = IOmem(c-str+1);
    for (c=str; *c!=EOS; c++) d[c-str] = *c;
    d[c-str] = EOS;
    return (d);
} /* end IOclone() */


#ifdef _IOmem_TEST
#define TRIALS  1000000
#define CHARS   IOMEM_UNIT
#define TIMES   IOMEM_SEGS/2
#include <Timer.h>
void main()
{
    int i, j;
    char *a[TIMES];

    printf ("\nIOmem/IOfree configuration:\n");
    printf ("static unit: IOMEM_UNIT = %d [chars] (%d bytes)\n",
	    IOMEM_UNIT, IOMEM_UNIT*sizeof(char));
    printf ("number of static units: IOMEM_SEGS = %d\n", IOMEM_SEGS);
    printf ("maximal static & dynamic clients: IOMEM_MAX_CLIENTS = %d\n\n",
	    IOMEM_MAX_CLIENTS);
    
    start_chronometer();
    for (i=0; i<TRIALS; i++)
    {
	for (j=0; j<TIMES; j++) a[j] = NULL;
	for (j=0; j<TIMES; j++) IOfree(a[j]);
    }
    stop_chronometer();
    printf ("empty loop: %d times %d [%d chars] took %s\n",
	    TRIALS, TIMES, CHARS, earthtime(stopped_usertime()));
    
    start_chronometer();
    for (i=0; i<TRIALS; i++)
    {
	for (j=0; j<TIMES; j++) a[j] = IOmem(CHARS);
	for (j=0; j<TIMES; j++) IOfree(a[j]);
    }
    stop_chronometer();
    printf ("IOmem/IOfree: %d times %d [%d chars] took %s\n",
	    TRIALS, TIMES, CHARS, earthtime(stopped_usertime()));
    
    start_chronometer();
    for (i=0; i<TRIALS; i++)
    {
	for (j=0; j<TIMES; j++) a[j] = IOmem(CHARS);
	IOfreeall();
    }
    stop_chronometer();
    printf ("IOmem/IOfreeall: %d times %d [%d chars] took %s\n",
	    TRIALS, TIMES, CHARS, earthtime(stopped_usertime()));
    
    start_chronometer();
    for (i=0; i<TRIALS; i++)
    {
	for (j=0; j<TIMES; j++) a[j] = (char *)malloc(CHARS*sizeof(char));
	for (j=0; j<TIMES; j++) free(a[j]);
    }
    stop_chronometer();
    printf ("malloc/free: %d times %d [%d chars] took %s\n\n",
	    TRIALS, TIMES, CHARS, earthtime(stopped_usertime()));
}
#endif  /* _IOmem_TEST */
