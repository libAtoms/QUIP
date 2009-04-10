/*********************************************/
/* libVecMat: -llapck -lblas -lm             */
/*            -lScalar -lIO                  */
/*                                           */
/* General Vector & Matrix Operation Library */
/*                                           */
/* Nov 19 1999 Ju Li <liju99@mit.edu>        */
/*********************************************/

#include "VecMat.h"


/* integer memory */

/* return pointer to an array of n ints */
int *Ialloc (int n)
{
    void *ptr;
    ptr = malloc(n*sizeof(int));
    if (ptr == NULL)
    {
	printf ("error: Ialloc: attempt to malloc %d "
		"ints, or %ld bytes, failed", n,
		LONG(n)*sizeof(int));
	exit(1);
    }
    return((int *)ptr);
} /* end Ialloc() */


/* return pointer to an array of n ints which are all cleared to 0 */
int *IALLOC (int n)
{
    void *ptr;
    ptr = calloc(n, sizeof(int));
    if (ptr == NULL)
    {
	printf ("error: IALLOC: attempt to calloc %d "
		"ints, or %ld bytes, failed", n,
		LONG(n)*sizeof(int));
	exit(1);
    }
    return((int *)ptr);
} /* end IALLOC() */


/* return pointer to an array of n ints which are all set to flag */
int *IAlloc (int n, int flag)
{
    int i, *ptr;
    ptr = (int *)malloc(n*sizeof(int));
    if (ptr == NULL)
    {
	printf ("error: IAlloc: attempt to malloc %d ints\n"
		"(set to %d), or %ld bytes, failed", n,
		flag, LONG(n)*sizeof(int));
	exit(1);
    }
    for (i=0; i<n; i++) ptr[i] = flag;
    return(ptr);
} /* end IAlloc() */
