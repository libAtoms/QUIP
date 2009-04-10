/**************************************/
/* libTimer:                          */
/*                                    */
/* Time and Schedule Manager          */
/*                                    */
/* Nov 19 1999 Ju Li <liju99@mit.edu> */
/**************************************/

#include "Timer.h"


/*************************************************/
/* Function prototypes to fool the compiler into */
/* believing that certain variable is useful --  */
/* for time benchmarking purposes.               */
/*************************************************/


void charfool (char c)
{
    return;
} /* end charfool() */


void ucharfool (unsigned char c)
{
    return;
} /* end ucharfool() */


void shortfool (short i)
{
    return;
} /* end shortfool() */


void ushortfool (unsigned short i)
{
    return;
} /* end ushortfool() */


void intfool (int i)
{
    return;
} /* end intfool() */


void uintfool (unsigned int i)
{
    return;
} /* end uintfool() */


void longfool (long l)
{
    return;
} /* end longfool() */


void ulongfool (unsigned long l)
{
    return;
} /* end ulongfool() */


void longlongfool (long long l)
{
    return;
} /* end longlongfool() */


void ulonglongfool (unsigned long long l)
{
    return;
} /* end ulonglongfool() */


void floatfool (float d)
{
    return;
} /* end floatfool() */


void doublefool (double d)
{
    return;
} /* end doublefool() */


void voidPfool (void *d)
{
    return;
} /* end voidPfool() */


void charPfool (char *c)
{
    return;
} /* end charPfool() */


void ucharPfool (unsigned char *c)
{
    return;
} /* end ucharPfool() */


void shortPfool (short *i)
{
    return;
} /* end shortPfool() */


void ushortPfool (unsigned short *i)
{
    return;
} /* end ushortPfool() */


void intPfool (int *i)
{
    return;
} /* end intPfool() */


void uintPfool (unsigned int *i)
{
    return;
} /* end uintPfool() */


void longPfool (long *l)
{
    return;
} /* end longPfool() */


void ulongPfool (unsigned long *l)
{
    return;
} /* end ulongPfool() */


void longlongPfool (long long *l)
{
    return;
} /* end longlongPfool() */


void ulonglongPfool (unsigned long long *l)
{
    return;
} /* end ulonglongPfool() */


void floatPfool (float *d)
{
    return;
} /* end floatPfool() */


void doublePfool (double *d)
{
    return;
} /* end doublePfool() */
