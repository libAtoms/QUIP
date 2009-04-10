/*  a-wada/qsortlib/ccg/qsortav.c                               Nov. 19, 1999 */
/*  Copyright (c) 1999    Akira Wada <a-wada@mars.dti.ne.jp>                  */

/*      <<< Example of QSORTBASE to generate sort subroutine >>>              */
/*          for "shelldrv.c"  modified from                                   */
/*              http://www.cs.princeton.edu/~rs/shell/driver.c                */
/*              comparing Shell-Sorts by R. Sedgewick                         */

/** define the object to be sorted **/
typedef long        typeA;

/** call the sorting algorithm (QSORTBASE) **/
#define QSNAME      qsortav
#define ARGS        typeA   A[], int ll, int rr
typedef int         PTRTYP;
#define DFWRK       STACK;   typeA t
#define VALID       ll < rr
#define HSKP        l = ll, r = rr
#define ESZ         1
#define COMP(i, j)  (A[i] - A[j])
#define SWAP(i, j)  t = A[i], A[i] = A[j], A[j] = t
#define ROTATE(i, j, k) t = A[i], A[i] = A[j], A[j] = A[k], A[k] = t
#define CLNUP       /* NOTHING */

#define STACK       PTRTYP s[STKSZ], *p = s
#define STKSZ       sizeof (int) * 8 * 2
#define PUSH(x, y)  *(p++) = x, *(p++) = y
#define EMPTY       (s >= p)
#define POP(y, x)   y = *(--p), x = *(--p)

#include "qsortbase.h"

    QSORTBASE (QSNAME)

/*  a-wada/qsortlib/ccg/qsortav.c end  */
/*============================================================================**
on web site:   http://www.mars.dti.ne.jp/~a-wada/qsortlib.html
   *shelldrv.c   http://www.mars.dti.ne.jp/~a-wada/qsortlib/ccn/shelldrv.c
    qsortbase.h  http://www.mars.dti.ne.jp/~a-wada/qsortlib/ccg/qsortbase.h
    qsortav.c    http://www.mars.dti.ne.jp/~a-wada/qsortlib/ccg/qsortav.c
   *random.c     http://www.mars.dti.ne.jp/~a-wada/qsortlib/ccn/random.c
   *shell-rs.c   http://www.mars.dti.ne.jp/~a-wada/qsortlib/ccn/shell-rs.c
**============================================================================*/
