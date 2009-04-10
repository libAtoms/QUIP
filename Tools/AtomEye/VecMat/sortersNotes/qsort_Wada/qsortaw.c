/*  a-wada/qsortlib/ccn/qsortaw.c                               Mar. 13, 1999 */
/*  Copyright (c) 1999    Akira Wada <a-wada@mars.dti.ne.jp>                  */

/*      <<< Example of QSORTBASE to generate sort subroutine >>>              */
/*          string sort for "radixdemo.c" modified from                       */
/*              http://www.ddj.com/ftp/1998/1998_11/aa118.zip                 */
/*              demo-program of "Radix Quick Sort"                            */
/*              by Jon Bentley and Robert Sedgewick                           */

/** define the object to be sorted **/
typedef char *      typeA;

/** call the sorting algorithm (QSORTBASE) **/
#define QSNAME      qsortaw
#define ARGS        typeA  P[], int N
typedef int         PTRTYP;
#define DFWRK       STACK;  typeA t
#define VALID       N > 1
#define HSKP        l = 0, r = N - 1
#define ESZ         1
#define ORDR(i, j)  (strcmp (P[i], P[j]) <= 0)
#define SWAP(i, j)  t = P[i], P[i] = P[j], P[j] = t
#define CLNUP       /* NOTHING */

#define STACK       PTRTYP s[STKSZ], *p = s
#define STKSZ       sizeof (int) * 8 * 2
#define PUSH(x, y)  *(p++) = x, *(p++) = y
#define EMPTY       (s >= p)
#define POP(y, x)   y = *(--p), x = *(--p)

#include <string.h>
#include "qsortbase.h"

    QSORTBASE (QSNAME)

/*  a-wada/qsortlib/ccn/qsortaw.c end  */
/*============================================================================**
** source codes on web site:
    index        http://www.mars.dti.ne.jp/~a-wada/qsortlib.html
   *radixdemo.c  http://www.mars.dti.ne.jp/~a-wada/qsortlib/ccn/radixdemo.c
    qsortbase.h  http://www.mars.dti.ne.jp/~a-wada/qsortlib/ccn/qsortbase.h
    qsortaw.c    http://www.mars.dti.ne.jp/~a-wada/qsortlib/ccn/qsortaw.c
   *random.c     http://www.mars.dti.ne.jp/~a-wada/qsortlib/ccn/random.c
    qsorts.c     http://www.mars.dti.ne.jp/~a-wada/qsortlib/ccn/qsorts.c
   *benchmark    http://www.mars.dti.ne.jp/~a-wada/qsortlib/ccn/qsaw0205.txt
                 http://www.mars.dti.ne.jp/~a-wada/qsortlib/ccn/qspt0225.txt
**============================================================================*/
