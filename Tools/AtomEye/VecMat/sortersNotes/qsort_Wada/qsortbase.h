/*  a-wada/qsortlib/ccn/qsortbase.h                             Mar. 12, 1999 */
/*  Copyright (c) 1999    Akira Wada <a-wada@mars.dti.ne.jp>                  */

#define QSORTBASE(qsortname)                                                   \
void qsortname (ARGS) {                                                  /*01*/\
    PTRTYP i, j, k, l, m, n, r; int psz, thr; DFWRK;              DFDGN; /*02*/\
    if (VALID) {                                                         /*03*/\
        HSKP; for ( ; ; ) {                                              /*04*/\
            m = l + ((psz = (r - l) / ESZ) / 2) * ESZ,                   /*05*/\
            i = l, j = r, thr = MDT, psz -= (MDA);                AVDGN; /*06*/\
            for ( ; ; ) {                                                /*07*/\
                if (ORDR (m, j)) {                                       /*08*/\
                    if (i < m  && !ORDR (i, m)) {                 RSTCS; /*09*/\
                        SWAP (i, m); FINDXC (MIN, m, j, r);}}            /*10*/\
                else {                                            RSTCS; /*11*/\
                    if (i >= m || !ORDR (i, m)) SWAP (i, j); else {      /*12*/\
                        SWAP (m, j); FINDXC (MAX, m, l, i);}}            /*13*/\
                i += ESZ, j -= ESZ;                                      /*14*/\
                if (psz > thr) thr *= MDM; else break;}           CKSRT; /*15*/\
            if (i < j) {                                                 /*16*/\
                for ( ; ; ) {                                            /*17*/\
                    while (i < m && ODNE (i, m)) i += ESZ;               /*18*/\
                    while (m < j && ODNE (m, j)) j -= ESZ;               /*19*/\
                    if    (i < j)   SWAP (i, j); else break;             /*20*/\
                    if (m == j) j -= ESZ, m = i; else                    /*21*/\
                    if (i == m) i += ESZ, m = j; else                    /*22*/\
                                j -= ESZ, i += ESZ;}              CKDGN; /*23*/\
                if ((i -= ESZ) - l < r - (j += ESZ)) {                   /*24*/\
                    if (l >= i)  l = j; else PUSH (j, r), r = i;}        /*25*/\
                else {                                                   /*26*/\
                    if (j >= r)  r = i; else PUSH (l, i), l = j;}}       /*27*/\
            else {                                                       /*28*/\
                if (!EMPTY) POP (r, l); else break;}}                    /*29*/\
        CLNUP;}}                                                         /*30*/

#define FINDXC(minmax, m, s, e)   n = m, k = s; do                             \
                if (!ORDR minmax) n = k; while ((k += ESZ) <= e);              \
                if (n != m) SWAP (m, n)
#define MIN               (n, k)
#define MAX               (k, n)
#define ODNE(i, j)  !ORDR (j, i)                    /* in order but not equal */

#define MDT     16  /* threshold for median of the 5 or the more,    MDT >= 8 */
#define MDA      2  /* adjusting threshold slightly,           MDT + MDA >= 4 */
#define MDM      2  /* multiplier for threshold,                     MDM >= 2 */

/*  for prevention from the degeneration caused by peculiar input  */
#define DGT    256  /* threshold of p-size for checking degeneration          */
#define DGV     16  /* degeneration assumed if the smaller < (p-size / DGV)   */
#define DGM      0  /* threshold for selecting no median                      */
#define DGS      2  /* threshold for samples distributed uniformly            */
#define DGU      4  /* threshold for sampling at random                       */
#define DFDGN   int dgn = 0, nsp, itv, cst
#define CKDGN   if (psz >= DGT && (m < l + (psz / DGV) * ESZ   ||              \
                                   m > r - (psz / DGV) * ESZ)) dgn++
#define AVDGN   cst = 1;                                                       \
                if (psz <= MDT || dgn <= DGM || dgn > DGS) {                   \
                    if (psz >  MDT && dgn > DGS) do {                          \
                        if (dgn <= DGU) {                                      \
                            nsp = 3; while (psz > thr) thr *= MDM, nsp += 2;   \
                            itv = (psz / (nsp - 1)) * ESZ, k = l, n = r;       \
                            for (thr = MDT; psz > thr; thr *= MDM) {           \
                                i += ESZ, k += itv; SWAP (i, k);               \
                                j -= ESZ, n -= itv; SWAP (n, j);}}             \
                        else {                                                 \
                            if (dgn == DGU + 1) dgn++, srandom (time (0));     \
                            for (thr /= MDM; psz > thr; thr *= MDM) {          \
                                SMPLRDM (i, i, j); i += ESZ;                   \
                                SMPLRDM (j, i, j); j -= ESZ;}                  \
                            SMPLRDM (m, i, j);}                                \
                        i = l, j = r, thr = MDT;} while (0)
#define SMPLRDM(m, s, e) if (m != (n = s + ESZ * (                             \
                    random () % ((e - (s)) / ESZ + 1)))) SWAP (m, n)
#define RSTCS   cst = 0
#define CKSRT } if (psz > MDT && cst != 0) do {                                \
                    k = l, n = l + ESZ;                                        \
                    while (k < r && ORDR (k, n)) k += ESZ, n += ESZ;           \
                    if (k >= r) i = j; else if (k >= m) i = m;} while (0)
#include <time.h>
#include <stdlib.h>

/*  the macros etc. below should be defined in the module calling this.
        ARGS, PTRTYP, DEFWRK, VALID, HSKP, ESZ,
        ORDR(i, j), SWAP(i, j), PUSH(x, y), EMPTY, POP(y, x),
        CLNUP, "qsortname", and the others required  <see e.g. qsorta.c etc.> */

/*  a-wada/qsortlib/ccn/qsortbase.h end  */
