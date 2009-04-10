#include "match.h"
#undef SWAP
#undef MIN
#undef MAX

#define QSNAME      jp
#define ARGS        QUICKSORT_INPUT_ARGS
#define VALID       N > 1
#define HSKP        {l=0; r=N-1; for (i=0; i<N; i++) idx[i] = i;}
#define ESZ         1
#define ORDR(i,j)   (x[idx[i]] <= x[idx[j]])
#define NORDR(i,j)  (x[idx[i]] > x[idx[j]])
#define ODNE(i,j)   (x[idx[i]] < x[idx[j]])  /* !ORDR(j,i) */
#define SWAP(i,j)   (t = idx[i], idx[i] = idx[j], idx[j] = t)
#define CLNUP
#define STKSZ       (8*sizeof(long int))
#define PUSH(x,y)   (*(p++)=x,*(p++)=y)
#define EMPTY       (s >= p)
#define NEMPTY       (s < p)
#define POP(y,x)    (y=*(--p),x=*(--p))

#define QSORTBASE(qsortname)                                                   \
void qsortname (ARGS) {                                                        \
    register IdxType i, j, k, l, t, r, m, n; \
    static IdxType s[STKSZ], *p = s; \
    static int dgn = 0, nsp, itv, cst, psz, thr; \
    if (VALID) {                                                         /*03*/\
        HSKP; for ( ; ; ) {                                              /*04*/\
            m = l + ((psz = (r - l) / ESZ) / 2) * ESZ,                   /*05*/\
            i = l, j = r, thr = MDT, psz -= (MDA);                AVDGN; /*06*/\
            for ( ; ; ) {                                                /*07*/\
                if (ORDR (m, j)) {                                       /*08*/\
                    if (i < m  && NORDR (i, m)) {                 RSTCS; /*09*/\
                        SWAP (i, m); FINDXC (MIN, m, j, r);}}            /*10*/\
                else {                                            RSTCS; /*11*/\
                    if (i >= m || NORDR (i, m)) SWAP (i, j); else {      /*12*/\
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
                if (NEMPTY) POP (r, l); else break;}}                    /*29*/\
        CLNUP;}}                                                         /*30*/

#define FINDXC(minmax, m, s, e)   n = m, k = s; do                             \
                if (NORDR minmax) n = k; while ((k += ESZ) <= e);              \
                if (n != m) SWAP (m, n)
#define MIN     (n, k)
#define MAX     (k, n)
#define MDT     16  /* threshold for median of the 5 or the more,    MDT >= 8 */
#define MDA      2  /* adjusting threshold slightly,           MDT + MDA >= 4 */
#define MDM      2  /* multiplier for threshold,                     MDM >= 2 */
/*  for prevention from the degeneration caused by peculiar input  */
#define DGT    256  /* threshold of p-size for checking degeneration          */
#define DGV     16  /* degeneration assumed if the smaller < (p-size / DGV)   */
#define DGM      0  /* threshold for selecting no median                      */
#define DGS      2  /* threshold for samples distributed uniformly            */
#define DGU      4  /* threshold for sampling at random                       */
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

			
QSORTBASE (QSNAME)
