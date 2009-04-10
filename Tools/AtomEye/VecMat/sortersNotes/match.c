/*****************************************/
/* Fastest Double Array Sorting Routines */
/*****************************************/

#include "match.h"

extern void merge(QUICKSORT_INPUT_ARGS);
extern void wright(QUICKSORT_INPUT_ARGS);
extern void nr(QUICKSORT_INPUT_ARGS);
extern void jp(QUICKSORT_INPUT_ARGS);
extern void anl(QUICKSORT_INPUT_ARGS);
extern void glibc(QUICKSORT_INPUT_ARGS);
struct List
{
    char *name;
    void (*fun)(QUICKSORT_INPUT_ARGS);
};
struct List funs[] =
{
    {"Glibc qsort", &glibc},
    {"Merge sort", &merge},
    {"Joe Wright", &wright},
    {"Numerical Recipes", &nr},
    {"Japanese (Akira Wada)", &jp},
    {"ANL Petsc", &anl},
};
#define NUM (sizeof(funs)/sizeof(struct List))

int main()
{
    static DataType x[DIM];
    static IdxType idx[DIM];
    IdxType i, j;
    printf ("%ld random numbers are generated...\n", (long)DIM);
    /* for (i=0; i<DIM; i++) x[i] = (DataType)FRANDOM()*10000.; */
    for (i=0; i<DIM; i++) x[i] = (DataType) i + 1.5 * FRANDOM();
    /* for (i=0; i<DIM; i++) x[i] = -i; */
    for (i=0; i<NUM; i++)
    {
	printf ("Testing %s..", funs[i].name);
	start_chronometer();
	funs[i].fun((IdxType)DIM, x, idx);
	stop_chronometer();
	/* check correctness */
	for (j=0; (j<DIM-1)&&(x[idx[j]]<=x[idx[j+1]]); j++);
	if (j==DIM-1) printf (" %lf seconds\n", stopped_usertime());
	else printf (" it gave wrong order\n");
    }
    return(0);
} /* end main() */


/*
  1. R.Sedgewick, Algorithms in C++, & the edition translated into japanese,
  2. D.E.Knuth,  The Art of Computer Programming, Vol.3 (the 2'nd printing),
  3. ftp://ring.aist.go.jp/pub/TeX/CTAN/indexing/makeindex/src/qsort.c  
  4. ftp://ftp.sra.co.jp/pub/cmd/postgres/6.2.1/postgresql-6.2.1/src/backend/lib/qsort.c
  5. ftp://ftp.jp.freebsd.org/pub/FreeBSD/FreeBSD-current/src/sys/libkern/qsort.c
  6. ftp://ftp.nizkor.vex.net/local/src/darcy/common/qsort.c
  7. http://www.tokuyama.ac.jp/home/~kawamura/ qs6.c & qs7.c
  8. ftp://ftp.groggs.group.cam.ac.uk/pub/software/qsort.c-1.12
  9. ftp://ftp.linux.org.uk/pub/linux/libc/libc-5.4.44.tar.gz
        /libc/stdlib/_quicksort.c  /libc/stdlib/qsort.c
  10. http://www.mars.dti.ne.jp/~a-wada/qsortlib.html
  */
