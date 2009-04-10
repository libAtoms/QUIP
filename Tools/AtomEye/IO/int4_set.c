/****************************************/
/* libIO:                               */
/*                                      */
/* Strings, parsers & files             */
/*                                      */
/* Dec.12, 1999  Ju Li <liju99@mit.edu> */
/****************************************/

#include "IO.h"
#include <Timer.h>

/********************************************************/
/* Duke it out with memset/bzero: and supporting longer */
/* period patterns like int2, int4 and int8.            */
/********************************************************/

#define REGS          4U
#define REGS_SHIFT    2
#define REGS_MASK    (REGS-1)
#define REGS_NASK    (~REGS_MASK)
#define REGTYPE       int4             /* DWORD */
#define BYTES         DWORD_BYTES
#undef BYTESHIFT
#define BYTESHIFT     DWORD_BYTESHIFT
#define BYTEMASK      DWORD_BYTEMASK
#define BYTENASK      DWORD_BYTENASK
#define INT2S        (BYTES/2)
#define INT2SHIFT    (BYTESHIFT-1)
#define INT2MASK     (BYTEMASK>>1)
#define INT2NASK     (~INT2MASK)
#define TESTBYTES    (1000000000)
#define TASKBYTES    (1000000/4)

/****************************************************************/
/* int4_aset() family requires destination pointer be aligned   */
/* with DWORD; however after adding the "count" it does not     */
/* have to, as it requires very little work to finish the tail. */
/* "int4" also suggests internal register implementation.       */
/****************************************************************/

/* Copy "char_count" x (char)"c" to "dest". "dest" must be aligned to DWORD */
void char_int4_aset (char *dest, char c, size_t char_count)
{
    register REGTYPE a0,a1,a2,a3;  /* REGS */
    register REGTYPE *q;
    register size_t i;

    a0  = (unsigned char)c;  a0 <<= 8 * sizeof(char);
    a0 |= (unsigned char)c;  a0 <<= 8 * sizeof(char);
    a0 |= (unsigned char)c;  a0 <<= 8 * sizeof(char);
    a0 |= (unsigned char)c;
    a1 = a2 = a3 = a0;
    
    /* there SHOULD be no misalignment */
    q = ((REGTYPE *) dest);
    /* number of register trains */
    i = char_count >> (BYTESHIFT + REGS_SHIFT);

#define PLAN_A() REPEAT64( if (i==0) goto exit1; q[0]=a0;q[1]=a1;q[2]=a2;q[3]=a3; if (i==1) goto exit1; q[4]=a0;q[5]=a1;q[6]=a2;q[7]=a3; if (i==2) goto exit1;q[8]=a0;q[9]=a1;q[10]=a2;q[11]=a3; if (i==3) goto exit1; q[12]=a0;q[13]=a1;q[14]=a2;q[15]=a3; if (i==4) goto exit1; q[16]=a0;q[17]=a1;q[18]=a2;q[19]=a3; if (i==5) goto exit1; q[20]=a0;q[21]=a1;q[22]=a2;q[23]=a3; if (i==6) goto exit1; \
q[24]=a0;q[25]=a1;q[26]=a2;q[27]=a3; if (i==7) goto exit1; q[28]=a0;q[29]=a1;q[30]=a2;q[31]=a3; i-=8; q+=32; )

#define PLAN_B() REPEAT64( if (i==0) goto exit1; q[0]=a0;q[1]=a1;q[2]=a2;q[3]=a3; if (i==1) goto exit1; q[4]=a0;q[5]=a1;q[6]=a2;q[7]=a3; if (i==2) goto exit1; q[8]=a0;q[9]=a1;q[10]=a2;q[11]=a3; if (i==3) goto exit1; q[12]=a0;q[13]=a1;q[14]=a2;q[15]=a3; i-=4; q+=16; )
    
  loop1:
    PLAN_B();
    goto loop1;
  exit1:

    /* what's left are still int4, but cannot fill a train */
    i = (char_count>>BYTESHIFT)&REGS_MASK;
    if (i==0) goto exit2;
    /* clean up the int4 tail: REGS-1 druple */
    q = ((REGTYPE *)dest) + ((char_count>>BYTESHIFT)&REGS_NASK);
    q[0] = a0;
    if (i==1) goto exit2; q[1] = a1;
    if (i==REGS-1) q[REGS-2] = a2;
  exit2:

    /* what's left are only chars */
    i = char_count & BYTEMASK;
    /* clean up the char tail: BYTES-1 druple */
    if (i==0) goto exit3; dest += (char_count & BYTENASK);
    dest[0] = c;
    if (i==1) goto exit3; dest[1] = c;
    if (i==BYTES-1) dest[BYTES-2] = c;
  exit3:
    return;
} /* end char_int4_aset() */


/* Copy "char_count" x (char)"c" to "dest" using four int4 registers */
void char_int4_set (char *dest, char c, size_t char_count)
{
    register unsigned long i;
    /* how many bytes are in front of DWORD alignment */
    i = ((unsigned long)dest) & BYTEMASK;
    if (i==0)
    {
        char_int4_aset (dest, c, char_count);
        return;
    }
    i = BYTES-i; if (i > char_count) i = char_count;
    if (i==0) return;    dest[0] = c;
    if (i==1) goto exit; dest[1] = c;
    if (i==BYTES-1) dest[BYTES-2] = c;
  exit:
    if (i == char_count) return;
    char_int4_aset (dest+i, c, char_count-i);
    return;
} /* end char_int4_set() */


#ifdef _char_int4_TEST
#define OP   (TESTBYTES/TASKBYTES)
#define N    (TASKBYTES/sizeof(char))
#define c    ((char)255)
int main (int argc, char *argv[])
{
    static char BUFFER[N+2]={0};
    char *a = BUFFER+1;
    register int i, j;

    memset (a, c^1, sizeof(char)*N); /* load it into L2 cache */

    start_chronometer();
    for (i=0; i<OP; i++)
    {
        for (j=0; j<N; j++) a[j] = c;
        charPfool(a);
    }
    stop_chronometer();
    printf ("simple approach: %g s -> %.3f MB/s\n", stopped_usertime(),
            1e-6*N*OP*sizeof(char)/stopped_usertime() );

    memset (a, c^1, sizeof(char)*N);
    start_chronometer();
    for (i=0; i<OP; i++)
    {
        char_int4_set (a, c, N);
        charPfool(a);
    }
    stop_chronometer();
    printf ("char_int4_set: %g s -> %.3f MB/s\n", stopped_usertime(),
            1e-6*N*OP*sizeof(char)/stopped_usertime() );

    start_chronometer();
    if (a[-1] != 0) pe ("a[%d] = %d wrong\n", -1, a[-1]);
    for (i=0; i<N; i++)
        if (a[i] != c) printf ("a[%d] = %d wrong\n", i, a[i]);
    if (a[i] != 0) pe ("a[%d] = %d wrong\n", i, a[i]);
    stop_chronometer();
    printf ("successful... checking took %g s -> %.3f MB/s\n",
            stopped_usertime(),
            1e-6*N*sizeof(unsigned char)/stopped_usertime());
    
    start_chronometer();
    for (i=0; i<OP; i++)
    {
        memset (a, c, N);
        charPfool(a);
    }
    stop_chronometer();
    printf ("memset: %g s -> %.3f MB/s\n", stopped_usertime(),
            1e-6*N*OP*sizeof(char)/stopped_usertime() );
    return (1);
}
#undef c
#undef N
#undef OP
#endif /* _char_int4_TEST */


/* Copy "uchar_count" x (unsigned char)"c" to */
/* "dest". "dest" must be aligned to DWORD.   */
void uchar_int4_aset (unsigned char *dest, unsigned char c, size_t uchar_count)
{
    register REGTYPE a0,a1,a2,a3;  /* REGS */
    register REGTYPE *q;
    register size_t i;

    a0  = c;  a0 <<= 8 * sizeof(unsigned char);
    a0 |= c;  a0 <<= 8 * sizeof(unsigned char);
    a0 |= c;  a0 <<= 8 * sizeof(unsigned char);
    a0 |= c;
    a1 = a2 = a3 = a0;
    
    /* there SHOULD be no misalignment */
    q = ((REGTYPE *) dest);
    /* number of register trains */
    i = uchar_count >> (BYTESHIFT + REGS_SHIFT);

  loop1:
    PLAN_B();
    goto loop1;
  exit1:

    /* what's left are still int4, but cannot fill a train */
    i = (uchar_count>>BYTESHIFT)&REGS_MASK;
    if (i==0) goto exit2;
    /* clean up the int4 tail: REGS-1 druple */
    q = ((REGTYPE *)dest) + ((uchar_count>>BYTESHIFT)&REGS_NASK);
    q[0] = a0;
    if (i==1) goto exit2; q[1] = a1;
    if (i==REGS-1) q[REGS-2] = a2;
  exit2:

    /* what's left are only uchars */
    i = uchar_count & BYTEMASK;
    /* clean up the uchar tail: BYTES-1 druple */
    if (i==0) goto exit3; dest += (uchar_count & BYTENASK);
    dest[0] = c;
    if (i==1) goto exit3; dest[1] = c;
    if (i==BYTES-1) dest[BYTES-2] = c;
  exit3:
    return;
} /* end uchar_int4_aset() */


/* Copy "uchar_count" x (unsigned char)"c" */
/* to "dest" using four int4 registers.    */
void uchar_int4_set (unsigned char *dest, unsigned char c, size_t uchar_count)
{
    register unsigned long i;
    /* how many bytes are in front of DWORD alignment */
    i = ((unsigned long)dest) & BYTEMASK;
    if (i==0)
    {
        uchar_int4_aset (dest, c, uchar_count);
        return;
    }
    i = BYTES-i; if (i > uchar_count) i = uchar_count;
    if (i==0) return;    dest[0] = c;
    if (i==1) goto exit; dest[1] = c;
    if (i==BYTES-1) dest[BYTES-2] = c;
  exit:
    if (i == uchar_count) return;
    uchar_int4_aset (dest+i, c, uchar_count-i);
    return;
} /* end uchar_int4_set() */


#ifdef _uchar_int4_TEST
#define OP   (TESTBYTES/TASKBYTES)
#define N    (TASKBYTES/sizeof(unsigned char))
#define c    255
int main (int argc, char *argv[])
{
    static unsigned char BUFFER[N+2]={0};
    unsigned char *a = BUFFER+1;
    register int i, j;

    memset (a, c^1, sizeof(unsigned char)*N); /* load it into L2 cache */

    start_chronometer();
    for (i=0; i<OP; i++)
    {
        for (j=0; j<N; j++) a[j] = c;
        ucharPfool(a);
    }
    stop_chronometer();
    printf ("simple approach: %g s -> %.3f MB/s\n", stopped_usertime(),
            1e-6*N*OP*sizeof(unsigned char)/stopped_usertime() );

    memset (a, c^1, sizeof(unsigned char)*N);
    start_chronometer();
    for (i=0; i<OP; i++)
    {
        uchar_int4_set (a, c, N);
        ucharPfool(a);
    }
    stop_chronometer();
    printf ("uchar_int4_set: %g s -> %.3f MB/s\n", stopped_usertime(),
            1e-6*N*OP*sizeof(unsigned char)/stopped_usertime() );

    start_chronometer();
    if (a[-1] != 0) pe ("a[%d] = %d wrong\n", -1, a[-1]);
    for (i=0; i<N; i++)
        if (a[i] != c) printf ("a[%d] = %d wrong\n", i, a[i]);
    if (a[i] != 0) pe ("a[%d] = %d wrong\n", i, a[i]);
    stop_chronometer();
    printf ("successful... checking took %g s -> %.3f MB/s\n",
            stopped_usertime(),
            1e-6*N*sizeof(unsigned char)/stopped_usertime());
    
    start_chronometer();
    for (i=0; i<OP; i++)
    {
        memset (a, c, N);
        ucharPfool(a);
    }
    stop_chronometer();
    printf ("memset: %g s -> %.3f MB/s\n", stopped_usertime(),
            1e-6*N*OP*sizeof(unsigned char)/stopped_usertime() );
    return (1);
}
#undef c
#undef N
#undef OP
#endif /* _uchar_int4_TEST */


/* Copy "int2_count" x (int2)"c" to "dest". "dest" must be */
/* aligned to DWORD, which is generally a good idea anyway */
void int2_int4_aset (int2 *dest, int2 c, size_t int2_count)
{
    register REGTYPE a0,a1,a2,a3;  /* REGS */
    register REGTYPE *q;
    register size_t i;

    a0  = c;  a0 <<= 8 * sizeof(int2);
    a0 |= c;
    a1 = a2 = a3 = a0;

    /* there SHOULD be no misalignment */
    q = ((REGTYPE *) dest);
    /* number of register trains */
    i = int2_count >> (INT2SHIFT + REGS_SHIFT);

  loop1:
    PLAN_B();
    goto loop1;
  exit1:

    /* what's left are still int4, but cannot fill a train */
    i = (int2_count>>INT2SHIFT)&REGS_MASK;
    if (i==0) goto exit2;
    /* clean up the int4 tail: REGS-1 druple */
    q = ((REGTYPE *)dest) + ((int2_count>>INT2SHIFT)&REGS_NASK);
    q[0] = a0;
    if (i==1) goto exit2; q[1] = a1;
    if (i==REGS-1) q[REGS-2] = a2;
  exit2:

    /* what's left are only int2s */
    if (int2_count & INT2MASK) dest[int2_count & INT2NASK] = c;
    return;
} /* end int2_int4_aset() */


/* Copy "int2_count" of (int2)"c" to "dest" using four int4   */
/* registers. "dest" must minimally be aligned to int2 (WORD) */
void int2_int4_set (int2 *dest, int2 c, size_t int2_count)
{
    register unsigned long i;
    /* how many int2 are in front of DWORD alignment */
    i = (((unsigned long)dest)>>1) & INT2MASK;
    if (i==0)
    {
        int2_int4_aset (dest, c, int2_count);
        return;
    }
    if (int2_count > 0)
    {
        dest[0] = c;
        int2_int4_aset (dest+1, c, int2_count-1);
    }
    return;
} /* end int2_int4_set() */


#ifdef _int2_int4_TEST
#define OP   (TESTBYTES/TASKBYTES)
#define N    (TASKBYTES/sizeof(int2))
#define c    255
int main (int argc, char *argv[])
{
    static int2 BUFFER[N+2]={0};
    int2 *a = BUFFER+1;
    register int i, j;
    
    memset (a, c^1, sizeof(int2)*N); /* load it into L2 cache */

    start_chronometer();
    for (i=0; i<OP; i++)
    {
        for (j=0; j<N; j++) a[j] = c;
        int2Pfool(a);
    }
    stop_chronometer();
    printf ("simple approach: %g s -> %.3f MB/s\n", stopped_usertime(),
            1e-6*N*OP*sizeof(int2)/stopped_usertime() );

    memset (a, c^1, sizeof(int2)*N);
    start_chronometer();
    for (i=0; i<OP; i++)
    {
        int2_int4_set (a, c, N);
        int2Pfool(a);
    }
    stop_chronometer();
    printf ("int2_int4_set: %g s -> %.3f MB/s\n", stopped_usertime(),
            1e-6*N*OP*sizeof(int2)/stopped_usertime() );

    start_chronometer();
    if (a[-1] != 0) pe ("a[%d] = %d wrong\n", -1, a[-1]);
    for (i=0; i<N; i++)
        if (a[i] != c) printf ("a[%d] = %d wrong\n", i, a[i]);
    if (a[i] != 0) pe ("a[%d] = %d wrong\n", i, a[i]);
    stop_chronometer();
    printf ("successful... checking took %g s -> %.3f MB/s\n",
            stopped_usertime(),
            1e-6*N*sizeof(int2)/stopped_usertime());
    
    start_chronometer();
    for (i=0; i<OP; i++)
    {
        bzero (a, N*sizeof(int2));
        int2Pfool(a);
    }
    stop_chronometer();
    printf ("bzero: %g s -> %.3f MB/s\n", stopped_usertime(),
            1e-6*N*OP*sizeof(int2)/stopped_usertime() );
    return (1);
}
#undef c
#undef N
#undef OP
#endif /* _int2_int4_TEST */


/* Copy "int4_count" x (int4)"c" to "dest" using four int4     */
/* registers. "dest" must minimally be aligned to int4 (DWORD) */
void int4_int4_set (int4 *dest, int4 c, size_t int4_count)
{
    register REGTYPE a0,a1,a2,a3;  /* REGS */
    register REGTYPE *q;
    register size_t i;

    a1 = a2 = a3 = a0 = c;
    /* there SHOULD be no misalignment */
    q = dest;
    /* number of register trains */
    i = int4_count >> REGS_SHIFT;
  loop1:
    PLAN_B();
    goto loop1;
  exit1:
    /* what's left are still int4, but cannot fill a train */
    i = int4_count & REGS_MASK;
    if (i==0) goto exit2;
    /* clean up the int4 tail: REGS-1 druple */
    dest[(int4_count&REGS_NASK)] = a0;
    if (i==1) goto exit2; dest[(int4_count&REGS_NASK)+1] = a1;
    if (i==REGS-1)   dest[(int4_count&REGS_NASK)+REGS-2] = a2;
  exit2:
    return;
} /* end int4_int4_set() */


#ifdef _int4_int4_TEST
#define OP   (TESTBYTES/TASKBYTES)
#define N    (TASKBYTES/sizeof(int4))
#define c    255
int main (int argc, char *argv[])
{
    static int4 BUFFER[N+2]={0};
    int4 *a = BUFFER+1;
    register int i, j;
    
    memset (a, c^1, sizeof(int4)*N); /* load it into L2 cache */

    start_chronometer();
    for (i=0; i<OP; i++)
    {
        for (j=0; j<N; j++) a[j] = c;
        int4Pfool(a);
    }
    stop_chronometer();
    printf ("simple approach: %g s -> %.3f MB/s\n", stopped_usertime(),
            1e-6*N*OP*sizeof(int4)/stopped_usertime() );

    memset (a, c^1, sizeof(int4)*N);
    start_chronometer();
    for (i=0; i<OP; i++)
    {
        int4_int4_set (a, c, N);
        int4Pfool(a);
    }
    stop_chronometer();
    printf ("int4_int4_set: %g s -> %.3f MB/s\n", stopped_usertime(),
            1e-6*N*OP*sizeof(int4)/stopped_usertime() );

    start_chronometer();
    if (a[-1] != 0) pe ("a[%d] = %d wrong\n", -1, a[-1]);
    for (i=0; i<N; i++)
        if (a[i] != c) printf ("a[%d] = %d wrong\n", i, a[i]);
    if (a[i] != 0) pe ("a[%d] = %d wrong\n", i, a[i]);
    stop_chronometer();
    printf ("successful... checking took %g s -> %.3f MB/s\n",
            stopped_usertime(),
            1e-6*N*sizeof(int4)/stopped_usertime());
    
    start_chronometer();
    for (i=0; i<OP; i++)
    {
        bzero (a, N*sizeof(int4));
        int4Pfool(a);
    }
    stop_chronometer();
    printf ("bzero: %g s -> %.3f MB/s\n", stopped_usertime(),
            1e-6*N*OP*sizeof(int4)/stopped_usertime() );
    return (1);
}
#undef c
#undef N
#undef OP
#endif /* _int4_int4_TEST */


/* Copy "int8_count" x (int8)"c" to "dest" using four int4     */
/* registers. "dest" must minimally be aligned to int4 (DWORD) */
void int8_int4_set (int8 *dest, int8 c, size_t int8_count)
{
    static union
    {
        int8 c;
        REGTYPE a[sizeof(int8)/sizeof(REGTYPE)];
    } BUFFER;
    register REGTYPE a0,a1,a2,a3;  /* REGS */
    register REGTYPE *q;
    register size_t i;

    BUFFER.c = c;
    a2 = a0 = BUFFER.a[0]; a3 = a1 = BUFFER.a[1];
    int8_count <<= 1;

    /* there SHOULD be no misalignment */
    q = ((REGTYPE *) dest);
    /* number of register trains */
    i = int8_count >> REGS_SHIFT;

  loop1:
    PLAN_B();
    goto loop1;
  exit1:

    /* what's left are still int8, but cannot fill a train */
    i = int8_count & REGS_MASK;
    if (i==0) goto exit2;
    q = ((REGTYPE *)dest) + (int8_count & REGS_NASK);
    q[0] = a0; q[1] = a1;
  exit2:
    return;
} /* end int8_int4_set() */


#ifdef _int8_int4_TEST
#define OP   (TESTBYTES/TASKBYTES)
#define N    (TASKBYTES/sizeof(int8))
#define c    255
int main (int argc, char *argv[])
{
    static int8 BUFFER[N+2]={0};
    int8 *a = BUFFER+1;
    register int i, j;
    
    memset (a, c^1, sizeof(int8)*N); /* load it into L2 cache */

    start_chronometer();
    for (i=0; i<OP; i++)
    {
        for (j=0; j<N; j++) a[j] = c;
        int8Pfool(a);
    }
    stop_chronometer();
    printf ("simple approach: %g s -> %.3f MB/s\n", stopped_usertime(),
            1e-6*N*OP*sizeof(int8)/stopped_usertime() );

    memset (a, c^1, sizeof(int8)*N);
    start_chronometer();
    for (i=0; i<OP; i++)
    {
        int8_int4_set (a, c, N);
        int8Pfool(a);
    }
    stop_chronometer();
    printf ("int8_int4_set: %g s -> %.3f MB/s\n", stopped_usertime(),
            1e-6*N*OP*sizeof(int8)/stopped_usertime() );

    start_chronometer();
    if (a[-1] != 0) pe ("a[%d] = %ld wrong\n", -1, a[-1]);
    for (i=0; i<N; i++)
        if (a[i] != c) printf ("a[%d] = %ld wrong\n", i, a[i]);
    if (a[i] != 0) pe ("a[%d] = %ld wrong\n", i, a[i]);
    stop_chronometer();
    printf ("successful... checking took %g s -> %.3f MB/s\n",
            stopped_usertime(),
            1e-6*N*sizeof(int8)/stopped_usertime());
    
    start_chronometer();
    for (i=0; i<OP; i++)
    {
        bzero (a, N*sizeof(int8));
        int8Pfool(a);
    }
    stop_chronometer();
    printf ("bzero: %g s -> %.3f MB/s\n", stopped_usertime(),
            1e-6*N*OP*sizeof(int8)/stopped_usertime() );
    return (1);
}
#undef c
#undef N
#undef OP
#endif /* _int8_int4_TEST */

#undef REGS
#undef REGS_SHIFT
#undef REGS_MASK
#undef REGS_NASK
#undef REGTYPE  
#undef BYTES    
#undef BYTESHIFT
#undef BYTEMASK 
#undef BYTENASK 
#undef INT2S
#undef INT2SHIFT
#undef INT2MASK
#undef INT2NASK
#undef TESTBYTES
#undef TASKBYTES
