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
/* period patterns like int2, int4, int8, real4, real8. */
/********************************************************/

#define REGS          8U
#define REGS_SHIFT    3
#define REGS_MASK    (REGS-1)
#define REGS_NASK    (~REGS_MASK)
#define REGTYPE       int8             /* QWORD */
#define BYTES         QWORD_BYTES
#undef BYTESHIFT
#define BYTESHIFT     QWORD_BYTESHIFT
#define BYTEMASK      QWORD_BYTEMASK
#define BYTENASK      QWORD_BYTENASK
#define INT2S        (BYTES/2)
#define INT2SHIFT    (BYTESHIFT-1)
#define INT2MASK     (BYTEMASK>>1)
#define INT2NASK     (~INT2MASK)
#define INT4S        (INT2S/2)
#define INT4SHIFT    (INT2SHIFT-1)
#define INT4MASK     (INT2MASK>>1)
#define INT4NASK     (~INT4MASK)
#define TESTBYTES    (1000000000)
#define TASKBYTES    (1000000/2)

/****************************************************************/
/* int8_aset() family requires destination pointer be aligned   */
/* with QWORD; however after adding the "count" it does not     */
/* have to, as it requires very little work to finish the tail. */
/* "int8" also suggests internal register implementation.       */
/****************************************************************/

/* Copy "char_count" x (char)"c" to "dest". "dest" must be aligned to QWORD */
void char_int8_aset (char *dest, char c, size_t char_count)
{
    register REGTYPE a0,a1,a2,a3,a4,a5,a6,a7;  /* REGS */
    register REGTYPE *q;
    register size_t i;

    a0  = (unsigned char)c;  a0 <<= 8 * sizeof(char);
    a0 |= (unsigned char)c;  a0 <<= 8 * sizeof(char);
    a0 |= (unsigned char)c;  a0 <<= 8 * sizeof(char);
    a0 |= (unsigned char)c;  a0 <<= 8 * sizeof(char);
    a0 |= (unsigned char)c;  a0 <<= 8 * sizeof(char);
    a0 |= (unsigned char)c;  a0 <<= 8 * sizeof(char);
    a0 |= (unsigned char)c;  a0 <<= 8 * sizeof(char);
    a0 |= (unsigned char)c;
    a1 = a2 = a3 = a4 = a5 = a6 = a7 = a0;
    
    /* there SHOULD be no misalignment */
    q = ((REGTYPE *) dest);
    /* number of register trains */
    i = char_count >> (BYTESHIFT + REGS_SHIFT);

  loop1: REPEAT32 ( if (i==0) goto exit1; \
q[0]=a0;q[1]=a1;q[2]=a2;q[3]=a3;     q[4]=a0;q[5]=a1;q[6]=a2;q[7]=a3; \
if (i==1) goto exit1; \
q[8]=a4;q[9]=a5;q[10]=a6;q[11]=a7;   q[12]=a4;q[13]=a5;q[14]=a6;q[15]=a7; \
if (i==2) goto exit1; \
q[16]=a0;q[17]=a1;q[18]=a2;q[19]=a3; q[20]=a0;q[21]=a1;q[22]=a2;q[23]=a3; \
if (i==3) goto exit1; \
q[24]=a4;q[25]=a5;q[26]=a6;q[27]=a7; q[28]=a4;q[29]=a5;q[30]=a6;q[31]=a7; \
i-=4; q+=32; );
    goto loop1;
  exit1:

    /* what's left are still int8, but cannot fill a train */
    i = (char_count>>BYTESHIFT)&REGS_MASK;
    if (i==0) goto exit2;
    /* clean up the int8 tail: REGS-1 druple */
    q = ((REGTYPE *)dest) + ((char_count>>BYTESHIFT)&REGS_NASK);
    q[0] = a0;
    if (i==1) goto exit2; q[1] = a1;
    if (i==2) goto exit2; q[2] = a2;
    if (i==3) goto exit2; q[3] = a3;
    if (i==4) goto exit2; q[4] = a4;
    if (i==5) goto exit2; q[5] = a5;
    if (i==REGS-1) q[REGS-2] = a6;
  exit2:

    /* what's left are only chars */
    i = char_count & BYTEMASK;
    /* clean up the char tail: BYTES-1 druple */
    if (i==0) goto exit3; dest += (char_count & BYTENASK);
    dest[0] = c;
    if (i==1) goto exit3; dest[1] = c;
    if (i==2) goto exit3; dest[2] = c;
    if (i==3) goto exit3; dest[3] = c;
    if (i==4) goto exit3; dest[4] = c;
    if (i==5) goto exit3; dest[5] = c;
    if (i==BYTES-1) dest[BYTES-2] = c;
  exit3:
    return;
} /* end char_int8_aset() */


/* Copy "char_count" x (char)"c" to "dest" using eight int8 registers */
void char_int8_set (char *dest, char c, size_t char_count)
{
    register unsigned long i;
    /* how many bytes are in front of QWORD alignment */
    i = ((unsigned long)dest) & BYTEMASK;
    if (i==0)
    {
        char_int8_aset (dest, c, char_count);
        return;
    }
    i = BYTES-i; if (i > char_count) i = char_count;
    if (i==0) return;    dest[0] = c;
    if (i==1) goto exit; dest[1] = c;
    if (i==2) goto exit; dest[2] = c;
    if (i==3) goto exit; dest[3] = c;
    if (i==4) goto exit; dest[4] = c;
    if (i==5) goto exit; dest[5] = c;
    if (i==BYTES-1) dest[BYTES-2] = c;
  exit:
    if (i == char_count) return;
    char_int8_aset (dest+i, c, char_count-i);
    return;
} /* end char_int8_set() */


#ifdef _char_int8_TEST
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
        char_int8_set (a, c, N);
        charPfool(a);
    }
    stop_chronometer();
    printf ("char_int8_set: %g s -> %.3f MB/s\n", stopped_usertime(),
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
#endif /* _char_int8_TEST */


/* Copy "uchar_count" x (unsigned char)"c" to */
/* "dest". "dest" must be aligned to QWORD.   */
void uchar_int8_aset
(unsigned char *dest, unsigned char c, size_t uchar_count)
{
    register REGTYPE a0,a1,a2,a3,a4,a5,a6,a7;  /* REGS */
    register REGTYPE *q;
    register size_t i;

    a0  = c;  a0 <<= 8 * sizeof(unsigned char);
    a0 |= c;  a0 <<= 8 * sizeof(unsigned char);
    a0 |= c;  a0 <<= 8 * sizeof(unsigned char);
    a0 |= c;  a0 <<= 8 * sizeof(unsigned char);
    a0 |= c;  a0 <<= 8 * sizeof(unsigned char);
    a0 |= c;  a0 <<= 8 * sizeof(unsigned char);
    a0 |= c;  a0 <<= 8 * sizeof(unsigned char);
    a0 |= c;
    a1 = a2 = a3 = a4 = a5 = a6 = a7 = a0;

    /* there SHOULD be no misalignment */
    q = ((REGTYPE *) dest);
    /* number of register trains */
    i = uchar_count >> (BYTESHIFT + REGS_SHIFT);

  loop1: REPEAT32 ( if (i==0) goto exit1; \
q[0]=a0;q[1]=a1;q[2]=a2;q[3]=a3;     q[4]=a0;q[5]=a1;q[6]=a2;q[7]=a3; \
if (i==1) goto exit1; \
q[8]=a4;q[9]=a5;q[10]=a6;q[11]=a7;   q[12]=a4;q[13]=a5;q[14]=a6;q[15]=a7; \
if (i==2) goto exit1; \
q[16]=a0;q[17]=a1;q[18]=a2;q[19]=a3; q[20]=a0;q[21]=a1;q[22]=a2;q[23]=a3; \
if (i==3) goto exit1; \
q[24]=a4;q[25]=a5;q[26]=a6;q[27]=a7; q[28]=a4;q[29]=a5;q[30]=a6;q[31]=a7; \
i-=4; q+=32; );
    goto loop1;
  exit1:

    /* what's left are still int8, but cannot fill a train */
    i = (uchar_count>>BYTESHIFT)&REGS_MASK;
    if (i==0) goto exit2;
    /* clean up the int8 tail: REGS-1 druple */
    q = ((REGTYPE *)dest) + ((uchar_count>>BYTESHIFT)&REGS_NASK);
    q[0] = a0;
    if (i==1) goto exit2; q[1] = a1;
    if (i==2) goto exit2; q[2] = a2;
    if (i==3) goto exit2; q[3] = a3;
    if (i==4) goto exit2; q[4] = a4;
    if (i==5) goto exit2; q[5] = a5;
    if (i==REGS-1) q[REGS-2] = a6;
  exit2:

    /* what's left are only uchars */
    i = uchar_count & BYTEMASK;
    /* clean up the uchar tail: BYTES-1 druple */
    if (i==0) goto exit3; dest += (uchar_count & BYTENASK);
    dest[0] = c;
    if (i==1) goto exit3; dest[1] = c;
    if (i==2) goto exit3; dest[2] = c;
    if (i==3) goto exit3; dest[3] = c;
    if (i==4) goto exit3; dest[4] = c;
    if (i==5) goto exit3; dest[5] = c;
    if (i==BYTES-1) dest[BYTES-2] = c;
  exit3:
    return;
} /* end uchar_int8_aset() */


/* Copy "uchar_count" x (unsigned char)"c" */
/* to "dest" using eight int8 registers.   */
void uchar_int8_set (unsigned char *dest, unsigned char c, size_t uchar_count)
{
    register unsigned long i;
    /* how many bytes are in front of QWORD alignment */
    i = ((unsigned long)dest) & BYTEMASK;
    if (i==0)
    {
        uchar_int8_aset (dest, c, uchar_count);
        return;
    }
    i = BYTES-i; if (i > uchar_count) i = uchar_count;
    if (i==0) return;    dest[0] = c;
    if (i==1) goto exit; dest[1] = c;
    if (i==2) goto exit; dest[2] = c;
    if (i==3) goto exit; dest[3] = c;
    if (i==4) goto exit; dest[4] = c;
    if (i==5) goto exit; dest[5] = c;
    if (i==BYTES-1) dest[BYTES-2] = c;
  exit:
    if (i == uchar_count) return;
    uchar_int8_aset (dest+i, c, uchar_count-i);
    return;
} /* end uchar_int8_set() */


#ifdef _uchar_int8_TEST
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
        uchar_int8_set (a, c, N);
        ucharPfool(a);
    }
    stop_chronometer();
    printf ("uchar_int8_set: %g s -> %.3f MB/s\n", stopped_usertime(),
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
#endif /* _uchar_int8_TEST */


/* Copy "int2_count" x (int2)"c" to "dest". "dest" must be */
/* aligned to QWORD, which is generally a good idea anyway */
void int2_int8_aset (int2 *dest, int2 c, size_t int2_count)
{
    register REGTYPE a0,a1,a2,a3,a4,a5,a6,a7;  /* REGS */
    register REGTYPE *q;
    register size_t i;

    a0  = c;  a0 <<= 8 * sizeof(int2);
    a0 |= c;  a0 <<= 8 * sizeof(int2);
    a0 |= c;  a0 <<= 8 * sizeof(int2);
    a0 |= c;
    a1 = a2 = a3 = a4 = a5 = a6 = a7 = a0;

    /* there SHOULD be no misalignment */
    q = ((REGTYPE *) dest);
    /* number of register trains */
    i = int2_count >> (INT2SHIFT + REGS_SHIFT);

  loop1: REPEAT32 ( if (i==0) goto exit1; \
q[0]=a0;q[1]=a1;q[2]=a2;q[3]=a3;     q[4]=a0;q[5]=a1;q[6]=a2;q[7]=a3; \
if (i==1) goto exit1; \
q[8]=a4;q[9]=a5;q[10]=a6;q[11]=a7;   q[12]=a4;q[13]=a5;q[14]=a6;q[15]=a7; \
if (i==2) goto exit1; \
q[16]=a0;q[17]=a1;q[18]=a2;q[19]=a3; q[20]=a0;q[21]=a1;q[22]=a2;q[23]=a3; \
if (i==3) goto exit1; \
q[24]=a4;q[25]=a5;q[26]=a6;q[27]=a7; q[28]=a4;q[29]=a5;q[30]=a6;q[31]=a7; \
i-=4; q+=32; );
    goto loop1;
  exit1:

    /* what's left are still int8, but cannot fill a train */
    i = (int2_count>>INT2SHIFT)&REGS_MASK;
    if (i==0) goto exit2;
    /* clean up the int8 tail: REGS-1 druple */
    q = ((REGTYPE *)dest) + ((int2_count>>INT2SHIFT)&REGS_NASK);
    q[0] = a0;
    if (i==1) goto exit2; q[1] = a1;
    if (i==2) goto exit2; q[2] = a2;
    if (i==3) goto exit2; q[3] = a3;
    if (i==4) goto exit2; q[4] = a4;
    if (i==5) goto exit2; q[5] = a5;
    if (i==REGS-1) q[REGS-2] = a6;
  exit2:

    /* what's left are only int2s */
    i = int2_count & INT2MASK;
    /* clean up the int2 tail: INT2S-1 druple */
    if (i==0) goto exit3; dest += (int2_count & INT2NASK);
    dest[0] = c;
    if (i==1) goto exit3; dest[1] = c;
    if (i==INT2S-1) dest[INT2S-2] = c;
  exit3:
    return;
} /* end int2_int8_aset() */


/* Copy "int2_count" of (int2)"c" to "dest" using eight int8  */
/* registers. "dest" must minimally be aligned to int2 (WORD) */
void int2_int8_set (int2 *dest, int2 c, size_t int2_count)
{
    register unsigned long i;
    /* how many int2 are in front of QWORD alignment */
    i = (((unsigned long)dest)>>1) & INT2MASK;
    if (i==0)
    {
        int2_int8_aset (dest, c, int2_count);
        return;
    }
    i = INT2S-i; if (i > int2_count) i = int2_count;
    if (i==0) return;    dest[0] = c;
    if (i==1) goto exit; dest[1] = c;
    if (i==INT2S-1) dest[INT2S-2] = c;
  exit:
    if (i == int2_count) return;
    int2_int8_aset (dest+i, c, int2_count-i);
    return;
} /* end int2_int8_set() */


#ifdef _int2_int8_TEST
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
        int2_int8_set (a, c, N);
        int2Pfool(a);
    }
    stop_chronometer();
    printf ("int2_int8_set: %g s -> %.3f MB/s\n", stopped_usertime(),
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
#endif /* _int2_int8_TEST */


/* Copy "int4_count" x (int4)"c" to "dest". "dest" must be */
/* aligned to QWORD, which is generally a good idea anyway */
void int4_int8_aset (int4 *dest, int4 c, size_t int4_count)
{
    register REGTYPE a0,a1,a2,a3,a4,a5,a6,a7;  /* REGS */
    register REGTYPE *q;
    register size_t i;

    a0  = c;  a0 <<= 8 * sizeof(int4);
    a0 |= c;
    a1 = a2 = a3 = a4 = a5 = a6 = a7 = a0;

    /* there SHOULD be no misalignment */
    q = ((REGTYPE *) dest);
    /* number of register trains */
    i = int4_count >> (INT4SHIFT + REGS_SHIFT);

  loop1: REPEAT32 ( if (i==0) goto exit1; \
q[0]=a0;q[1]=a1;q[2]=a2;q[3]=a3;     q[4]=a0;q[5]=a1;q[6]=a2;q[7]=a3; \
if (i==1) goto exit1; \
q[8]=a4;q[9]=a5;q[10]=a6;q[11]=a7;   q[12]=a4;q[13]=a5;q[14]=a6;q[15]=a7; \
if (i==2) goto exit1; \
q[16]=a0;q[17]=a1;q[18]=a2;q[19]=a3; q[20]=a0;q[21]=a1;q[22]=a2;q[23]=a3; \
if (i==3) goto exit1; \
q[24]=a4;q[25]=a5;q[26]=a6;q[27]=a7; q[28]=a4;q[29]=a5;q[30]=a6;q[31]=a7; \
i-=4; q+=32; );
    goto loop1;
  exit1:

    /* what's left are still int8, but cannot fill a train */
    i = (int4_count>>INT4SHIFT)&REGS_MASK;
    if (i==0) goto exit2;
    /* clean up the int8 tail: REGS-1 druple */
    q = ((REGTYPE *)dest) + ((int4_count>>INT4SHIFT)&REGS_NASK);
    q[0] = a0;
    if (i==1) goto exit2; q[1] = a1;
    if (i==2) goto exit2; q[2] = a2;
    if (i==3) goto exit2; q[3] = a3;
    if (i==4) goto exit2; q[4] = a4;
    if (i==5) goto exit2; q[5] = a5;
    if (i==REGS-1)   q[REGS-2] = a6;
  exit2:

    /* what's left are only int4s */
    if (int4_count & INT4MASK) dest[int4_count & INT4NASK] = c;
    return;
} /* end int4_int8_aset() */


/* Copy "int4_count" of (int4)"c" to "dest" using eight int8   */
/* registers. "dest" must minimally be aligned to int4 (DWORD) */
void int4_int8_set (int4 *dest, int4 c, size_t int4_count)
{
    register unsigned long i;
    /* how many int4 are in front of QWORD alignment */
    i = (((unsigned long)dest)>>2) & INT4MASK;
    if (i==0)
    {
        int4_int8_aset (dest, c, int4_count);
        return;
    }
    if (int4_count > 0)
    {
        dest[0] = c;
        int4_int8_aset (dest+1, c, int4_count-1);
    }
    return;
} /* end int4_int8_set() */


#ifdef _int4_int8_TEST
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
        int4_int8_set (a, c, N);
        int4Pfool(a);
    }
    stop_chronometer();
    printf ("int4_int8_set: %g s -> %.3f MB/s\n", stopped_usertime(),
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
#endif /* _int4_int8_TEST */


/* Copy "int8_count" x (int8)"c" to "dest" using eight int8    */
/* registers. "dest" must minimally be aligned to int8 (QWORD) */
void int8_int8_set (int8 *dest, int8 c, size_t int8_count)
{
    register REGTYPE a0,a1,a2,a3,a4,a5,a6,a7;  /* REGS */
    register REGTYPE *q;
    register size_t i;

    a1 = a2 = a3 = a4 = a5 = a6 = a7 = a0 = c;
    /* there SHOULD be no misalignment */
    q = ((REGTYPE *) dest);
    /* number of register trains */
    i = int8_count >> REGS_SHIFT;

  loop1: REPEAT32 ( if (i==0) goto exit1; \
q[0]=a0;q[1]=a1;q[2]=a2;q[3]=a3;     q[4]=a4;q[5]=a5;q[6]=a6;q[7]=a7; \
if (i==1) goto exit1; \
q[8]=a0;q[9]=a1;q[10]=a2;q[11]=a3;   q[12]=a4;q[13]=a5;q[14]=a6;q[15]=a7; \
if (i==2) goto exit1; \
q[16]=a0;q[17]=a1;q[18]=a2;q[19]=a3; q[20]=a4;q[21]=a5;q[22]=a6;q[23]=a7; \
if (i==3) goto exit1; \
q[24]=a0;q[25]=a1;q[26]=a2;q[27]=a3; q[28]=a4;q[29]=a5;q[30]=a6;q[31]=a7; \
i-=4; q+=32; );
    goto loop1;
  exit1:

    /* what's left are still int8, but cannot fill a train */
    i = int8_count & REGS_MASK;
    if (i==0) goto exit2;
    /* clean up the int8 tail: REGS-1 druple */
    dest[(int8_count&REGS_NASK)] = a0;
    if (i==1) goto exit2; dest[(int8_count&REGS_NASK)+1] = a1;
    if (i==2) goto exit2; dest[(int8_count&REGS_NASK)+2] = a2;
    if (i==3) goto exit2; dest[(int8_count&REGS_NASK)+3] = a3;
    if (i==4) goto exit2; dest[(int8_count&REGS_NASK)+4] = a4;
    if (i==5) goto exit2; dest[(int8_count&REGS_NASK)+5] = a5;
    if (i==REGS-1)   dest[(int8_count&REGS_NASK)+REGS-2] = a6;
  exit2:
    return;
} /* end int8_int8_set() */


#ifdef _int8_int8_TEST
#define OP   (TESTBYTES/TASKBYTES)
#define N    (TASKBYTES/sizeof(int8))
#define c    255
int main (int argc, char *argv[])
{
    static int8 BUFFER[N+2]={0};
    int8 *a = BUFFER+1;
    register int i, j;
    
    memset (a, c^1, sizeof(int8)*N);  /* load it into L2 cache */
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
        int8_int8_set (a, c, N);
        int8Pfool(a);
    }
    stop_chronometer();
    printf ("int8_int8_set: %g s -> %.3f MB/s\n", stopped_usertime(),
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
#endif /* _int8_int8_TEST */



/** REGS isomorphic booty: **/

/* Copy "real8_count" x (real8)"c" to "dest" using eight real8  */
/* registers. "dest" must minimally be aligned to real8 (QWORD) */
void real8_real8_set (real8 *dest, real8 c, size_t real8_count)
{
    register real8 a0,a1,a2,a3,a4,a5,a6,a7;  /* REGS */
    register real8 *q;
    register size_t i;

    a1 = a2 = a3 = a4 = a5 = a6 = a7 = a0 = c;
    /* there SHOULD be no misalignment */
    q = ((real8 *) dest);
    /* number of register trains */
    i = real8_count >> REGS_SHIFT;

  loop1: REPEAT32 ( if (i==0) goto exit1; \
q[0]=a0;q[1]=a1;q[2]=a2;q[3]=a3;     q[4]=a4;q[5]=a5;q[6]=a6;q[7]=a7; \
if (i==1) goto exit1; \
q[8]=a0;q[9]=a1;q[10]=a2;q[11]=a3;   q[12]=a4;q[13]=a5;q[14]=a6;q[15]=a7; \
if (i==2) goto exit1; \
q[16]=a0;q[17]=a1;q[18]=a2;q[19]=a3; q[20]=a4;q[21]=a5;q[22]=a6;q[23]=a7; \
if (i==3) goto exit1; \
q[24]=a0;q[25]=a1;q[26]=a2;q[27]=a3; q[28]=a4;q[29]=a5;q[30]=a6;q[31]=a7; \
i-=4; q+=32; );
    goto loop1;
  exit1:

    /* what's left are still real8, but cannot fill a train */
    i = real8_count & REGS_MASK;
    if (i==0) goto exit2;
    /* clean up the real8 tail: REGS-1 druple */
    dest[(real8_count&REGS_NASK)] = a0;
    if (i==1) goto exit2; dest[(real8_count&REGS_NASK)+1] = a1;
    if (i==2) goto exit2; dest[(real8_count&REGS_NASK)+2] = a2;
    if (i==3) goto exit2; dest[(real8_count&REGS_NASK)+3] = a3;
    if (i==4) goto exit2; dest[(real8_count&REGS_NASK)+4] = a4;
    if (i==5) goto exit2; dest[(real8_count&REGS_NASK)+5] = a5;
    if (i==REGS-1)   dest[(real8_count&REGS_NASK)+REGS-2] = a6;
  exit2:
    return;
} /* end real8_real8_set() */

/**************************************************************************/
/* Could not do better because if we use real8, we have to union it with  */
/* two real4's first as in int8_int4_set(). But it may lead to floating   */
/* point exception as we try to store it to a register, since the joint   */
/* or decomposition of an integer is always integer but could be NaN for  */
/* floating points. This hardware exception occurs on Pentium even when   */
/* no arithmetic ops are carried out and we are just loading the meory to */
/* register. On Origin 2000 this does not occur, still it's bad practice. */
/**************************************************************************/

#ifdef _real8_real8_TEST
#define OP   (TESTBYTES/TASKBYTES)
#define N    (TASKBYTES/sizeof(real8))
#define c    ((real8)0.2)
int main (int argc, char *argv[])
{
    static real8 BUFFER[N+2]={0};
    real8 *a = BUFFER+1;
    register int i, j;
    
    memset (a, 0, sizeof(real8)*N); /* load it into L2 cache */

    start_chronometer();
    for (i=0; i<OP; i++)
    {
        for (j=0; j<N; j++) a[j] = c;
        real8Pfool(a);
    }
    stop_chronometer();
    printf ("simple approach: %g s -> %.3f MB/s\n", stopped_usertime(),
            1e-6*N*OP*sizeof(real8)/stopped_usertime() );

    memset (a, 0, sizeof(real8)*N);
    start_chronometer();
    for (i=0; i<OP; i++)
    {
        real8_real8_set (a, c, N);
        real8Pfool(a);
    }
    stop_chronometer();
    printf ("real8_real8_set: %g s -> %.3f MB/s\n", stopped_usertime(),
            1e-6*N*OP*sizeof(real8)/stopped_usertime() );

    start_chronometer();
    if (a[-1] != 0) pe ("a[%d] = %f wrong\n", -1, a[-1]);
    for (i=0; i<N; i++)
        if (a[i] != c) printf ("a[%d] = %f wrong\n", i, a[i]);
    if (a[i] != 0) pe ("a[%d] = %f wrong\n", i, a[i]);
    stop_chronometer();
    printf ("successful... checking took %g s -> %.3f MB/s\n",
            stopped_usertime(),
            1e-6*N*sizeof(real8)/stopped_usertime());
    
    start_chronometer();
    for (i=0; i<OP; i++)
    {
        bzero (a, N*sizeof(real8));
        real8Pfool(a);
    }
    stop_chronometer();
    printf ("bzero: %g s -> %.3f MB/s\n", stopped_usertime(),
            1e-6*N*OP*sizeof(real8)/stopped_usertime() );
    return (1);
}
#undef c
#undef N
#undef OP
#endif /* _real8_real8_TEST */


/* Copy "real4_count" x (real4)"c" to "dest" using eight real4  */
/* registers. "dest" must minimally be aligned to real4 (DWORD) */
void real4_real4_set (real4 *dest, real4 c, size_t real4_count)
{
    register real4 a0,a1,a2,a3,a4,a5,a6,a7;  /* REGS */
    register real4 *q;
    register size_t i;

    a1 = a2 = a3 = a4 = a5 = a6 = a7 = a0 = c;
    /* there SHOULD be no misalignment */
    q = dest;
    /* number of register trains */
    i = real4_count >> REGS_SHIFT;

  loop1: REPEAT32 ( if (i==0) goto exit1; \
q[0]=a0;q[1]=a1;q[2]=a2;q[3]=a3;     q[4]=a4;q[5]=a5;q[6]=a6;q[7]=a7; \
if (i==1) goto exit1; \
q[8]=a0;q[9]=a1;q[10]=a2;q[11]=a3;   q[12]=a4;q[13]=a5;q[14]=a6;q[15]=a7; \
if (i==2) goto exit1; \
q[16]=a0;q[17]=a1;q[18]=a2;q[19]=a3; q[20]=a4;q[21]=a5;q[22]=a6;q[23]=a7; \
if (i==3) goto exit1; \
q[24]=a0;q[25]=a1;q[26]=a2;q[27]=a3; q[28]=a4;q[29]=a5;q[30]=a6;q[31]=a7; \
i-=4; q+=32; );
    goto loop1;
  exit1:

    /* what's left are still real4, but cannot fill a train */
    i = real4_count & REGS_MASK;
    if (i==0) goto exit2;
    /* clean up the real4 tail: REGS-1 druple */
    dest[(real4_count&REGS_NASK)] = a0;
    if (i==1) goto exit2; dest[(real4_count&REGS_NASK)+1] = a1;
    if (i==2) goto exit2; dest[(real4_count&REGS_NASK)+2] = a2;
    if (i==3) goto exit2; dest[(real4_count&REGS_NASK)+3] = a3;
    if (i==4) goto exit2; dest[(real4_count&REGS_NASK)+4] = a4;
    if (i==5) goto exit2; dest[(real4_count&REGS_NASK)+5] = a5;
    if (i==REGS-1)   dest[(real4_count&REGS_NASK)+REGS-2] = a6;
  exit2:
    return;
} /* end real4_real4_set() */


#ifdef _real4_real4_TEST
#define OP   (TESTBYTES/TASKBYTES)
#define N    (TASKBYTES/sizeof(real4))
#define c    ((real4)0.2)
int main (int argc, char *argv[])
{
    static real4 BUFFER[N+2]={0};
    real4 *a = BUFFER+1;
    register int i, j;

    memset (a, 0, sizeof(real4)*N);  /* load it into L2 cache */
    start_chronometer();
    for (i=0; i<OP; i++)
    {
        for (j=0; j<N; j++) a[j] = c;
        real4Pfool(a);
    }
    stop_chronometer();
    printf ("simple approach: %g s -> %.3f MB/s\n", stopped_usertime(),
            1e-6*N*OP*sizeof(real4)/stopped_usertime() );

    memset (a, 0, sizeof(real4)*N);
    start_chronometer();
    for (i=0; i<OP; i++)
    {
        real4_real4_set (a, c, N);
        real4Pfool(a);
    }
    stop_chronometer();
    printf ("real4_real4_set: %g s -> %.3f MB/s\n", stopped_usertime(),
            1e-6*N*OP*sizeof(real4)/stopped_usertime() );

    start_chronometer();
    if (a[-1] != 0) pe ("a[%d] = %f wrong\n", -1, a[-1]);
    for (i=0; i<N; i++)
        if (a[i] != c) printf ("a[%d] = %f wrong\n", i, a[i]);
    if (a[i] != 0) pe ("a[%d] = %f wrong\n", i, a[i]);
    stop_chronometer();
    printf ("successful... checking took %g s -> %.3f MB/s\n",
            stopped_usertime(),
            1e-6*N*sizeof(real4)/stopped_usertime());
    
    start_chronometer();
    for (i=0; i<OP; i++)
    {
        bzero (a, N*sizeof(real4));
        real4Pfool(a);
    }
    stop_chronometer();
    printf ("bzero: %g s -> %.3f MB/s\n", stopped_usertime(),
            1e-6*N*OP*sizeof(real4)/stopped_usertime() );
    return (1);
}
#undef c
#undef N
#undef OP
#endif /* _real4_real4_TEST */


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
#undef INT4S   
#undef INT4SHIFT
#undef INT4MASK 
#undef INT4NASK 
#undef TESTBYTES
#undef TASKBYTES
