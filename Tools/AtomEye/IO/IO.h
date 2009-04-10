/***************************************/
/* libIO:                              */
/*                                     */
/* Strings, parsers & files            */
/*                                     */
/* Dec.12, 1999 Ju Li <liju99@mit.edu> */
/***************************************/

#ifndef _IO_h
#define _IO_h

#include <math.h>
#include <fcntl.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <termios.h>
#include <unistd.h>
#include <signal.h>
#include <limits.h>
#include <ctype.h>
#include <sys/wait.h>
#include <sys/types.h>
#include <sys/file.h>
#include <sys/stat.h>
#include <sys/errno.h>
#include <sys/param.h>
#include <sys/resource.h>
#include <readline/readline.h>
#include <readline/history.h>
#include <glob.h>
extern char *xmalloc ();

/* IO is the base of all Ju Li's C performance libraries */
#define JLCPL_TITLE "C performance libraries"
#define JLCPL_VERSION "1.0"
#define JLCPL_YR  "2004"
#define JLCPL_YRS "1999-2004"
#define JLCPL_DATE "Mar. 2, 2004"
#define JLCPL_AUTHOR_FN "Ju"
#define JLCPL_AUTHOR_LN "Li"
#define JLCPL_AUTHOR JLCPL_AUTHOR_FN " " JLCPL_AUTHOR_LN
#define JLCPL_EMAIL "liju99@alum.mit.edu"
#define JLCPL_URL "http://alum.mit.edu/www/liju99/"

#define WORD_BYTES       2U
#define WORD_BYTESHIFT   1
#define WORD_BYTEMASK   (WORD_BYTES-1)
#define WORD_BYTENASK   (~WORD_BYTEMASK)
#define DWORD_BYTES      4U
#define DWORD_BYTESHIFT  2
#define DWORD_BYTEMASK  (DWORD_BYTES-1)
#define DWORD_BYTENASK  (~DWORD_BYTEMASK)
#define QWORD_BYTES      8U
#define QWORD_BYTESHIFT  3
#define QWORD_BYTEMASK  (QWORD_BYTES-1)
#define QWORD_BYTENASK  (~QWORD_BYTEMASK)

/************************/
/* Machine Architecture */
/************************/

#if defined(_IRIX64) || defined(_IRIX)

#define INT_REGISTER_BYTES      (8 * 4)
typedef unsigned short int int2;
typedef unsigned int int4;
typedef unsigned long int8;
#define REAL_REGISTER_BYTES     (8 * 8)
typedef float real4;
typedef double real8;
/* special feature, non-transportable */
typedef long double real16;
/* L1 data cache line size */
#define CACHELINE_BYTES  32
/* instruction cache size: 32 Kbytes, data cache size: 32 Kbytes */
#define L1_CODE_BYTES    (32*1024)
#define L1_DATA_BYTES    (32*1024)
/* page size: 4 kB */
#define PAGE_BYTES       (4*1024)
/* secondary unified instruction/data cache size: 4 Mbytes */
#define L2_BYTES         (4*1024*1024)
/* main memory: 2240 Mbytes */
#define MEMORY_BYTES     (2048L*1024*1024)

/* see Timer.h */
#define int2fool(a) ushortfool(a)
#define int4fool(a) uintfool(a)
#define int8fool(a) ulongfool(a)
#define real4fool(a) floatfool(a)
#define real8fool(a) doublefool(a)
#define int2Pfool(a) ushortPfool(a)
#define int4Pfool(a) uintPfool(a)
#define int8Pfool(a) ulongPfool(a)
#define real4Pfool(a) floatPfool(a)
#define real8Pfool(a) doublePfool(a)
/* at O2000 int8 movers are much faster than the simple "for-loop" write */
#define charset(dest, c, char_count) \
  memset((void *)(dest), (int)(c), (size_t)(char_count))
/* char_int8_set((char *)(dest), (char)(c), (size_t)(char_count)) */
#define ucharset(dest, c, uchar_count) \
  memset((void *)(dest), (int)(c), (size_t)(uchar_count))
/* uchar_int8_set((unsigned char *)(dest), (unsigned char)(c), \ */
/* (size_t)(uchar_count)) */
#define int2set(dest, i2, int2_count) \
  int2_int8_set((int2 *)(dest), (int2)(i2), (size_t)(int2_count))
#define int4set(dest, i4, int4_count) \
  int4_int8_set((int4 *)(dest), (int4)(i4), (size_t)(int4_count))
#define int8set(dest, i8, int8_count) \
  int8_int8_set((int8 *)(dest), (int8)(i8), (size_t)(int8_count))
/* almost 3 times the speed when both compiled under the best flags */
/* the above all match the throughput of system memset/bzero */

#else                                      /** Pentium II **/

/* EAX EBX ECX EDX EDI ESI EBP(frame pointer) ESP(stack pointer) */
#define INT_REGISTER_BYTES       (8 * 4)
typedef unsigned short int int2;
typedef unsigned int int4;
typedef unsigned long long int8;
/* 16 DWORD float registers = 8 QWORD registers */
#define REAL_REGISTER_BYTES      (8 * 8)
typedef float real4;
typedef double real8;
#define CACHELINE_BYTES  32
/* Level 1 cache: data 16 kB, code 16 kB */
#define L1_DATA_BYTES    (16*1024)
#define L1_CODE_BYTES    (16*1024)
/* page size: 4 kB */
#if defined(_CYGWIN)
#define PAGE_SHIFT        12
#define PAGE_SIZE        (1UL << PAGE_SHIFT)
#else
#include <sys/user.h>
#endif
#define PAGE_BYTES       (PAGE_SIZE)
/* Level 2 cache: 512 kB */
#define L2_BYTES         (512*1024)
/* main memory: 256 MB */
#define MEMORY_BYTES     (256*1024*1024)

/* see Timer.h */
#define int2fool(a) ushortfool(a)
#define int4fool(a) uintfool(a)
#define int8fool(a) ulonglongfool(a)
#define real4fool(a) floatfool(a)
#define real8fool(a) doublefool(a)
#define int2Pfool(a) ushortPfool(a)
#define int4Pfool(a) uintPfool(a)
#define int8Pfool(a) ulonglongPfool(a)
#define real4Pfool(a) floatPfool(a)
#define real8Pfool(a) doublePfool(a)

/* on PII int4 movers are fast and stable */
#define charset(dest, c, char_count) \
  memset((void *)(dest), (int)(c), (size_t)(char_count))
/* char_int4_set((char *)(dest), (char)(c), (size_t)(char_count)) */
#define ucharset(dest, c, uchar_count) \
  memset((void *)(dest), (int)(c), (size_t)(uchar_count))
/* uchar_int4_set((unsigned char *)(dest), (unsigned char)(c), \ */
/* (size_t)(uchar_count)) */
#define int2set(dest, i2, int2_count) \
  int2_int4_set((int2 *)(dest), (int2)(i2), (size_t)(int2_count))
#define int4set(dest, i4, int4_count) \
  int4_int4_set((int4 *)(dest), (int4)(i4), (size_t)(int4_count))
#define int8set(dest, i8, int8_count) \
  int8_int8_set((int8 *)(dest), (int8)(i8), (size_t)(int8_count))
/* int8_int4_set((int8 *)(dest), (int8)(i8), (size_t)(int8_count)) */
/* the above all match the throughput of system memset/bzero */

#endif  /* Machine Architecture */

#define ISNULL(p)  (VOIDP(p)==NULL)
#define NOTNULL(p) (VOIDP(p)!=NULL)

/**************/
/* Convention */
/**************/

/* for single use, char really destroys alignment in a stack */
typedef int bool;
/* when used as array, char or Bmap can be considered */
/* typedef char boolean; */
#define TRUE          1
#define FALSE         0
#define ZERO_ONE(a)  ((a)?1:0)
#define bword(a)  ((a)?"yes":"no")
#define Bword(a)  ((a)?"Yes":"No")
#define BWORD(a)  ((a)?"YES":"NO")

/* Flag for "Not Available": */
/* standard flag for nonnegative int and double */
#define NVL       (-999999999)
/* for (-inf,+inf) doubles: because it is so big, risk of data loss by */
/* encoding conflict is minimal if the parameters stored are properly  */
/* non-dimensionalized so they are most likely to be within 1e-3 to    */
/* 1e3, a good idea anyway for storing multi-dimensional floating-pts. */
#define DOUBLENVL (-999999999999999.125)
/* default terminal width */
#define TERMCHAR      79
#define TERMSIZE     (TERMCHAR+1)
/* End Of String */
#define EOS         ((char)0)
/* maximal length of an integer written in decimal */
#define INTEGER_DECIMALCHAR  (31)

typedef char TermString[TERMSIZE];

/* CHAR means effective content, SIZE means memory allocation */

/* forced conversion of all transportable types */
#define INT2(x)      ((int2)(x))
#define INT4(x)      ((int4)(x))
#define INT8(x)      ((int8)(x))
#define CHAR(x)      ((char)(x))
#define UCHAR(x)     ((unsigned char)(x))
#define SHORT(x)     ((short)(x))
#define USHORT(x)    ((unsigned short)(x))
#define INT(x)       ((int)(x))
#define UINT(x)      ((unsigned int)(x))
#define LONG(x)      ((long)(x))
#define ULONG(x)     ((unsigned long)(x))
#define REAL4(x)     ((real4)(x))
#define REAL8(x)     ((real8)(x))
#define FLOAT(x)     ((float)(x))
#define DOUBLE(x)    ((double)(x))
#define VOIDP(x)     ((void *)(x))
#define INT2P(x)     ((int2 *)(x))
#define INT4P(x)     ((int4 *)(x))
#define INT8P(x)     ((int8 *)(x))
#define CHARP(x)     ((char *)(x))
#define UCHARP(x)    ((unsigned char *)(x))
#define SHORTP(x)    ((short *)(x))
#define USHORTP(x)   ((unsigned short *)(x))
#define INTP(x)      ((int *)(x))
#define UINTP(x)     ((unsigned int *)(x))
#define LONGP(x)     ((long *)(x))
#define ULONGP(x)    ((unsigned long *)(x))
#define REAL4P(x)    ((real4 *)(x))
#define REAL8P(x)    ((real8 *)(x))
#define FLOATP(x)    ((float *)(x))
#define DOUBLEP(x)   ((double *)(x))

/* most types bear interpretation as numerical value: */
#define IOVAL_CHAR    0
#define IOVAL_UCHAR   1
#define IOVAL_SHORT   2
#define IOVAL_USHORT  3
#define IOVAL_INT     4
#define IOVAL_UINT    5
#define IOVAL_LONG    6
#define IOVAL_ULONG   7
#define IOVAL_FLOAT   8
#define IOVAL_DOUBLE  9
#define IOVAL(ptr,valtype) ( \
  (valtype==IOVAL_CHAR)   ? (*CHARP(ptr))   : \
  (valtype==IOVAL_UCHAR)  ? (*UCHARP(ptr))  : \
  (valtype==IOVAL_SHORT)  ? (*SHORTP(ptr))  : \
  (valtype==IOVAL_USHORT) ? (*USHORTP(ptr)) : \
  (valtype==IOVAL_INT)    ? (*INTP(ptr))    : \
  (valtype==IOVAL_UINT)   ? (*UINTP(ptr))   : \
  (valtype==IOVAL_LONG)   ? (*LONGP(ptr))   : \
  (valtype==IOVAL_ULONG)  ? (*ULONGP(ptr))  : \
  (valtype==IOVAL_FLOAT)  ? (*FLOATP(ptr))  : \
  (valtype==IOVAL_DOUBLE) ? (*DOUBLEP(ptr)) : 0 )

#define ISINT(i)     ((i)==INT(i))
#define ISDIGIT(c)   ((((char)(c))>='0')&&(((char)(c))<='9'))
#define ISLOWER(c)   ((((char)(c))>='a')&&(((char)(c))<='z'))
#define ISUPPER(c)   ((((char)(c))>='A')&&(((char)(c))<='Z'))
#define ISALPHA(c)   (ISLOWER(c)||ISUPPER(c))
#define ISALNUM(c)   (ISDIGIT(c)||ISALPHA(c))
#define ISBLANK(c)   ((((char)(c))==' ')||(((char)(c))=='\t'))

#define ISNOTINT(i)   ((i)!=INT(i))
#define ISNOTDIGIT(c) ((((char)(c))<'0')||(((char)(c))>'9'))
#define ISNOTLOWER(c) ((((char)(c))<'a')||(((char)(c))>'z'))
#define ISNOTUPPER(c) ((((char)(c))<'A')||(((char)(c))>'Z'))
#define ISNOTALPHA(c) (ISNOTLOWER(c)&&ISNOTUPPER(c))
#define ISNOTALNUM(c) (ISNOTDIGIT(c)&&ISNOTALPHA(c))
#define ISNOTBLANK(c) ((((char)(c))!=' ')&&(((char)(c))!='\t'))

#define SRTRTRT(x)    #x
#define STR(x)        SRTRTRT(x)
#define CAN(x)        do {x;} while(0)

#define ONE_OR_ZERO(x) ((x)&&1)

/* Application Binary Interface with fortran */
#ifndef _HPUX
#define FORTRAN_SYMBOL(a) a##_
#else
#define FORTRAN_SYMBOL(a) a
#endif

/* not using CAN to package because "break" would not work in "x" */
#define REPEAT2(x)    x;x;
#define REPEAT4(x)    x;x;x;x;
#define REPEAT8(x)    x;x;x;x; x;x;x;x;
#define REPEAT16(x)   x;x;x;x; x;x;x;x; x;x;x;x; x;x;x;x;
#define REPEAT32(x)   REPEAT16(x); REPEAT16(x);
#define REPEAT64(x)   REPEAT32(x); REPEAT32(x);
#define REPEAT128(x)  REPEAT64(x); REPEAT64(x);
#define REPEAT256(x)  REPEAT128(x); REPEAT128(x);
/* in "if/else" statements please explicitly use { REPEAT8(x); } */

#define fcr(fp)            fputc('\n', fp)
#define Fcr(fp)            { if (fp) fputc('\n', fp); }
#define cr()               fcr(stdout)
#define fsp(fp)            fputc(' ', fp)
#define Fsp(fp)            { if (fp) fputc(' ', fp); }
#define sp()               fsp(stdout)
#define fqt(fp)            fputc('"', fp)
#define Fqt(fp)            { if (fp) fputc('"', fp); }
#define qt()               fqt(stdout)
#define VALID_FILENO(x)    (((x) >= 0) && ((x) <  NOFILE))
#define INVALID_FILENO(x)  (((x) <  0) || ((x) >= NOFILE))

#define fscanf_skipline(fp)   while (getc(fp)!='\n');
#define scanf_skipline()      fscanf_skipline(stdin)

/** IO.c: miscellaneous I/O subroutines **/

void press_return_to (FILE *in, FILE *out, char *do_what);
#define Press_return() press_return_to(stdin,stdout,"continue")

/* Same as fprintf() except it does not print if "stream" is NULL */
void Fprintf (FILE *stream, char *format, ...);

/* same as fgets() except newline is treated as EOF (not stored) */
/* and what is left in that line is read without storing.        */
char *Fgets (char *s, int size, FILE *stream);

#define FGETS_MAXCHAR 1023
#define FGETS_MAXSIZE 1024
#define FGets(s,stream) Fgets(s,FGETS_MAXSIZE,stream)
#define Gets(s) FGets(s,stdin)

/* clear input stream buffer */
#define clear_input_stream_buffer(stream) tcflush(fileno(stream),TCIFLUSH)
#define clear_stdin_buffer() clear_input_stream_buffer(stdin)

/* Same as fgets() except we replace ^M in the stream by '\n' */
char *FGETS(char *s, int size, FILE *stream);

/********************************************************/
/* Scan a line from "in", and see if "keys" strings are */
/* at the front. If the line is longer than "linesize", */
/* claim that the file might be corrupted and exit. If  */
/* the line is a simple '\n', return -1. If the line is */
/* a simple EOF, return -4. If multiple "keys" match    */
/* the line, use the first match (therefore if you want */
/* to match "book" AND anything else starting with "bo" */
/* but not "book", you put "book" ahead of "bo" in      */
/* "keys"). Ending EOF or '\n' will always be stripped  */
/* from "linebuffer". If there is a match, the "key"    */
/* string will also be stripped. In the event of '\n'   */
/* termination, return keyword index; in the event of   */
/* EOF termination, return -(index+5). If the line has  */
/* content but does not match, return -2 if '\n'-       */
/* terminated, or -3 if EOF-terminated.                 */
/* -------------------- Summary ----------------------- */
/*  >=0: match by '\n'-terminated lines;                */
/*  -1:  empty '\n' line                                */
/*  -2: '\n'-terminated with content but does not match */
/*  -3:  EOF-terminated with content but does not match */
/*  -4:  empty EOF line                                 */
/*  <=-5: match by EOF-terminated lines                 */
/* "linebuffer" will be stripped of EOF or '\n', and    */
/* if there is a match, the matching string itself.     */
/********************************************************/
int freadline_matchheader
(FILE *in, int n_key, char *keys[], int linesize, char *linebuffer);

/* return string "1st" for i=1, "2nd" for i=2, "100th" for i=100, etc. */
char *word_for_order (int i);

#define LETTER_LOWERCASE(c) (ISUPPER(c)?CHAR(c+32):CHAR(c))
#define LETTER_UPPERCASE(c) (ISLOWER(c)?CHAR(c-32):CHAR(c))

/* Convert a string to all uppercase and return back its pointer */
char *convert_to_uppercase (char *s);

/* Convert a string to all lowercase and return back its pointer */
char *convert_to_lowercase (char *s);

/* Similar as UNIX tr: translate characters in s[] from set1 to set2 */
char *convert_tr (char *s, char *set1, char *set2);

/*  "2.13339347864174        9.99200229300000D-002"
 -> "2.13339347864174        9.99200229300000E-002" */
#define STANDARDIZE_NUMERICAL_STRING(s) (convert_tr(s,"Dd","Ee"))

/* check if filename is valid (not NULL or empty string) */
#define Fvalid(filename)  ( ((filename)!=NULL) && ((filename)[0]!=EOS) )
#define Finvalid(filename)  ( !Fvalid(filename) )

#define FSTAT_LSTAT_RUN_ERROR       0
#define FSTAT_IS_REGULAR_FILE       1
#define FSTAT_IS_DIRECTORY          2
#define FSTAT_IS_CHARACTER_SPACIAL  3
#define FSTAT_IS_BLOCK_SPACIAL      4
#define FSTAT_IS_FIFO               5
#define FSTAT_IS_SYMBOLIC_LINK      6
#define FSTAT_IS_SOCKET             7
#define FSTAT_IS_UNKNOWN            8
/* Return nature of the file */
int Fstat (char filename[]);

/* check if file exists */
#define Fexists(filename) (access(filename,F_OK)>=0)
#define Freadable(filename) (access(filename,R_OK)>=0)
#define Fwritable(filename) (access(filename,W_OK)>=0)
#define Fexecutable(filename) (access(filename,X_OK)>=0)

/* Return current umask */
mode_t Umask();

/* For UNIX systems */
#define DIRECTORY_SEPARATOR                 '/'

#define DOPEN_DIR_NAME_INVALID               0
#define DOPEN_DIR_EXISTS                     1
#define DOPEN_NAME_EXISTS_BUT_IS_NOT_A_DIR   2
#define DOPEN_CANNOT_CREATE_DIR              3
/* Create a directory if necessary. Return success or error code */
int dopen (char dirname[]);
/* Create a directory if necessary; print error message upon failure */
void Dopen (char dirname[]);

#define ropen(a)  fopen(CHARP(a), "r")
#define wopen(a)  fopen(CHARP(a), "w")
#define rwopen(a) fopen(CHARP(a), "rw")
#define Wopen(a)  wopen(BLITZKRIEG_FMASK(a))

/* obtain file handle for specified operations, else print error */
FILE *Fopen (char filename[], char operations[]);
#define rOpen(a)  Fopen(CHARP(a), "r")
#define wOpen(a)  Fopen(CHARP(a), "w")
#define rwOpen(a) Fopen(CHARP(a), "rw")
#define WOpen(a)  wOpen(BLITZKRIEG_FMASK(a))
/* safely close a file handle */
#define Fclose(x) CAN( if ((x)!=NULL) fclose(x); )

#define MALLOC(fun,ptr,n,type) { \
  if ((((ptr)=(type *)malloc((n)*sizeof(type)))==NULL)&&((n)>0)) \
  pe (__FILE__ ": line " STR(__LINE__) ":\n" \
  STR(fun) ":MALLOC: error occurred when allocating %d bytes\n" \
  "for pointer \"" STR(ptr) "\" of type \"" STR(type) "\".\n", \
  (n)*sizeof(type)); }
#define CALLOC(fun,ptr,n,type) { \
  if ((((ptr)=(type *)calloc((n)*sizeof(type),1))==NULL)&&((n)>0)) \
  pe (__FILE__ ": line " STR(__LINE__) ":\n" \
  STR(fun) ":CALLOC: error occurred when allocating %d bytes\n" \
  "for pointer \"" STR(ptr) "\" of type \"" STR(type) "\".\n", \
  (n)*sizeof(type)); }
#define REALLOC(fun,ptr,n,type) { \
  if ((((ptr)=(type *)realloc((ptr),(n)*sizeof(type)))==NULL)&&((n)>0)) \
  pe (__FILE__ ": line " STR(__LINE__) ":\n" \
  STR(fun) ":REALLOC: error occurred when allocating %d bytes\n" \
  "for pointer \"" STR(ptr) "\" of type \"" STR(type) "\".\n", \
  (n)*sizeof(type)); }

/* initialize allocated memory to zero at first: */
#define RECALLOC(fun,ptr,n,type) { \
  if ((ptr)==NULL) CALLOC(fun,ptr,n,type) else REALLOC(fun,ptr,n,type); }
/* as accumulator, shrinkable but not expandable */

/* reallocate memory and set all zero: */
#define REALLOC_ZERO(fun,ptr,n,type) { \
  REALLOC(fun,ptr,n,type); bzero(VOIDP(ptr),(n)*sizeof(type)); }
/* as accumulator, clear each time reallocated */

#define CLONE(ptr_original,n,type,ptr_cloned) { \
  MALLOC(CLONE,ptr_cloned,n,type); \
  memcpy((void *)(ptr_cloned), (void *)(ptr_original), (n)*sizeof(type)); }
#define RECLONE(ptr_original,n,type,ptr_cloned) { \
  REALLOC(RECLONE,ptr_cloned,n,type); \
  memcpy((void *)(ptr_cloned), (void *)(ptr_original), (n)*sizeof(type)); }

/* safely free a piece of malloc/calloc/realloc memory */
#define Free(a) { if ((void *)(a)!=NULL) {free((void *)(a)); a=NULL;} }

#define MEMCPY(dest,src,n,type) \
  memcpy( (void *)(dest), (void *)(src), (n)*sizeof(type) )
#define MEMCpy(dest,src,type) \
  memcpy( (void *)(dest), (void *)(src), sizeof(type) )

#define BEQV(n,src,dest) \
  memcpy( (void *)(dest), (void *)(src), (n)*sizeof(char) )

/* Location of EOS in string q */
#define strend(q) ((q)+strlen(q))

/* Same as strstr() except there is no distinction in case */
char *strcasestr (const char *haystack, const char *needle);

/* If "haystack" contains "needle" as its first word, return */
/* pointer at that first occurrence; otherwise return NULL.  */
#define str_begin_with(haystack,needle) \
  ( (strstr(CHARP(haystack),CHARP(needle)) == CHARP(haystack)) ? \
  CHARP(haystack) : NULL )

/* Same as str_begin_with() except there is no distinction in case */
#define str_casebegin_with(haystack,needle) \
  ( (strcasestr(CHARP(haystack),CHARP(needle)) == CHARP(haystack)) ? \
  CHARP(haystack) : NULL )

/* If "haystack" contains "needle" as its last word, return */
/* pointer at that last occurrence; otherwise return NULL.  */
char *str_end_with (char *haystack, char *needle);

/* Same as str_end_with() except there is no distinction in case */
char *str_caseend_with (char *haystack, char *needle);

/* COMPATIBILITY: strsep - extract token from string          */
/* The  strsep()  function  returns  the  next token from the */
/* string stringp which is delimited by delim.  The token  is */
/* terminated with a `\0' character and stringp is updated to */
/* point past the token.   RETURN VALUE:                      */
/* The strsep() function returns a pointer to the  token,  or */
/* NULL if delim is not found in stringp.                     */
char *STRSEP (char **stringp, char *delim);

/* Get the absolute pathname of a shell path */
char *absolute_pathname (char *pathname);

/* Return the basename of a file (without suffix): "/tmp/a.b.c"->"a.b"  */
/* same as echo /tmp/a.b.c | sed -e "s/.*\///g" | sed -e 's/\.[^.]*$//' */
char *file_basename1 (char *pathname);

/* Check whether an executable "command" exists in current PATH */
int command_exists (char *command);
int commands_exists (int n_commands, char *commands[]);

#define NUMBER_POSTSCRIPT_VIEWERS 6
extern char *postscript_viewers[NUMBER_POSTSCRIPT_VIEWERS];
#define NUMBER_RASTER_VIEWERS 4
extern char *raster_viewers[NUMBER_RASTER_VIEWERS];

/* run command in foreground; return -1 if command not found */
int try_to_runfg (int n_commands, char *commands[], char *arg);
int TRY_to_runfg (char *arg, ...);
#define TRY_TO_RUNFG(command,arg) TRY_to_runfg(arg,command,NULL);

/* run command in background; return -1 if command not found */
int try_to_runbg (int n_commands, char *commands[], char *arg);
int TRY_to_runbg (char *arg, ...);
#define TRY_TO_RUNBG(command,arg) try_to_runbg(arg,command,NULL);

#ifndef _HPUX
#define TMPNAM(a) tmpnam(a)
#else
#define TMPNAM(a) P_tmpdir"IO"
#endif

/****************************************************************/
/* Check if the file is in GZIP format: if it is, decompress it */
/* using environmental 'gzip' and return a temporary filename   */
/* that points to the decompressed file. Otherwise return NULL. */
/****************************************************************/
char *RGZIP_fname (char *original_fname);

/*****************************************************************/
/* Check if the file is in BZIP2 format: if it is, decompress it */
/* using environmental 'bzip2' and return a temporary filename   */
/* that points to the decompressed file. Otherwise return NULL.  */
/*****************************************************************/
char *RBZIP2_fname (char *original_fname);

/******************************************************************/
/* Check if the file is in UNIX .z format: if it is, decompress   */
/* using environmental 'compress' and return a temporary filename */
/* that points to the decompressed file. Otherwise return NULL.   */
/******************************************************************/
char *RUNIXZ_fname (char *original_fname);

/***************************************************************/
/* If the file is .bz2, .gz or .z compressed, decompress using */
/* the environment 'bzip2', 'gzip' or 'compress' and return a  */
/* handle to file which will be automatically removed when the */
/* handle is closed. Otherwise return an ordinary read handle. */
/***************************************************************/
FILE *ROpen (char *original_fname);

/***********************************************************************/
/* Any fname that ends with '.blitz' is going to have that taken down  */
/* by 'rename', so the actual fname appears instantaneously and ready. */
/***********************************************************************/
#define BLITZKRIEG_SUFFIX ".blitz"
#define BLITZKRIEG_FMASK(fname) str2(CHARP(fname),BLITZKRIEG_SUFFIX)

/*************************************************************/
/* Realize the compression of a file whose name ends with    */
/* ".bz2" or ".gz". Return after-processed filesize in bytes */
/* Also, any original_fname that ends with ".blitz" is going */
/* to have ".blitz" taken down in a blitzkrieg fashion.      */
/*************************************************************/
long Zrealize (char *original_fname);
#define Zclose(fp,fname) { Fclose(fp); Zrealize(fname); }

/* Test if a filename is writable by opening and closing it */
int tested_to_be_writable (char filename[]);

/* check if a file exists, and if so ask permission to overwrite it */
int Freetowrite (char filename[]);

/* obtain size in bytes of a file with filename "fname" */
long Fsize (char *fname);

/* obtain size in bytes of an opened file with file number "fid" */
long fsize (int fid);
#define fpsize(fp) fsize(fileno(fp))

/* end of string */
char *eos (char *s);

/* copy char to string at cell_numberth static cell and return its pointer */
#define C2S_MAX_CELLS  8
char *c2s (int cell_number, char c);
#define cs0(c)  c2s(0,(char)(c))
#define cs1(c)  c2s(1,(char)(c))

/* forward to char position in "s" that is not a white space, */
/* i.e., form-feed '\f', newline '\n', carriage return '\r',  */
/* horizontal tab '\t' and vertical tab ('\v').               */
char *space_advance(char *s);

/* forward to char position in "s" that is not a blank, i.e., space or tab */
char *blank_advance(char *s);

/* forward to char position in "s" that IS a space or tab, or stop at EOS */
char *nonblank_advance(char *s);

/* Remove trailing blanks: "Ju Li  " -> "Ju Li" */
char *remove_trailing_blanks(char *p);

#define str_append(dest,src,dest_max_size) ( \
  strncpy((dest)+strlen(dest),src,dest_max_size-strlen(dest)), \
  dest[dest_max_size-1]=EOS )
#define STR_append(dest,src) str_append(dest,src,sizeof(dest))


/***************************************************************/
/* Concatenate several strings together: at most STRS_MAX_CHAR */
/* characters in all. strs() is reusable STRS_CARTRIDGE times. */
/***************************************************************/
#define STRS_CARTRIDGE  4
#define STRS_MAX_CHAR   1023
#define STRS_MAX_SIZE   1024
char *strs (int number_of_strings, ...);

#define str2(a,b)  strs(2,CHARP(a),CHARP(b))
#define str3(a,b,c)  strs(3,CHARP(a),CHARP(b),CHARP(c))
#define str4(a,b,c,d)  strs(4,CHARP(a),CHARP(b),CHARP(c),CHARP(d))
#define str5(a,b,c,d,e) \
  strs(5,CHARP(a),CHARP(b),CHARP(c),CHARP(d),CHARP(e))

/* static memory driver of vsnprintf() */
#define STRF_MAX_CHAR  1023
#define STRF_MAX_SIZE  1024
char *vstrf (char *format, va_list ap);
char *strf (char *format, ...);

/* change stderr to a different stream */
void redirect_stderr_to (FILE *fp);

/* perror() driver using strf() */
void pe (char *format, ...);

/* print to stderr stream */
void pr (char *format, ...);

#ifndef _HPUX
/* Find the last occurrence of "needle" in "haystack"     */
/* and return pointer to that location; else return NULL. */
char *strrstr (char *haystack, char *needle);
#endif

/***********************************************************/
/* assume c[0] is the first byte of char rendering device  */
/* crd[rows][cols+1]. We snprintf in fmt to &crd[i][0] the */
/* ith row of A for at most cols characters (not including */
/* EOS), i=0..rows-1, then return the largest width (not   */
/* including EOS). crd[i][?..width-1] is padded with ' '.  */
/***********************************************************/
int Asprintf (char c[], int cols, int rows, char *fmt, void *A);
/* fprintf to fp rows x format from A; return largest row strlen */
int Afprintf (FILE *fp, int cols, int rows, char *fmt, void *A);
#define Aprint(rows,fmt,A) Afprintf(stdout,TERMCHAR,rows,fmt,(void *)(A))

/* fork out a child process running slave_work_till_die() */
void spinoff( void (*slave_work_till_die)() );

/* base64 encoding map */
extern char Base64Alphabet[64];

/* Economically print floating-pt numbers by "0.123" -> ".123" reduction. */
void vzfprintf (FILE *stream, char *format, va_list ap);
void zfprintf (FILE *stream, char *format, ...);
void zprintf (char *format, ...);

typedef struct
{
    int N;
    int idx_min;
    double min;
    int idx_max;
    double max;
    double average;
    double variance;
    double standard_deviation;
} SimpleStatistics;

/*************************************************************/
/* Starting at "series", there are "N" elements of "valtype" */
/* and separated in memory by regular "bytes_separation",    */
/* calculate the statistical properties of common interest.  */
/*************************************************************/
void CalculateSimpleStatistics
(int N, char *series, int bytes_separation, int valtype, SimpleStatistics *s);


/* IOmem.c: */

/* return pointer to an array of n chars */
char *IOalloc (int n);

/* return pointer to an array of n chars which are all cleared to 0 */
char *IOALLOC (int n);

/* clone a string using malloc() memory. Remember to free it with free()! */
char *IOClone (char *str);

/*********************************************************************/
/* fast IO (static.. dynamic) memory manager for string manipulation */
/*********************************************************************/

/* chunk is big enough for most string manipulation tasks */
#define IOMEM_UNIT  (2048)
/* "small" number of chunks for fast search and allocation */
#define IOMEM_SEGS  (8)
/* maximum static and dynamic clients managed by IOmem() */
#define IOMEM_MAX_CLIENTS (8)
/* also "small" for fast search and allocation */

char *IOmem (size_t size_in_chars);
void IOfree (char *ptr);
void IOfreeall ();
/* clone a string using IOmem() memory. Remember to free it with IOfree()! */
char *IOclone (char *str);

/**********************************************************************/
/* Cmem() works in companion with an existent char array: if that     */
/* array is big enough, there is no need for other, possibly dynamic, */
/* allocation; otherwise function backup_alloc() is called to do the  */
/* job. Remember to call FREE(ptr,char_array,backup_free) afterwards. */
/**********************************************************************/
#define Cmem(size_in_chars,char_array,backup_alloc) (((size_in_chars)<= \
sizeof(char_array)/sizeof(char))?char_array : backup_alloc(size_in_chars))


/* scratch.c: */
    
/***********************************************************/
/* indefinite temporary scratch space: provides fp_scratch */
/***********************************************************/

#define SCRATCH_BUFFER_SIZE (16*1024)
extern FILE *fp_scratch;
/* stream fp_scratch is fully buffered: unless fflush()ed or */
extern char scratch_buffer [SCRATCH_BUFFER_SIZE];
/* SCRATCH_BUFFER_SIZE reached, all fprintf()ed stuff is there */

/* At first call, set up fp_scratch which is fully buffered */
/* at scratch_buffer with SCRATCH_BUFFER_SIZE chars. In     */
/* subsequent calls, previous data are discarded and the    */
/* stream is back to SEEK_SET. When kill_scratch() called   */
/* or process terminates, there would be nothing on disk.   */
void reset_scratch();

/* free fp_scratch: as if reset_scratch() was never called */
void kill_scratch();

/* shorthand for fprintf (fp_scratch, ..) */
void scratch (char *format, ...);

/* Dump what has been scratched down to another stream "fp", but adding */
/* string "esc" before each line. esc==NULL would be interpreted as "". */
void dump_scratched_to (FILE *fp, char *esc);


/* dump.c: */

/******************************************************************/
/* Dump indefinite amount of data from a single object to stream, */
/* formatted as linefmt x minline\maxline. One can add labeling   */
/* info at top and bottom if separated by '\n' in linefmt. Use %N */
/* in linefmt to denote line index, %-2N to denote line index-2;  */
/* %I to denote token index, %2I to denote token index +2. Lastly */
/* only minchar\maxchar of the line will be printed to stream.    */
/******************************************************************/

#define DUMP_LINENO  'N'
#define DUMP_TOKENO  'I'

void fDump (FILE *fp, char *linefmt, int minline, int maxline, void *A, 
	    int minchar, int maxchar);
#define Dump(linefmt,minline,maxline,A,minchar,maxchar) \
fDump(stdout,linefmt,minline,maxline,(void *)(A),minchar,maxchar)

#define fdump(fp,linefmt,maxline,A) \
  fDump(fp,linefmt,0,maxline,(void *)(A),0,TERMCHAR)
#define dump(linefmt,maxline,A) fdump(stdout,linefmt,maxline,A)

#define fdumpI2(fp,maxline,A) \
  fdump(fp,strs(7,"\n",STR(A),"[%",cs0(DUMP_TOKENO),\
  "]:  < %d %d >  \t# %1", cs1(DUMP_LINENO)," #\n\n"), maxline, A)
#define dumpI2(maxline,A) fdumpI2(stdout,maxline,A)

#define fdumpI3(fp,maxline,A) \
  fdump(fp,strs(7,"\n",STR(A),"[%",cs0(DUMP_TOKENO),\
  "]:  < %d %d %d >  \t# %1", cs1(DUMP_LINENO)," #\n\n"), maxline, A)
#define dumpI3(maxline,A) fdumpI3(stdout,maxline,A)

#define fdump2(fp,maxline,A) \
  fdump(fp,strs(7,"\n",STR(A),"[%",cs0(DUMP_TOKENO),\
  "]:  (%14.8le %15.8le)\t# %1", cs1(DUMP_LINENO)," #\n\n"), maxline, A)
#define dump2(maxline,A) fdump2(stdout,maxline,A)

#define fdump3(fp,maxline,A) \
  fdump(fp,strs(7,"\n",STR(A),"[%",cs0(DUMP_TOKENO),\
  "]:  (%14.8le %15.8le %15.8le)\t# %1", cs1(DUMP_LINENO)," #\n\n"), \
  maxline, A)
#define dump3(maxline,A) fdump3(stdout,maxline,A)

#define fdumph(fp,maxbytes,A) \
 fdump(fp,strs(8,"\n",STR(A),"[%",cs0(DUMP_TOKENO),"]:   \t", \
"%02x%02x %02x%02x %02x%02x %02x%02x %02x%02x %02x%02x %02x%02x %02x%02x"\
  "  \t# %1",cs1(DUMP_LINENO)," #\n\n"),(maxbytes)/16,CHARP(A))
#define dumph(maxbytes,A) fdumph(stdout,maxbytes,A)

/* Very crude printout of a double precision real matrix */
#define FMUMP_COLUMNS  3
#define FMUMP_FORMAT   "%15.8e"
#define FMUMP_SPACING  "  "
void fMump (FILE *fp, int rows, int cols, double *A);

#define Mump(rows,cols,A) fMump(stdout,rows,cols,A)
#define fmump(fp,rank,A) fMump(fp,rank,rank,A)
#define mump(rank,A) fmump(stdout,rank,A)

/* print.c: */
    
/********************************************************************/
/* Line/matrix multi-object formatter using printf escape sequences */
/* "%[0-9.l]*[a-zA-Z]" with alignment directives. Works for finite  */
/* & pretty-printable multi-object data. For example, object width  */
/* is limited by terminal width, height is limited by 9, etc.       */
/* For large or indefinite amount of data, use stream single-object */
/* fDump() that has brute-force screening and cutoffs.              */
/********************************************************************/

/* Maximal number of objects pieced together and respective maximal size */
#define MPRINT_MAX_OBJS  16
#define MPRINT_MAX_ROWS  9
#define MPRINT_MAX_COLS  TERMCHAR

/***********************************************************/
/* Add "%3M][ %lf ..]" control to printf format, where "3" */
/* denotes number of rows; first ] is termination symbol,  */
/* [ %lf ..] denotes how you would print out one row.      */
/***********************************************************/
#define MATRIX_ESC_CHAR      'M'    /* "%3M][ %lf ..]" */
/* alignment directives: finite state machine style */
#define ALIGN_TOP_ESC_CHAR   't'    /* "%t" */
#define ALIGN_MID_ESC_CHAR   'm'    /* "%m" */
#define ALIGN_BOT_ESC_CHAR   'q'    /* "%q" */
/* means starting from the fourth line from the top */
#define ALIGN_USR_ESC_CHAR   'a'    /* "%3a" */

void Mprintf (char *fmt, ... );
void Mfprintf (FILE *fp, char *fmt, ... );
void Mvprintf (FILE *fp, char *fmt, va_list ap);

/* Sole double matrix Mfprintf() driver: catches the lone "%M" in */
/* sole_fmt and automatically justifies elements column by column */

#define SPRINT_MAX_FMT   512    /* expanded fmt size */
#define SPRINT_MAX_BUF   64     /* double element string size */
#define SPRINT_MAX_ELM   16     /* max # of elements in one row */

void Sprint (FILE *fp, char *sole_fmt, double *A, int rows, int cols,
	     char open_char, char close_char, int open_space,
	     int close_space, int space_between);

#define Sfpr(fp,sole,A,rows,cols) Sprint(fp,sole,A,rows,cols,'|','|',1,1,2)
#define Spr(sole,A,rows,cols)     Sfpr(stdout,sole,A,rows,cols)
#define S2fpr(fp,sole,A)          Sfpr(fp,sole,A,2,2)
#define S2pr(sole,A)              S2fpr(stdout,sole,A)
#define S3fpr(fp,sole,A)          Sfpr(fp,sole,A,3,3)
#define S3pr(sole,A)              S3fpr(stdout,sole,A)
#define S6fpr(fp,sole,A)          Sfpr(fp,sole,A,6,6)
#define S6pr(sole,A)              S6fpr(stdout,sole,A)

#define SfPR(fp,sole,A,rows,cols) Sprint(fp,sole,A[0],rows,cols,'|','|',1,1,2)
#define SPR(sole,A,rows,cols)     Sfpr(stdout,sole,A[0],rows,cols)
#define S2fPR(fp,sole,A)          Sfpr(fp,sole,A[0],2,2)
#define S2PR(sole,A)              S2fpr(stdout,sole,A[0])
#define S3fPR(fp,sole,A)          Sfpr(fp,sole,A[0],3,3)
#define S3PR(sole,A)              S3fpr(stdout,sole,A[0])
#define S6fPR(fp,sole,A)          Sfpr(fp,sole,A[0],6,6)
#define S6PR(sole,A)              S6fpr(stdout,sole,A[0])

#define Vfpr(fp,sole,A,cols)      Sprint(fp,sole,A,1,cols,'(',')',1,1,2)
#define Vpr(sole,A,cols)          Vfpr(stdout,sole,A,cols)
#define V2fpr(fp,sole,A)          Vfpr(fp,sole,A,2)
#define V2pr(sole,A)              V2fpr(stdout,sole,A)
#define V3fpr(fp,sole,A)          Vfpr(fp,sole,A,3)
#define V3pr(sole,A)              V3fpr(stdout,sole,A)
#define V6fpr(fp,sole,A)          Vfpr(fp,sole,A,6)
#define V6pr(sole,A)              V6fpr(stdout,sole,A)

/* Driver of Sprint() to print A[][]*factor */
void Sprintmul (FILE *fp, char *sole_fmt, double *A, double factor,
                int rows, int cols, char open_char, char close_char,
                int open_space, int close_space, int space_between);

#define Sfprmul(fp,sole,A,factor,rows,cols) \
  Sprintmul(fp,sole,A,factor,rows,cols,'|','|',1,1,2)
#define Sprmul(sole,A,factor,rows,cols) Sfprmul(stdout,sole,A,factor,rows,cols)
#define S2fprmul(fp,sole,A,factor)      Sfprmul(fp,sole,A,factor,2,2)
#define S2prmul(sole,A,factor)          S2fprmul(stdout,sole,A,factor)
#define S3fprmul(fp,sole,A,factor)      Sfprmul(fp,sole,A,factor,3,3)
#define S3prmul(sole,A,factor)          S3fprmul(stdout,sole,A,factor)
#define SfPRmul(fp,sole,A,factor,rows,cols) \
  Sprintmul(fp,sole,A[0],factor,rows,cols,'|','|',1,1,2)
#define SPRmul(sole,A,factor,rows,cols) \
  Sfprmul(stdout,sole,A[0],factor,rows,cols)
#define S2fPRmul(fp,sole,A,factor)      Sfprmul(fp,sole,A[0],factor,2,2)
#define S2PRmul(sole,A,factor)          S2fprmul(stdout,sole,A[0],factor)
#define S3fPRmul(fp,sole,A,factor)      Sfprmul(fp,sole,A[0],factor,3,3)
#define S3PRmul(sole,A,factor)          S3fprmul(stdout,sole,A[0],factor)
#define S6fPRmul(fp,sole,A,factor)      Sfprmul(fp,sole,A[0],factor,6,6)
#define S6PRmul(sole,A,factor)          S6fprmul(stdout,sole,A[0],factor)
#define Vfprmul(fp,sole,A,factor,cols) \
  Sprintmul(fp,sole,A,factor,1,cols,'(',')',1,1,2)
#define Vprmul(sole,A,factor,cols)      Vfprmul(stdout,sole,A,factor,cols)
#define V2fprmul(fp,sole,A,factor)      Vfprmul(fp,sole,A,factor,2)
#define V2prmul(sole,A,factor)          V2fprmul(stdout,sole,A,factor)
#define V3fprmul(fp,sole,A,factor)      Vfprmul(fp,sole,A,factor,3)
#define V3prmul(sole,A,factor)          V3fprmul(stdout,sole,A,factor)

/* Sole int matrix Mfprintf() driver: catches the lone "%M" in    */
/* sole_fmt and automatically justifies elements column by column */

#define IPRINT_MAX_FMT   512    /* expanded fmt size */
#define IPRINT_MAX_BUF   64     /* double element string size */
#define IPRINT_MAX_ELM   16     /* max # of elements in one row */

void Iprint (FILE *fp, char *sole_fmt, int *A, int rows, int cols,
	     char open_char, char close_char, int open_space,
	     int close_space, int space_between);

#define Ifpr(fp,sole,A,rows,cols) Iprint(fp,sole,A,rows,cols,'|','|',1,1,2)
#define Ipr(sole,A,rows,cols)     Ifpr(stdout,sole,A,rows,cols)
#define IfPR(fp,sole,A,rows,cols) Iprint(fp,sole,A[0],rows,cols,'|','|',1,1,2)
#define IPR(sole,A,rows,cols)     Ifpr(stdout,sole,A[0],rows,cols)
#define ifpr(fp,sole,A,cols)      Iprint(fp,sole,A,1,cols,'(',')',1,1,2)
#define ipr(sole,A,cols)          ifpr(stdout,sole,A,cols)

/* fjoin.c: */

/* create an output stream that feeds to multiple open streams */
FILE *fjoin (int number_of_streams, ...);
/* delete the pipe and bury the dead child */
void fbreakup (FILE *fp);

/***************************************************************************/
/* 1. The above, like fopen(), is re-usable, i.e., many joins can co-exist */
/* 2. One can still write to an individual component after it was joined   */
/* but it is not guaranteed to be in the right "apparent order", To ensure */
/* that, use sync() before printing to the individual component.           */
/***************************************************************************/

/* system-wise global like errno, wonder whether that is wise... */
extern FILE *ft;
/* single pipe object handlers */
#define fbrk()        CAN( if (ft!=NULL) fbreakup(ft); ft = NULL; )
#define ftie(fp1,fp2) CAN( if (ft!=NULL) fbreakup(ft); ft = fjoin(2,fp1,fp2); )
#define ftie3(fp1,fp2,fp3) \
CAN( if (ft!=NULL) fbreakup(ft); ft=fjoin(3,fp1,fp2,fp3); )
#define ftie4(fp1,fp2,fp3,fp4) \
CAN( if (ft!=NULL) fbreakup(ft); ft=fjoin(4,fp1,fp2,fp3,fp4); )

    
/* jterm.c: */

/**************************************************************************/
/* open a dumb terminal at "DISPLAY" (NULL means default $DISPLAY), piped */
/* to open streams "in" (NULL means no) and "out" (NULL means no); return */
/* a unique identification number which will be needed in jclose().       */
/**************************************************************************/
pid_t jopen (FILE *in, FILE *out, char *DISPLAY);

/* close the pipe(s) and bury the child. If pid==0 we close all pipes */
void jclose (pid_t pid);

/*****************************************************************/
/* Dumb terminal emulation based on a tiny hack of XFree86 3.3.5 */
/* xterm, that pipes your "out" stream to X Window and keyboard  */
/* inputs to your "in" stream, at X display "DISPLAY". You need  */
/* to get the xterm source code and insert the following:        */
/*                                                               */
/* main.c/main() front:  argv[0] = "xterm";                      */
/*                                                               */
/* main.c/#ifndef AMOEBA/spawn() front:                          */
/*                                                               */
/*   int Ju_STDIN_FILENO, Ju_STDOUT_FILENO, Ju_pid, Ju_status;   */
/*   char Ju_c;                                                  */
/*   Ju_STDIN_FILENO = dup(STDIN_FILENO);                        */
/*   Ju_STDOUT_FILENO = dup(STDOUT_FILENO);                      */
/*                                                               */
/* before  if (command_to_exec) {                                */
/* 	        execvp(*command_to_exec, command_to_exec); .. :  */
/*                                                               */
/*   if ((Ju_pid=fork()) == 0)                                   */
/*   {                                                           */
/*     close (Ju_STDOUT_FILENO);                                 */
/*     while ( read(Ju_STDIN_FILENO, &Ju_c, 1) > 0 )             */
/*       while ( write(STDOUT_FILENO, &Ju_c, 1) == 0 );          */
/*     exit(0);                                                  */
/*   }                                                           */
/*   close (Ju_STDIN_FILENO);                                    */
/*   while ( read(STDIN_FILENO, &Ju_c, 1) > 0 )                  */
/*     while ( write(Ju_STDOUT_FILENO, &Ju_c, 1) == 0 );         */
/*   kill(Ju_pid,SIGKILL); waitpid (Ju_pid,&Ju_status,0);        */
/*   return(0);                                                  */
/*                                                               */
/* then compile and install as "jterm" in your run path.         */
/*                                                               */
/* jterms share all behaviors of xterm including ~/.Xdefaults    */
/*****************************************************************/

/* single window version */
extern pid_t JTERM_STDID;
#define Jopen(DISPLAY) CAN(Jclose(); JTERM_STDID=jopen(stdin,stdout,DISPLAY))
#define Jclose() CAN(if (JTERM_STDID>=0) jclose(JTERM_STDID); JTERM_STDID=-1)


/* bitmap.c */

/* Bitmap object: as if b[], which is char array, */
/* is a bit array indexed by i and has n members. */
typedef unsigned char Bmap;
#define BMAP_UNIT_IN_BITS  8
#define BMAP_UNIT_MASK     255
#define bunit(i)        ((i) >> 3)  /* which unit */
#define bslot(i)        ((i) &  7)  /* which bit of unit */
/* value of the ith bit of certain elementary type */
#define BIT(a,j)        ((a)&(1<<(j)))
#define CHAR_IN_BITS    (sizeof(char)*8)
#define INT_IN_BITS     (sizeof(int)*8)
#define DOUBLE_IN_BITS  (sizeof(double)*8)
#define LONG_IN_BITS    (sizeof(long)*8)
#define BITS_TO_BMAP_UNITS(n)  (bunit((n)-1)+1)  /* n should be > 0 */
/* n should be > 0, and the BYTES is rounded to be commensurate with Bmap */
#define BITS_TO_BYTES(n)       (BITS_TO_BMAP_UNITS(n)*sizeof(Bmap))
/* declare a Bitmap object of n bits */
#define BMAP(b,n)  Bmap b[BITS_TO_BMAP_UNITS(n)]
#define bmask(i)   (1 << bslot(i))  /* ..000100.. */
#define bnask(i)   (~bmask(i))       /* ..111011.. */
/* operations on a single Bitmap unit left value */
/* index j is wraparounded by BMAP_UNIT_IN_BITS  */
#define Bset(B,j)  (B|=bmask(j))   /* set bit to 1 */
#define Bclr(B,j)  (B&=bnask(j))   /* set bit to 0 */
#define Bval(B,j)  (B&bmask(j))
#define Bassign(B,j,val) if (val) Bset(B,j); else Bclr(B,j)
/* accessing certain unit of a Bitmap array */
#define Bhost(b,i)  (((Bmap *)(b))+bunit(i))  /* pointer */
#define BHOST(b,i)  (*Bhost(b,i))             /* l&r value */
/* operations on a Bitmap array */
#define BSET(b,i)  Bset(BHOST(b,i),i)
#define BCLR(b,i)  Bclr(BHOST(b,i),i)
#define BVAL(b,i)  Bval(BHOST(b,i),i)
#define BASSIGN(b,i,val)  Bassign(BHOST(b,i),i,val)

/* copy >=n bits from Bitmap array src to dest: rounded to last Bitmap unit */
#define BCPY(dest,src,n) \
((Bmap *)memcpy((void *)(dest),(void *)(src),BITS_TO_BYTES(n))

/* set >=n bits of Bitmap array b to 0: rounded to last Bitmap unit */
#define BZERO(b,n) (memset((void *)(b),0,BITS_TO_BYTES(n)))

/* set >=n bits of Bitmap array b to 1: rounded to last Bitmap unit */
#define B_ONE(b,n) \
(memset((void *)(b),BMAP_UNIT_MASK,BITS_TO_BYTES(n)))

/* set >=n bits of Bitmap array b to val: rounded to last Bitmap unit */
#define BASSIGNALL(b,n,val)  if (val) B_ONE(b,n); else BZERO(b,n)

/* Above are mainly for Bitmap array initialization, */
/* since they could contaminate the last Bitmap unit */

/* byte-wise assign units of Bitmap array b to bytevalue */
#define BBYTEASSIGN(b,units,bytevalue) \
 (memset((void *)(b),(int)(bytevalue),(size_t)(units)*sizeof(Bmap)))
/* equivalent to Bperiod() with period = 8 */

/* return pointer to a bitmap of n bits (rounded) */
Bmap *Balloc (int n);

/* return pointer to a bitmap of n bits (rounded) which are all set to 0 */
Bmap *BALLOC (int n);

/* return pointer to a bitmap of n bits (rounded) which are */
/* all set to 0 or 1, depending on whether val is 0 or not. */
Bmap *BAlloc (int n, int val);

/* "..." means set of (Bmap *b, int offset, int num_bits), escape seq. is */
/* "%B". Print out num_bits from Bitmap array b, starting from bit offset */
void Bfprintf (FILE *out, char *fmt, ...);
void Bprintf (char *fmt, ...);
void Bvfprintf (FILE *out, char *fmt, va_list ap);
#define Bfdump(out,b,offset,num_bits) \
  Bfprintf(out,"\n%B\n\n",(Bmap *)(b),(int)(offset),(int)(num_bits))
#define Bdump(b,offset,num_bits) Bfdump(stdout,b,offset,num_bits)

/* Bfwrite: reverse bitstream direction to emulate human handwriting */
void Bfwrite (FILE *out, char *c, int length_in_bytes);
#define Bfwritechar(out,s) Bfwrite(out,CHARP(s),sizeof(char))
#define Bwritechar(s) Bfwritechar(stdout,s)
#define BfwriteChar(out,S) Bfwritechar(out,&(S))
#define BwriteChar(S) BfwriteChar(stdout,S)
#define Bfwriteshort(out,s) Bfwrite(out,CHARP(s),sizeof(short))
#define Bwriteshort(s) Bfwriteshort(stdout,s)
#define BfwriteShort(out,S) Bfwriteshort(out,&(S))
#define BwriteShort(S) BfwriteShort(stdout,S)
#define Bfwriteint(out,s) Bfwrite(out,CHARP(s),sizeof(int))
#define Bwriteint(s) Bfwriteint(stdout,s)
#define BfwriteInt(out,S) Bfwriteint(out,&(S))
#define BwriteInt(S) BfwriteInt(stdout,S)
#define Bfwritelong(out,s) Bfwrite(out,CHARP(s),sizeof(long int))
#define Bwritelong(s) Bfwritelong(stdout,s)
#define BfwriteLong(out,S) Bfwritelong(out,&(S))
#define BwriteLong(S) BfwriteLong(stdout,S)
#define Bfwritedouble(out,s) Bfwrite(out,CHARP(s),sizeof(double))
#define Bwritedouble(s) Bfwritedouble(stdout,s)
#define BfwriteDouble(out,S) Bfwritedouble(out,&(S))
#define BwriteDouble(S) BfwriteDouble(stdout,S)
/*********************************************************/
/* When we run "i=0x020103e0; BwriteInt(i)"; we get,     */
/* "00000010 00000001 00000011 11100000".                */
/* In the printout, rightmost 8-bits are the first byte, */
/* leftmost are the 4th byte. Bwrite accommodates human  */
/* convention, for example, when we write 0xabff, ff is  */
/* the first byte, ab is the second byte. Left / right   */
/* shifts are with respect to this rather weird human    */
/* convention. On the other hand, text stream usually    */
/* flows from left to right in increasing index. Thus    */
/* this invariably will cause some confusion. Bdump()    */
/* printout is the exact bitwise mirror of Bwrite(),     */
/* and because its origin is on the left, Bdump() is     */
/* more suitable for indefinite bitstream printout.      */
/*********************************************************/


/* memusage.c: */

/* Memory usage snapshot */

#if \
  defined(_Linux) || \
  defined(_SunOS) || \
  defined(_alpha) || \
  defined(_Darwin)

/* calling Linux procps command "ps m": */
long memusage();

#else

extern struct rusage mem_usage;
#define memusage() \
  (getrusage(RUSAGE_SELF,&mem_usage),mem_usage.ru_maxrss*1000)

#endif

/* quick view of memory */
char *strmem (long bytes);
#define fpmemusage(fp) { if (fp) fprintf(fp, \
  "This process has used up to %s.\n", strmem(memusage())); }
#define pmemusage() fpmemusage(stdout)


/* int8_set.c: */

/********************************************************/
/* Duke it out with memset/bzero: and supporting longer */
/* period patterns like int2, int4, int8, real4, real8. */
/********************************************************/

/****************************************************************/
/* int8_aset() family requires destination pointer be aligned   */
/* with QWORD; however after adding the "count" it does not     */
/* have to, as it requires very little work to finish the tail. */
/* "int8" also suggests internal register implementation.       */
/****************************************************************/

/* Copy "char_count" x (char)"c" to "dest". "dest" must be aligned to QWORD */
void char_int8_aset (char *dest, char c, size_t char_count);

/* Copy "char_count" x (char)"c" to "dest" using eight int8 registers */
void char_int8_set (char *dest, char c, size_t char_count);

/* Copy "uchar_count" x (unsigned char)"c" to */
/* "dest". "dest" must be aligned to QWORD.   */
void uchar_int8_aset
(unsigned char *dest, unsigned char c, size_t uchar_count);

/* Copy "uchar_count" x (unsigned char)"c" */
/* to "dest" using eight int8 registers.   */
void uchar_int8_set (unsigned char *dest, unsigned char c, size_t uchar_count);

/* Copy "int2_count" x (int2)"c" to "dest". "dest" must be */
/* aligned to QWORD, which is generally a good idea anyway */
void int2_int8_aset (int2 *dest, int2 c, size_t int2_count);

/* Copy "int2_count" of (int2)"c" to "dest" using eight int8  */
/* registers. "dest" must minimally be aligned to int2 (WORD) */
void int2_int8_set (int2 *dest, int2 c, size_t int2_count);

/* Copy "int4_count" x (int4)"c" to "dest". "dest" must be */
/* aligned to QWORD, which is generally a good idea anyway */
void int4_int8_aset (int4 *dest, int4 c, size_t int4_count);

/* Copy "int4_count" of (int4)"c" to "dest" using eight int8   */
/* registers. "dest" must minimally be aligned to int4 (DWORD) */
void int4_int8_set (int4 *dest, int4 c, size_t int4_count);

/* Copy "int8_count" x (int8)"c" to "dest" using eight int8    */
/* registers. "dest" must minimally be aligned to int8 (QWORD) */
void int8_int8_set (int8 *dest, int8 c, size_t int8_count);

/** REGS isomorphic booty: **/

/* Copy "real8_count" x (real8)"c" to "dest" using eight real8  */
/* registers. "dest" must minimally be aligned to real8 (QWORD) */
void real8_real8_set (real8 *dest, real8 c, size_t real8_count);
#define real8set(dest, r8, real8_count) \
  real8_real8_set((real8 *)(dest), (real8)(r8), (size_t)(real8_count))

/* Copy "real4_count" x (real4)"c" to "dest" using eight real4  */
/* registers. "dest" must minimally be aligned to real4 (DWORD) */
void real4_real4_set (real4 *dest, real4 c, size_t real4_count);
/**************************************************************************/
/* Could not do better because if we use real8, we have to union it with  */
/* two real4's first as in int8_int4_set(). But it may lead to floating   */
/* point exception as we try to store it to a register, since the joint   */
/* or decomposition of an integer is always integer but could be NaN for  */
/* floating points. This hardware exception occurs on Pentium even when   */
/* no arithmetic ops are carried out and we are just loading the meory to */
/* register. On Origin 2000 this does not occur, still it's bad practice. */
/**************************************************************************/
#define real4set(dest, r4, real4_count) \
  real4_real4_set((real4 *)(dest), (real4)(r4), (size_t)(real4_count))


/* int4_set.c: */

/********************************************************/
/* Duke it out with memset/bzero: and supporting longer */
/* period patterns like int2, int4 and int8.            */
/********************************************************/

/****************************************************************/
/* int4_aset() family requires destination pointer be aligned   */
/* with DWORD; however after adding the "count" it does not     */
/* have to, as it requires very little work to finish the tail. */
/* "int4" also suggests internal register implementation.       */
/****************************************************************/

/* Copy "char_count" x (char)"c" to "dest". "dest" must be aligned to DWORD */
void char_int4_aset (char *dest, char c, size_t char_count);

/* Copy "char_count" x (char)"c" to "dest" using four int4 registers */
void char_int4_set (char *dest, char c, size_t char_count);

/* Copy "uchar_count" x (unsigned char)"c" to */
/* "dest". "dest" must be aligned to DWORD.   */
void uchar_int4_aset(unsigned char *dest, unsigned char c, size_t uchar_count);

/* Copy "uchar_count" x (unsigned char)"c" */
/* to "dest" using four int4 registers.    */
void uchar_int4_set(unsigned char *dest, unsigned char c, size_t uchar_count);

/* Copy "int2_count" x (int2)"c" to "dest". "dest" must be */
/* aligned to DWORD, which is generally a good idea anyway */
void int2_int4_aset (int2 *dest, int2 c, size_t int2_count);

/* Copy "int2_count" of (int2)"c" to "dest" using four int4   */
/* registers. "dest" must minimally be aligned to int2 (WORD) */
void int2_int4_set (int2 *dest, int2 c, size_t int2_count);

/* Copy "int4_count" x (int4)"c" to "dest" using four int4     */
/* registers. "dest" must minimally be aligned to int4 (DWORD) */
void int4_int4_set (int4 *dest, int4 c, size_t int4_count);

/* Copy "int8_count" x (int8)"c" to "dest" using four int4     */
/* registers. "dest" must minimally be aligned to int4 (DWORD) */
void int8_int4_set (int8 *dest, int8 c, size_t int8_count);


/** redirect transportable prototypes **/

#ifdef _IRIX64                        /** SGI Origin 2000 **/

#define SHORTset(dest, i2, int2_count)   int2set(dest, i2, int2_count)
#define INTset(dest, i4, int4_count)     int4set(dest, i4, int4_count)
#define LONGset(dest, i8, int8_count)    int8set(dest, i8, int8_count)
#define FLOATset(dest, r4, real4_count)  real4set(dest, r4, real4_count)
#define DOUBLEset(dest, r8, real8_count) real8set(dest, r8, real8_count)

#else                                 /** Pentium II **/

#define SHORTset(dest, i2, int2_count)   int2set(dest, i2, int2_count)
#define INTset(dest, i4, int4_count)     int4set(dest, i4, int4_count)
#define LONGset(dest, i4, int4_count)    int4set(dest, i4, int4_count) /* ! */
#define FLOATset(dest, r4, real4_count)  real4set(dest, r4, real4_count)
#define DOUBLEset(dest, r8, real8_count) real8set(dest, r8, real8_count)

#endif


/* readline.c: */

/* Read a string, and return a pointer to it. Returns NULL on EOF. */
char *readline_gets (char *prompt, char *default_answer);
#define readline_Gets(prompt) readline_gets(prompt,NULL);
#define readline_GETS() readline_gets(NULL,NULL);


/* filelist.c: */

/* Return a list of "numerically similar" filenames to "seed_fname" */
/* and sort them; return 0 if no digits are found in "seed_fname".  */
/* Remember to free dynamical memory later using globfree(glob_t*). */
int numerically_sorted_glob (char *seed_fname, glob_t *globbuf);

/* Find the fname that is numerically similar to "start_fname" */
/* but is "glob_advance" down the list */
char *numerically_sorted_glob_advance
(char *start_fname, glob_t *globbuf, int glob_advance);

/* memory efficient implementation */
char *Numerically_sorted_glob_advance (char *fname, int glob_advance);

/* seek the first in the similarity family */
char *Numerically_sorted_glob_first (char *fname);

/* seek the last in the similarity family */
char *Numerically_sorted_glob_last (char *fname);


/* gnuplot.c */

/****************************************************/
/* pipe graphics scripts to powerful gnuplot engine */
/****************************************************/

typedef struct
{
    char *gnuplotoptions;
    FILE *pipestream;
    char *epsfilename;
    char *epsoptions;
} Gnuplot;

#define GNUPLOT_NOT_OK              0
#define GNUPLOT_OK                  1
#define GNUPLOT_DEFAULT_OPTIONS     "-noraise"
#define GNUPLOT_DEFAULT_EPSFILENAME "gnuplot.eps"
/* Note that if more than one Gnuplot pipes are open and both are */
/* using the same default eps filename, this could lead to racing */
#define GNUPLOT_DEFAULT_EPSOPTIONS  "enhanced color"

/* Close pipe, kill gnuplot, free memory, reset pointers */
void gnuplot_free (Gnuplot *G);

/* Set eps filename: if NULL the eps filename default will be presumed */
void gnuplot_epsfilename (Gnuplot *G, char *epsfilename);
#define GNUPLOT_EPSFILENAME_NOW(G)  \
  (NOTNULL((G)->epsfilename)?(G)->epsfilename:GNUPLOT_DEFAULT_EPSFILENAME)
/* Set the eps filename to be default */
#define gnuplot_default_epsfilename(G) gnuplot_epsfilename(G,NULL)

/* Open gnuplot window. Returns GNUPLOT_NOT_OK or GNUPLOT_OK */
int gnuplot_recreate (Gnuplot *G, char *epsfilename);
/* Open gnuplot window with default eps filename */
#define gnuplot_Recreate(G) gnuplot_recreate(G,NULL)

/* Ask gnuplot to run a script file */
void gnuplot_run_file (Gnuplot *G, char *script_fname);

/* Ask gnuplot to run a script character-buffer */
void gnuplot_run_buffer (Gnuplot *G, char *buffer);

/* Ask gnuplot to run print-outs */
void gnuplot_run_printf (Gnuplot *G, char *format, ...);

/* Set gnuplot output to X11 */
void gnuplot_term_x11 (Gnuplot *G);

/* Set gnuplot output to eps: if epsfilename=NULL, the old */
/* G->epsfilename is used; else the new name will be used. */
void gnuplot_term_eps (Gnuplot *G, char *epsfilename);
/* Set gnuplot output to eps with G->epsfilename unchanged */
#define gnuplot_term_Eps(G) gnuplot_term_eps(G,NULL)

/**************************************************************/
/* Synchronize with the child gnuplot process by asking it to */
/* save the drawing into an eps file and make sure of getting */
/* it. If epsfilename=NULL, the old G->epsfilename will be    */
/* used; else the new epsfilename will be used.               */
/**************************************************************/
void gnuplot_eps_sync (Gnuplot *G, char *epsfilename);
/* Synchronize using the old unchanged G->epsfilename */
#define gnuplot_eps_Sync(G) gnuplot_eps_sync(G,NULL)

#endif  /* _IO_h */
