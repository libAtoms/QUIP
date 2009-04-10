/*********************************************/
/* libVecMat: -llapck -lblas -lm             */
/*            -lScalar -lIO                  */
/*                                           */
/* General Vector & Matrix Operation Library */
/*                                           */
/* Nov 19 1999 Ju Li <liju99@mit.edu>        */
/*********************************************/

#ifndef _VecMat_h
#define _VecMat_h

#include <IO.h>
#include <Scalar.h>

/* local static memory / global static.. dynamic memory hierarchy */
#define FREE(ptr,array,backup_free) \
 {if (((void *)(ptr))!=(void *)array) backup_free(ptr);}

#define Imem(size_in_ints,int_array,backup_alloc) ((size_in_ints <= \
sizeof(int_array)/sizeof(int))?int_array : backup_alloc(size_in_ints))

#define Dmem(size_in_dbles,dble_array,backup_alloc) ((size_in_dbles <= \
sizeof(dble_array)/sizeof(double))?dble_array : backup_alloc(size_in_dbles))

#define Bmem(size_in_bits,bit_array,backup_alloc) \
 ((BITS_TO_BYTES(size_in_bits) <= sizeof(bit_array))? \
  bit_array : backup_alloc(size_in_bits))


/* Ialloc.c: */

/* integer memory allocation and access */

/* return pointer to an array of n ints */
int *Ialloc (int n);

/* return pointer to an array of n ints which are all cleared to 0 */
int *IALLOC (int n);

/* return pointer to an array of n ints which are all set to flag */
int *IAlloc (int n, int flag);

/* a[] := 0 (int); */
#define IZERO(n,a) bzero((void *)(a),(n)*sizeof(int))

/* b[] := a[] (int); then return b[] */
#define IEQV(n,b,a) ((int *)memcpy((void *)(b),(void *)(a),(n)*sizeof(int)))


/* ordering.c: */

/********************************************/
/* elementary order, swapping & permutation */
/********************************************/

/* assign index array idx[] := a..b (1..10, 10..1); return idx[] */
int *sequentially_index (int idx[], int a, int b);

/* assign index array idx[] := a:b:c; return idx[] */
int *sequentially_Index (int idx[], int a, int b, int c);

/* assign index array idx[] := 0..N-1 */
void Sequentially_index (int N, int idx[]);

/* assign index array idx[] := a, a+b, a+2*b, .., a+(N-1)*b; a=0..b-1 */
void Sequentially_Index (int N, int a, int b, int idx[]);

/* swap contents of *a and *b, each of size_in_bytes bytes */
void swapb (char *a, char *b, size_t size_in_bytes);
#define SWAPB(a,b,c) swapb((char *)(a),(char *)(b),(size_t)(c))

/* swap contents of *a and *b, each of size_in_ints integers */
void swapi (int *a, int *b, size_t size_in_ints);
#define SWAPI(a,b,c) swapi((int *)(a),(int *)(b),(size_t)(c))

/* swap contents of *a and *b, each of size_in_doubles doubles */
void swapd (double *a, double *b, size_t size_in_doubles);
#define SWAPD(a,b,c) swapd((double *)(a),(double *)(b),(size_t)(c))

/* randomly permute array a[min..max] (e.g. a[0..n-1]),  */
/* each element, like a[max], is of size_in_bytes bytes. */
char *r_permuteb (char *a, int min, int max, size_t size_in_bytes);
#define R_PERMUTEB(a,min,max,c) \
((char *)r_permuteb((char *)(a),(int)(min),(int)(max),(size_t)(c)))

/* randomly permute array a[min..max] (e.g. a[0..n-1]), */
/* each element, like a[max], is of size_in_ints ints.  */
int *r_permutei (int *a, int min, int max, size_t size_in_ints);

/* randomly permute array a[min..max] (e.g. a[0..n-1]), */
/* each element, like a[max], is of size_in_doubles doubles. */
double *r_permuted (double *a, int min, int max, size_t size_in_doubles);

/* randomly permute index array a[min..max] (e.g. a[0..n-1]) */
int *r_permute (int *a, int min, int max);

/* determine if values in idx[] is a permutation of min..max */
#define is_permutation_STATIC_BITS (1024L*32)
int is_permutation (int idx[], int min, int max);

/* b[i] := a[idx[i]], i=0..N-1, each a[i] is object of size_in_chars chars */
void rearrange (int N, size_t size_in_chars, char a[], int idx[], char b[]);
#define REARRANGE(N,size_in_chars,a,idx,b) \
rearrange(N,size_in_chars,(char *)(a),idx,(char *)(b))

/* a[] := a[idx[]], each a[i] is an object of size_in_chars chars */
void Rearrange (int N, size_t size_in_chars, char a[], int idx[]);

/* b[i] := a[idx[i]], i=0..N-1, each a[i] is object of size_in_ints ints */
void rearrangei (int N, size_t size_in_ints, int a[], int idx[], int b[]);

/* a[] := a[idx[]], each a[i] is an object of size_in_ints ints */
void Rearrangei (int N, size_t size_in_ints, int a[], int idx[]);

/* b[i] := a[idx[i]], i=0..N-1, each a[i] is object of size_in_dbles doubles */
void rearranged
(int N, size_t size_in_dbles, double a[], int idx[], double b[]);

/* a[] := a[idx[]], each a[i] is an object of size_in_dbles doubles */
void Rearranged (int N, size_t size_in_dbles, double a[], int idx[]);

/* returns max(a[]) */
int IMAX (int n, int a[]);

/* returns index imax, such that a[imax] >= a[] */
int Imax (int n, int a[]);

/* returns min(a[]) */
int IMIN (int n, int a[]);

/* returns index imin, such that a[imin] <= a[] */
int Imin (int n, int a[]);

/* summing all elements in a[] */
int Isum (int n, int a[]);

/* see sorters.c: qsort_glibc */
void Iqsort_glibc (int N,  int x[], int idx[], int option);
#define IqSort_glibc(N,x,idx,a,b) ( Sequentially_Index(N,a,b,idx), \
  Iqsort_glibc(N,x,idx,USE_OLD_IDX) )

/* see sorters.c: qsort_numerical_recipes */
void Iqsort_numerical_recipes (int N,  int x[], int idx[], int option);
#define IqSort_numerical_recipes(N,x,idx,a,b) ( Sequentially_Index(N,a,b,idx),\
  Iqsort_numerical_recipes(N,x,idx,USE_OLD_IDX) )

/* determine whether x[idx[i]], i=0..N-1, is non-decreasing */
int Is_nondecreasing (int N, int x[], int idx[]);

/* Find out if there are identical elements in an integer array */
bool Iall_different (int N,  int x[]);


/* sorters.c: */

/**********************************************/
/* fast sorting routines for array of doubles */
/**********************************************/

/* determine whether x[idx[i]], i=0..N-1, is non-decreasing */
int is_nondecreasing (int N, double x[], int idx[]);

/* rearrange array of doubles: y[i] := x[idx[i]], i=0..N-1 */
void Vrearrange (int N, double x[], int idx[], double y[]);

/* rearrange array of doubles: x[] := x[idx[]], then set idx[] := 0..N-1 */
void VRearrange (int N, double x[], int idx[]);

#define USE_OLD_IDX  0  /* incoming index permutation but not sequential */
#define USE_NEW_IDX  1  /* incoming index undefined: must initialize */
#define USE_SEQ_IDX  2  /* incoming index sequential */

/* non-recursive quicksort category:  */
void qsort_glibc (int N,  double x[], int idx[], int option);
#define qSort_glibc(N,x,idx,a,b) ( Sequentially_Index(N,a,b,idx), \
  qsort_glibc(N,x,idx,USE_OLD_IDX) )

void qsort_numerical_recipes (int N,  double x[], int idx[], int option);
#define qSort_numerical_recipes(N,x,idx,a,b) ( Sequentially_Index(N,a,b,idx), \
  qsort_numerical_recipes(N,x,idx,USE_OLD_IDX) )

/* recursive mergesort category: */
void msort_Lee (int N,  double x[], int idx[], int option);

/* non-recursive mergesort category: */
void msort_Ju (int N,  double x[], int idx[], int option);

/* Instead of sorting x[0..N-1], sort x[m*(0..N-1)+n] with n=0..m-1 */


/* neighborlist.c: */

/******************************************************************/
/* 1D memory tectonics to represent 2D structure of neighbor list */
/******************************************************************/

/* clear all records in an uncompressed list */
void Iclearlist (int *idx, int n);

/* return a one-dimensional list[] evenly indexed by (*idx)[] */
/* where list[(*idx)[2*i]<=k<(*idx)[2*i+1]] belong to owner i */
/* total n owners, each "maximally" keeps likely_size entries */
int *Icreatelist(int **idx, int n, int likely_size);

/* Same as Icreatelist() except using realloc() instead of malloc() */
int *Irecreatelist(int **idx, int n, int likely_size);

/* Kill memory holes, free extra space, and change */
/* [2*i]-[2*i+1] notation to [i]-[i+1] notation.   */
int *Icompresslist (int **idx, int n);

/* Insert memory holes of size "pad_each" to each entry and  */
/* change [i]-[i+1] notation back to [2*i]-[2*i+1] notation. */
int *Idecompresslist (int **idx, int n, int pad_each);

/* free memory allocated by Icreatelist() / Irecreatelist() and set NULL */
void Ifreelist(int **idx, int **list);


/***************************************************************/
/* Above are for combined idx/list allocation models which can */
/* be regarded as a subclass. Next are more general operators. */
/***************************************************************/

/* print allocation and occupation statistics */
void Ilist_statistics
(char *name, int idx[], int n, bool compressed, FILE *out);

/* Without changing idx[2*min] and idx[2*max], solve the memory */
/* conflict at idx[2*i+1], idx[2*i+2]. You must be sure there   */
/* is conflict idx[2*i+1]==idx[2*i+2] before calling Imakespace */
/* Imakespace (faster than IMAKESPACE) does NOT preserve order. */
void Imakespace (int idx[], int list[], int i, int min, int max);

/* Safely append value to the end of owner i's (unordered) list; */
/* original entry order of this and other owners' may be altered */
void Iappend (int idx[], int list[], int value, int i, int min, int max);

/* Without changing idx[2*min] and idx[2*max], solve the memory */
/* conflict at idx[2*i+1], idx[2*i+2]. You must be sure there   */
/* is conflict idx[2*i+1]==idx[2*i+2] before calling IMAKESPACE */
/* IMAKESPACE (slower than Imakespace) preserves list[] order.  */
void IMAKESPACE (int idx[], int list[], int i, int min, int max);

/* safely append value to the end of owner i's (ordered) list */
void IAPPEND (int idx[], int list[], int value, int i, int min, int max);



/* Mmem.c: */

/*******************************************/
/* Matrix scale memory, a.k.a. mega-memory */
/*                                         */
/* scenarios expected:                     */
/*                                         */
/* a) Hermitian matrix diagonalization.    */
/* b) Complex matrix inversion.            */
/* c) Huge array merge sort indexing.      */
/* d) Huge (I mean huge) bitmap.           */
/*******************************************/

/* fast matrix (static.. dynamic) memory manager */

/* max complex double square matrices w/o dynamic allocation */
#define MAX_READY_CMATRIX_RANK   (32)
#define MAX_READY_CMATRIX_DBLES  (2*SQUARE(MAX_READY_CMATRIX_RANK))
/* zhpev_ main memory for Hermitian matrix and eigenvectors */
#define MMEM_SEGS         2
#define MMEM_UNIT         MAX_READY_CMATRIX_DBLES
#define MMEM_MAX_CLIENTS  6

double *Mmem (size_t size_in_dbles);
void Mfree (double *ptr);
void Mfreeall();
/* clone an array of size_in_dbles using Mmem() memory */
double *Mclone (double *ptr, size_t size_in_dbles);
/* remember to free it later with Mfree() */

/* mega-memory support for other types */
#define MImem(size_in_ints) \
((int *)Mmem(SEIL((size_in_ints)*sizeof(int), sizeof(double))))
#define MCmem(size_in_chars) ((char *) \
Mmem(SEIL((size_in_chars)*sizeof(char), sizeof(double))))
#define MBmem(size_in_bits) ((Bmap *) \
Mmem(SEIL(BITS_TO_BYTES(size_in_bits), sizeof(double))))
#define MFREE(ptr) Mfree((double *)(ptr))

/* Vmem.c: */

/* double type vector allocation and memory manager */

/* return pointer to an array of n doubles */
double *Valloc (int n);

/* return pointer to an array of n doubles which are all cleared to 0 */
double *VALLOC (int n);

/* return pointer to an array of n doubles which are all set to "value" */
double *VAlloc (int n, double value);

/* return pointer to an array of n doubles (reallocated) */
double *Vrealloc (double *a, int n);

/* a[] := 0. */
#define VZERO(n,a) bzero((void *)(a),(n)*sizeof(double))

/****************************************************************/
/* (eigen)Vector scale memory, or, sub-mega memory.             */
/*                                                              */
/* scenarios expected:                                          */
/*                                                              */
/* a) single eigenvector from Hermitian matrix diagonalization. */
/* b) small matrices, say, k-space dynamical matrix of quartz.  */
/* c) medium sized bitmap for ordinary tagging purposes.        */
/****************************************************************/

/* fast vector (static.. dynamic) memory manager */

/* max complex double vector w/o dynamic allocation */
#define MAX_READY_CVECTOR_RANK   (16*MAX_READY_CMATRIX_RANK)
#define MAX_READY_CVECTOR_DBLES  (2*MAX_READY_CVECTOR_RANK)
/* zhpev_: subsidiary memory for vectors and pivot indices */
#define VMEM_SEGS         8
#define VMEM_UNIT         MAX_READY_CVECTOR_DBLES
#define VMEM_MAX_CLIENTS  8

double *Vmem (size_t size_in_dbles);
void Vfree (double *ptr);
void Vfreeall();
/* clone an array of size_in_dbles using Vmem() memory */
double *Vclone (double *ptr, size_t size_in_dbles);
/* remember to free it later with Vfree() */

/* mega-memory support for other types */
#define VImem(size_in_ints) \
((int *)Vmem(SEIL((size_in_ints)*sizeof(int), sizeof(double))))
#define VCmem(size_in_chars) ((char *) \
Vmem(SEIL((size_in_chars)*sizeof(char), sizeof(double))))
#define VBmem(size_in_bits) ((Bmap *) \
Vmem(SEIL(BITS_TO_BYTES(size_in_bits), sizeof(double))))
#define VFREE(ptr) Vfree((double *)(ptr))

/* Vector.c: */

/* double type vector operations */

/* x[0:n-1] := Frandom() [uniform on (0,1)]; then return x[] */
double *Vfrandom (int n, double x[]);
#define VFrandom(n,x,i) for ((i)=(n); (i)--;) x[i]=Frandom();

/* x[0:n-1] := FRANDOM() [uniform on (-0.5,0.5)]; then return x[] */
double *VFRANDOM (int n, double x[]);
#define VFrandoM(n,x,i) for ((i)=(n); (i)--;) x[i]=FRANDOM();

/* x[0:n-1] := Frandnorm(E,sigma2); then return x[] */
double *Vfrandnorm (int n, double E, double sigma2, double x[]);

/* x[0:n-1] := Frandnorm(0,1); then return x[] */
double *VFRANDNORM (int n, double x[]);

/* b[] := a[]; then return b[] */
#define VEQV(n,a,b) \
  ((double *)memcpy((void *)(b),(void *)(a),(n)*sizeof(double)))

/* b[] := a[] */
#define VCLONE(n,a,b)    ((b)=Valloc(n), VEQV(n,a,b))
#define VRECLONE(n,a,b)  ((b)=Vrealloc(b,n), VEQV(n,a,b))

/* return the Euclidean norm squared := |x|^2 */
double Vlength2 (int n, double x[]);
/* length2 := |a|^2 */
#define VLENGTH2(n,a,i,length2) \
  for ((length2)=0,(i)=(n); (i)--;) (length2) += SQUARE((a)[i]);

/* length := |a| */
#define Vlength(n,x) sqrt(Vlength2((n),(x)))
#define VLENGTH(n,a,i,length) { \
  for ((length)=0,(i)=(n); (i)--;) (length) += SQUARE((a)[i]); \
  length = sqrt(length); }

/* return the Euclidean distance squared := |a-b|^2 */
double Vdistance2 (int n, double a[], double b[]);
#define Vdistance(n,a,b) sqrt(Vdistance2((n),(a),(b)))

/* normalize x[] to unit length; then return x[] */
double *Vnormalize (int n, double x[]);
#define VNORMALIZE(n,x,length) ((length)=Vlength(n,x),VDIV(n,x,(length)))

/* diff2 := |a-b|^2 */
#define VDIFF2(n,a,b,i,diff2) \
  for ((diff2)=0,(i)=(n); (i)--;) (diff2) += SQUARE((a)[i]-(b)[i]);

/* diff := |a-b| */
#define VDIFF(n,a,b,i,diff) { \
  for ((diff)=0,(i)=(n); (i)--;) (diff) += SQUARE((a)[i]-(b)[i]); \
  diff = sqrt(diff); }

/* b[] := -a[]; then return b[] */
double *Vneg (int n, double *a, double *b);

/* a[] := -a[]; then return a[] */
double *VNeg (int n, double *a);

/* make b[] an image of a[] in [0,1)^3; then return b[] */
double *Vtrim (int n, double a[], double b[]);

/* change a[] to its own image in [0,1)^3; then return a[] */
double *VTRIM (int n, double a[]);
#define VTriM(n,a,i) for ((i)=(n); (i)--;) Trim((a)[i]);

/* make b[] image of a[]'s in [-0.5,0.5)^3; then return b[] */
double *Vimage (int n, double a[], double b[]);

/* change a[] to its own image in [-0.5,0.5)^3; then return a[] */
double *VIMAGE (int n, double a[]);

/* returns max(a[]) */
double VMAX (int n, double a[]);
/* a[] must have at least one element */
#define VMaX(n,a,i,max) { for ((max)=(a)[0],(i)=1; (i)<(n); (i)++) \
  if ((a)[i]>(max)) (max)=(a)[i]; }

/* returns index imax, such that a[imax] >= a[] */
int Vmax (int n, double a[]);

/* *vmax := max(a[]); and returns the index imax, s.t. a[imax] = max(a[]) */
int VMax (int n, double a[], double *vmax);

/* returns min(a[]) */
#define VMiN(n,a,i,min) { for ((min)=(a)[0],(i)=1; (i)<(n); (i)++) \
  if ((a)[i]<(min)) (min)=(a)[i]; }

/* returns index imin, such that a[imin] <= a[] */
int Vmin (int n, double a[]);

/* *vmin := min(a[]); and returns the index imin, s.t. a[imin] = min(a[]) */
int VMin (int n, double a[], double *vmin);

/* summing all elements in a[] */
double Vsum (int n, double a[]);

/* the algebraic average of a[] */
#define Vavg(n,a) (Vsum(n,a)/(n))

/* scalar & vector operations */

/* b[] := multiplier * a[]; then return b[] */
double *Vmul (int n, double multiplier, double a[], double b[]);
#define VMul(n,multiplier,a,b,i) \
  for ((i)=(n); (i)--;) (b)[i] = (multiplier) * (a)[i];
#define VCLONEMul(n,multiplier,a,b,i) \
  { (b)=Valloc(n); VMul(n,multiplier,a,b,i); }
#define VRECLONEMul(n,multiplier,a,b,i) \
  { (b)=Vrealloc(b,n); VMul(n,multiplier,a,b,i); }

/* a[] := multiplier * a[]; then return a[] */
double *VMUL (int n, double multiplier, double a[]);
#define VMuL(n,multiplier,a,i) for ((i)=(n); (i)--;) (a)[i]*=(multiplier);

/* b[] := a[] / divisor; then return b[] */
double *Vdiv (int n, double a[], double divisor, double b[]);

/* a[] := a[] / divisor; then return a[] */
double *VDIV (int n, double a[], double divisor);
#define VDiV(n,a,divisor,i) for ((i)=(n); (i)--;) (a)[i]/=(divisor);


/*******************************/
/* vector & vector with scalar */
/*******************************/

/* c[] := a[] + b[]; then return c[] */
double *Vadd (int n, double a[], double b[], double c[]);

/* b[] := a[] + b[]; then return b[] */
double *VADD (int n, double a[], double b[]);

/* c[] := a[] - b[]; then return c[] */
double *Vsub (int n, double a[], double b[], double c[]);

/* a[] := a[] - b[]; then return a[] */
double *VSUB (int n, double a[], double b[]);
#define VSuB(n,a,b,i) for ((i)=(n); (i)--;) (a)[i] -= (b)[i];

/* dot product of a[] and b[] */
double Vdot (int n, double a[], double b[]);

/* c[] := a[] + multiplier * b[]; then return c[] */
double *Vaddmul (int n, double a[], double multiplier, double b[], double c[]);

/* c[] := a[] + b[] / divisor; then return c[] */
double *Vadddiv (int n, double a[], double b[], double divisor, double c[]);

/* b[] := b[] + multiplier * a[]; then return b[] */
double *VADDmul (int n, double multiplier, double a[], double b[]);

/* c[] := aa * a[] + bb * b[]; then return c[] */
double *Vaddmulmul
(int n, double aa, double a[], double bb, double b[], double c[]);

/* b[] := b[] + a[] / divisor; then return b[] */
double *VADDdiv (int n, double a[], double divisor, double b[]);

/* c[] := a[] - multiplier * b[]; then return c[] */
double *Vsubmul (int n, double a[], double multiplier, double b[], double c[]);

/* c[] := a[] - b[] / divisor; then return c[] */
double *Vsubdiv (int n, double a[], double b[], double divisor, double c[]);

/* a[] := a[] - multiplier * b[]; then return a[] */
double *VSUBmul (int n, double a[], double multiplier, double b[]);

/* a[] := a[] - b[] / divisor; then return a[] */
double *VSUBdiv (int n, double a[], double b[], double divisor);


/* Matrix.c: */

#define KroneckerDelta(i,j) (((i)==(j))?1:0)

/* LapackDrivers.c: */

int FORTRAN_SYMBOL(dgetrf)
    (int *m, int *n, double *a, int *lda, int *ipiv, int *info);

int FORTRAN_SYMBOL(dgetri)
    (int *n, double *a, int *lda, int *ipiv, double *work,
     int *lwork, int *info);

int FORTRAN_SYMBOL(zhpev)
    (char *jobz, char *uplo, int *n, double 
     *ap, double *w, double *z__, int *ldz, double *
     work, double *rwork, int *info);
 
/***************************************************************/
/* BLAS and LAPACK library drivers with simplified interface   */
/* and temporary memory hierarchy (static..dynamic) management */
/* invisible to the end-user.                                  */
/***************************************************************/

/* B := inv(A), both A and B are real matrices of rank x rank */
void Minv (int rank, double A[], double B[]);

/*******************************************************************/
/* Diagonalize a complex Hermitian matrix described by ap[], order */
/* its (real) eigenvalues in ascending order, and compute all      */
/* (complex) eigenvectors associated with them.                    */
/*                                                                 */
/* ap=double[rank*(rank+1)]: upper triangle of A packed columnwise */
/* ap[j*(j+1)+2*i] = Re(A_{ij}), ap[j*(j+1)+2*i+1] = Img(A_{ij})   */
/* i, j = 0 .. rank-1, i <= j.                                     */
/*                                                                 */
/* eigenvalues = double[rank], eigenvectors = double[2*rank*rank]. */
/*******************************************************************/
void MDiag(int rank, double ap[], double eigenvalues[], double eigenvectors[]);

#endif  /* _VecMat_h */
