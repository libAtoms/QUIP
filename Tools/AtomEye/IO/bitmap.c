/***************************************/
/* libIO:                              */
/*                                     */
/* Strings, parsers & files            */
/*                                     */
/* Dec.12, 1999 Ju Li <liju99@mit.edu> */
/***************************************/

#include "IO.h"

/**************************************************/
/* Bitmap object: as if b[], which is char array, */
/* is a bit array indexed by i and has n members. */
/**************************************************/

/* return pointer to a bitmap of n bits (rounded) */
Bmap *Balloc (int n)
{
    void *ptr;
    ptr = malloc((bunit(n-1)+1)*sizeof(Bmap));
    if (ptr == NULL)
    {
	printf ("error: Balloc: attempt to malloc %d\n"
		"bitmap units, or %ld bytes, failed", bunit(n-1)+1,
		(long)((bunit(n-1)+1)*sizeof(Bmap)) );
	exit(1);
    }
    return((Bmap *)ptr);
} /* end Balloc() */


/* return pointer to a bitmap of n bits (rounded) which are all set to 0 */
Bmap *BALLOC (int n)
{
    void *ptr;
    ptr = calloc(bunit(n-1)+1,sizeof(Bmap));
    if (ptr == NULL)
    {
	printf ("error: BALLOC: attempt to calloc %d\n"
		"bitmap units, or %ld bytes, failed", bunit(n-1)+1,
		(long)((bunit(n-1)+1)*sizeof(Bmap)) );
	exit(1);
    }
    return((Bmap *)ptr);
} /* end BALLOC() */


/* return pointer to a bitmap of n bits (rounded) which are */
/* all set to 0 or 1, depending on whether val is 0 or not. */
Bmap *BAlloc (int n, int val)
{
    void *ptr;
    ptr = malloc((bunit(n-1)+1)*sizeof(Bmap));
    if (ptr == NULL)
    {
	printf ("error: BAlloc: attempt to malloc %d bitmap\n"
		"units, or %ld bytes, failed", bunit(n-1)+1,
		(long)((bunit(n-1)+1)*sizeof(Bmap)) );
	exit(1);
    }
    BASSIGNALL(ptr,n,val);
    return((Bmap *)ptr);
} /* end BAlloc() */


/* "..." means set of (Bmap *b, int offset, int num_bits), escape seq. is */
/* "%B". Print out num_bits from Bitmap array b, starting from bit offset */
void Bfprintf (FILE *out, char *fmt, ...)
{
    va_list ap;
    va_start (ap, fmt);
    Bvfprintf (out, fmt, ap);
    va_end (ap);
    return;
} /* end Bfprintf() */

void Bprintf (char *fmt, ...)
{
    va_list ap;
    va_start (ap, fmt);
    Bvfprintf (stdout, fmt, ap);
    va_end (ap);
    return;
} /* end Bprintf() */

void Bvfprintf (FILE *out, char *fmt, va_list ap)
{
    register int i, j, offset, num_bits;
    char *start;
    Bmap *b;
    
    for ( start=fmt; *start!=EOS; ++start )
    {
        if (*start == '%')
	{
	    if (*(start+1) == '%')
	    {  /* %% stands for % */
                fputc ('%', out);
		start++;
		continue;
	    }
            else if (*(start+1) == 'B')
            {
                b = va_arg (ap, Bmap *);
                offset = va_arg (ap, int);
                num_bits = va_arg (ap, int);
                j = (num_bits>0)? 1:-1;
                for (i=0; i!=num_bits; i+=j)
                    fprintf(out, "%1d", ZERO_ONE(BVAL(b,offset+i)));
                start++;
                continue;
            }
            pe ("Bvfprintf: \"%s\" \nhas escape sequences "
                "other than %%B\n", fmt);
        }
        else fputc (*start, out);
    }
    return;
} /* end Bvfprintf() */

/* Bfwrite: reverse bitstream direction to emulate human handwriting */
void Bfwrite (FILE *out, char *c, int length_in_bytes)
{
    register int i;
    for (i=length_in_bytes-1; i>=0; i--)
    {
        Bfprintf(out,"%B",(Bmap *)(c+i),7,-8);
        if (i != 0) fsp(out);
    }
    return;
} /* end Bfwrite() */

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

#ifdef _Bwrite_TEST
int main (int argc, char *argv[])
{
    char C = 0x8c;
    unsigned short S = 0xfe12;
    int I = 0xfe1234;
    long L = ((long)I) << 6;
    double D = 1.0;
    printf ("long int = %d bytes\n", sizeof(long int));
    qt(); BwriteChar(C); qt(); cr();
    qt(); BwriteShort(S); qt(); cr();
    qt(); BwriteInt(I); qt(); cr();
    qt(); BwriteLong(L); qt(); cr();
    qt(); BwriteDouble(D); qt(); cr();
    return(0);
}
#endif  /* _Bwrite_TEST */
