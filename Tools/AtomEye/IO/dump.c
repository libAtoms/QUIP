/****************************************/
/* libIO:                               */
/*                                      */
/* Strings, parsers & files             */
/*                                      */
/* Dec.12, 1999  Ju Li <liju99@mit.edu> */
/****************************************/

#include "IO.h"

/******************************************************************/
/* Dump indefinite amount of data from a single object to stream, */
/* formatted as line x minline\maxline. One can add labeling info */
/* at top and bottom if separated by '\n' in *line. Use %n in     */
/* *line to denote line index, %-2n to denote line index-2; %i to */
/* denote token index, %2i to denote token index+2. Finally only  */
/* buffer[minchar\maxchar] will be printed to the *fp stream.     */
/******************************************************************/

#define printtoken(TYPE) do { if ( (i>=minline) && unprinted ) \
{ snprintf(r, &buf[maxchar]-r+1, inert, *((TYPE *)B)); \
 while (*r != EOS) r++; if (r==&buf[maxchar]) \
 { fprintf(fp, "%s\n", buf+minchar); unprinted=FALSE; } }; \
 B += sizeof(TYPE); tokens++; } while(0)

void fDump (FILE *fp, char *linefmt, int minline, int maxline, void *A, 
	    int minchar, int maxchar)
{
    int i, j, tokens, line_bytes=0, line_tokens=0;
    char *obj_start, *obj_end, *inert, *start, *end=NULL, *r;
    char *format, *buf, *s, islong, unprinted, c, d, e, f;
    unsigned char *B;

    if (!fp) return;
    /* perhaps calling from some quirky functions */
    if (maxline < minline) return;
    if (maxchar < minchar) return;
    
    for (obj_start=linefmt; *obj_start!=EOS; obj_start++)
	if ( *obj_start == '%' )
        {
	    if ( *(obj_start+1) != '%' ) break;
	    else obj_start++;
        }
    if ( *obj_start == EOS )
    { /* there is no linefmt operator */
	fprintf (fp, "%s\n", linefmt);
	return;
    }
    /* backtrack to first '\n' */
    for (; (obj_start>=linefmt)&&(*obj_start!='\n'); obj_start--);
    if (obj_start>=linefmt)
    {
	obj_start++;
	buf = IOmem(obj_start-linefmt+1);
	for (r=linefmt; r<obj_start; r++) buf[r-linefmt] = *r;
	buf[r-linefmt] = EOS;
	fprintf (fp, "%s", buf);
	IOfree (buf);
    }
    else obj_start++;
    
    for (obj_end=eos(linefmt)-1; obj_end>linefmt; obj_end--)
	if ( *obj_end == '%' )
        {
	    if ( *(obj_end-1) != '%' ) break;
	    else obj_end--;
        }
    for (; (obj_end<eos(linefmt)) && (*obj_end!='\n'); obj_end++);
    
    format = IOmem(obj_end-obj_start+1);
    for (r=obj_start; r<obj_end; r++) format[r-obj_start] = *r;
    format[r-obj_start] = EOS;
    
    buf = IOmem(maxchar+1);
    buf[maxchar] = EOS;
    tokens = 0;
    B = (unsigned char *)A;
    
    for (i=0; i<maxline; i++)
    {
	if (i==1)
	{
	    line_bytes = B - (unsigned char *)A;
	    line_tokens = tokens;
	    if (minline > 1)
	    {
		B += (minline-1) * line_bytes;
		tokens += (minline-1) * line_tokens;
		i += minline-2;
		continue;
	    }
	}
	r = buf;
	buf[minchar] = EOS;
	unprinted = TRUE;
	for ( inert=format; ; inert = end+1 )
	{
	    for ( start=inert; *start!=EOS; start++ )
	    {
		if ( *start == '%' )
		{
		    if ( *(start+1) == '%' )
		    { 
			start++;
			continue;
		    }
		    else
		    {
			for ( islong=FALSE, end=start+1;
			      ISDIGIT(*end)
				  || (*end=='.')
				  || (*end=='#')
				  || (*end=='+')
				  || (*end=='-')
				  || (*end=='l');
			      end++ )
			    if (*end == 'l') islong = TRUE;
			if ( !ISALPHA(*end) )
			{
			    printf ("error: fDump: \"%s\" "
				    "contains single %%\n", linefmt);
			    exit(1);
			}
			c = *(end+1);
			*(end+1) = EOS;
			switch ( *end )
			{
			case 'o':
			case 'u':
			case 'x':
			case 'X':
			    f = *end;
			    *end = EOS;
			    if (end == start+1) j = 0;
			    else j = atoi(start+1);
			    *end = f;
			    if (j == 2)
			    {
                                s = (char *) B;
                                j = *B;
                                B = (unsigned char *)(&j);
                                printtoken (int);
                                B = (unsigned char *)s + sizeof(char);
			    }
			    else if ( islong ) printtoken(unsigned int);
			    else printtoken(unsigned short);
			    break;
			case 'd':
			case 'i':
			  if ( islong ) printtoken(long int);
			    else printtoken(int);
			    break;
			case 'f':
			case 'g':
			case 'G':
			case 'e':
			case 'E':
			    if ( islong ) printtoken(double);
			    else printtoken(float);
			    break;
			case 'c':
			    printtoken(char);
			    break;
			case 's':
			    printtoken(char *);
			    break;
			case DUMP_LINENO:
			case DUMP_TOKENO:
			    f = *end;
			    *end = EOS;
			    if (end == start+1) j = 0;
			    else j = atoi(start+1);
			    *end = f;
			    d = *(start+1);
			    *(start+1) = 'd';
			    e = *(start+2);
			    *(start+2) = EOS;
			    if ( (i >= minline) && unprinted )
			    {
				snprintf(r, &buf[maxchar]-r+1, inert,
					 ((f==DUMP_LINENO)?i:tokens)+j);
				while (*r != EOS) r++;
				if (r == &buf[maxchar])
				{
				    fprintf(fp, "%s\n", buf+minchar);
				    unprinted = FALSE;
				}
			    }
			    *(start+2) = e;
			    *(start+1) = d;
			    break;
			default:
			    printf ("error: fDump: \"%s\" descriptor in "
				    "\"%s\" unsupported\n", start, linefmt);
			    exit(1);
			}
			*(end+1) = c;
			break;
		    }
		}
	    }  /* find escape sequence start and treat it */
	    if (*start == EOS) break;  /* run out of commands */
	    if ( (i>1) && (!unprinted) )
	    {
		B = ((unsigned char *)A) + (i+1) * line_bytes;
		tokens = (i+1) * line_tokens;
		break;
	    }
	}  /* loop over inert beginnings */
	if ( (i >= minline) && unprinted )
	{
	    snprintf(r, &buf[maxchar]-r+1, inert);
	    fprintf(fp, "%s\n", buf+minchar);
	}
    }  /* loop over line number i */
    IOfree (buf);
    IOfree (format);
    if (obj_end+1<eos(linefmt)) fprintf (fp,"%s",obj_end+1);
    return;
} /* end fDump() */


#ifdef _dump_TEST
#define MAXLINE 10
#define MAXBYTE 256
void main()
{
    int i;
    int x[MAXLINE*3];
    /* char x[MAXBYTE]; */
    /* double x[MAXLINE*3]; */
    for (i=0; i<MAXLINE*3; i++) x[i] = 0xffffffff;
    /* FILE *file = ropen("dump.c"); */
    /* fread(x,1,MAXBYTE,file); */
    /* Fclose(file); */
    dumph(MAXBYTE,x);
    /* dump3(3,x); */
    /* dumpI3(MAXLINE,x); */
}
#endif


/* Very crude printout of a double precision real matrix */
void fMump (FILE *fp, int rows, int cols, double *A)
{
    int i, j;
    static char sep[TERMCHAR+1]="---------------------------------------"
	"----------------------";
    if (!fp) return;
    fprintf(fp, "this matrix is %d x %d (Row by Column)",
	    rows, cols);
    fcr(fp);
    fprintf(fp, sep); 
    for (i=0; i<rows; i++)
	for (j=0; j<cols; j++)
	{
	    if ( j % FMUMP_COLUMNS == 0 )
	    {
		fcr(fp);
		if ( (j == 0) && ((rows > 3)||(cols > 2*FMUMP_COLUMNS)) )
		    fprintf(fp, "<R%d>:\t", i);
		if (j != 0)
		    fprintf(fp, "  [C%d]\t", j);
	    }
	    fprintf(fp, FMUMP_SPACING FMUMP_FORMAT, A[i*cols+j]);
	}
    fcr(fp);
    fprintf(fp, sep);
    fcr(fp);
    return;
} /* end fMump() */

#ifdef _Mump_TEST
#define ROWS 10
#define COLS 100
void main()
{
    int i, j;
    double A[ROWS*COLS];
    for (i=0,j=1; i<ROWS*COLS; i++,j=-j) A[i] = i*j;
    Mump (ROWS, COLS, A);
}
#endif
