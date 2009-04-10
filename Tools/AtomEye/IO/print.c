/****************************************/
/* libIO:                               */
/*                                      */
/* Strings, parsers & files             */
/*                                      */
/* Dec.12, 1999  Ju Li <liju99@mit.edu> */
/****************************************/

#include "IO.h"

/********************************************************************/
/* Line/matrix multi-object formatter using printf escape sequences */
/* "%[0-9.l]*[a-zA-Z]" with alignment directives. Works for finite  */
/* & pretty-printable multi-object data. For example, object width  */
/* is limited by terminal width, height is limited by 9, etc.       */
/* For large or indefinite amount of data, use stream single-object */
/* fDump() that has brute-force screening and cutoffs.              */
/********************************************************************/

/* escape sequence type */
#define ESC_ALIGN_TOP  -1
#define ESC_ALIGN_MID  -2
#define ESC_ALIGN_BOT  -3
#define ESC_ALIGN_USR   0
#define ESC_MATRIX      1
#define ESC_PRINTF      2   /* all other printf sequences */

/* escape sequence value */
#define ALIGN_TOP  -1 /* alignment directives */
#define ALIGN_MID  -2
#define ALIGN_BOT  -3

#define PRECISION_SHORT  0
#define PRECISION_LONG   1

static int nobj;
static struct Obj
{
    int width;
    int rows;
    int alignment;
    /* character rendering device */
    char crd[MPRINT_MAX_ROWS][MPRINT_MAX_COLS+1];
    int start_row;
} o[MPRINT_MAX_OBJS];

/* parsed escape sequence from stream: redundancy is welcome */
struct Esc
{
    char *start; /* start of the escape sequence in format[] */
    char *end;   /* end of the escape sequence in format[] */
    int type;  /* nature of the escape sequence: ESC_something */
    int value; /* value for that type: for example, ALIGN_something */
    char *eff;  /* effective content of the sequence */
    char *inert;  /* inert content before start */
};


/******************************************************************/
/* Locate the first occurrence of %[0-9.l]*[a-zA-Z] in format[];  */
/* if find one, assign start,end,type,inert[], and possibly value */
/* and eff[]; else assign type,inert[]; return end+1.             */
/******************************************************************/
static char *parse_printf_escape (char *format, struct Esc *s)
{
    int before_escape;
    char *r, c;
    for ( s->start=format; *(s->start)!=EOS; ++s->start )
	if (*s->start == '%')
	{
	    if (*(s->start+1) == '%')
	    {  /* %% stands for % */
		s->start++;
		continue;
	    }
	    else
	    {
		for ( s->end=s->start+1;
		      ISDIGIT(*s->end)
			  || (*s->end=='.')
			  || (*s->end=='#')
			  || (*s->end=='+')
			  || (*s->end=='-')
			  || (*s->end=='l');
		      s->end++ );
		if ( !ISALPHA(*s->end) )
		{
		    printf ("error: parse_printf_escape: \"%s\" "
			    "contains single %%\n", format);
		    exit(1);
		}
		switch ( *s->end )
		{
		case ALIGN_TOP_ESC_CHAR:
		    s->type = ESC_ALIGN_TOP; /* use */
		    s->value = ALIGN_TOP;    /* use */
		    goto exit;
		case ALIGN_MID_ESC_CHAR:
		    s->type = ESC_ALIGN_MID; /* use */
		    s->value = ALIGN_MID;    /* use */
		    goto exit;
		case ALIGN_BOT_ESC_CHAR:
		    s->type = ESC_ALIGN_BOT; /* use */
		    s->value = ALIGN_BOT;    /* use */
		    goto exit;
		case ALIGN_USR_ESC_CHAR:
		    s->type = ESC_ALIGN_USR; /* use */
		    for (r=s->start; r<s->end; r++)
			s->eff[r-s->start] = *r;
		    s->eff[r-s->start] = EOS;
		    s->value = atoi(s->eff+1); /* use */
		    goto exit;
		case MATRIX_ESC_CHAR:
		    s->type = ESC_MATRIX; /* use */
		    for (r=s->start; r<s->end; r++)
			s->eff[r-s->start] = *r;
		    s->eff[r-s->start] = EOS;
		    s->value = atoi(s->eff+1); /* use */
		    c = *(s->end+1); /* sequence termination symbol */
		    for ( before_escape=TRUE, r=s->end+2;
			  (before_escape||(*r!=c)) && (*r!=EOS); r++)
		    {
			s->eff[r-s->end-2] = *r;
			/* search for termination symbol after % */
			if  (*r == '%')
			    before_escape = FALSE;
		    }
		    if ( *r == EOS )
		    {
			printf ("error: parse_printf_escape: unclosed "
				"matrix formatter\n\"%s\".\n", s->start);
			exit(1);
		    }
		    s->eff[r-s->end-2] = c;
		    s->eff[r-s->end-1] = EOS; /* use */
		    s->end = r;               /* use */
		    goto exit;
		default:
		    s->type = ESC_PRINTF;       /* use */
		    s->value = PRECISION_SHORT; /* use */
		    for (r=s->start; r<=s->end; r++)
		    {
			s->eff[r-s->start] = *r;
			if (*r == 'l') s->value = PRECISION_LONG;
		    }
		    s->eff[r-s->start] = EOS;   /* use */
		    goto exit;
		}
	    }
	}
    if (*(s->start) == EOS)
    {
	s->type = ESC_PRINTF;
	s->end = s->start - 1;
	s->eff[0] = EOS;
    }
exit:
    for (r=format; r<s->start; r++)
	s->inert[r-format] = *r;
    s->inert[r-format] = EOS;  /* use */
    return (s->end + 1);
} /* end parse_printf_escape() */


#define inc_nobj() if (++nobj==MPRINT_MAX_OBJS) break;
#define line_obj(o)  (o.rows==1)
#define virgin_obj(o)  (line_obj(o)&&(o.width==0))
#define make_virgin_obj(o) \
(o.rows=1,o.width=0,o.crd[0][0]=EOS,o.alignment=alignment)
#define sticky_line(o) (virgin_obj(o)||(line_obj(o)&&(alignment==o.alignment)))
#define islong (s.value == PRECISION_LONG)
#define printtoken(TYPE) \
snprintf(r, &(o[nobj].crd[0][MPRINT_MAX_COLS])-r+1, s.inert, va_arg(ap, TYPE))

void Mvprintf (FILE *fp, char *fmt, va_list ap)
{
    struct Esc s;
    void *matrix;
    int alignment = ALIGN_MID;
    char *format, *fmt_start, *fmt_end, *p, *q, *r, f;
    int i, j, k, maxrow;

    if (!fp) return;
    format  = IOclone(fmt);
    s.inert = IOclone(fmt);
    s.eff   = IOclone(fmt);
    
    for (fmt_start=format; *fmt_start!=EOS; fmt_start++)
	if ( *fmt_start == '%' )
        {
	    if ( *(fmt_start+1) != '%' ) break;
	    else fmt_start++;
        }
    if ( *fmt_start == EOS )
    { 
	fprintf (fp, format);
	fcr(fp);
	return;
    }
    /* backtrack to first '\n' */
    for (; (fmt_start>=format)&&(*fmt_start!='\n'); fmt_start--);
    if (fmt_start>=format)
    {
	*fmt_start = EOS;
	fprintf (fp, format);
	fcr(fp);
	fmt_start++;
    }
    else fmt_start++;
    
    for (fmt_end=eos(fmt_start)-1; fmt_end>fmt_start; fmt_end--)
	if ( *fmt_end == '%' )
        {
	    if ( *(fmt_end-1) != '%' ) break;
	    else fmt_end--;
        }
    for (; (fmt_end<eos(fmt_start)) && (*fmt_end!='\n'); fmt_end++);
    if (*fmt_end == '\n') *(fmt_end++) = EOS;
    
    nobj = 0;
    make_virgin_obj(o[nobj]);
    p = fmt_start;
    while (TRUE)
    {
	q = parse_printf_escape (p, &s);
	switch ( s.type )
	{
            case ESC_ALIGN_TOP:
            case ESC_ALIGN_MID:
            case ESC_ALIGN_BOT:
            case ESC_ALIGN_USR:
                if (*s.inert != EOS)
                { /* there are some inert content before align directive */
                    if  (!sticky_line(o[nobj]))
                    {
                        inc_nobj();
                        make_virgin_obj(o[nobj]);
                    }
                    o[nobj].alignment=alignment;
                    for (r=&(o[nobj].crd[0][0]); *r!=EOS; r++);
                    snprintf(r, &(o[nobj].crd[0][MPRINT_MAX_COLS])-r+1,
                             s.inert);
                    for (; *r!=EOS; r++);
                    o[nobj].width = r - &(o[nobj].crd[0][0]);
                }
                alignment = s.value;
                break;
            case ESC_MATRIX:
                if (*s.inert != EOS)
                { /* there are some inert content before %M */
                    if  (!sticky_line(o[nobj]))
                    {
                        inc_nobj();
                        make_virgin_obj(o[nobj]);
                    }
                    o[nobj].alignment=alignment;
                    for (r=&(o[nobj].crd[0][0]); *r!=EOS; r++);
                    snprintf(r, &(o[nobj].crd[0][MPRINT_MAX_COLS])-r+1,
                             s.inert);
                    for (; *r!=EOS; r++);
                    o[nobj].width = r - &(o[nobj].crd[0][0]);
                }
                if (!virgin_obj(o[nobj])) inc_nobj();
                matrix = va_arg(ap, void *);
                o[nobj].alignment = alignment;
                o[nobj].rows = s.value;
                if (o[nobj].rows > MPRINT_MAX_ROWS)
                    o[nobj].rows = MPRINT_MAX_ROWS;
                o[nobj].width = Asprintf(&(o[nobj].crd[0][0]), MPRINT_MAX_COLS,
                                         o[nobj].rows, s.eff, matrix);
                inc_nobj();
                make_virgin_obj(o[nobj]);
                break;
            case ESC_PRINTF:
                /* combine inert and active parts */
                strcat(s.inert, s.eff);
                if  (!sticky_line(o[nobj]))
                { /* only when necessary do we start a new one */
                    inc_nobj();
                    make_virgin_obj(o[nobj]);
                }
                o[nobj].alignment = alignment;
                for (r=&(o[nobj].crd[0][0]); *r!=EOS; r++);
                switch ( *s.end )
                {
                    case 'o':
                    case 'u':
                    case 'x':
                    case 'X':
                        f = *s.end;
                        *s.end = EOS;
                        if (s.end == s.start+1) j = 0;
                        else j = atoi(s.start+1);
                        *s.end = f;
                        if (j == 2)
                            snprintf(r, &(o[nobj].crd[0][MPRINT_MAX_COLS])-r+1,
                                     s.inert, (int)(va_arg(ap, int)));
                        else if ( islong ) printtoken(unsigned int);
                        else printtoken(unsigned int);
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
                        else printtoken(double);
                        break;
                    case 'c':
                        printtoken(int);
                        break;
                    case 's':
                        printtoken(char *);
                        break;
                    default:
                        snprintf(r, &(o[nobj].crd[0][MPRINT_MAX_COLS])-r+1,
                                 s.inert);
                }
                for (; *r!=EOS; r++);
                o[nobj].width = r - &(o[nobj].crd[0][0]);
                break;
            default:
                printf ("error: Mvprintf: unrecognized "
                        "escape sequence type number %d\n", s.type);
                exit(1);
	}
	if (*q == EOS)
	{
	    if (!virgin_obj(o[nobj])) nobj++;
	    break;
	}
	p = q;
    }
    /** render line and matrix objects according to their alignments **/
    maxrow = 0;
    /* first pass: calculate maximum (total) number of rows */
    for (i=0; i<nobj; i++)
	switch ( o[i].alignment )
	{
            case ALIGN_TOP:
            case ALIGN_MID:
            case ALIGN_BOT:
                if (o[i].rows > maxrow) maxrow = o[i].rows;
                continue;
            default: /* value in ESC_ALIGN_USR */
                if (o[i].rows + o[i].alignment > maxrow)
                    maxrow = o[i].rows + o[i].alignment;
                continue;
	}
    /* second pass: calculate start row of each object */
    for (i=0; i<nobj; i++)
	switch ( o[i].alignment )
	{
            case ALIGN_TOP:
                o[i].start_row = 0;
                continue;
            case ALIGN_MID:
                o[i].start_row = (maxrow - o[i].rows) >> 1;
                continue;
            case ALIGN_BOT:
                o[i].start_row = maxrow - o[i].rows;
                continue;
            default: 
                o[i].start_row = o[i].alignment;
                continue;
	}
    /* last pass: render each object */
    for (j=0; j<maxrow; j++)
    {
	for (i=0; i<nobj; i++)
	    for (k=0; k<o[i].width; k++)
		if ( ( j >= o[i].start_row ) &&
		     ( j < o[i].start_row + o[i].rows ) )
		    fputc ((int)o[i].crd[j-o[i].start_row][k], fp);
		else fputc (' ', fp);
	fcr(fp);
    }
    if (*fmt_end != EOS)
    {
	fprintf (fp, fmt_end);
	fcr(fp);
    }
    IOfree(s.eff);
    IOfree(s.inert);
    IOfree(format);
    return;
} /* end Mvprintf() */


void Mfprintf (FILE *fp, char *fmt, ... )
{
    va_list ap;
    va_start(ap, fmt);
    Mvprintf (fp, fmt, ap);
    return;
} /* end Mfprintf() */


void Mprintf (char *fmt, ... )
{
    va_list ap;
    va_start(ap, fmt);
    Mvprintf (stdout, fmt, ap);
    return;
} /* end Mprintf() */

#ifdef _Mprintf_TEST
void main()
{
    double A[3][3] = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    Mprintf ("\n\tSo ...\n\n%t%s%m = %3M|| %.1lf %.1lf %.1lf | %%%q%s%m"
	     "%1M><%.0lf %.0lf %.0lf>\n\n\t     That's sick!\n",
	     "Buddha", A, "Buddha", &A[1][0]);
    /* Mprintf ("\n%3M|| %.2lf %.2lf %.2lf |\n ", A); */
}
#endif


/* Sole double matrix Mfprintf() driver: catches the lone "%M" in */
/* sole_fmt and automatically justifies elements column by column */

void Sprint (FILE *fp, char *sole_fmt, double *A, int rows, int cols,
	     char open_char, char close_char, int open_space,
	     int close_space, int space_between)
{
    static char fmt[SPRINT_MAX_FMT];
    static char buf[SPRINT_MAX_BUF];
    static int pre_dot[SPRINT_MAX_ELM];
    static int aft_dot[SPRINT_MAX_ELM];
    int i, j, k;
    char *p, *q, *r;
    
    if (!fp) return;
    if (cols > SPRINT_MAX_ELM)
    {
	printf ("error: Sprint: matrix cols = %d > SPRINT_MAX_ELM = %d\n",
		cols, SPRINT_MAX_ELM);
	exit(1);
    }

    /* ..%2M][ %16.12lf  %16.12lf ].. */
    if ( strlen(sole_fmt) + 4 +
	 (3+open_space+close_space+(cols-1)*space_between+cols*8)
	 >= SPRINT_MAX_FMT )
    {
	printf ("error: Sprint: matrix cols = %d makes SPRINT_MAX_FMT\n"
		"= %d overflow possible\n", cols, SPRINT_MAX_FMT);
	exit(1);
    }
    
    p = strstr (sole_fmt, "%M");
    if ( p == NULL )
    {
	printf ("error: Sprint: no %%M escape sequence in\n"
		"\"%s\"\n", sole_fmt);
	exit(1);
    }

    for (r=sole_fmt; r<p; r++) fmt[r-sole_fmt] = *r;
    *(r=&fmt[r-sole_fmt]) = EOS;

    for (j=0; j<cols; j++)
    {
	pre_dot[j] = 0;
	aft_dot[j] = 0;
	for (i=0; i<rows; i++)
	{
	    sprintf (buf, "%.6f", A[i*cols+j]);
	    /* cosmetics */
	    for (k=FALSE,q=buf; *q!=EOS; q++)
		if (*q=='.') k=TRUE;
	    if (k)
		for (q--;
		     ((*q==' ')||(*q=='0')||(*q=='.'));
		     q--)
		{
		    if (*q=='.')
		    {
			*q = EOS;
			break;
		    }
		    else *q = EOS;
		}
	    if ( (q=strchr(buf, '.')) != NULL )
	    {
		if (q-buf > pre_dot[j]) pre_dot[j] = q-buf;
		if ((k=eos(buf)-q-1) > aft_dot[j]) aft_dot[j] = k;
	    }
	    else if ((k=strlen(buf)) > pre_dot[j]) pre_dot[j] = k;
	}
    }
    
    sprintf (r, "%%%dM%c%c", rows, close_char, open_char);
    while (*r != EOS) r++;
    for (i=0; i<open_space; i++) *(r++) = ' ';
    for (j=0; j<cols; j++)
    {
	r += sprintf (r, "%%%d.%dlf",
		      pre_dot[j]+aft_dot[j]+((aft_dot[j]>0)?1:0),
		      aft_dot[j]);
	if (j != cols-1)
	    for (i=0; i<space_between; i++) *(r++) = ' ';
    }
    for (i=0; i<close_space; i++) *(r++) = ' ';
    *(r++) = close_char;
    sprintf (r, p+2);
    Mfprintf (fp, fmt, A);
    return;
} /* end Sprint() */

#ifdef _Sprint_TEST
#define ROWS 3
void main()
{
    double A[ROWS][ROWS] = {{0.2, 2, 3}, {4, -25.789, 6}, {-7, 8, 19}};
    S3pr("\nA = %M\n ", A[0]);
    V3pr("A(1,:) = %M parsec.\n ", A[1]);
}
#endif


/* Driver of Sprint() to print A[][]*factor. */
void Sprintmul (FILE *fp, char *sole_fmt, double *A, double factor,
                int rows, int cols, char open_char, char close_char,
                int open_space, int close_space, int space_between)
{
    int i;
    double *B;
    MALLOC ( Sprintmul, B, rows*cols, double );
    for (i=0; i<rows*cols; i++) B[i] = A[i] * factor;
    Sprint (fp, sole_fmt, B, rows, cols, open_char, close_char,
            open_space, close_space, space_between);
    free (B);
    return;
} /* end Sprintmul() */


/* Sole int matrix Mfprintf() driver: catches the lone "%M" in */
/* sole_fmt and automatically justifies elements column by column */

void Iprint (FILE *fp, char *sole_fmt, int *A, int rows, int cols,
	     char open_char, char close_char, int open_space,
	     int close_space, int space_between)
{
    static char fmt[IPRINT_MAX_FMT];
    static char buf[IPRINT_MAX_BUF];
    static int pre_dot[IPRINT_MAX_ELM];
    static int aft_dot[IPRINT_MAX_ELM];
    int i, j, k;
    char *p, *q, *r;
    
    if (!fp) return;
    if (cols > IPRINT_MAX_ELM)
    {
	printf ("error: Iprint: matrix cols = %d > IPRINT_MAX_ELM = %d\n",
		cols, IPRINT_MAX_ELM);
	exit(1);
    }

    /* ..%2M][ %16.12lf  %16.12lf ].. */
    if ( strlen(sole_fmt) + 4 +
	 (3+open_space+close_space+(cols-1)*space_between+cols*8)
	 >= IPRINT_MAX_FMT )
    {
	printf ("error: Iprint: matrix cols = %d makes IPRINT_MAX_FMT\n"
		"= %d overflow possible\n", cols, IPRINT_MAX_FMT);
	exit(1);
    }
    
    p = strstr (sole_fmt, "%M");
    if ( p == NULL )
    {
	printf ("error: Iprint: no %%M escape sequence in\n"
		"\"%s\"\n", sole_fmt);
	exit(1);
    }

    for (r=sole_fmt; r<p; r++) fmt[r-sole_fmt] = *r;
    *(r=&fmt[r-sole_fmt]) = EOS;

    for (j=0; j<cols; j++)
    {
	pre_dot[j] = 0;
	aft_dot[j] = 0;
	for (i=0; i<rows; i++)
	{
	    sprintf (buf, "%d", A[i*cols+j]);
	    /* cosmetics */
	    for (k=FALSE,q=buf; *q!=EOS; q++)
		if (*q=='.') k=TRUE;
	    if (k)
		for (q--;
		     ((*q==' ')||(*q=='.'));
		     q--)
		{
		    if (*q=='.')
		    {
			*q = EOS;
			break;
		    }
		    else *q = EOS;
		}
	    if ( (q=strchr(buf, '.')) != NULL )
	    {
		if (q-buf > pre_dot[j]) pre_dot[j] = q-buf;
		if ((k=eos(buf)-q-1) > aft_dot[j]) aft_dot[j] = k;
	    }
	    else if ((k=strlen(buf)) > pre_dot[j]) pre_dot[j] = k;
	}
    }
    
    sprintf (r, "%%%dM%c%c", rows, close_char, open_char);
    while (*r != EOS) r++;
    for (i=0; i<open_space; i++) *(r++) = ' ';
    for (j=0; j<cols; j++)
    {
	r += sprintf (r, "%%%dd",
                      pre_dot[j]+aft_dot[j]+((aft_dot[j]>0)?1:0));
	if (j != cols-1)
	    for (i=0; i<space_between; i++) *(r++) = ' ';
    }
    for (i=0; i<close_space; i++) *(r++) = ' ';
    *(r++) = close_char;
    sprintf (r, p+2);
    Mfprintf (fp, fmt, A);
    return;
} /* end Iprint() */

#ifdef _Iprint_TEST
#define ROWS 3
void main()
{
    int A[ROWS][ROWS] = {{0, 2, 3}, {4, 25, 6}, {-7189, 8, 19}};
    Ipr("\nA = %M\n ", A[0], 3, 3);
    ipr("A(1,:) = %M parsec.\n ", A[1], 3);
}
#endif
