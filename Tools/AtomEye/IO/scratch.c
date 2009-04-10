/****************************************/
/* libIO:                               */
/*                                      */
/* Strings, parsers & files             */
/*                                      */
/* Dec.12, 1999  Ju Li <liju99@mit.edu> */
/****************************************/

#include "IO.h"

/***********************************************************/
/* indefinite temporary scratch space: provides fp_scratch */
/***********************************************************/

FILE *fp_scratch = NULL;
/* stream fp_scratch is fully buffered: unless fflush()ed or */
char scratch_buffer [SCRATCH_BUFFER_SIZE];
/* SCRATCH_BUFFER_SIZE reached, all fprintf()ed stuff are still there */


/* At first call, set up fp_scratch which is fully buffered */
/* at scratch_buffer with SCRATCH_BUFFER_SIZE chars. In     */
/* subsequent calls, previous data are discarded and the    */
/* stream is back to SEEK_SET. When kill_scratch() called   */
/* or process terminates, there would be nothing on disk.   */
void reset_scratch()
{
    char *path;
    if (fp_scratch == NULL)
    {
	fp_scratch = fopen (path=TMPNAM(NULL), "w+");
	if (fp_scratch == NULL)
	    pe ("reset_scratch: cannot open temporary file\n");
	/* the file is deleted at the termination of the */
	if (unlink(path) < 0)
	    pe ("reset_scratch: cannot unlink temporary file\n");
	/* process because its link count would be 0.    */
	setvbuf (fp_scratch, scratch_buffer, _IOFBF, SCRATCH_BUFFER_SIZE);
	/* fully buffered: unless fflush() or page full, it is there */
	return;
    }
    rewind (fp_scratch);
    return;
} /* end reset_scratch() */


/* free fp_scratch: as if reset_scratch() was never called */
void kill_scratch()
{
    if (fp_scratch == NULL) return;
    fclose (fp_scratch);
    fp_scratch = NULL;
    return;
} /* end kill_scratch() */


/* shorthand for fprintf (fp_scratch, ..) */
void scratch (char *format, ...)
{
    va_list ap;
    va_start (ap, format);
    vfprintf (fp_scratch, format, ap);
    va_end (ap);
    return;
} /* end scratch() */


#define BUFFER_SIZE (TERMCHAR+1)
/* Dump what has been scratched down to another stream "fp", but adding */
/* string "esc" before each line. esc==NULL would be interpreted as "". */
void dump_scratched_to (FILE *fp, char *esc)
{
    static char BUFFER [BUFFER_SIZE];
    char *end;
    int should_print_esc = TRUE;
    long position = 0, end_position;
    end_position = ftell (fp_scratch);
    if (end_position == 0) return;
    rewind (fp_scratch);
    while ( fgets (BUFFER, BUFFER_SIZE, fp_scratch) != NULL )
    {  /* there are stuff read */
	if (should_print_esc && (esc != NULL)) fputs (esc, fp);
	end = eos(BUFFER);
	position += end - BUFFER;
	if (position >= end_position)
	{
	    *(end + end_position - position) = EOS;
	    fputs (BUFFER, fp);
	    return;
	}
	else fputs (BUFFER, fp);
	should_print_esc = (*(end-1)=='\n');
    }
    clearerr (fp_scratch);
    return;
} /* end dump_scratched_to() */
#undef BUFFER_SIZE


#ifdef _scratch_TEST
void main()
{
    reset_scratch();
    scratch ("This is bad ..\n"); 
    dump_scratched_to (stdout, "# ");
    /* reset_scratch(); */
    scratch ("This is worse.\n");
    dump_scratched_to (stdout, "# ");
    /* kill_scratch(); */
    reset_scratch();
    scratch ("This is better\n");
    dump_scratched_to (stdout, "# ");
    reset_scratch();
    dump_scratched_to (stdout, "# ");
}
#endif  /* _scratch_TEST */
