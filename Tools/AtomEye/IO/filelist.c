/***************************************/
/* libIO:                              */
/*                                     */
/* Strings, parsers & files            */
/*                                     */
/* Dec.12, 1999 Ju Li <liju99@mit.edu> */
/***************************************/

#include "IO.h"

/* compile a list of similar filenames */

#define MAXDIGITS 16
/* Return a list of "numerically similar" filenames to "seed_fname" */
/* and sort them; return 0 if no digits are found in "seed_fname".  */
/* Remember to free dynamical memory later using globfree(glob_t*). */
int numerically_sorted_glob (char *seed_fname, glob_t *globbuf)
{
    int i,j;
    char *p,*q,*suffix,*prefix,*pattern;
    for (p=eos(seed_fname); p>=seed_fname; p--)
        if (ISDIGIT(*p) &&
            (! ( (p-2>=seed_fname) &&
                 (*p    =='2') &&
                 (*(p-1)=='z') &&
                 (*(p-2)=='b') ) )) break;
    if (p==seed_fname-1) return(0);
    suffix = IOClone(p+1);
    for (q=p; q>=seed_fname; q--)
        if (ISNOTDIGIT(*q)) break;
    prefix = IOClone(seed_fname);
    prefix[q+1-seed_fname] = EOS;
    pattern = IOalloc(strlen(prefix)+5*MAXDIGITS+strlen(suffix)+1);
    for (i=1; i<MAXDIGITS; i++)
    {
        strcpy (pattern, prefix);
        for (j=0; j<i; j++) strcat (pattern, "[0-9]");
        strcat (pattern, suffix);
        glob(pattern, (i==1) ? 0 : GLOB_APPEND, NULL, globbuf);
    }
    free(pattern);
    free(prefix);
    free(suffix);
    return (globbuf->gl_pathc);
} /* end numerically_sorted_glob() */
#undef MAXDIGITS


#ifdef _numerically_sorted_glob_TEST
int main (int argc, char *argv[])
{
    int i;
    glob_t globbuf;
    numerically_sorted_glob ("/tmp/pg1.bz2", &globbuf);
    for (i=0; i<globbuf.gl_pathc; i++)
        printf ("%s\n", globbuf.gl_pathv[i]);
    globfree (&globbuf);
    return (0);
}
#endif /* _numerically_sorted_glob_TEST */


/* Find the fname that is numerically similar to "start_fname" */
/* but is "glob_advance" down the list */
char *numerically_sorted_glob_advance
(char *start_fname, glob_t *globbuf, int glob_advance)
{
    int i,k;
    if (globbuf->gl_pathc<=0)
        pe ("numerically_sorted_glob_advance:\n"
            "glob list is not well defined.\n");
    for (i=0; i<globbuf->gl_pathc; i++)
        if (!strcmp(start_fname,globbuf->gl_pathv[i])) break;
    if (i == globbuf->gl_pathc)
        pe ("numerically_sorted_glob_advance: file\n"
            "%s\ndoes not exist in the glob list.\n", start_fname);
    k = (i + glob_advance) % INT(globbuf->gl_pathc);
    if (k < 0) k += globbuf->gl_pathc;
    return (globbuf->gl_pathv[k]);
} /* end numerically_sorted_glob_advance() */


/* memory efficient implementation */
char *Numerically_sorted_glob_advance (char *fname, int glob_advance)
{
    int i;
    glob_t globbuf;
    i = numerically_sorted_glob (fname, &globbuf);
    if (i <= 0)
    {
        /* pe ("Numerically_sorted_glob_advance: \n" */
        /* "glob list is not well defined.\n"); */
        return(fname);
    }
    strcpy(fname, numerically_sorted_glob_advance
           (fname, &globbuf, glob_advance));
    globfree (&globbuf);
    return(fname);
} /* end Numerically_sorted_glob_advance() */


/* seek the first in the similarity family */
char *Numerically_sorted_glob_first (char *fname)
{
    int i;
    glob_t globbuf;
    i = numerically_sorted_glob(fname, &globbuf);
    if (i <= 0)
    {
        /* pe ("Numerically_sorted_glob_first: \n" */
        /* "glob list is not well defined.\n"); */
        return(fname);
    }
    strcpy(fname, globbuf.gl_pathv[0]);
    globfree (&globbuf);
    return(fname);
} /* end Numerically_sorted_glob_first() */


/* seek the last in the similarity family */
char *Numerically_sorted_glob_last (char *fname)
{
    int i;
    glob_t globbuf;
    i = numerically_sorted_glob(fname, &globbuf);
    if (i <= 0)
    {
        /* pe ("Numerically_sorted_glob_last: \n" */
        /* "glob list is not well defined.\n"); */
        return(fname);
    }
    strcpy(fname, globbuf.gl_pathv[globbuf.gl_pathc-1]);
    globfree (&globbuf);
    return(fname);
} /* end Numerically_sorted_glob_last() */
