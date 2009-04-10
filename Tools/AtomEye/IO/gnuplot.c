/****************************************/
/* libIO:                               */
/*                                      */
/* Strings, parsers & files             */
/*                                      */
/* Dec.12, 1999  Ju Li <liju99@mit.edu> */
/****************************************/

#include "IO.h"


/****************************************************/
/* pipe graphics scripts to powerful gnuplot engine */
/****************************************************/


/* Close pipe, kill gnuplot, free memory, reset pointers */
void gnuplot_free (Gnuplot *G)
{
    if (NOTNULL(G))
    {
        G->gnuplotoptions = NULL;
        if (NOTNULL(G->pipestream)) pclose(G->pipestream);
        G->pipestream = NULL;
        Free(G->epsfilename);
        G->epsoptions = NULL;
    }
    return;
} /* end gnuplot_free() */


/* Set eps filename: if NULL the eps filename default will be presumed */
void gnuplot_epsfilename (Gnuplot *G, char *epsfilename)
{
    if (NOTNULL(epsfilename))
    {
        REALLOC( gnuplot_epsfilename,
                 G->epsfilename,
                 strlen(epsfilename)+1,
                 char );
        strcpy( G->epsfilename, epsfilename );
    }
    else Free(G->epsfilename);
    return;
} /* end gnuplot_epsfilename() */


#define GNUPLOT_PATH "gnuplot"  /* name of the executable in sh */
/* Open gnuplot window. Returns GNUPLOT_NOT_OK or GNUPLOT_OK */
int gnuplot_recreate (Gnuplot *G, char *epsfilename)
{
    if (ISNULL(G)) return(GNUPLOT_NOT_OK);
    if (ISNULL(G->pipestream))
        G->pipestream =
            popen( str3(GNUPLOT_PATH, " ", 
                        (G->gnuplotoptions!=NULL) ?
                        G->gnuplotoptions :
                        GNUPLOT_DEFAULT_OPTIONS ),
                   "w");
    if (ISNULL(G->pipestream))
    {
        pr ("gnuplot_recreate: \"%s\" is not under your run path.\n",
            GNUPLOT_PATH);
        return (GNUPLOT_NOT_OK);
    }
    gnuplot_epsfilename (G, epsfilename);
    return(GNUPLOT_OK);
} /* end gnuplot_recreate() */
#undef GNUPLOT_PATH


/* Ask gnuplot to run a script file */
void gnuplot_run_file (Gnuplot *G, char *fname)
{
    register int c;
    FILE *script = ropen(fname);
    while ( (c=fgetc(script)) != EOF ) fputc( c, G->pipestream );
    fclose (script);
    fflush (G->pipestream);
    return;
} /* end gnuplot_run_file() */


/* Ask gnuplot to run a script character-buffer */
void gnuplot_run_buffer (Gnuplot *G, char *buffer)
{
    while ( *buffer != EOS ) fputc( *(buffer++), G->pipestream );
    fflush (G->pipestream);
    return;
} /* end gnuplot_run_buffer() */


/* Ask gnuplot to run print-outs */
void gnuplot_run_printf (Gnuplot *G, char *format, ...)
{
    va_list ap;
    va_start (ap, format);
    vfprintf (G->pipestream, format, ap);
    va_end (ap);
    fflush (G->pipestream);
    return;
} /* end gnuplot_run_printf() */


/* Set gnuplot output to X11 */
void gnuplot_term_x11 (Gnuplot *G)
{
    gnuplot_run_buffer (G, "set term x11\n");
} /* end gnuplot_term_x11() */


/* Set gnuplot output to eps: if epsfilename=NULL, the old */
/* G->epsfilename is used; else the new name will be used. */
void gnuplot_term_eps (Gnuplot *G, char *epsfilename)
{
    if (NOTNULL(epsfilename)) gnuplot_epsfilename (G, epsfilename);
    gnuplot_run_printf (G,
                        "set term postscript eps %s\n"
                        "set output '%s'\n",
                        (G->epsoptions!=NULL)?
                        G->epsoptions :
                        GNUPLOT_DEFAULT_EPSOPTIONS, 
                        GNUPLOT_EPSFILENAME_NOW(G));
    return;
} /* end gnuplot_term_eps() */


/**************************************************************/
/* Synchronize with the child gnuplot process by asking it to */
/* save the drawing into an eps file and make sure of getting */
/* it. If epsfilename=NULL, the old G->epsfilename will be    */
/* used; else the new epsfilename will be used.               */
/**************************************************************/
void gnuplot_eps_sync (Gnuplot *G, char *epsfilename)
{
    remove (GNUPLOT_EPSFILENAME_NOW(G));
    gnuplot_term_eps (G, epsfilename);
    gnuplot_run_buffer (G, "replot\n");
    while (!Fexists(GNUPLOT_EPSFILENAME_NOW(G)));
    return;
} /* end gnuplot_term_eps() */


#ifdef _gnuplot_anim_TEST
#define STEPS 1000
#define FNAME "gnuplot_anim.eps"
int main (int argc, char *argv[])
{
    int i;
    double phase;
    Gnuplot G[1] = {{0}};
    gnuplot_recreate (G,FNAME);
    gnuplot_run_buffer (G, "set xrange [0:2*pi]\n");
    for (i=0,phase=0; i<STEPS; i++,phase+=0.05)
    {
        gnuplot_term_x11 (G);
        gnuplot_run_printf (G,"plot cos(x+%g)\n", phase);
        gnuplot_eps_Sync (G);
        usleep(100000);
    }
    gnuplot_free(G);
    return (0);
}
#endif /* _gnuplot_anim_TEST */
