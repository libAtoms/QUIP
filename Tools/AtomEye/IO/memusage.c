/***************************************/
/* libIO:                              */
/*                                     */
/* Strings, parsers & files            */
/*                                     */
/* Dec.12, 1999 Ju Li <liju99@mit.edu> */
/***************************************/

#include "IO.h"

/* Memory usage snapshot */

#if \
  defined(_Linux) || \
  defined(_SunOS) || \
  defined(_alpha) || \
  defined(_Darwin)
#define BUFFERSIZE 256

/* Linux /proc inquiry (faster? safer?, than ps) */
static long linux_proc_mem_inquiry()
{
    FILE *fp;
    char buffer[BUFFERSIZE], *p;
    int memsize;
    sprintf(buffer, "/proc/%d/status", getpid());
    if (!Fexists(buffer)) return(-1);
    fp = ropen(buffer);
    while (fgets(buffer,BUFFERSIZE,fp))
        if (strcasestr(buffer,"rss"))
        {
            while (1)
            {
                for (p=buffer; *p!=EOS; p++)
                    if (ISDIGIT(*p))
                    {
                        sscanf(p, "%d", &memsize);
                        fclose(fp);
                        return (memsize*1000);
                    }
                if (*p==EOS) fgets(buffer,BUFFERSIZE,fp);
            }
        }
    fclose(fp);
    return(-1);
} /* linux_proc_mem_inquiry() */

long memusage()
{
    FILE *fp;
    char buffer[BUFFERSIZE];
    int memsize;
#if \
  defined(_Linux) || \
  defined(_alpha)
    if ( (memsize=linux_proc_mem_inquiry()) >= 0 )
        return (memsize);
#endif
#define PSPATH "ps"
    sprintf(buffer, PSPATH" -p %d -o rss", getpid());
    fp = popen(buffer, "r");
    if (ISNULL(fp))
    {
        pr ("memusage: \"%s\" is not under your run path.\n", PSPATH);
        return (0);
    }
    fscanf(fp, "%s\n%d", buffer, &memsize);
    pclose(fp);
    return (memsize*1000);
} /* memusage() */
#undef PSPATH
#undef BUFFERSIZE

#else

struct rusage mem_usage;

#endif

#define BUFFERSIZE 16
/* quick view of memory */
char *strmem (long bytes)
{
    static char BUFFER [BUFFERSIZE];
    if (bytes >= 1E12)
        snprintf (BUFFER, BUFFERSIZE, "%.3g TB", bytes/1E12);
    else if (bytes >= 1E9)
        snprintf (BUFFER, BUFFERSIZE, "%.3g GB",  bytes/1E9);
    else if (bytes >= 1E6)
        snprintf (BUFFER, BUFFERSIZE, "%.3g MB",  bytes/1E6);
    else if (bytes >= 1E3)
        snprintf (BUFFER, BUFFERSIZE, "%.3g KB",  bytes/1E3);
    else snprintf (BUFFER, BUFFERSIZE, "%ld Bytes", bytes);
    BUFFER [BUFFERSIZE-1] = 0;
    return (BUFFER);
} /* end strmem() */
#undef BUFFERSIZE

#ifdef _pmemusage_TEST
#define MATSIZE (1024 * 256)
int main (int argc, char *argv[])
{
    double x[MATSIZE];
    bzero (x, MATSIZE*sizeof(double));
    printf ("pid = %d\n", getpid());
    pmemusage();
    /* printf ("%d\n", mem_usage.ru_maxrss); */
    /* printf ("%d\n", mem_usage.ru_idrss); */
    printf ("Press Ctrl-C to stop.\n");
    while(1);
    return (0);
}
#endif /* _pmemusage_TEST */
