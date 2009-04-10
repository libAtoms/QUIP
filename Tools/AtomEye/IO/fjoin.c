/****************************************/
/* libIO:                               */
/*                                      */
/* Strings, parsers & files             */
/*                                      */
/* Dec.12, 1999  Ju Li <liju99@mit.edu> */
/****************************************/

#include "IO.h"

/* each open pipe of the parent is sucked by one live child */
static pid_t childpid [NOFILE] = {0};

/* create an output stream that feeds to multiple open streams */
FILE *fjoin (int number_of_streams, ...)
{
    pid_t pid;
    int i, j, streams_fd[NOFILE], fd[2];
    FILE *streams[NOFILE], *fp_parent;
    struct stat status;
    va_list ap;
    char c;
    
    va_start (ap, number_of_streams);
    for (i=0; i<number_of_streams; i++)
    {
        streams[i] = va_arg (ap, FILE *);
        if (streams[i] == NULL)
        { /* seg-fault protection */
            printf ("error: fjoin: stupid! the %s stream is NULL\n",
                    word_for_order(i+1));
            exit(1);
        }
        streams_fd[i] = fileno (streams[i]);
        if (INVALID_FILENO(streams_fd[i]))
        {
            printf ("error: fjoin: the %s stream has fileno %d\n",
                    word_for_order(i+1), streams_fd[i]);
            perror ("fjoin: ");
            exit(1);
        }
        if ( fstat(streams_fd[i], &status) < 0 )
        {
            printf ("error: fjoin: the %s stream is invalid (errno = %d)\n",
                    word_for_order(i+1), errno);
            perror ("fjoin: ");
            exit(1);
        }
        if ( ! ( S_ISREG(status.st_mode) || S_ISCHR(status.st_mode) ) )
        {  /* only joins regular files: this of course excludes pipes */
            printf ("error: fjoin: the %s stream is not a regular file\n",
                    word_for_order(i+1));
            exit(1);
        }
        /* parent clears stdio buffer because child inherits everything */
        fflush (streams[i]);
        /* stdio, after all, is written by somebody of my stand */
    }
    va_end (ap);
    
    if (pipe(fd) < 0)
    {
        printf ("error: fjoin: cannot create pipe\n");
        perror ("fjoin: ");
        exit(1);
    }
    
    if ( (pid=fork()) < 0 )
    {
        printf ("error: fjoin: cannot create child\n");
        perror ("fjoin: ");
        exit(1);
    }
    else if (pid == 0)  /* child process */
    {   
        close (fd[1]);
        /***************************************************************/
        /* If I am a second child, I inherit all open file descriptors */
        /* from my parent (as explicitly dup()ed), which of course     */
        /* include the write end of the first child's pipe. For kernel */
        /* to kill that pipe when it should, I have to close all open  */
        /* ends so the parent has the sole control over the pipe life. */
        /***************************************************************/
        for (i=0; i<NOFILE;)
        {  /* this child should only have what it needs to know */
            if ( i == fd[0] ) goto keep_it;
            for (j=0; j<number_of_streams; j++)
                if ( i == streams_fd[j] ) goto keep_it;
/*********************************************************************/
/* the parent still has the streams' descriptors and can still write */
/* but I hope he uses buffered fprintf().. since I use atom pipe and */
/* blocked read/write, I should be (very probably) ahead of him.     */
/*********************************************************************/
            close(i);
          keep_it:
            i++;           /* to satisfy the IRIX compiler */
        }
        /***************************************************************/
        while ( read(fd[0],&c,1) > 0 )
            for (i=0; i<number_of_streams; i++)
            {
                while ( (j=write(streams_fd[i],&c,1)) == 0 );
                if (j < 0)
                {
                    printf ("error: fjoin: errno %d writing to the %s "
                            "stream\n", errno, word_for_order(i+1));
                    perror ("fjoin: ");
                    exit(1);
                }
            }
        exit(0);
    }
    /* parent process */
    close (fd[0]);
    if ( (fp_parent=fdopen(fd[1],"w")) == NULL )
    {
        printf ("error: fjoin: cannot attach output fileno to a stream\n");
        perror ("fjoin: ");
        exit(1);
    }
    /* this stdio stream is buffer-less like stderr */
    setbuf (fp_parent, NULL);
    childpid [fd[1]] = pid;
    return (fp_parent);
} /* end fjoin() */


/* delete the pipe and bury the dead child */
void fbreakup (FILE *fp)
{
    int fd, process_status;
    struct stat fd_status;
    pid_t pid;

    if (fp == NULL)
    { /* seg-fault protection */
        printf ("error: fbreakup: stupid! the stream is NULL\n");
        exit(1);
    }
    fd = fileno (fp);
    if ( INVALID_FILENO(fd) )
    {
        printf ("error: fbreakup: illegal file number %d\n", fd);
        perror ("fbreakup: ");
        exit(1);
    }
    if ( fstat(fd, &fd_status) < 0 )
    {
        printf ("error: fbreakup: the stream is invalid (errno = %d)\n",
                errno);
        perror ("fbreakup: ");
        exit(1);
    }
    if ( (pid=childpid[fd]) == 0 )
    {
        printf ("error: fbreakup: this stream was not created by fjoin()!\n");
        exit(1);
    }
    childpid [fd] = 0;
    fclose (fp);  /* logic confirmed, I wait for the consequences ... */
    while (waitpid (pid, &process_status, 0) < 0)
        if (errno != EINTR)
        {
            printf ("error: fbreakup: non-interrupt errno %d while "
                    "waiting for child %d\n", errno, pid);
            perror ("fbreakup: ");
            exit(1);
        }
    return;
} /* end fbreakup() */

/************************************************************************/
/* 1. only open, writable, and regular files (not pipes) can be joined  */
/* 2. fjoin(), like fopen(), is re-usable, i.e. many joins can co-exist */
/* 3. one can still fprintf to an individual stream while it is joined  */
/*    but it is not guaranteed to be in the apparent order: to ensure   */
/*    that, use sync() before printing to that individual component     */
/************************************************************************/

#ifdef _fjoin_TEST
void main()
{
    FILE *a, *b, *c, *combination1, *combination2;
    a = wopen("/tmp/a");
    b = wopen("/tmp/b");
    c = wopen("/tmp/c");

    fprintf (a, "a\n"); fprintf (b, "b\n"); fprintf (c, "c\n");
    
    combination1 = fjoin(1,a);
    combination2 = fjoin(2,b,c);
    
    fprintf (combination1, "damn you\n");
    fprintf (combination2, "what the hell\n");
    
    fbreakup (combination1);
    fbreakup (combination2);
    
    fprintf (a, "a\n"); fprintf (b, "b\n"); fprintf (c, "c\n");
    
    fclose(a);
    fclose(b);
    fclose(c);
}
#endif  /* _fjoin_TEST */


/* for use with single pipe handlers ftie()/fbrk() macros */
FILE *ft = NULL;

#ifdef _ft_TEST
void main()
{
    FILE *a, *b, *c;
    a = wopen("/tmp/a");
    b = wopen("/tmp/b");
    c = wopen("/tmp/c");

    fprintf (a, "a\n"); fprintf (b, "b\n"); fprintf (c, "c\n");
    ftie(a,b);
    
    fprintf (ft, "writing a,b\n");
    fprintf (a, "a\n"); fprintf (b, "b\n"); fprintf (c, "c\n");

    ftie(b,c);
    fprintf (ft, "writing b,c\n");

    ftie4(stdout,a,b,c);
    fprintf (ft, "bye-bye!\n");
    
    fbrk();  /* free a file handler for the parent and kill a child */
}
#endif  /* _ft_TEST */
