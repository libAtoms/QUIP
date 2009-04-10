/****************************************/
/* libIO:                               */
/*                                      */
/* Strings, parsers & files             */
/*                                      */
/* Dec.12, 1999  Ju Li <liju99@mit.edu> */
/****************************************/

#include "IO.h"

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

/* a simple mechanism for keeping track of multiple jterms */

static pid_t in_pid [NOFILE] = {0};
static pid_t out_pid [NOFILE] = {0};
static pid_t shifted_in_pid [NOFILE] = {0};
static pid_t shifted_out_pid [NOFILE] = {0};

/**************************************************************************/
/* open a dumb terminal at "DISPLAY" (NULL means default $DISPLAY), piped */
/* to open streams "in" (NULL means no) and "out" (NULL means no); return */
/* a unique identification number which will be needed in jclose().       */
/**************************************************************************/
pid_t jopen (FILE *in, FILE *out, char *DISPLAY)
{
    int i, in_fd, out_fd, fa[2], fb[2];
    char BUFFER [INTEGER_DECIMALCHAR+1];
    pid_t pid;

    if (in == NULL)   /* he doesn't want an in-stream */
	in_fd = -1;
    else
    {
	in_fd = fileno (in);
	if (VALID_FILENO(in_fd)) fflush (in);
	else in_fd = NOFILE;
	if (pipe(fa) < 0)
	{
	    printf ("error: jopen: cannot create pipe to stdin\n");
	    perror ("jopen: ");
	    exit(1);
	}
    }

    if (out == NULL)   /* he doesn't want an out-stream */
 	out_fd = -1;
    else
    {
	out_fd = fileno (out);
	if (VALID_FILENO(out_fd)) fflush (out);
	else out_fd = NOFILE;
	setbuf (out, NULL);
	if (pipe(fb) < 0)
	{
	    printf ("error: jopen: cannot create pipe to stdout\n");
	    perror ("jopen: ");
	    exit(1);
	}
    }

    if ( (in_fd==-1) && (out_fd==-1) )
    {
	printf ("error: jopen: are you playing with me?\n"
		"both \"in\" and \"out\" streams are bad.\n");
	exit(1);
    }
    
    if ( (pid=fork()) < 0 )
    {
        printf ("error: jopen: cannot create child\n");
        perror ("jopen: ");
        exit(1);
    }
    else if (pid == 0)   /* child process */
    {   
        if ( VALID_FILENO(in_fd) )
	{
	    close (fa[0]);
	    /* shift the pipe end to STDOUT_FILENO */
	    if (fa[1] != STDOUT_FILENO)
	    {
		dup2 (fa[1], STDOUT_FILENO);
		close (fa[1]);
	    }
	}
	if ( VALID_FILENO(out_fd) )
	{
	    close (fb[1]);
	    /* shift the pipe end to STDIN_FILENO */
	    if (fb[0] != STDIN_FILENO)
	    {
		dup2 (fb[0], STDIN_FILENO);
		close (fb[0]);
	    }
	}
	/* application context now has spiked STDIN_FILENO/STDOUT_FILENO */
	for (i=0; i<NOFILE; i++)
	{ /* this child should only have what it needs to know */
	    if (VALID_FILENO(in_fd)&&(i==STDOUT_FILENO)) continue;
	    if (VALID_FILENO(out_fd)&&(i==STDIN_FILENO)) continue;
	    close(i);
	}
	snprintf(BUFFER, INTEGER_DECIMALCHAR+1, "jterm %d", (int)getpid());
	if (DISPLAY != NULL)
	    execlp ( "jterm", "jterm",
		     "-display", DISPLAY,
		     "-title", BUFFER,
		     "-n", BUFFER,
		     "-bg", "gray40",
		     "-fg", "white",
		     "-cr", "yellow",
		     "-fn", "7x13",
		     (char *)0 );
	else
	    execlp ( "jterm", "jterm",
		     "-title", BUFFER,
		     "-n", BUFFER,
		     "-bg", "gray40",
		     "-fg", "white",
		     "-cr", "yellow",
		     "-fn", "7x13",
		     (char *)0 );
	exit(0);
    }
    
    /* parent process */
    if (VALID_FILENO(in_fd))
    {
	close (fa[1]);
	/* save the old in_fd in an efficient manner */
	i = dup(in_fd);
	in_pid [in_fd] = pid;
	shifted_in_pid [i] = pid;
	/* shift the pipe end to in_fd */
	dup2 (fa[0], in_fd);
	close (fa[0]);
    }
    if (VALID_FILENO(out_fd))
    {
	close (fb[0]);
	/* save the old out_fd in an efficient manner */
	i = dup(out_fd);
	out_pid [out_fd] = pid;
	shifted_out_pid [i] = pid;
	/* shift the pipe end to out_fd */
	dup2 (fb[1], out_fd);
	close (fb[1]);
    }
    return (pid);
} /* end jopen() */


/* close the pipe(s) and bury the child. If pid==0 we close all pipes */
void jclose (pid_t pid)
{
    int i, j, process_status;
    if (pid <= 0)
    { /* close all pipes */
	for (i=0; i<NOFILE; i++)
	{
	    if (out_pid[i] > 0)
	    { /* shift the stream fileno from pipe end back to out_pid[i] */
		for (j=0; j<NOFILE; j++)
		    if (shifted_out_pid[j] == out_pid[i])
			break;
		dup2 (j, i); /* old i (pipe end) will be closed by dup2() */
		close (j);
		/* by closing the write end the child needn't quit */
		shifted_out_pid[j] = out_pid[i] = 0;
	    }
	    if (in_pid[i] > 0)
	    { /* shift the stream fileno from pipe end back to in_pid[i] */
		for (j=0; j<NOFILE; j++)
		    if (shifted_in_pid[j] == in_pid[i])
			break;
		dup2 (j, i); /* old i (pipe end) will be closed by dup2() */
		close (j);
		/*****************************************************/
		/* Only when ctrl-D in jterm shuts down the inflow   */
		/* or when SIGPIPE encountered (by our action above) */
		/* does the forked child terminate. In the first     */
		/* case the child is already dead now and we just    */
		/* fetch its cold corpse from the kernel. In the     */
		/* latter case we'll wait for a warm one.            */
		/*****************************************************/
		kill (in_pid[i], SIGKILL);  /* make sure it's dead */
		waitpid (in_pid[i], &process_status, 0);
		shifted_in_pid[j] = in_pid[i] = 0;
	    }
	}
    }
    else
    {
	for (i=0; i<NOFILE; i++)
	{
	    if (out_pid[i] == pid)
	    { /* shift the stream fileno from pipe end back to out_pid[i] */
		for (j=0; j<NOFILE; j++)
		    if (shifted_out_pid[j] == out_pid[i])
			break;
		dup2 (j, i); /* old i (pipe end) will be closed by dup2() */
		close (j);
		/* by closing the write end the child needn't quit */
		shifted_out_pid[j] = out_pid[i] = 0;
	    }
	    if (in_pid[i] == pid)
	    { /* shift the stream fileno from pipe end back to in_pid[i] */
		for (j=0; j<NOFILE; j++)
		    if (shifted_in_pid[j] == in_pid[i])
			break;
		dup2 (j, i); /* old i (pipe end) will be closed by dup2() */
		close (j);
		/*****************************************************/
		/* Only when ctrl-D in jterm shuts down the inflow   */
		/* or when SIGPIPE encountered (by our action above) */
		/* does the forked child terminate. In the first     */
		/* case the child is already dead now and we just    */
		/* fetch its cold corpse from the kernel. In the     */
		/* latter case we'll wait for a warm one.            */
		/*****************************************************/
		kill (in_pid[i], SIGKILL);  /* make sure it's dead */
		waitpid (in_pid[i], &process_status, 0);
		shifted_in_pid[j] = in_pid[i] = 0;
	    }
	}
    }
    return;
} /* end jclose() */


#ifdef _jterm_TEST
#define BUFFER_SIZE 20
void main()
{
    int jid1, age;
    double height;
    char c, name[BUFFER_SIZE];
    jid1 = jopen (stdin, stdout, NULL);
    printf ("Name please?\n");
    Gets (name);
    printf ("\nPlease type your age and height:\n");
    scanf ("%d %lf", &age, &height); getchar();
    printf ("\nYou have lived for %d years...\nway too long.\n",
	    age);
    sleep(3);
    jclose (jid1);
    printf ("back alive, %s (%d %g)!\n", name, age, height);
    jid1 = jopen (stdin, stdout, NULL);
    printf ("back alive again, %s (%d %g)!\n", name, age, height);
    sleep(3);
    printf ("back alive again, %s (%d %g)!\n", name, age, height);
    sleep(5);
    jclose (jid1);
}
#endif  /* _jterm_TEST */


pid_t JTERM_STDID = -1;

#ifdef _Jterm_TEST
void main()
{
    int c;
    FILE *a = ropen("/etc/passwd");
    Jopen(NULL);
    while ((c=fgetc(a)) != EOF) putchar(c);
    printf ("\npress a key to exit ... ");
    getchar();
    Jclose();
    printf ("you've get your stdout back\n");
}
#endif  /* _Jterm_TEST */
