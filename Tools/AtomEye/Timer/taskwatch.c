/**************************************/
/* libTimer:                          */
/*                                    */
/* Time and Schedule Manager          */
/*                                    */
/* Nov 19 1999 Ju Li <liju99@mit.edu> */
/**************************************/

#include "Timer.h"


/******************************************************************/
/* Real time watches of multiple (< TASKWATCH_MAX_RECORDS) tasks, */
/* each id. by a token[] string (< TASKWATCH_MAX_TOKENSIZE bytes) */
/* Commands: 0: init, 1: check, 2: reset time origin, 3: destroy. */
/* Returns real time elapsed in seconds.                          */
/******************************************************************/

double taskwatch (char token[], int command)
{
    struct List
    {
	char token[TASKWATCH_MAX_TOKENSIZE];
	double init_time;
	struct List *next;
    };
    static struct List mempool[TASKWATCH_MAX_RECORDS]
	= {{{(char)0}, 0., (struct List *)NULL}};
    static struct timeval taskwatch_t;
#define TASKWATCH_TIME ((1E-6)*taskwatch_t.tv_usec+taskwatch_t.tv_sec)
    struct List *p, *q;
    double cc, dd;
    switch (command)
    { 
        case TASKWATCH_INIT:  /* start a new task watch */
            for ( p = mempool; /* root node is mempool[0] */
                  (p->next != NULL) &&
                      strncmp(p->next->token,
                              token, TASKWATCH_MAX_TOKENSIZE-1);
                  p = p->next );
            if ( p->next == NULL )
            { /* look for free space */
                for (q=mempool+1; q<mempool+TASKWATCH_MAX_RECORDS; q++)
                    if ( *(q->token)==(char)0 )
                    {
                        p->next = q;
                        q->next = NULL;
                        strncpy (q->token, token, TASKWATCH_MAX_TOKENSIZE-1);
                        gettimeofday (&taskwatch_t, NULL); /* last step */
                        return(q->init_time=TASKWATCH_TIME);
                    }
                printf("error: taskwatch: mempool %d records all used up.\n",
                       TASKWATCH_MAX_RECORDS-1);
                exit(1);
            }
            else
            {
                printf("error: taskwatch: token name \"%s\" already used.\n",
                       token);
                exit(1);
            }
            break;
        case TASKWATCH_CHECK: /* return time elapsed in seconds */
            gettimeofday (&taskwatch_t, NULL);  /* first step */
            for ( p = mempool;
                  (p->next != NULL) &&
                      strncmp(p->next->token,
                              token, TASKWATCH_MAX_TOKENSIZE-1);
                  p = p->next);
            if ( p->next == NULL )
            {
                printf("error: taskwatch: token name \"%s\" not found.\n",
                       token);
                exit(1);
            }
            else return (TASKWATCH_TIME - p->next->init_time);
            break;
        case TASKWATCH_RESET: /* return time elapsed in seconds */
            gettimeofday (&taskwatch_t, NULL); /* first step */
            dd = TASKWATCH_TIME;
            for ( p = mempool; /* and refresh the watch's start-point */
                  (p->next != NULL) &&
                      strncmp(p->next->token,
                              token, TASKWATCH_MAX_TOKENSIZE-1);
                  p = p->next);
            if ( p->next == NULL )
            {
                printf("error: taskwatch: token name \"%s\" not found.\n",
                       token);
                exit(1);
            }
            else
            {
                cc = p->next->init_time;
                gettimeofday (&taskwatch_t, NULL); /* last step */
                p->next->init_time = TASKWATCH_TIME;
                return (dd - cc);
            }
            break;
        case TASKWATCH_KILL: /* destroy the record and free resources */
            gettimeofday (&taskwatch_t, NULL); /* first step */
            for ( p = mempool;
                  (p->next != NULL) &&
                      strncmp(p->next->token,
                              token, TASKWATCH_MAX_TOKENSIZE-1);
                  p = p->next);
            if ( p->next == NULL )
            {
                printf("error: taskwatch: token name \"%s\" not found.\n",
                       token);
                exit(1);
            }
            else
            {
                cc = p->next->init_time;
                *(p->next->token) = (char)0;
                p->next = p->next->next;
                return (TASKWATCH_TIME - cc);
            }
            break;
        default:
            printf("error: taskwatch: invalid command %d.\n", command);
            exit(1);
    }
    return (SECONDS_SINCE_EPOCH);
} /* end taskwatch() */


#ifdef _taskwatch_TEST
int main(int argc, char *argv[])
{
    double a, b, c;
    printf ("%f %f\n", SECONDS_SINCE_EPOCH,
	    taskwatch("task1", TASKWATCH_INIT));
    waitfor(1.);
    printf ("%f\n", taskwatch("task1", TASKWATCH_RESET));
    waitfor(1.);
    printf ("%f\n", taskwatch("task1", TASKWATCH_RESET));
    waitfor(1.);
    printf ("%f\n", taskwatch("task1", TASKWATCH_KILL));
    return(0);
}
#endif

/* cc -D_Taskwatch_TEST taskwatch.c Timer.c -lm; time a.out */
