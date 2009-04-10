/**************************************/
/* libTimer: -lm                      */
/*                                    */
/* Time and Schedule Manager          */
/*                                    */
/* Nov 19 1999 Ju Li <liju99@mit.edu> */
/**************************************/

#include "Timer.h"


/************************************************************************/
/* Multiple task, piecewise linear schedule manager:                    */
/*                                                                      */
/* No more than SCHEDULE_MAX_TASKS-1 tasks;                             */
/* Each task id. by a unique token[] string (< SCHEDULE_MAX_TOKENSIZE); */
/*                                                                      */
/* A task has a time domain [tmin, tmax], a value domain x, and a set   */
/* of control points (t_i, x_i) which are ordered by t_i. Inquiries of  */
/* x using schedule_value(t) produces linear interpolation if t is      */
/* between two control points; otherwise, it is set to the value of the */
/* nearest control point.                                               */
/*                                                                      */
/* schedule_create supports simplified OR unordered specification of    */
/* the control points (these two features don't work well together) via */
/* schedule_str[] = "cntl_pt_str_1 cntl_pt_str_2 ..."                   */
/*                                                                      */
/* Format of cntl_pt_str:                                               */
/*                                                                      */
/* "[Ttime] value"  or  "[Ffraction] value"  or  simply "value"         */
/*                                                                      */
/* Processing Flowchart:                                                */
/*                                                                      */
/* 1. Sequentially parse schedule_str[] to cntl_pt_str's                */
/*                                                                      */
/* 2. Plug in tmin, tmax to evaluate the times of those cntl_pt_str's   */
/* with [Ttime] or [Ffraction].                                         */
/*                                                                      */
/* 3. For simple cntl_pt_str's, assign it tmin if it is the first       */
/* cntl_pt_str in schedule_str[], tmax if it is the last cntl_pt_str in */
/* schedule_str[]. Then, we mesh t regularly between two already        */
/* assigned cntl_pt_str's, to those yet unassigned simple cntl_pt_str's */
/* in between, respecting the original order in schedule_str[].         */
/*                                                                      */
/* 4. Re-order the control points in time.                              */
/*                                                                      */
/* For example, (assume tmin=0, tmax=1),                                */
/*                                                                      */
/* "0.  1.  2." regularly mesh [0,1] to (0, 0.), (0.5, 1.), (1, 2.)     */
/*                                                                      */
/* "F0.1 -1.  T0.4 0.2" will be understood as having control points     */
/* (f0.1, -1.) and (t0.4, 0.2), in which f0.1 could be smaller or       */
/* larger than t0.4.                                                    */
/*                                                                      */
/* "-2.  F0.2 -1.  -0.5 -0.1 F0.5 0.2 1." would generate 6 control      */
/* points, (tmin, -2.), (f0.2, -1.), (f0.3, -0.5), (f0.4, -0.1),        */
/* (f0.5, 0.2), (tmax, 1.)                                              */
/************************************************************************/

static ScheduleList pool[SCHEDULE_MAX_TASKS]={0};

/* for private use in stdlib qsort() */
static int schedule_control_compare (const void *a, const void *b)
{
    return( *((double *)a) > *((double *)b) ? 1 : -1 );
} /* end schedule_control_compare() */


/* Add a schedule of name token[]; return pointer to the */
/* parsed, interpreted, and time-ordered control points. */
ScheduleList * schedule_create
(char token[], char schedule_str[], double tmin, double tmax)
{
    int i, j, k, l, first_of_a_pair;
    char *start, *end;
    ScheduleList *p, *q;
    static char buffer[SCHEDULE_MAX_STRLEN];
    static double control[SCHEDULE_MAX_PIECES][2];
#define SCHEDULE_NONSENSE (-176943285.21301e37)
    /* copy schedule_str[] to buffer[] */
    strncpy (buffer, schedule_str, SCHEDULE_MAX_STRLEN-1);
    buffer[SCHEDULE_MAX_STRLEN-1] = '\0';
    for ( p = pool; (p->next != NULL) && /* root node is pool[0] */
              strncmp(p->next->token, token, SCHEDULE_MAX_TOKENSIZE-1);
          p = p->next );
    if ( p->next != NULL )
    {
        printf("error: schedule_create: token name \"%s\" "
               "already used.\n", token);
        exit(1);
    }
    for ( q = pool+1; q < pool+SCHEDULE_MAX_TASKS; q++)
        if (*(q->token) == '\0')
        { /* look for free space in pool[] */
            p->next = q;
            q->next = NULL;
            strncpy (q->token, token, SCHEDULE_MAX_TOKENSIZE-1);
            q->token[SCHEDULE_MAX_TOKENSIZE-1] = '\0';
            first_of_a_pair = 1; /* beginning of a control record */
            for (start=buffer; (*start=='\t') || (*start==' '); start++);
            for (j=0; (j<SCHEDULE_MAX_PIECES)&&(*start!='\0'); start=end)
            {  /* j is the total number of control points */
                for (end=start; (*end!='\t')&&(*end!=' ')&&(*end!='\0');
                     end++);
                if (*end != '\0')
                { /* eat up ensuing white space */
                    while ((*end=='\t') || (*end==' ')) end++;
                    *(end-1) = '\0';
                }
                if ( ((*start>='0') && (*start<='9')) ||
                     (*start=='-') || (*start=='+') ) /* numeric string */
                {
                    if (first_of_a_pair) control[j][0] = SCHEDULE_NONSENSE; 
                    control[j++][1] = atof(start); /* control value    */
                    first_of_a_pair = 1;             /* start a new pair */
                }
                else if ( ((*start=='T') || (*start=='t')) && first_of_a_pair )
                { /* explicit control time */
                    control[j][0] = atof(start+1);
                    first_of_a_pair = 0;  /* wait for value input */
                }
                else if ( ((*start=='F') || (*start=='f')) && first_of_a_pair )
                { /* fractional control time */
                    control[j][0] = tmin + (tmax-tmin)*atof(start+1);
                    first_of_a_pair = 0;  /* wait for value input */
                }
                else
                {
                    printf ("error: schedule_create: cannot parse "
                            "substring \"%s\".\n", start);
                    exit(1);
                }
            }
            if (j >= SCHEDULE_MAX_PIECES)
            {
                printf("error: schedule_create: SCHEDULE_MAX_PIECES =\n"
                       "%d overflow for\n"
                       "\"%s\".\n", SCHEDULE_MAX_PIECES, buffer);
                exit(1);
            }
            q->s = (double *) malloc(2*(q->n=j)*sizeof(double));
            if (control[0][0] == SCHEDULE_NONSENSE)
                control[0][0] = tmin;
            if (control[j-1][0] == SCHEDULE_NONSENSE)
                control[j-1][0] = tmax;
            for (i=1; i<j-1; i++)
                if (control[i][0] == SCHEDULE_NONSENSE)
                { /* equal mesh in t between defined control points */
                    for (k=i-1; control[k][0]==SCHEDULE_NONSENSE; k--);
                    for (l=i+1; control[l][0]==SCHEDULE_NONSENSE; l++);
                    control[i][0] = control[k][0] +
                        (control[l][0]-control[k][0]) * (i-k) / (l-k);
                }
            /* sort control times in ascending order */
            qsort (&control[0][0], q->n, 2*sizeof(double),
                   schedule_control_compare);
            memcpy(q->s, &control[0][0], 2*q->n*sizeof(double));
            return (q);
        }
    printf("error: schedule_create: SCHEDULE_MAX_TASKS = "
           "%d overflow.\n", SCHEDULE_MAX_TASKS);
    exit(1);
    return(NULL);
} /* end schedule_create() */


#define schedule_fetch_for(client_function) \
for (q = pool; (q->next != NULL) && \
strncmp(q->next->token, token, SCHEDULE_MAX_TOKENSIZE-1); q = q->next); \
if (q->next == NULL) \
{ printf("error: client_function: token name \"%s\" " \
         "not found.\n", token); exit(1);}


/* print out this schedule's description <fast> */
void schedulelist_print (FILE *out, ScheduleList *q)
{
    int i;
    if (!out) return;
    fprintf(out, "Schedule token = \"%s\", ", q->token);
    fprintf(out, "%d control points of (time, value):\n", q->n);
    for (i=0; i<q->n; i++)
        fprintf(out, "(%g, %g)\n", q->s[2*i], q->s[2*i+1]);
    return;
} /* end schedulelist_print() */


/* print out this schedule's description <slow> */
void schedule_print (FILE *out, char token[])
{
    ScheduleList *q;
    schedule_fetch_for(schedule_print);
    schedulelist_print(out, q->next);
    return;
} /* end schedule_print() */


/* Evaluate a schedule at time t using a set of properly */
/* arranged control points of ScheduleList structure. <fast> */
double schedulelist_value (ScheduleList *q, double t)
{
    int i;
    if (t <= q->s[0]) return (q->s[1]);
    if (t >= q->s[2*q->n-2]) return (q->s[2*q->n-1]);
    for (i=0; i<q->n-1; i++)
        if ( (t > q->s[2*i]) && (t <= q->s[2*i+2]) ) break;
    if (fabs(q->s[2*i] - q->s[2*i+2]) < SCHEDULE_TIME_TINY)
        return( (q->s[2*i+1] + q->s[2*i+3]) / 2 );
    else return (q->s[2*i+1] + (q->s[2*i+3] - q->s[2*i+1]) *
                 (t - q->s[2*i]) / (q->s[2*i+2] - q->s[2*i]));
} /* end schedulelist_value() */


/* Evaluate a schedule of a registered token at time t. <slow> */
double schedule_value (char token[], double t)
{
    ScheduleList *q;
    schedule_fetch_for(schedule_value);
    return (schedulelist_value(q->next, t));
}  /* end schedule_value() */


/* Multiply the values of all control points by factor <fast> */
void schedulelist_scale_value (ScheduleList *q, double factor)
{
    int i;
    for (i=0; i<q->n; i++) q->s[2*i+1] *= factor;
    return;
} /* end schedulelist_scale_value() */


/* Multiply the values of all control points by factor <slow> */
void schedule_scale_value (char token[], double factor)
{
    ScheduleList *q;
    schedule_fetch_for (schedule_scale_value);
    schedulelist_scale_value (q->next, factor);
    return;
} /* end schedule_scale_value() */


/* See if the schedule value remains constant during period */
/* [t1,t2]: if it does, then return the upper and lower bound */
/* that contains [t1,t2]; otherwise return NULL. <fast>     */
double * schedulelist_const_period
(ScheduleList *q, double t1, double t2)
{
    int i, j, k;
    double ti, tf;
    static double bound[2];
    ti = (t1 < t2)? t1:t2;
    tf = (t1 < t2)? t2:t1;
    for (i=0; i<q->n; i++)
        if (fabs(ti - q->s[2*i]) < SCHEDULE_TIME_TINY) break;
        else if (ti < q->s[2*i])
            if (i==0) break;
            else {i--; break;}
    if (ti < q->s[2*0]) bound[0] = ti;
    else /* search left for lower bounds */
        for (k=i; k>=0; k--)
            if (fabs(q->s[2*k+1]-q->s[2*i+1])<SCHEDULE_TIME_TINY)
                bound[0] = q->s[2*k];
    for (j=q->n-1; j>=0; j--)
        if (fabs(tf - q->s[2*j]) < SCHEDULE_TIME_TINY) break;
        else if (tf > q->s[2*j])
            if (j==q->n-1) break;
            else {j++; break;}
    if (tf > q->s[2*(q->n-1)]) bound[1] = tf;
    else /* search right for upper bounds */
        for (k=j; k<q->n; k++)
            if (fabs(q->s[2*k+1]-q->s[2*j+1])<SCHEDULE_TIME_TINY)
                bound[1] = q->s[2*k];
    for (; i<j; i++)
        if ((fabs(q->s[2*i+1]-q->s[2*j+1]) > SCHEDULE_TIME_TINY))
            return(NULL);
    return(bound);
} /* end schedulelist_const_period() */


/* See if the schedule value remains constant during period */
/* [t1,t2]: if it does, then return the upper and lower bound */
/* that contains [t1,t2]; otherwise return NULL. <slow> */
double *schedule_const_period (char token[], double t1, double t2)
{
    ScheduleList *q;
    schedule_fetch_for (schedule_const_period);
    return (schedulelist_const_period(q->next, t1, t2));
} /* end schedule_const_period() */


/* delete a schedule and free resources */
void schedule_kill (char token[])
{
    ScheduleList *q;
    schedule_fetch_for (schedule_kill);
    free (q->next->s);
    q->next->n = 0;
    *(q->next->token) = '\0';
    q->next = q->next->next;
    return;
} /* end schedule_kill() */


/* kill all currently registered schedules */
void schedule_kill_all()
{
    ScheduleList *q;
    for (q=pool; q->next != NULL;)
    {
        free (q->next->s);
        q->next->n = 0;
        *(q->next->token) = '\0';
        q->next = q->next->next;
    }
    return;
} /* end schedule_kill_all() */


#ifdef _schedule_TEST
int main(int argc, char *argv[])
{
    ScheduleList *a =
        schedule_create ("a", "  F0.1 10 20 t0.9 30  ", 0., 1.);
    ScheduleList *b =
        schedule_create ("b", "  F0.1 10 F0.11 20 t0.9 15  ", 0., 1.);
    ScheduleList *c;
    double *x, *u, *y, *z;
    schedule_print(stdout,"a"); printf("\n");
    schedule_print(stdout,"b"); printf("\n");
    printf ("%g %g %g\n\n", schedule_value("a", -0.1),
            schedule_value("a", 0.2),
            schedule_value("a", 0.7));
    schedule_kill ("a");
    printf ("%g %g %g %g\n\n", schedule_value("b", -0.1),
            schedule_value("b", 0.2),
            schedule_value("b", 1.7),
            schedule_value("b", 0.103));
    schedule_kill_all();
    schedule_create("a", "1 2 3 4 5", 0., 1.);
    schedule_print(stdout,"a"); printf("\n");
    c = schedule_create
        ("c", "  F0.1 10   10 10 10 10   F0.9 10 F0.95 20 ", 0., 1.);
    schedule_print(stdout,"c"); printf("\n");
    x = schedule_const_period("c", -10, 1.5);
    u = schedule_const_period("c", 0.89, 0.901);
    y = schedule_const_period("c", 0.4, 0.5);
    z = schedule_const_period("c", 0.9, 0.4);
    printf ("c = %ld %ld %g %g %g %g\n", (long)x, (long)u,
            y[0], y[1], z[0], z[1]);
    return(0);
}
#endif


/*******************************/
/* Multiple-schedule functions */
/*******************************/

/* Grab a schedule string (ended by INTER_SCHEDULE_SEPARATOR) from stream */
char *grab_schedule_str (FILE *fp)
{
    static char *buffer = NULL;
    int i, c;
    if (!fp)
        if (buffer)
        { /* pass in fp = NULL, then free the allocation */
            free (buffer);
            return (buffer=NULL);
        }
    for (i=0; (c=getc(fp))!=EOF;)
    { /* ignore line-breaks */
        if (c == '\n') continue;
        buffer = realloc (buffer, (i+1)*sizeof(char));
        buffer[i] = c;
        if (buffer[i] == INTER_SCHEDULE_SEPARATOR)
        {
            buffer[i] = 0;
            return (buffer);
        }
        i++;
    }
    /* if EOF encountered without INTER_SCHEDULE_SEPARATOR */
    return (NULL);
} /* end grab_schedule_str() */


/* Download n schedules from stream, set them up according */
/* to tokens[], and save pointers to schedules[]. Return 0 */
/* if ALL schedules are constant during (tmin, tmax).      */
int grab_schedules
(FILE *fp, double tmin, double tmax, int n, char *tokens[],
 ScheduleList *schedules[], FILE *info)
{
    int i, not_constant=0;
    char *schedule_str;
    for (i=0; i<n; i++)
    {
        schedule_str = grab_schedule_str(fp);
        if (!schedule_str)
        {
            fprintf(stderr, "error: grab_schedules: stream ends \n"
                    "without separator '%c' for token[%d] = \"%s\"\n",
                    INTER_SCHEDULE_SEPARATOR, i, tokens[i]);
            exit(1);
        }
        schedules[i] = schedule_create
            (tokens[i], schedule_str, tmin, tmax);
        schedulelist_print (info, schedules[i]);
        if ( schedulelist_const_period(schedules[i], tmin, tmax)
             == NULL ) not_constant = 1;
    }
    grab_schedule_str (NULL);  /* free space */
    return (not_constant);
} /* end grab_schedules() */


#ifdef _grab_schedules_TEST
#define tmin 0
#define tmax 1
int main (int argc, char *argv[])
{
    int n;
    char *tokens[] = {"H11","H12","H13","H22","H23","H33"};
    ScheduleList *schedules[SCHEDULE_MAX_TASKS];
    printf ("Input number of schedules: ");
    fscanf (stdin, "%d", &n);
    printf ("Input %d schedules, separated by '%c' :\n", n,
            INTER_SCHEDULE_SEPARATOR);
    n = grab_schedules (stdin, tmin, tmax, n, tokens, schedules, stdout);
    printf ("schedules are %s.\n", n?"not all constants":"all constants");
    return (0);
}
#undef tmax
#undef tmin
#endif /* _grab_schedules_TEST */


/******************************************************/
/* Schedules: (dimension, ScheduleList *) repeat, 0.  */
/* First determine if the Schedules are all constants */
/* in [t1,t2]. If so, return the earliest time they   */
/* still all remain so. Otherwise return t2.          */
/******************************************************/
double schedules_calm (double t1, double t2, ...)
{
    int i,n;
    va_list ap;
    ScheduleList **q;
    double tcalm, *this;
    /* IEEE 754 64-bit standard: 2.225074e-308 to 1.797693e308 */
    tcalm = -1e307;
    va_start (ap, t2);
    while ((n=va_arg(ap,int)) > 0)
        for (q=va_arg(ap,ScheduleList **),i=0; i<n; i++)
            if (this = schedulelist_const_period (q[i], t1, t2))
            {
                if (*this > tcalm) tcalm = *this;
            }
            else
            {
                va_end(ap);
                return(t2);
            }
    va_end(ap);
    return(tcalm);
} /* end schedules_calm() */
