/**************************************/
/* libTimer: -lm                      */
/*                                    */
/* Time and Schedule Manager          */
/*                                    */
/* Nov 19 1999 Ju Li <liju99@mit.edu> */
/**************************************/

#ifndef _Timer_h
#define _Timer_h

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/resource.h>

/** Timer.c: **/

/* Return a string like "Wed Jun 30 21:49:08 1993" */
char *date_string();

extern struct timeval TIMER;
/* seconds (double) since Jan 1, 1970. More accurate than time() */
#define SECONDS_SINCE_EPOCH \
(gettimeofday(&TIMER,NULL)?(1E-6*TIMER.tv_usec+TIMER.tv_sec):\
 (1E-6*TIMER.tv_usec+TIMER.tv_sec))

/* sleep for a prescribed time, occupying no CPU resources */
void waitfor (double period_in_seconds);
#ifdef _IRIX64
#define WAITFOR_OVERHEAD (1./CLOCKS_PER_SEC) /* subtract off certain overhead */
#else
#define WAITFOR_OVERHEAD (1./CLOCKS_PER_SEC) /* subtract off certain overhead */
#endif

/* Sufficient time to see a warning message or notice */
#define SUFFICIENT_TIME_TO_NOTICE_IN_SECONDS  3.0

/* convert seconds to time string appreciable by human beings */
char *earthtime (double seconds);


/* chronometer.c: */

/* accurate usertime benchmark of a single task */
extern struct rusage chronometer_start, chronometer_stop;
#define start_chronometer() getrusage(RUSAGE_SELF, &chronometer_start)
#define stop_chronometer()  getrusage(RUSAGE_SELF, &chronometer_stop)
#define stopped_usertime() ((chronometer_stop.ru_utime.tv_usec - \
			     chronometer_start.ru_utime.tv_usec)*1E-6 + \
			     chronometer_stop.ru_utime.tv_sec - \
			     chronometer_start.ru_utime.tv_sec)

/* taskwatch.c */
    
/******************************************************************/
/* Real time watches of multiple (< TASKWATCH_MAX_RECORDS) tasks, */
/* each id. by a token[] string (< TASKWATCH_MAX_TOKENSIZE bytes) */
/* Commands: 0: init, 1: check, 2: reset time origin, 3: destroy. */
/* Returns real time elapsed in seconds.                          */
/******************************************************************/
#define TASKWATCH_MAX_RECORDS  256
#define TASKWATCH_MAX_TOKENSIZE 16
#define TASKWATCH_INIT  0
#define TASKWATCH_CHECK 1
#define TASKWATCH_RESET 2
#define TASKWATCH_KILL  3
double taskwatch (char token[], int command);



/* linear_scheduler.c: */

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

#define SCHEDULE_MAX_TASKS     16
#define SCHEDULE_MAX_TOKENSIZE 16
#define SCHEDULE_MAX_PIECES    128
#define SCHEDULE_MAX_STRLEN   (16*SCHEDULE_MAX_PIECES)
#define SCHEDULE_TIME_TINY    (1000*2.220446e-16)
/* recursive typedef */
typedef struct ScheduleList
{
    char token[SCHEDULE_MAX_TOKENSIZE]; /* schedule identifier */
    int n;                         /* number of control points */
    double *s;              /* (time, value) of control points */
    struct ScheduleList *next;
} ScheduleList;

/* Add a schedule of name token[]; return pointer to the */
/* parsed, interpreted, and time-ordered control points. */
ScheduleList * schedule_create
(char token[], char schedule_str[], double tmin, double tmax);

/* print out this schedule's description <fast> */
void schedulelist_print (FILE *out, ScheduleList *q);

/* print out this schedule's description <slow> */
void schedule_print (FILE *out, char token[]);

/* Evaluate a schedule at time t using a set of properly */
/* arranged control points of ScheduleList structure. <fast> */
double schedulelist_value (ScheduleList *q, double t);

/* Evaluate a schedule of a registered token at time t. <slow> */
double schedule_value (char token[], double t);

/* Multiply the values of all control points by factor <fast> */
void schedulelist_scale_value (ScheduleList *q, double factor);

/* Multiply the values of all control points by factor <slow> */
void schedule_scale_value (char token[], double factor);

/* See if the schedule value remains constant during period */
/* [t1,t2]: if it does, then return the upper and lower bound */
/* that contains [t1,t2]; otherwise return NULL. <fast>     */
double * schedulelist_const_period
(ScheduleList *q, double t1, double t2);

/* See if the schedule value remains constant during period */
/* [t1,t2]: if it does, then return the upper and lower bound */
/* that contains [t1,t2]; otherwise return NULL. <slow> */
double * schedule_const_period (char token[], double t1, double t2);

/* delete a schedule and free resources */
void schedule_kill (char token[]);

/* kill all currently registered schedules */
void schedule_kill_all();


/*******************************/
/* Multiple-schedule functions */
/*******************************/

#define INTER_SCHEDULE_SEPARATOR  ';'
/* Grab a schedule string (ended by INTER_SCHEDULE_SEPARATOR) from stream */
char *grab_schedule_str (FILE *fp);

#define SCHEDULES_ALL_CONSTANTS    0
/* Download n schedules from stream, set them up according */
/* to tokens[], and save pointers to schedules[]. Return 0 */
/* if ALL schedules are constant during (tmin, tmax).      */
int grab_schedules
(FILE *fp, double tmin, double tmax, int n, char *tokens[],
 ScheduleList *schedules[], FILE *info);

/******************************************************/
/* Schedules: (dimension, ScheduleList *) repeat, 0.  */
/* First determine if the Schedules are all constants */
/* in [t1,t2]. If so, return the earliest time they   */
/* still all remain so. Otherwise return t2.          */
/******************************************************/
double schedules_calm (double t1, double t2, ...);


/* Gregorian.c: */

/**********************/
/* Gregorian Calendar */
/**********************/

typedef struct
{
    int year;
    int month;
    int day;
    int hour;
    int minute;
    int second;
} GregorianDate;

/* Matlab serial number for date (1 corresponds to 1-Jan-0000) */
double datenummx (GregorianDate *G);

/* Matlab date serial number from ISO 8601 input string like "1970-01-01" */
double datenum29 (char *str);

/* Calculate date components from serial date number */
void datevecmx (double datenum, GregorianDate *G);


/* fools.c: */

/*************************************************/
/* Function prototypes to fool the compiler into */
/* believing that certain variable is useful --  */
/* for time benchmarking purposes.               */
/*************************************************/

void charfool (char c);
void ucharfool (unsigned char c);
void shortfool (short i);
void ushortfool (unsigned short i);
void intfool (int i);
void uintfool (unsigned int i);
void longfool (long l);
void ulongfool (unsigned long l);
void longlongfool (long long l);
void ulonglongfool (unsigned long long l);
void floatfool (float d);
void doublefool (double d);
void voidPfool (void *d);
void charPfool (char *c);
void ucharPfool (unsigned char *c);
void shortPfool (short *i);
void ushortPfool (unsigned short *i);
void intPfool (int *i);
void uintPfool (unsigned int *i);
void longPfool (long *l);
void ulongPfool (unsigned long *l);
void longlongPfool (long long *l);
void ulonglongPfool (unsigned long long *l);
void floatPfool (float *d);
void doublePfool (double *d);

#endif /* _Timer_h */
