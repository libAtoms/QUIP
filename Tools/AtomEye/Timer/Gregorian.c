/**************************************/
/* libTimer:                          */
/*                                    */
/* Time and Schedule Manager          */
/*                                    */
/* Nov 19 1999 Ju Li <liju99@mit.edu> */
/**************************************/

#include "Timer.h"

/**********************/
/* Gregorian Calendar */
/**********************/


/* Cumulative days per month in a nonleap year. */
static double cdm[] = {0,31,59,90,120,151,181,212,243,273,304,334};


/* Matlab serial number for date (1 corresponds to 1-Jan-0000) */
double datenummx (GregorianDate *G)
{
    double result;
    if ((G->month < 1) || (G->month > 12))
        printf ("datenummx: month = %d\n", G->month);
    result = 365.*G->year + ceil(G->year/4.) -
        ceil(G->year/100.) + ceil(G->year/400.) + 
        cdm[G->month-1] + G->day;
    if (G->month > 2)
    {
        if ((G->year%4 == 0) && (G->year%100 != 0) || (G->year%400 == 0))
        {
            result += 1.;
        }
    }
    result += (G->hour*3600. + G->minute*60. + G->second) / 86400.;
    return (result);
} /* end datenummx() */


/* Matlab date serial number from ISO 8601 input string like "1970-01-01" */
double datenum29 (char *str)
{
    GregorianDate G[1]={{0,1,1,0,0,0}};
    char *newstr = IOClone(str);
    newstr[4] = 0;
    newstr[7] = 0;
    newstr[10] = 0;
    G->year   = atoi(newstr);
    G->month  = atoi(newstr+5);
    G->day    = atoi(newstr+8);
    free (newstr);
    return(datenummx(G));
} /* end datenum29() */


#ifdef _datenum
int main (int argc, char *argv[])
{
    double result;
    GregorianDate G[1]={{0,1,1,0,0,0}};
    if (argc == 1)
    {
        printf ("\nPurpose: Matlab serial number for date "
                "(1 corresponds to 1-Jan-0000).\n\n");
        printf ("Usage: %s year month day\n", argv[0]);
        printf ("   or, %s year month day hour minute second\n\n",
                argv[0]);
        printf ("Matlab function datenum(1970,01,01) = 719529\n\n");
        return (1);
    }
    if (argc > 1) G->year   = atoi(argv[1]);
    if (argc > 2) G->month  = atoi(argv[2]);
    if (argc > 3) G->day    = atoi(argv[3]);
    if (argc > 4) G->hour   = atoi(argv[4]);
    if (argc > 5) G->minute = atoi(argv[5]);
    if (argc > 6) G->second = atoi(argv[6]);
    result = datenummx(G);
    if (result==(int)result) printf ("%d\n", (int)result);
    else printf ("%.15g\n", result);
    return (0);
}
#endif /* _datenum */


/* Cumulative days per month in both nonleap and leap years. */
static double cdm0[] = {0,31,59,90,120,151,181,212,243,273,304,334,365};
static double cdml[] = {0,31,60,91,121,152,182,213,244,274,305,335,366};


/* Calculate date components from serial date number */
void datevecmx (double datenum, GregorianDate *G)
{ 
    double *cdm; 
    double ts, y;
    int iy, leap, mon;
    
    datenum = 86400*datenum;
    if (1)
    {
        datenum = floor(datenum+0.5);
    }
    ts = datenum;
    datenum = floor(datenum/60.);
    G->second = ts - 60.*datenum;
    ts = datenum;
    datenum = floor(datenum/60.);
    G->minute = ts - 60.*datenum;
    ts = datenum;
    datenum = floor(datenum/24.);
    G->hour = ts - 24.*datenum;

    datenum = floor(datenum);
    y = floor(datenum/365.2425);
    ts = datenum - (365.0*y + ceil(0.25*y)-ceil(0.01*y)+ceil(0.0025*y));
    if (ts <= 0)
    {
        y = y - 1.;
        datenum = datenum -
            (365.0*y + ceil(0.25*y)-ceil(0.01*y)+ceil(0.0025*y));
    }
    else datenum = ts;

    G->year = y;
    iy = (int) y;
    leap = ((iy%4 == 0) && (iy%100 != 0) || (iy%400 == 0));

    cdm = (leap ? cdml : cdm0);

    mon = (int) datenum/29.-1;
    if (datenum > cdm[mon+1]) mon++;
    G->month = mon+1;

    datenum = datenum - cdm[mon];
    G->day = datenum;

    return;
} /* end datevecmx() */


#ifdef _datevec
int main (int argc, char *argv[])
{
    double datenum;
    GregorianDate G[1]={{0,1,1,0,0,0}};
    if (argc == 1)
    {
        printf ("\nPurpose: inverse function of datenum.\n\n");
        printf ("Usage: %s 719529\n", argv[0]);
        printf ("   or, %s 719529.23\n\n",
                argv[0]);
        printf ("Matlab function datevec(719529.23) = 1970 1 1 5 31 12\n\n");
        return (1);
    }
    else datenum = atof(argv[1]);
    datevecmx(datenum, G);
    printf ("%d\n", G->year);
    printf ("%02d\n", G->month);
    printf ("%02d\n", G->day);
    if (datenum!=(int)datenum)
    {
        printf ("%02d\n", G->hour);
        printf ("%02d\n", G->minute);
        printf ("%02d\n", G->second);
    }
    return (0);
}
#endif /* _datevec */
