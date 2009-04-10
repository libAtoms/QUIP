/**************************************/
/* libTimer:                          */
/*                                    */
/* Time and Schedule Manager          */
/*                                    */
/* Nov 19 1999 Ju Li <liju99@mit.edu> */
/**************************************/

#include "Timer.h"


/* Return a string like "Wed Jun 30 21:49:08 1993" */
char *date_string()
{
    char *p, *q;
    time_t elapsed;
    time (&elapsed);
    p = ctime (&elapsed);
    for (q=p; (*q!='\n')&&(*q!=0); q++);
    *q=0;
    return(p);
} /* end date_string() */


#ifdef _date_string_TEST
int main (int argc, char *argv[])
{
    printf ("%s\n", date_string());
    return (0);
}
#endif /* _date_string_TEST */



struct timeval TIMER;


/* sleep for a prescribed time, occupying no CPU resources */
void waitfor (double period_in_seconds)
{
    static struct timeval time;
#ifdef WAITFOR_OVERHEAD
    period_in_seconds -= WAITFOR_OVERHEAD;
#endif
    if (period_in_seconds < 0.) period_in_seconds=0.;
    time.tv_sec  = floor(period_in_seconds);
    time.tv_usec = (period_in_seconds-time.tv_sec)*1E6;
    select (0, NULL, NULL, NULL, &time);
    /* expect microseconds of context switching time */
    return;
} /* end waitfor() */

#ifdef _waitfor_TEST
#define COUNTDOWN 10
int main(int argc, char *argv[])
{
    int i = COUNTDOWN;
    printf("Countdown to %d seconds ...\n", COUNTDOWN);
    do {
	printf(" %d", i--);
	fflush(stdout);
	waitfor(1.);
    } while (i>0);
    printf("\nNow goto hell.\n");
    return(0);
}
/* cc -D_waitfor_TEST Timer.c -lm; time a.out */
#endif

/* 1999, long-march.MIT.EDU: */
/* 100MHz motherboard freq. x 4.5 Pentium II CPU (32bit)  */
/* 512K level-2 cache, 256MB main memory, 25 GB hard-disk */
/* l00Mb/s Ethernet, 128b/8MB video card, 17' 1024x768 CRT */
/* Athena / RedHat 5.2 / Linux 2.0.36 / 24 bit X-server */


/* convert seconds to time string appreciable by human beings */
char *earthtime (double seconds)
{
#define SOLAR_YEAR_IN_S       (31556925.974592)
#define SIDEREAL_YEAR_IN_S    (31558149.540288)
#define ANOMALISTIC_YEAR_IN_S (31558433.011776)
    const struct time_unit
    {
	char *name;
	double scale;
    } earthman[] =
      { {"s",  1},
	{"min", 60},
	{"hr", 60},
	{"day", 24},
	{"mon", 30},
	{"yr", SOLAR_YEAR_IN_S/2592000.} };
    /* use solar year because it's the most relevant to season */
#define MAXUNITS (sizeof(earthman)/sizeof(struct time_unit))
#define MAXTIMESTRLEN 128
    static char time[MAXTIMESTRLEN], *p;
    int i, denom;
    double scale;
    seconds /= earthman[0].scale;
    for (scale=1,i=1; i<MAXUNITS; i++)
	scale*=earthman[i].scale;
    for (p=time,i=MAXUNITS-1; i>0; i--)
    {
	denom = floor(seconds/scale);
	if (denom > 0)
	{
	    snprintf(p, time+MAXTIMESTRLEN-p,
		     "%d %s, ", denom, earthman[i].name);
	    while((*p)!=(char)0) p++;
	}
	seconds -= denom*scale;
	scale /= earthman[i].scale;
    }
    snprintf(p, time+MAXTIMESTRLEN-p,
	     "%.3f %s", seconds, earthman[0].name);
    return(time);
} /* end earthtime() */
/* ... living on earth */

#ifdef _earthtime_TEST
int main(int argc, char *argv[])
{
    time_t elapsed;
    time(&elapsed);
    printf("Current time is %s", ctime(&elapsed));
    printf("nearly %s\n"
	   "has passed since Epoch (January 1, 1970), Earthlings.\n",
	   earthtime((double)elapsed));
    return(0);
}
#endif


/*
  Time Notes:

  The best frame to understand earth is to co-move with earth's center
  of mass but not rotating. In that frame, the sun moves around us in
  nearly closed & elliptical, planar path; all other stars are fixed;
  and the earth is spinning like crazy. We call the earth's spinning J
  the north-pole, which is the normal to the equatorial plane. On the
  equator we pick Greenwich to be 0 degree called meridian, and the
  plane of north-pole to meridian to south-pole is called the meridian
  plane. The nearly elliptical path sun traverses in our frame is
  called the ecliptic, its angular version called zodiac. A section of
  the zodiac points to a certain constellation (which corresponds to
  an animal). A zodiac section also corresponds to
  month/season. Viewed on earth, certain constellation near the sun
  co-spins with the sun about the north-pole (Polaris), but the sun's
  partner changes by the month.

  Season is not mainly due to the closeness to the sun but by the
  sun's angle.  When the sun is spinning near the north-pole, people
  in northern sphere feel warm because 1) day (defined as when the sun
  is above horizon) is longer, 2) sunlight angle is more plumb,
  whereas when the sun is spinning closer to the south-pole, 1) day is
  shorter, 2) sunlight angle is more oblique. The point in space that
  sun crosses the equatorial plane from south to north last time is
  called vernal equinox. 2000 years ago vernal equinox points to
  Aries, now vernal equinox points to Pisces. The reason is because J
  itself, and therefore the equatorial plane, is slowly changing, due
  to a precession cycle of 26,000 yrs (wobbling).

  Our standing assumption would be earth's radius is 0 compared to
  ecliptic, and the ecliptic is 0 compared to galaxy.  The "small
  numbers" we are going to account for are a) the gravitational torque
  on earth as one side of it is closer to the sun at a time, which has
  non-zero time average and causes the 26,000 yr precession (similar
  to a precessing top, I guess the physical point J pierces earth does
  not change, so Amundsen's north pole does not shift on earth even if
  where it points to has changed). And b) other planets' attractions
  that makes earth's trajectory, and therefore the ecliptic,
  non-closed. This is called orbital precession. Earth's orbital
  precession is entirely normal (Newtonian), whereby Mercury's
  abnormal, or excess orbital precession can only be explained by
  general relativity.

  The basic SI time unit is the second, which has atomic clock
  definition, and so by it are hour and day. A solar year is defined
  to be the interval between two vernal equinoxes, which turns out to
  be (365.25636042 + 1.1E-7 TE) days, where TE is one Temporal Epoch
  (36525. days). A sidereal year is the time interval between two
  ecliptic points of the same orientation in our previously defined
  frame, irrespective of earth's spinning and precessing. It is longer
  than the solar year because earth's equatorial plane precesses
  slowly "to meet" the sun. A sidereal year is theoretically cleaner
  because it does not depend on earth's J at all. But practically
  solar year is more important because season drift is more important
  to us than star drift (sun changes spring-break partner from Aries
  to Pisces in the last two millennium, so what?).

  Note that the primary unit of the sidereal time system is the the
  sidereal day, which is the time interval between the crossings of
  the current (this year's) vernal equinox (Pisces, of course) to the
  earth's meridian plane. Sidereal day, hour and second are mainly for
  the use of earth's star-watching buff; it depends on the magnitude
  of J, not its direction (though you could argue that the start of a
  sidereal day does depend on the vernal equinox and therefore on the
  precession). Sidereal year is mainly for the comparing of galactic
  time with one canonical definition of the earth's orbital
  pseudo-period. There is no logical connection between the magnitudes
  of the sidereal year and the sidereal day: one could be defined
  "without" the earth, and the other could be defined "without" the
  sun.

  Same happens for the solar year and the "atomic" day, or for that
  matter, the sidereal day with "atomic" day - which is what the TV,
  and therefore we, are using. That is, there is no longer any logical
  connection between 7:15 on our watch (presumably adjusted to a clock
  in Paris) and the sun's position, although scientists, when they
  were making the standards, had put up the best effort to make 86,400
  atomic seconds equal to the "present" sidereal day, or roughly the
  current sunrise-sunrise cycle. There is always some error in that,
  because of 1) measurement errors in both the sidereal day AND the
  atomic frequency then, which screwed up the conversion factor and
  unfortunately set permanent; and more fundamentally, 2), the
  sidereal day is NOT a constant and change with Epochs, so there CAN
  be no rigorous conversion whatsoever. The scientists' solution to
  that linearly accumulating error is: whenever "significant" drift
  has happened between our watch and say, sunrise-sunrise cycle, we
  just shift the clock in Paris by one second, without much ado. This
  is a drastic departure from the attitude of the first mechanical
  clock inventor: though he could not implement it in full, he
  nevertheless defined one second to be an integer partition of a
  sidereal day (which he held to be a constant because he didn't know
  better) and then tried to mimic it, i.e., discretizing the day into
  86,400 parts. We, now have already defined the second but it turns
  out not to fit exactly, so we say we don't care anymore and just cut
  out the fringe once for a while. Therefore we see units based on
  less constant periods are secondary to more constant periods, and
  humans always like the most constant ones.
  
  Discretization of the solar year into sidereal or atomic days (which
  is also called a calendar) involves a much larger fringe of error,
  and subsequent drifts in season is more dramatic, and demands
  constant pruning. Julius Caesar and Gregory XIII decreed adding to
  the calendar a leap year every 4 (except centennial) and 400 years,
  which is a nice little trick of discretizating 365.2422-day solar
  year into 365 97/400. This will work well for a long time in keeping
  the vernal equinox within one day's error from each year's March 21.

  So-called "anomalistic year" is another definition of the earth's
  orbital pseudo-period, using successive perihelions. All three
  definitions of the year are NOT constants and change very slowly
  with TE. Human tried to trust in the periodicity of celestial bodies
  including sun/earth orbit and earth's spinning with respect to
  Galaxy. None of them turn out to be usable.

  Thu Nov 25 12:37:58 1999  Ju Li
  
  */
