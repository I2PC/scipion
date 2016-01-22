/*****************************************************************************
 *	System includes
 ****************************************************************************/
#include	<stdio.h>
#include	<string.h>
#include	<time.h>

/*****************************************************************************
 *	Toolbox defines
 ****************************************************************************/
#include	"configs.h"
#include	"debug.h"
#include	"error.h"
#include	"ttimestamp.h"

/*****************************************************************************
 *	Toolbox includes
 ****************************************************************************/
#include	"messagedisplay.h"
#include	"timestamp.h"

/*****************************************************************************
 *	Conditional includes
 ****************************************************************************/
/* None */

/*****************************************************************************
 *	Declaration of static procedures
 ****************************************************************************/
 /* None */

/*****************************************************************************
 *	Definition of static procedures
 ****************************************************************************/
/* None */

/*****************************************************************************
 *	Definition of extern procedures
 ****************************************************************************/
/*--------------------------------------------------------------------------*/
extern int		GetTimeStamp
				(
					struct TTimeStamp
							*TimeStamp			/* output time stamp */
				)

/* get a time stamp */
/* system time is the number of seconds elapsed since
	the last system boot (MacOS) or the start of the current process (UNIX) */
/* calendar time is the number of seconds since Midnight, January 1, 1904 A.D. */
/* UTC time is Coordinated Universal Time */
/* local time is UTC time after correction for time zone and Daylight Saving Time */

{ /* begin GetTimeStamp */

	struct tm
			*p;
	time_t	Now;
	clock_t	ticks;
	int		Status = !ERROR;

/**/DEBUG_WRITE_ENTERING(GetTimeStamp,
/**/	"About to get a time stamp")

	ticks = clock();
	Status = (ticks == (clock_t)(-1));
	if (Status == ERROR) {
		WRITE_ERROR(GetTimeStamp, "System processing time is not available")
/**/	DEBUG_WRITE_LEAVING(GetTimeStamp, "Done")
		return(Status);
	}
	Now = time(&Now);
	Status = (Now == (time_t)(-1));
	if (Status == ERROR) {
		WRITE_ERROR(GetTimeStamp, "Calendar time is not available")
/**/	DEBUG_WRITE_LEAVING(GetTimeStamp, "Done")
		return(Status);
	}
	TimeStamp->SecondsSinceProcessBirth = (double)ticks / (double)CLOCKS_PER_SEC;
	TimeStamp->SecondsSince1904_01_01_00H00 = (double)Now;
	p = gmtime(&Now);
	Status = (p == (struct tm *)NULL);
	if (Status == ERROR) {
		WRITE_ERROR(GetTimeStamp, "UTC time is not available")
/**/	DEBUG_WRITE_LEAVING(GetTimeStamp, "Done")
		return(Status);
	}
	p = (struct tm *) memcpy(&(TimeStamp->UTCTime), p, sizeof(struct tm));
	p = localtime(&Now);
	p = (struct tm *) memcpy(&(TimeStamp->LocalTime), p, sizeof(struct tm));
/**/DEBUG_WRITE_LEAVING(GetTimeStamp, "Done")
	return(Status);
} /* end GetTimeStamp */

/*--------------------------------------------------------------------------*/
extern int		GetElapsedTime
				(
					struct TTimeStamp
							OldTimeStamp,		/* oldest time stamp */
					struct TTimeStamp
							NewTimeStamp,		/* newest time stamp */
					double	*SystemTimeDifference,
												/* elapsed system time in seconds */
					double	*RealTimeDifference	/* true elapsed time in seconds */
				)

/* get the time elapsed between two time stamps */
/* system time is the fractional number of seconds since the last system boot */
/* calendar time is the integer number of seconds since Midnight, January 1, 1904 A.D. */
/* UTC time is Coordinated Universal Time */
/* local time is UTC time after correction for time zone and Daylight Saving Time (if available) */
/* note that, if it is possible to compute the elapsed time between two historical dates,
	the elapsed time between two dates that follow the release of the compiler
	may be inaccurate because of the time drift of UTC time with respect to
	International Atomic Time (TAI). This drift is +/- 1 leap second about once a year,
	and is corrected by the International Earth Rotation Service (IERS),
	which impedes any algorithmic approach. Unfortunately, TAI is not available in ANSI-C.
	The list of leap seconds (so far, all positive) is available at
	ftp://maia.usno.navy.mil/ser7/tai-utc.dat
	The first correction was +10 seconds on January 1, 1972 A.D.

	1972 JAN  1 =JD 2441317.5  TAI-UTC=  10.0       S + (MJD - 41317.) X 0.0      S
	1972 JUL  1 =JD 2441499.5  TAI-UTC=  11.0       S + (MJD - 41317.) X 0.0      S
	1973 JAN  1 =JD 2441683.5  TAI-UTC=  12.0       S + (MJD - 41317.) X 0.0      S
	1974 JAN  1 =JD 2442048.5  TAI-UTC=  13.0       S + (MJD - 41317.) X 0.0      S
	1975 JAN  1 =JD 2442413.5  TAI-UTC=  14.0       S + (MJD - 41317.) X 0.0      S
	1976 JAN  1 =JD 2442778.5  TAI-UTC=  15.0       S + (MJD - 41317.) X 0.0      S
	1977 JAN  1 =JD 2443144.5  TAI-UTC=  16.0       S + (MJD - 41317.) X 0.0      S
	1978 JAN  1 =JD 2443509.5  TAI-UTC=  17.0       S + (MJD - 41317.) X 0.0      S
	1979 JAN  1 =JD 2443874.5  TAI-UTC=  18.0       S + (MJD - 41317.) X 0.0      S
	1980 JAN  1 =JD 2444239.5  TAI-UTC=  19.0       S + (MJD - 41317.) X 0.0      S
	1981 JUL  1 =JD 2444786.5  TAI-UTC=  20.0       S + (MJD - 41317.) X 0.0      S
	1982 JUL  1 =JD 2445151.5  TAI-UTC=  21.0       S + (MJD - 41317.) X 0.0      S
	1983 JUL  1 =JD 2445516.5  TAI-UTC=  22.0       S + (MJD - 41317.) X 0.0      S
	1985 JUL  1 =JD 2446247.5  TAI-UTC=  23.0       S + (MJD - 41317.) X 0.0      S
	1988 JAN  1 =JD 2447161.5  TAI-UTC=  24.0       S + (MJD - 41317.) X 0.0      S
	1990 JAN  1 =JD 2447892.5  TAI-UTC=  25.0       S + (MJD - 41317.) X 0.0      S
	1991 JAN  1 =JD 2448257.5  TAI-UTC=  26.0       S + (MJD - 41317.) X 0.0      S
	1992 JUL  1 =JD 2448804.5  TAI-UTC=  27.0       S + (MJD - 41317.) X 0.0      S
	1993 JUL  1 =JD 2449169.5  TAI-UTC=  28.0       S + (MJD - 41317.) X 0.0      S
	1994 JUL  1 =JD 2449534.5  TAI-UTC=  29.0       S + (MJD - 41317.) X 0.0      S
	1996 JAN  1 =JD 2450083.5  TAI-UTC=  30.0       S + (MJD - 41317.) X 0.0      S
	1997 JUL  1 =JD 2450630.5  TAI-UTC=  31.0       S + (MJD - 41317.) X 0.0      S
	1999 JAN  1 =JD 2451179.5  TAI-UTC=  32.0       S + (MJD - 41317.) X 0.0      S
*/

{ /* begin GetElapsedTime */

	time_t	Old, New;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(GetElapsedTime, SystemTimeDifference, Status,
/**/	"No SystemTimeDifference")
/**/DEBUG_CHECK_NULL_POINTER(GetElapsedTime, RealTimeDifference, Status,
/**/	"No RealTimeDifference")
/**/DEBUG_RETURN_ON_ERROR(GetElapsedTime, Status)
/**/DEBUG_WRITE_ENTERING(GetElapsedTime,
/**/	"About to compute time difference")

	Old = mktime(&(OldTimeStamp.UTCTime));
	Status = (Old == (time_t)(-1));
	if (Status == ERROR) {
		WRITE_ERROR(GetElapsedTime, "Invalid OldTimeStamp")
/**/	DEBUG_WRITE_LEAVING(GetElapsedTime, "Done")
		return(Status);
	}
	New = mktime(&(NewTimeStamp.UTCTime));
	Status = (New == (time_t)(-1));
	if (Status == ERROR) {
		WRITE_ERROR(GetElapsedTime, "Invalid NewTimeStamp")
/**/	DEBUG_WRITE_LEAVING(GetElapsedTime, "Done")
		return(Status);
	}
	*RealTimeDifference = difftime(New, Old);
	Status = (*RealTimeDifference < 0.0);
	if (Status == ERROR) {
		WRITE_ERROR(GetElapsedTime, "Loss of causality")
/**/	DEBUG_WRITE_LEAVING(GetElapsedTime, "Done")
		return(Status);
	}
	*SystemTimeDifference = NewTimeStamp.SecondsSinceProcessBirth
		- OldTimeStamp.SecondsSinceProcessBirth;
	if (*SystemTimeDifference < 0.0) {
		WRITE_ERROR(GetElapsedTime, "Wrap-around of the system clock occurred")
/**/	DEBUG_WRITE_LEAVING(GetElapsedTime, "Done")
		return(Status);
	}
/**/DEBUG_WRITE_LEAVING(GetElapsedTime, "Done")
	return(Status);
} /* end GetElapsedTime */

