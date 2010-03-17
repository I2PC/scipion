/*****************************************************************************
 *	System includes
 ****************************************************************************/
#include	<math.h>

/*****************************************************************************
 *	Toolbox defines
 ****************************************************************************/
#include	"configs.h"
#include	"debug.h"

/*****************************************************************************
 *	Toolbox includes
 ****************************************************************************/
#include	"positivepower.h"

/*****************************************************************************
 *	Conditional includes
 ****************************************************************************/
#ifdef		DEBUG
#include	<limits.h>
#include	<stdio.h>
#include	<string.h>
#include	"error.h"
#include	"messagedisplay.h"
#endif

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
extern double	PositiveIntPower
				(
					double	Argument,			/* argument */
					long	Exponent			/* exponent */
				)

/* recursive computation of a positive integer power of the argument */
/* negative arguments are allowed */

{ /* begin PositiveIntPower */

	double	PseudoRoot;
#ifdef DEBUG
/**/int		Status = !ERROR;
#endif

/**/DEBUG_CHECK_RANGE_LONG(PositiveIntPower, Exponent, 0L, LONG_MAX, Status,
/**/	"Invalid exponent (should be strictly positive)")
/**/DEBUG_WRITE_ENTERING(PositiveIntPower,
/**/	"About to raise a (double)argument to a positive (int)exponent")
#ifdef DEBUG
/**/if (Status == ERROR) {
/**/	WRITE_ERROR(PositiveIntPower, "Unable to perform PositiveIntPower")
/**/	DEBUG_WRITE_LEAVING(PositiveIntPower, "Done")
/**/	return(HUGE_VAL);
/**/}
#endif

	if ((Exponent & 1L) != 0L) {
		if (Exponent == 1L) {
/**/		DEBUG_WRITE_LEAVING(PositiveIntPower, "Done")
			return(Argument);
		}
		if (Exponent == 3L) {
/**/		DEBUG_WRITE_LEAVING(PositiveIntPower, "Done")
			return(Argument * Argument * Argument);
		}
		if (Exponent == 5L) {
			PseudoRoot = Argument * Argument;
/**/		DEBUG_WRITE_LEAVING(PositiveIntPower, "Done")
			return(Argument * PseudoRoot * PseudoRoot);
		}
		if (Exponent == 7L) {
			PseudoRoot = Argument * Argument * Argument;
/**/		DEBUG_WRITE_LEAVING(PositiveIntPower, "Done")
			return(Argument * PseudoRoot * PseudoRoot);
		}
		if (Exponent <= POSITIVE_POWER_TRADE_OFF) {
			PseudoRoot = PositiveIntPower(Argument, Exponent / 2L);
/**/		DEBUG_WRITE_LEAVING(PositiveIntPower, "Done")
			return(Argument * PseudoRoot * PseudoRoot);
		}
/**/	DEBUG_WRITE_LEAVING(PositiveIntPower, "Done")
		return(pow(Argument, (double)Exponent));
	}
	else {
		if (Exponent == 0L) {
			if (Argument == 0.0) {
/**/			DEBUG_WRITE_LEAVING(PositiveIntPower, "Done")
				return(0.5);
			}
			else {
/**/			DEBUG_WRITE_LEAVING(PositiveIntPower, "Done")
				return(1.0);
			}
		}
		if (Exponent == 2L) {
/**/		DEBUG_WRITE_LEAVING(PositiveIntPower, "Done")
			return(Argument * Argument);
		}
		if (Exponent == 4L) {
			Argument *= Argument;
/**/		DEBUG_WRITE_LEAVING(PositiveIntPower, "Done")
			return(Argument * Argument);
		}
		if (Exponent == 6L) {
			Argument *= Argument * Argument;
/**/		DEBUG_WRITE_LEAVING(PositiveIntPower, "Done")
			return(Argument * Argument);
		}
		if (Exponent == 8L) {
			Argument *= Argument;
			Argument *= Argument;
/**/		DEBUG_WRITE_LEAVING(PositiveIntPower, "Done")
			return(Argument * Argument);
		}
		if (Exponent <= POSITIVE_POWER_TRADE_OFF) {
			PseudoRoot = PositiveIntPower(Argument, Exponent / 2L);
/**/		DEBUG_WRITE_LEAVING(PositiveIntPower, "Done")
			return(PseudoRoot * PseudoRoot);
		}
/**/	DEBUG_WRITE_LEAVING(PositiveIntPower, "Done")
		return(pow(Argument, (double)Exponent));
	}
} /* end PositiveIntPower */

