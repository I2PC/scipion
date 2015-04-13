/*****************************************************************************
 *	System includes
 ****************************************************************************/
#include	<math.h>
#include	<stdio.h>

/*****************************************************************************
 *	Toolbox defines
 ****************************************************************************/
#include	"configs.h"
#include	"debug.h"
#include	"error.h"

/*****************************************************************************
 *	Toolbox includes
 ****************************************************************************/
#include	"kernel.h"
#include	"kerneldiff2.h"
#include	"messagedisplay.h"

/*****************************************************************************
 *	Conditional includes
 ****************************************************************************/
#ifdef		DEBUG
#include	<limits.h>
#include	<string.h>
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
extern int		BlipDiff2
				(
					long	Degree,				/* degree */
					double	Argument,			/* input */
					double	*Result				/* output */
				)

/* computes the second derivative of a Blu interpolant function */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin BlipDiff2 */

	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(BlipDiff2, Result, Status,
/**/	"No output data")
/**/DEBUG_RETURN_ON_ERROR(BlipDiff2, Status)
/**/DEBUG_WRITE_ENTERING(BlipDiff2,
/**/	"About to compute the second derivative of a generic blip")

	switch (Degree) {
		case 0L:
			*Result = Blip00Diff2(Argument);
			break;
		case 1L:
			*Result = Blip01Diff2(Argument);
			break;
		case 3L:
			*Result = Blip03Diff2(Argument);
			break;
		default:
			Status = ERROR;
			WRITE_ERROR(BlipDiff2, "Invalid blip degree")
			break;
	}
/**/DEBUG_WRITE_LEAVING(BlipDiff2, "Done")
	return(Status);
} /* end BlipDiff2 */

/*--------------------------------------------------------------------------*/
extern double	Blip00Diff2
				(
					double	Argument			/* input */
				)

/* returns the value of the second derivative of a Blu interpolant function
	of degree 0 (order 1) evaluated at Argument */

{ /* begin Blip00Diff2 */

/**/DEBUG_WRITE_ENTERING(Blip00Diff2,
/**/	"About to compute the second derivative of a blip of degree 0")

	if (fabs(Argument) == 0.5) {
		WRITE_ERROR(Blip00Diff2, "Undefined value at |Argument| == 0.5")
/**/	DEBUG_WRITE_LEAVING(Blip00Diff2, "Done")
		return(HUGE_VAL);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Blip00Diff2, "Done")
		return(0.0);
	}
} /* end Blip00Diff2 */

/*--------------------------------------------------------------------------*/
extern double	Blip01Diff2
				(
					double	Argument			/* input */
				)

/* returns the value of the second derivative of a Blu interpolant function
	of degree 1 (order 2) evaluated at Argument */

{ /* begin Blip01Diff2 */

/**/DEBUG_WRITE_ENTERING(Blip01Diff2,
/**/	"About to compute the second derivative of a blip of degree 1")

	if (Argument == 0.0) {
		WRITE_ERROR(Blip01Diff2, "Undefined value at Argument == 0.0")
/**/	DEBUG_WRITE_LEAVING(Blip01Diff2, "Done")
		return(-HUGE_VAL);
	}
	else if (fabs(Argument) == 1.0) {
		WRITE_ERROR(Blip01Diff2, "Undefined value at |Argument| == 1.0")
/**/	DEBUG_WRITE_LEAVING(Blip01Diff2, "Done")
		return(HUGE_VAL);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Blip01Diff2, "Done")
		return(0.0);
	}
} /* end Blip01Diff2 */

/*--------------------------------------------------------------------------*/
extern double	Blip03Diff2
				(
					double	Argument			/* input */
				)

/* returns the value of the second derivative of a Blu interpolant function
	of degree 3 (order 4) evaluated at Argument */

{ /* begin Blip03Diff2 */

/**/DEBUG_WRITE_ENTERING(Blip03Diff2,
/**/	"About to compute the second derivative of a blip of degree 3")

	Argument = fabs(Argument);
	if (Argument == 0.0) {
		WRITE_ERROR(Blip03Diff2, "Undefined value at |Argument| == 0.0")
/**/	DEBUG_WRITE_LEAVING(Blip03Diff2, "Done")
		return(HUGE_VAL);
	}
	else if (Argument < 1.0) {
/**/	DEBUG_WRITE_LEAVING(Blip03Diff2, "Done")
		return(Argument * 3.0 - 2.0);
	}
	else if (Argument == 1.0) {
		WRITE_ERROR(Blip03Diff2, "Undefined value at |Argument| == 1.0")
/**/	DEBUG_WRITE_LEAVING(Blip03Diff2, "Done")
		return(-HUGE_VAL);
	}
	else if (Argument < 2.0) {
/**/	DEBUG_WRITE_LEAVING(Blip03Diff2, "Done")
		return(2.0 - Argument);
	}
	else if (Argument == 2.0) {
		WRITE_ERROR(Blip03Diff2, "Undefined value at |Argument| == 2.0")
/**/	DEBUG_WRITE_LEAVING(Blip03Diff2, "Done")
		return(HUGE_VAL);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Blip03Diff2, "Done")
		return(0.0);
	}
} /* end Blip03Diff2 */

/*--------------------------------------------------------------------------*/
extern int		BsplineDiff2
				(
					long	Degree,				/* degree */
					double	Argument,			/* input */
					double	*Result				/* output */
				)

/* computes the second derivative of a Basic spline function */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin BsplineDiff2 */

	double	Left, Center, Right;
	int		Status = !ERROR;

/**/DEBUG_CHECK_RANGE_LONG(BsplineDiff2, Degree, 0L, LONG_MAX, Status,
/**/	"Invalid degree (should be positive)")
/**/DEBUG_CHECK_NULL_POINTER(BsplineDiff2, Result, Status,
/**/	"No output data")
/**/DEBUG_RETURN_ON_ERROR(BsplineDiff2, Status)
/**/DEBUG_WRITE_ENTERING(BsplineDiff2,
/**/	"About to compute the second derivative of a generic B-spline")

	switch (Degree) {
		case 0L:
			*Result = Bspline00Diff2(Argument);
			break;
		case 1L:
			*Result = Bspline01Diff2(Argument);
			break;
		case 2L:
			*Result = Bspline02Diff2(Argument);
			break;
		case 3L:
			*Result = Bspline03Diff2(Argument);
			break;
		case 4L:
			*Result = Bspline04Diff2(Argument);
			break;
		case 5L:
			*Result = Bspline05Diff2(Argument);
			break;
		case 6L:
			*Result = Bspline06Diff2(Argument);
			break;
		case 7L:
			*Result = Bspline07Diff2(Argument);
			break;
		case 8L:
			*Result = Bspline08Diff2(Argument);
			break;
		case 9L:
			*Result = Bspline09Diff2(Argument);
			break;
		case 10L:
			*Result = Bspline10Diff2(Argument);
			break;
		case 11L:
			*Result = Bspline11Diff2(Argument);
			break;
		default:
			Status = Bspline(Degree - 2L, Argument - 1.0, &Left);
			if (Status == ERROR) {
/**/			DEBUG_WRITE_LEAVING(BsplineDiff2, "Done")
				return(Status);
			}
			Status = Bspline(Degree - 2L, Argument, &Center);
			if (Status == ERROR) {
/**/			DEBUG_WRITE_LEAVING(BsplineDiff2, "Done")
				return(Status);
			}
			Status = Bspline(Degree - 2L, Argument + 1.0, &Right);
			if (Status == ERROR) {
/**/			DEBUG_WRITE_LEAVING(BsplineDiff2, "Done")
				return(Status);
			}
			*Result = Right - 2.0 * Center + Left;
			break;
	}
/**/DEBUG_WRITE_LEAVING(BsplineDiff2, "Done")
	return(Status);
} /* end BsplineDiff2 */

/*--------------------------------------------------------------------------*/
extern double	Bspline00Diff2
				(
					double	Argument			/* input */
				)

/* returns the value of the second derivative of a Basic spline function
	of degree 0 (order 1) evaluated at Argument */

{ /* begin Bspline00Diff2 */

/**/DEBUG_WRITE_ENTERING(Bspline00Diff2,
/**/	"About to compute the second derivative of a B-spline of degree 0")

	if (fabs(Argument) == 0.5) {
		WRITE_ERROR(Bspline00Diff2, "Undefined value at |Argument| == 0.5")
/**/	DEBUG_WRITE_LEAVING(Bspline00Diff2, "Done")
		return(HUGE_VAL);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Bspline00Diff2, "Done")
		return(0.0);
	}
} /* end Bspline00Diff2 */

/*--------------------------------------------------------------------------*/
extern double	Bspline01Diff2
				(
					double	Argument			/* input */
				)

/* returns the value of the second derivative of a Basic spline function
	of degree 1 (order 2) evaluated at Argument */

{ /* begin Bspline01Diff2 */

/**/DEBUG_WRITE_ENTERING(Bspline01Diff2,
/**/	"About to compute the second derivative of a B-spline of degree 1")

	if (Argument == 0.0) {
		WRITE_ERROR(Bspline01Diff2, "Undefined value at Argument == 0.0")
/**/	DEBUG_WRITE_LEAVING(Bspline01Diff2, "Done")
		return(-HUGE_VAL);
	}
	else if (fabs(Argument) == 1.0) {
		WRITE_ERROR(Bspline01Diff2, "Undefined value at |Argument| == 1.0")
/**/	DEBUG_WRITE_LEAVING(Bspline01Diff2, "Done")
		return(HUGE_VAL);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Bspline01Diff2, "Done")
		return(0.0);
	}
} /* end Bspline01Diff2 */

/*--------------------------------------------------------------------------*/
extern double	Bspline02Diff2
				(
					double	Argument			/* input */
				)

/* returns the value of the second derivative of a Basic spline function
	of degree 2 (order 3) evaluated at Argument */

{ /* begin Bspline02Diff2 */

/**/DEBUG_WRITE_ENTERING(Bspline02Diff2,
/**/	"About to compute the second derivative of a B-spline of degree 2")

	Argument = fabs(Argument);
	if (Argument < 0.5) {
/**/	DEBUG_WRITE_LEAVING(Bspline02Diff2, "Done")
		return(-2.0);
	}
	else if (Argument == 0.5) {
/**/	DEBUG_WRITE_LEAVING(Bspline02Diff2, "Done")
		return(-1.0 / 2.0);
	}
	else if (Argument < 1.5) {
/**/	DEBUG_WRITE_LEAVING(Bspline02Diff2, "Done")
		return(1.0);
	}
	else if (Argument == 1.5) {
/**/	DEBUG_WRITE_LEAVING(Bspline02Diff2, "Done")
		return(1.0 / 2.0);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Bspline02Diff2, "Done")
		return(0.0);
	}
} /* end Bspline02Diff2 */

/*--------------------------------------------------------------------------*/
extern double	Bspline03Diff2
				(
					double	Argument			/* input */
				)

/* returns the value of the second derivative of a Basic spline function
	of degree 3 (order 4) evaluated at Argument */

{ /* begin Bspline03Diff2 */

/**/DEBUG_WRITE_ENTERING(Bspline03Diff2,
/**/	"About to compute the second derivative of a B-spline of degree 3")

	Argument = fabs(Argument);
	if (Argument < 1.0) {
/**/	DEBUG_WRITE_LEAVING(Bspline03Diff2, "Done")
		return(Argument * 3.0 - 2.0);
	}
	else if (Argument < 2.0) {
/**/	DEBUG_WRITE_LEAVING(Bspline03Diff2, "Done")
		return(2.0 - Argument);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Bspline03Diff2, "Done")
		return(0.0);
	}
} /* end Bspline03Diff2 */

/*--------------------------------------------------------------------------*/
extern double	Bspline04Diff2
				(
					double	Argument			/* input */
				)

/* returns the value of the second derivative of a Basic spline function
	of degree 4 (order 5) evaluated at Argument */

{ /* begin Bspline04Diff2 */

/**/DEBUG_WRITE_ENTERING(Bspline04Diff2,
/**/	"About to compute the second derivative of a B-spline of degree 4")

	Argument = fabs(Argument);
	if (Argument < 0.5) {
/**/	DEBUG_WRITE_LEAVING(Bspline04Diff2, "Done")
		return(Argument * Argument * 3.0 - 5.0 / 4.0);
	}
	else if (Argument < 1.5) {
/**/	DEBUG_WRITE_LEAVING(Bspline04Diff2, "Done")
		return(Argument * (5.0 - Argument * 2.0) - 5.0 / 2.0);
	}
	else if (Argument < 2.5) {
		Argument -= 5.0 / 2.0;
/**/	DEBUG_WRITE_LEAVING(Bspline04Diff2, "Done")
		return(Argument * Argument * (1.0 / 2.0));
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Bspline04Diff2, "Done")
		return(0.0);
	}
} /* end Bspline04Diff2 */

/*--------------------------------------------------------------------------*/
extern double	Bspline05Diff2
				(
					double	Argument			/* input */
				)

/* returns the value of the second derivative of a Basic spline function
	of degree 5 (order 6) evaluated at Argument */

{ /* begin Bspline05Diff2 */

/**/DEBUG_WRITE_ENTERING(Bspline05Diff2,
/**/	"About to compute the second derivative of a B-spline of degree 5")

	Argument = fabs(Argument);
	if (Argument < 1.0) {
/**/	DEBUG_WRITE_LEAVING(Bspline05Diff2, "Done")
		return(Argument * Argument * (3.0 - Argument * (5.0 / 3.0)) - 1.0);
	}
	else if (Argument < 2.0) {
/**/	DEBUG_WRITE_LEAVING(Bspline05Diff2, "Done")
		return(Argument * (Argument * (Argument * (5.0 / 6.0) - 9.0 / 2.0)
			+ 15.0 / 2.0) - 7.0 / 2.0);
	}
	else if (Argument < 3.0) {
		Argument -= 3.0;
/**/	DEBUG_WRITE_LEAVING(Bspline05Diff2, "Done")
		return(Argument * Argument * Argument * (-1.0 / 6.0));
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Bspline05Diff2, "Done")
		return(0.0);
	}
} /* end Bspline05Diff2 */

/*--------------------------------------------------------------------------*/
extern double	Bspline06Diff2
				(
					double	Argument			/* input */
				)

/* returns the value of the second derivative of a Basic spline function
	of degree 6 (order 7) evaluated at Argument */

{ /* begin Bspline06Diff2 */

/**/DEBUG_WRITE_ENTERING(Bspline06Diff2,
/**/	"About to compute the second derivative of a B-spline of degree 6")

	Argument = fabs(Argument);
	if (Argument < 0.5) {
		Argument *= Argument;
/**/	DEBUG_WRITE_LEAVING(Bspline06Diff2, "Done")
		return(Argument * (7.0 / 4.0 - Argument * (5.0 / 6.0)) - 77.0 / 96.0);
	}
	else if (Argument < 1.5) {
/**/	DEBUG_WRITE_LEAVING(Bspline06Diff2, "Done")
		return(Argument * (Argument * (Argument * (Argument * 5.0 / 8.0
			- 35.0 / 12.0) + 63.0 / 16.0) - 35.0 / 48.0) - 91.0 / 128.0);
	}
	else if (Argument < 2.5) {
/**/	DEBUG_WRITE_LEAVING(Bspline06Diff2, "Done")
		return(Argument * (Argument * (Argument * (7.0 / 3.0 - Argument * (1.0 / 4.0))
			- 63.0 / 8.0) + 133.0 / 12.0) - 329.0 / 64.0);
	}
	else if (Argument < 3.5) {
		Argument -= 7.0 / 2.0;
		Argument *= Argument;
/**/	DEBUG_WRITE_LEAVING(Bspline06Diff2, "Done")
		return(Argument * Argument * (1.0 / 24.0));
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Bspline06Diff2, "Done")
		return(0.0);
	}
} /* end Bspline06Diff2 */

/*--------------------------------------------------------------------------*/
extern double	Bspline07Diff2
				(
					double	Argument			/* input */
				)

/* returns the value of the second derivative of a Basic spline function
	of degree 7 (order 8) evaluated at Argument */

{ /* begin Bspline07Diff2 */

	double	a2;

/**/DEBUG_WRITE_ENTERING(Bspline07Diff2,
/**/	"About to compute the second derivative of a B-spline of degree 7")

	Argument = fabs(Argument);
	if (Argument < 1.0) {
		a2 = Argument * Argument;
/**/	DEBUG_WRITE_LEAVING(Bspline07Diff2, "Done")
		return(a2 * (a2 * (Argument * (7.0 / 24.0) - 5.0 / 6.0) + 4.0 / 3.0)
			- 2.0 / 3.0);
	}
	else if (Argument < 2.0) {
/**/	DEBUG_WRITE_LEAVING(Bspline07Diff2, "Done")
		return(Argument * (Argument * (Argument * (Argument * (3.0 / 2.0
			- Argument * (7.0 / 40.0)) - 14.0 / 3.0) + 6.0) - 7.0 / 3.0) - 1.0 / 5.0);
	}
	else if (Argument < 3.0) {
/**/	DEBUG_WRITE_LEAVING(Bspline07Diff2, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (7.0 / 120.0)
			- 5.0 / 6.0) + 14.0 / 3.0) - 38.0 / 3.0) + 49.0 / 3.0) - 23.0 / 3.0);
	}
	else if (Argument < 4.0) {
		Argument -= 4.0;
		a2 = Argument * Argument;
/**/	DEBUG_WRITE_LEAVING(Bspline07Diff2, "Done")
		return(Argument * a2 * a2 * (-1.0 / 120.0));
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Bspline07Diff2, "Done")
		return(0.0);
	}
} /* end Bspline07Diff2 */

/*--------------------------------------------------------------------------*/
extern double	Bspline08Diff2
				(
					double	Argument			/* input */
				)

/* returns the value of the second derivative of a Basic spline function
	of degree 8 (order 9) evaluated at Argument */

{ /* begin Bspline08Diff2 */

/**/DEBUG_WRITE_ENTERING(Bspline08Diff2,
/**/	"About to compute the second derivative of a B-spline of degree 8")

	Argument = fabs(Argument);
	if (Argument < 0.5) {
		Argument *= Argument;
/**/	DEBUG_WRITE_LEAVING(Bspline08Diff2, "Done")
		return(Argument * (Argument * (Argument * (7.0 / 72.0) - 15.0 / 32.0)
			+ 129.0 / 128.0) - 289.0 / 512.0);
	}
	else if (Argument < 1.5) {
/**/	DEBUG_WRITE_LEAVING(Bspline08Diff2, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (21.0 / 40.0
			- Argument * (7.0 / 90.0)) - 9.0 / 8.0) + 7.0 / 16.0) + 27.0 / 32.0)
			+ 21.0 / 640.0) - 363.0 / 640.0);
	}
	else if (Argument < 2.5) {
/**/	DEBUG_WRITE_LEAVING(Bspline08Diff2, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (7.0 / 180.0) - 21.0 / 40.0) + 45.0 / 16.0) - 119.0 / 16.0)
			+ 621.0 / 64.0) - 3381.0 / 640.0) + 195.0 / 256.0);
	}
	else if (Argument < 3.5) {
/**/	DEBUG_WRITE_LEAVING(Bspline08Diff2, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (9.0 / 40.0
			- Argument * (1.0 / 90.0)) - 15.0 / 8.0) + 131.0 / 16.0) - 627.0 / 32.0)
			+ 15369.0 / 640.0) - 1465.0 / 128.0);
	}
	else if (Argument < 4.5) {
		Argument -= 9.0 / 2.0;
		Argument *= Argument * Argument;
/**/	DEBUG_WRITE_LEAVING(Bspline08Diff2, "Done")
		return(Argument * Argument * (1.0 / 720.0));
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Bspline08Diff2, "Done")
		return(0.0);
	}
} /* end Bspline08Diff2 */

/*--------------------------------------------------------------------------*/
extern double	Bspline09Diff2
				(
					double	Argument			/* input */
				)

/* returns the value of the second derivative of a Basic spline function
	of degree 9 (order 10) evaluated at Argument */

{ /* begin Bspline09Diff2 */

	double	a2;

/**/DEBUG_WRITE_ENTERING(Bspline09Diff2,
/**/	"About to compute the second derivative of a B-spline of degree 9")

	Argument = fabs(Argument);
	if (Argument < 1.0) {
		a2 = Argument * Argument;
/**/	DEBUG_WRITE_LEAVING(Bspline09Diff2, "Done")
		return(a2 * (a2 * (a2 * (7.0 / 72.0 - Argument * (1.0 / 40.0))
			- 25.0 / 72.0) + 19.0 / 24.0) - 35.0 / 72.0);
	}
	else if (Argument < 2.0) {
/**/	DEBUG_WRITE_LEAVING(Bspline09Diff2, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (Argument * (1.0 / 60.0) - 7.0 / 36.0) + 7.0 / 8.0) - 65.0 / 36.0)
			+ 35.0 / 24.0) - 1.0 / 12.0) + 7.0 / 24.0) - 19.0 / 36.0);
	}
	else if (Argument < 3.0) {
/**/	DEBUG_WRITE_LEAVING(Bspline09Diff2, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (5.0 / 36.0 - Argument * (1.0 / 140.0)) - 9.0 / 8.0) + 175.0 / 36.0)
			- 95.0 / 8.0) + 191.0 / 12.0) - 83.0 / 8.0) + 635.0 / 252.0);
	}
	else if (Argument < 4.0) {
/**/	DEBUG_WRITE_LEAVING(Bspline09Diff2, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (Argument * (1.0 / 560.0) - 7.0 / 144.0) + 9.0 / 16.0) - 515.0 / 144.0)
			+ 215.0 / 16.0) - 1423.0 / 48.0) + 563.0 / 16.0) - 2449.0 / 144.0);
	}
	else if (Argument < 5.0) {
		Argument -= 5.0;
		a2 = Argument * Argument * Argument;
/**/	DEBUG_WRITE_LEAVING(Bspline09Diff2, "Done")
		return(Argument * a2 * a2 * (-1.0 / 5040.0));
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Bspline09Diff2, "Done")
		return(0.0);
	}
} /* end Bspline09Diff2 */

/*--------------------------------------------------------------------------*/
extern double	Bspline10Diff2
				(
					double	Argument			/* input */
				)

/* returns the value of the second derivative of a Basic spline function
	of degree 10 (order 11) evaluated at Argument */

{ /* begin Bspline10Diff2 */

/**/DEBUG_WRITE_ENTERING(Bspline10Diff2,
/**/	"About to compute the second derivative of a B-spline of degree 10")

	Argument = fabs(Argument);
	if (Argument < 0.5) {
		Argument *= Argument;
/**/	DEBUG_WRITE_LEAVING(Bspline10Diff2, "Done")
		return(Argument * (Argument * (Argument * (77.0 / 1440.0 - Argument
			* (1.0 / 160.0)) - 583.0 / 2304.0) + 14597.0 / 23040.0)
			- 156409.0 / 368640.0);
	}
	else if (Argument < 1.5) {
/**/	DEBUG_WRITE_LEAVING(Bspline10Diff2, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (Argument * (Argument * (1.0 / 192.0) - 11.0 / 240.0) + 77.0 / 576.0)
			- 77.0 / 960.0) - 935.0 / 4608.0) - 77.0 / 3840.0) + 5885.0 / 9216.0)
			- 11.0 / 15360.0) - 62557.0 / 147456.0);
	}
	else if (Argument < 2.5) {
/**/	DEBUG_WRITE_LEAVING(Bspline10Diff2, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (Argument * (11.0 / 210.0 - Argument * (1.0 / 336.0)) - 55.0 / 144.0)
			+ 22.0 / 15.0) - 3575.0 / 1152.0) + 1661.0 / 480.0) - 4543.0 / 2304.0)
			+ 7513.0 / 6720.0) - 163603.0 / 258048.0);
	}
	else if (Argument < 3.5) {
/**/	DEBUG_WRITE_LEAVING(Bspline10Diff2, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (Argument * (Argument * (1.0 / 896.0) - 33.0 / 1120.0)
			+ 385.0 / 1152.0) - 1353.0 / 640.0) + 74525.0 / 9216.0) - 48433.0 / 2560.0)
			+ 479281.0 / 18432.0) - 1352153.0 / 71680.0) + 1654543.0 / 294912.0);
	}
	else if (Argument < 4.5) {
/**/	DEBUG_WRITE_LEAVING(Bspline10Diff2, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (Argument * (11.0 / 1260.0 - Argument * (1.0 / 4032.0)) - 77.0 / 576.0)
			+ 209.0 / 180.0) - 28765.0 / 4608.0) + 61061.0 / 2880.0) - 407429.0 / 9216.0)
			+ 2070343.0 / 40320.0) - 3702215.0 / 147456.0);
	}
	else if (Argument < 5.5) {
		Argument -= 11.0 / 2.0;
		Argument *= Argument;
		Argument *= Argument;
/**/	DEBUG_WRITE_LEAVING(Bspline10Diff2, "Done")
		return(Argument * Argument * (1.0 / 40320.0));
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Bspline10Diff2, "Done")
		return(0.0);
	}
} /* end Bspline10Diff2 */

/*--------------------------------------------------------------------------*/
extern double	Bspline11Diff2
				(
					double	Argument			/* input */
				)

/* returns the value of the second derivative of a Basic spline function
	of degree 11 (order 12) evaluated at Argument */

{ /* begin Bspline11Diff2 */

	double	a2;

/**/DEBUG_WRITE_ENTERING(Bspline11Diff2,
/**/	"About to compute the second derivative of a B-spline of degree 11")

	Argument = fabs(Argument);
	if (Argument < 1.0) {
		a2 = Argument * Argument;
/**/	DEBUG_WRITE_LEAVING(Bspline11Diff2, "Done")
		return(a2 * (a2 * (a2 * (a2 * (Argument * (11.0 / 8640.0) - 1.0 / 160.0)
			+ 7.0 / 180.0) - 23.0 / 120.0) + 31.0 / 60.0) - 809.0 / 2160.0);
	}
	else if (Argument < 2.0) {
/**/	DEBUG_WRITE_LEAVING(Bspline11Diff2, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (Argument * (Argument * (3.0 / 224.0 - Argument * (11.0 / 12096.0))
			- 11.0 / 140.0) + 2.0 / 9.0) - 11.0 / 40.0) + 1.0 / 12.0) - 11.0 / 60.0)
			+ 25.0 / 42.0) - 11.0 / 560.0) - 563.0 / 1512.0);
	}
	else if (Argument < 3.0) {
/**/	DEBUG_WRITE_LEAVING(Bspline11Diff2, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (Argument * (Argument * (Argument * (11.0 / 24192.0) - 5.0 / 448.0)
			+ 33.0 / 280.0) - 25.0 / 36.0) + 99.0 / 40.0) - 65.0 / 12.0) + 143.0 / 20.0)
			- 239.0 / 42.0) + 1749.0 / 560.0) - 1619.0 / 1512.0);
	}
	else if (Argument < 4.0) {
/**/	DEBUG_WRITE_LEAVING(Bspline11Diff2, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (Argument * (Argument * (1.0 / 192.0 - Argument * (11.0 / 72576.0))
			- 11.0 / 140.0) + 49.0 / 72.0) - 297.0 / 80.0) + 631.0 / 48.0)
			- 1199.0 / 40.0) + 1009.0 / 24.0) - 36597.0 / 1120.0) + 9385.0 / 864.0);
	}
	else if (Argument < 5.0) {
/**/	DEBUG_WRITE_LEAVING(Bspline11Diff2, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (Argument * (Argument * (Argument * (11.0 / 362880.0) - 3.0 / 2240.0)
			+ 11.0 / 420.0) - 107.0 / 360.0) + 517.0 / 240.0) - 2477.0 / 240.0)
			+ 11737.0 / 360.0) - 54797.0 / 840.0) + 250657.0 / 3360.0) - 1113317.0 / 30240.0);
	}
	else if (Argument < 6.0) {
		Argument -= 6.0;
		a2 = Argument * Argument;
		a2 *= a2;
		a2 *= a2;
/**/	DEBUG_WRITE_LEAVING(Bspline11Diff2, "Done")
		return(Argument * a2 * (-1.0 / 362880.0));
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Bspline11Diff2, "Done")
		return(0.0);
	}
} /* end Bspline11Diff2 */

/*--------------------------------------------------------------------------*/
extern double	DodgsonDiff2
				(
					double	Argument			/* input */
				)

/* returns the value of the second derivative of a Dodgson kernel
	(order 2) evaluated at Argument */

{ /* begin DodgsonDiff2 */

/**/DEBUG_WRITE_ENTERING(DodgsonDiff2,
/**/	"About to compute the second derivative of a Dodgson kernel")

	Argument = fabs(Argument);
	if (Argument < 0.5) {
/**/	DEBUG_WRITE_LEAVING(DodgsonDiff2, "Done")
		return(-4.0);
	}
	else if (Argument == 0.5) {
		WRITE_ERROR(DodgsonDiff2, "Undefined value at Argument == 0.5")
/**/	DEBUG_WRITE_LEAVING(DodgsonDiff2, "Done")
		return(-HUGE_VAL);
	}
	else if (Argument < 1.5) {
/**/	DEBUG_WRITE_LEAVING(DodgsonDiff2, "Done")
		return(2.0);
	}
	else if (fabs(Argument) == 1.5) {
		WRITE_ERROR(DodgsonDiff2, "Undefined value at |Argument| == 1.5")
/**/	DEBUG_WRITE_LEAVING(DodgsonDiff2, "Done")
		return(HUGE_VAL);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(DodgsonDiff2, "Done")
		return(0.0);
	}
} /* end DodgsonDiff2 */

/*--------------------------------------------------------------------------*/
extern double	German04Diff2
				(
					double	Argument			/* input */
				)

/* returns the value of the second derivative of a German function
	of degree 4 (order 5) evaluated at Argument */

{ /* begin German04Diff2 */

/**/DEBUG_WRITE_ENTERING(German04Diff2,
/**/	"About to compute the second derivative of a German of degree 4")

	Argument = fabs(Argument);
	if (Argument < 1.0) {
/**/	DEBUG_WRITE_LEAVING(German04Diff2, "Done")
		return(Argument * (Argument * (1.0 / 2.0) + 185.0 / 24.0) - 335.0 / 72.0);
	}
	else if (Argument == 1.0) {
/**/	DEBUG_WRITE_LEAVING(German04Diff2, "Done")
		return(421.0 / 144.0);
	}
	else if (Argument < 2.0) {
/**/	DEBUG_WRITE_LEAVING(German04Diff2, "Done")
		return(Argument * (-29.0 / 24.0 - Argument * (3.0 / 4.0)) + 17.0 / 4.0);
	}
	else if (Argument == 2.0) {
/**/	DEBUG_WRITE_LEAVING(German04Diff2, "Done")
		return(-49.0 / 72.0);
	}
	else if (Argument < 3.0) {
/**/	DEBUG_WRITE_LEAVING(German04Diff2, "Done")
		return(Argument * (Argument * (1.0 / 4.0) - 7.0 / 8.0) + 5.0 / 9.0);
	}
	else if (Argument == 3.0) {
/**/	DEBUG_WRITE_LEAVING(German04Diff2, "Done")
		return(11.0 / 144.0);
	}
	else if (Argument < 4.0) {
/**/	DEBUG_WRITE_LEAVING(German04Diff2, "Done")
		return(Argument * (1.0 / 24.0) - 11.0 / 72.0);
	}
	else if (Argument == 4.0) {
/**/	DEBUG_WRITE_LEAVING(German04Diff2, "Done")
		return(1.0 / 144.0);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(German04Diff2, "Done")
		return(0.0);
	}
} /* end German04Diff2 */

/*--------------------------------------------------------------------------*/
extern double	KeysOptimalDiff2
				(
					double	Argument			/* input */
				)

/* returns the value of the second derivative of the cubic Keys optimal kernel
	(order 3) evaluated at Argument */

{ /* begin KeysOptimalDiff2 */

/**/DEBUG_WRITE_ENTERING(KeysOptimalDiff2,
/**/	"About to compute the second derivative of a Keys kernel")

	Argument = fabs(Argument);
	if (Argument < 1.0) {
/**/	DEBUG_WRITE_LEAVING(KeysOptimalDiff2, "Done")
		return(Argument * 9.0 - 5.0);
	}
	else if (Argument == 1.0) {
/**/	DEBUG_WRITE_LEAVING(KeysOptimalDiff2, "Done")
		return(3.0);
	}
	else if (Argument < 2.0) {
/**/	DEBUG_WRITE_LEAVING(KeysOptimalDiff2, "Done")
		return(5.0 - Argument * 3.0);
	}
	else if (Argument == 2.0) {
/**/	DEBUG_WRITE_LEAVING(KeysOptimalDiff2, "Done")
		return(-1.0 / 2.0);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(KeysOptimalDiff2, "Done")
		return(0.0);
	}
} /* end KeysOptimalDiff2 */

/*--------------------------------------------------------------------------*/
extern int		OmomsDiff2
				(
					long	Degree,				/* degree */
					double	Argument,			/* input */
					double	*Result				/* output */
				)

/* computes the second derivative of a Blu optimum function */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin OmomsDiff2 */

	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(OmomsDiff2, Result, Status,
/**/	"No output data")
/**/DEBUG_RETURN_ON_ERROR(OmomsDiff2, Status)
/**/DEBUG_WRITE_ENTERING(OmomsDiff2,
/**/	"About to compute the second derivative of a generic oMoms")

	switch (Degree) {
		case 0L:
			*Result = Omoms00Diff2(Argument);
			break;
		case 1L:
			*Result = Omoms01Diff2(Argument);
			break;
		case 3L:
			*Result = Omoms03Diff2(Argument);
			break;
		default:
			Status = ERROR;
			WRITE_ERROR(OmomsDiff2, "Invalid oMoms degree")
			break;
	}
/**/DEBUG_WRITE_LEAVING(OmomsDiff2, "Done")
	return(Status);
} /* end OmomsDiff2 */

/*--------------------------------------------------------------------------*/
extern double	Omoms00Diff2
				(
					double	Argument			/* input */
				)

/* returns the value of the second derivative of a Blu optimum function
	of degree 0 (order 1) evaluated at Argument */

{ /* begin Omoms00Diff2 */

/**/DEBUG_WRITE_ENTERING(Omoms00Diff2,
/**/	"About to compute the second derivative of an oMoms of degree 0")

	if (fabs(Argument) == 0.5) {
		WRITE_ERROR(Omoms00Diff2, "Undefined value at |Argument| == 0.5")
/**/	DEBUG_WRITE_LEAVING(Omoms00Diff2, "Done")
		return(HUGE_VAL);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Omoms00Diff2, "Done")
		return(0.0);
	}
} /* end Omoms00Diff2 */

/*--------------------------------------------------------------------------*/
extern double	Omoms01Diff2
				(
					double	Argument			/* input */
				)

/* returns the value of the second derivative of a Blu optimum function
	of degree 1 (order 2) evaluated at Argument */

{ /* begin Omoms01Diff2 */

/**/DEBUG_WRITE_ENTERING(Omoms01Diff2,
/**/	"About to compute the second derivative of an oMoms of degree 1")

	if (Argument == 0.0) {
		WRITE_ERROR(Omoms01Diff2, "Undefined value at Argument == 0.0")
/**/	DEBUG_WRITE_LEAVING(Omoms01Diff2, "Done")
		return(-HUGE_VAL);
	}
	else if (fabs(Argument) == 1.0) {
		WRITE_ERROR(Omoms01Diff2, "Undefined value at |Argument| == 1.0")
/**/	DEBUG_WRITE_LEAVING(Omoms01Diff2, "Done")
		return(HUGE_VAL);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Omoms01Diff2, "Done")
		return(0.0);
	}
} /* end Omoms01Diff2 */

/*--------------------------------------------------------------------------*/
extern double	Omoms02Diff2
				(
					double	Argument			/* input */
				)

/* returns the value of the second derivative of a Blu optimum function
	of degree 2 (order 3) evaluated at Argument */

{ /* begin Omoms02Diff2 */

/**/DEBUG_WRITE_ENTERING(Omoms02Diff2,
/**/	"About to compute the second derivative of an oMoms of degree 2")

	Argument = fabs(Argument);
	if (Argument < 0.5) {
/**/	DEBUG_WRITE_LEAVING(Omoms02Diff2, "Done")
		return(-2.0);
	}
	else if (Argument == 0.5) {
		WRITE_ERROR(Omoms02Diff2, "Undefined value at |Argument| == 0.5")
/**/	DEBUG_WRITE_LEAVING(Omoms02Diff2, "Done")
		return(HUGE_VAL);
	}
	else if (Argument < 1.5) {
/**/	DEBUG_WRITE_LEAVING(Omoms02Diff2, "Done")
		return(1.0);
	}
	else if (Argument == 1.5) {
		WRITE_ERROR(Omoms02Diff2, "Undefined value at |Argument| == 1.5")
/**/	DEBUG_WRITE_LEAVING(Omoms02Diff2, "Done")
		return(HUGE_VAL);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Omoms02Diff2, "Done")
		return(0.0);
	}
} /* end Omoms02Diff2 */

/*--------------------------------------------------------------------------*/
extern double	Omoms03Diff2
				(
					double	Argument			/* input */
				)

/* returns the value of the second derivative of a Blu optimum function
	of degree 3 (order 4) evaluated at Argument */

{ /* begin Omoms03Diff2 */

/**/DEBUG_WRITE_ENTERING(Omoms03Diff2,
/**/	"About to compute the second derivative of an oMoms of degree 3")

	Argument = fabs(Argument);
	if (Argument == 0.0) {
		WRITE_ERROR(Omoms03Diff2, "Undefined value at |Argument| == 0.0")
/**/	DEBUG_WRITE_LEAVING(Omoms03Diff2, "Done")
		return(HUGE_VAL);
	}
	else if (Argument < 1.0) {
/**/	DEBUG_WRITE_LEAVING(Omoms03Diff2, "Done")
		return(Argument * 3.0 - 2.0);
	}
	else if (Argument == 1.0) {
		WRITE_ERROR(Omoms03Diff2, "Undefined value at |Argument| == 1.0")
/**/	DEBUG_WRITE_LEAVING(Omoms03Diff2, "Done")
		return(-HUGE_VAL);
	}
	else if (Argument < 2.0) {
/**/	DEBUG_WRITE_LEAVING(Omoms03Diff2, "Done")
		return(2.0 - Argument);
	}
	else if (Argument == 2.0) {
		WRITE_ERROR(Omoms03Diff2, "Undefined value at |Argument| == 2.0")
/**/	DEBUG_WRITE_LEAVING(Omoms03Diff2, "Done")
		return(HUGE_VAL);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Omoms03Diff2, "Done")
		return(0.0);
	}
} /* end Omoms03Diff2 */

/*--------------------------------------------------------------------------*/
extern double	PositiveDiff2
				(
					double	Argument			/* input */
				)

/* returns the value of the second derivative of a Positive function
	of degree 3 (order 1) evaluated at Argument */

{ /* begin PositiveDiff2 */

/**/DEBUG_WRITE_ENTERING(PositiveDiff2,
/**/	"About to compute the second derivative of Positive of degree 3")

	Argument = fabs(Argument);
	if (Argument < 0.5) {
/**/	DEBUG_WRITE_LEAVING(PositiveDiff2, "Done")
		return(Argument * 24.0 - 9.0);
	}
	else if (Argument < 1.0) {
/**/	DEBUG_WRITE_LEAVING(PositiveDiff2, "Done")
		return(3.0);
	}
	else if (Argument < 1.5) {
/**/	DEBUG_WRITE_LEAVING(PositiveDiff2, "Done")
		return(15.0 - Argument * 12.0);
	}
	else if (Argument < 2.0) {
/**/	DEBUG_WRITE_LEAVING(PositiveDiff2, "Done")
		return(Argument * 12.0 - 21.0);
	}
	else if (Argument == 2.0) {
/**/	DEBUG_WRITE_LEAVING(PositiveDiff2, "Done")
		return(3.0 / 2.0);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(PositiveDiff2, "Done")
		return(0.0);
	}
} /* end PositiveDiff2 */

/*--------------------------------------------------------------------------*/
extern double	SincDiff2
				(
					double	Argument			/* input */
				)

/* returns the value of the second derivative of the sinc kernel evaluated at Argument */

{ /* begin SincDiff2 */

/**/DEBUG_WRITE_ENTERING(SincDiff2,
/**/	"About to compute the second derivative of a sinc kernel")

	if (Argument == 0.0) {
		return((-1.0 / 3.0) * PI * PI);
/**/	DEBUG_WRITE_LEAVING(SincDiff2, "Done")
	}
	else {
/**/	DEBUG_WRITE_LEAVING(SincDiff2, "Done")
		return(2.0 * (Sinc(Argument) * (1.0 - (1.0 / 2.0) * PI * PI * Argument * Argument)
			- cos(PI * Argument)) / (Argument * Argument));
	}
} /* end SincDiff2 */

