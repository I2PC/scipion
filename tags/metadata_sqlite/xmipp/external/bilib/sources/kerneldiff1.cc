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
#include	"kerneldiff1.h"
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
extern int		BlipDiff1
				(
					long	Degree,				/* degree */
					double	Argument,			/* input */
					double	*Result				/* output */
				)

/* computes the first derivative of a Blu interpolant function */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin BlipDiff1 */

	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(BlipDiff1, Result, Status,
/**/	"No output data")
/**/DEBUG_RETURN_ON_ERROR(BlipDiff1, Status)
/**/DEBUG_WRITE_ENTERING(BlipDiff1,
/**/	"About to compute the first order derivative of a generic blip")

	switch (Degree) {
		case 0L:
			*Result = Blip00Diff1(Argument);
			break;
		case 1L:
			*Result = Blip01Diff1(Argument);
			break;
		case 3L:
			*Result = Blip03Diff1(Argument);
			break;
		default:
			Status = ERROR;
			WRITE_ERROR(BlipDiff1, "Invalid blip degree")
			break;
	}
/**/DEBUG_WRITE_LEAVING(BlipDiff1, "Done")
	return(Status);
} /* end BlipDiff1 */

/*--------------------------------------------------------------------------*/
extern double	Blip00Diff1
				(
					double	Argument			/* input */
				)

/* returns the value of the first derivative of a Blu interpolant function
	of degree 0 (order 1) evaluated at Argument */

{ /* begin Blip00Diff1 */

/**/DEBUG_WRITE_ENTERING(Blip00Diff1,
/**/	"About to compute the first derivative of a blip of degree 0")

	if (Argument == -0.5) {
		WRITE_ERROR(Blip00Diff1, "Undefined value at Argument == -0.5")
/**/	DEBUG_WRITE_LEAVING(Blip00Diff1, "Done")
		return(HUGE_VAL);
	}
	else if (Argument == 0.5) {
		WRITE_ERROR(Blip00Diff1, "Undefined value at Argument == 0.5")
/**/	DEBUG_WRITE_LEAVING(Blip00Diff1, "Done")
		return(-HUGE_VAL);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Blip00Diff1, "Done")
		return(0.0);
	}
} /* end Blip00Diff1 */

/*--------------------------------------------------------------------------*/
extern double	Blip01Diff1
				(
					double	Argument			/* input */
				)

/* returns the value of the first derivative of a Blu interpolant function
	of degree 1 (order 2) evaluated at Argument */

{ /* begin Blip01Diff1 */

	double	a;

/**/DEBUG_WRITE_ENTERING(Blip01Diff1,
/**/	"About to compute the first derivative of a blip of degree 1")

	a = fabs(Argument);
	if (a == 0.0) {
	}
	else if (a < 1.0) {
		a = -1.0;
	}
	else if (a == 1.0) {
		a = -1.0 / 2.0;
	}
	else {
		a = 0.0;
	}
/**/DEBUG_WRITE_LEAVING(Blip01Diff1, "Done")
	return((0.0 < Argument) ? (a) : (-a));
} /* end Blip01Diff1 */

/*--------------------------------------------------------------------------*/
extern double	Blip03Diff1
				(
					double	Argument			/* input */
				)

/* returns the value of the first derivative of a Blu interpolant function
	of degree 3 (order 4) evaluated at Argument */

{ /* begin Blip03Diff1 */

	double	a;

/**/DEBUG_WRITE_ENTERING(Blip03Diff1,
/**/	"About to compute the first derivative of a blip of degree 3")

	a = fabs(Argument);
	if (a == 0.0) {
	}
	else if (a < 1.0) {
		a = a * (a * (3.0 / 2.0) - 2.0) + 3.0 / 56.0;
	}
	else if (a == 1.0) {
		a = -27.0 / 56.0;
	}
	else if (a < 2.0) {
		a = a * (2.0 - a * (1.0 / 2.0)) - 113.0 / 56.0;
	}
	else if (a == 2.0) {
		a = -1.0 / 112.0;
	}
	else {
		a = 0.0;
	}
/**/DEBUG_WRITE_LEAVING(Blip03Diff1, "Done")
	return((0.0 < Argument) ? (a) : (-a));
} /* end Blip03Diff1 */

/*--------------------------------------------------------------------------*/
extern int		BsplineDiff1
				(
					long	Degree,				/* degree */
					double	Argument,			/* input */
					double	*Result				/* output */
				)

/* computes the first derivative of a Basic spline function */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin BsplineDiff1 */

	double	Left, Right;
	int		Status = !ERROR;

/**/DEBUG_CHECK_RANGE_LONG(BsplineDiff1, Degree, 0L, LONG_MAX, Status,
/**/	"Invalid degree (should be positive)")
/**/DEBUG_CHECK_NULL_POINTER(BsplineDiff1, Result, Status,
/**/	"No output data")
/**/DEBUG_RETURN_ON_ERROR(BsplineDiff1, Status)
/**/DEBUG_WRITE_ENTERING(BsplineDiff1,
/**/	"About to compute the first order derivative of a generic B-spline")

	switch (Degree) {
		case 0L:
			*Result = Bspline00Diff1(Argument);
			break;
		case 1L:
			*Result = Bspline01Diff1(Argument);
			break;
		case 2L:
			*Result = Bspline02Diff1(Argument);
			break;
		case 3L:
			*Result = Bspline03Diff1(Argument);
			break;
		case 4L:
			*Result = Bspline04Diff1(Argument);
			break;
		case 5L:
			*Result = Bspline05Diff1(Argument);
			break;
		case 6L:
			*Result = Bspline06Diff1(Argument);
			break;
		case 7L:
			*Result = Bspline07Diff1(Argument);
			break;
		case 8L:
			*Result = Bspline08Diff1(Argument);
			break;
		case 9L:
			*Result = Bspline09Diff1(Argument);
			break;
		case 10L:
			*Result = Bspline10Diff1(Argument);
			break;
		case 11L:
			*Result = Bspline11Diff1(Argument);
			break;
		default:
			Status = Bspline(Degree - 1L, Argument - 0.5, &Left);
			if (Status == ERROR) {
/**/			DEBUG_WRITE_LEAVING(BsplineDiff1, "Done")
				return(Status);
			}
			Status = Bspline(Degree - 1L, Argument + 0.5, &Right);
			if (Status == ERROR) {
/**/			DEBUG_WRITE_LEAVING(BsplineDiff1, "Done")
				return(Status);
			}
			*Result = Right - Left;
			break;
	}
/**/DEBUG_WRITE_LEAVING(BsplineDiff1, "Done")
	return(Status);
} /* end BsplineDiff1 */

/*--------------------------------------------------------------------------*/
extern double	Bspline00Diff1
				(
					double	Argument			/* input */
				)

/* returns the value of the first derivative of a Basic spline function
	of degree 0 (order 1) evaluated at Argument */

{ /* begin Bspline00Diff1 */

/**/DEBUG_WRITE_ENTERING(Bspline00Diff1,
/**/	"About to compute the first derivative of a B-spline of degree 0")

	if (Argument == -0.5) {
		WRITE_ERROR(Bspline00Diff1, "Undefined value at Argument == -0.5")
/**/	DEBUG_WRITE_LEAVING(Bspline00Diff1, "Done")
		return(HUGE_VAL);
	}
	else if (Argument == 0.5) {
		WRITE_ERROR(Bspline00Diff1, "Undefined value at Argument == 0.5")
/**/	DEBUG_WRITE_LEAVING(Bspline00Diff1, "Done")
		return(-HUGE_VAL);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Bspline00Diff1, "Done")
		return(0.0);
	}
} /* end Bspline00Diff1 */

/*--------------------------------------------------------------------------*/
extern double	Bspline01Diff1
				(
					double	Argument			/* input */
				)

/* returns the value of the first derivative of a Basic spline function
	of degree 1 (order 2) evaluated at Argument */

{ /* begin Bspline01Diff1 */

	double	a;

/**/DEBUG_WRITE_ENTERING(Bspline01Diff1,
/**/	"About to compute the first derivative of a B-spline of degree 1")

	a = fabs(Argument);
	if (a == 0.0) {
	}
	else if (a < 1.0) {
		a = -1.0;
	}
	else if (a == 1.0) {
		a = -1.0 / 2.0;
	}
	else {
		a = 0.0;
	}
/**/DEBUG_WRITE_LEAVING(Bspline01Diff1, "Done")
	return((0.0 < Argument) ? (a) : (-a));
} /* end Bspline01Diff1 */

/*--------------------------------------------------------------------------*/
extern double	Bspline02Diff1
				(
					double	Argument			/* input */
				)

/* returns the value of the first derivative of a Basic spline function
	of degree 2 (order 3) evaluated at Argument */

{ /* begin Bspline02Diff1 */

	double	a;

/**/DEBUG_WRITE_ENTERING(Bspline02Diff1,
/**/	"About to compute the first derivative of a B-spline of degree 2")

	a = fabs(Argument);
	if (a < 0.5) {
		a *= -2.0;
	}
	else if (a < 1.5) {
		a -= 3.0 / 2.0;
	}
	else {
		a = 0.0;
	}
/**/DEBUG_WRITE_LEAVING(Bspline02Diff1, "Done")
	return((0.0 < Argument) ? (a) : (-a));
} /* end Bspline02Diff1 */

/*--------------------------------------------------------------------------*/
extern double	Bspline03Diff1
				(
					double	Argument			/* input */
				)

/* returns the value of the first derivative of a Basic spline function
	of degree 3 (order 4) evaluated at Argument */

{ /* begin Bspline03Diff1 */

	double	a;

/**/DEBUG_WRITE_ENTERING(Bspline03Diff1,
/**/	"About to compute the first derivative of a B-spline of degree 3")

	a = fabs(Argument);
	if (a < 1.0) {
		a *= a * (3.0 / 2.0) - 2.0;
	}
	else if (a < 2.0) {
		a = 2.0 - a;
		a *= a * (-1.0 / 2.0);
	}
	else {
		a = 0.0;
	}
/**/DEBUG_WRITE_LEAVING(Bspline03Diff1, "Done")
	return((0.0 < Argument) ? (a) : (-a));
} /* end Bspline03Diff1 */

/*--------------------------------------------------------------------------*/
extern double	Bspline04Diff1
				(
					double	Argument			/* input */
				)

/* returns the value of the first derivative of a Basic spline function
	of degree 4 (order 5) evaluated at Argument */

{ /* begin Bspline04Diff1 */

	double	a;

/**/DEBUG_WRITE_ENTERING(Bspline04Diff1,
/**/	"About to compute the first derivative of a B-spline of degree 4")

	a = fabs(Argument);
	if (a < 0.5) {
		a *= a * a - 5.0 / 4.0;
	}
	else if (a < 1.5) {
		a = a * (a * (5.0 / 2.0 - a * (2.0 / 3.0)) - 5.0 / 2.0) + 5.0 / 24.0;
	}
	else if (a < 2.5) {
		a -= 5.0 / 2.0;
		a *= a * a * (1.0 / 6.0);
	}
	else {
		a = 0.0;
	}
/**/DEBUG_WRITE_LEAVING(Bspline04Diff1, "Done")
	return((0.0 < Argument) ? (a) : (-a));
} /* end Bspline04Diff1 */

/*--------------------------------------------------------------------------*/
extern double	Bspline05Diff1
				(
					double	Argument			/* input */
				)

/* returns the value of the first derivative of a Basic spline function
	of degree 5 (order 6) evaluated at Argument */

{ /* begin Bspline05Diff1 */

	double	a;

/**/DEBUG_WRITE_ENTERING(Bspline05Diff1,
/**/	"About to compute the first derivative of a B-spline of degree 5")

	a = fabs(Argument);
	if (a < 1.0) {
		a *= a * a * (1.0 - a * (5.0 / 12.0)) - 1.0;
	}
	else if (a < 2.0) {
		a = a * (a * (a * (a * (5.0 / 24.0) - 3.0 / 2.0) + 15.0 / 4.0)
			- 7.0 / 2.0) + 5.0 / 8.0;
	}
	else if (a < 3.0) {
		a -= 3.0;
		a *= a;
		a *= a * (-1.0 / 24.0);
	}
	else {
		a = 0.0;
	}
/**/DEBUG_WRITE_LEAVING(Bspline05Diff1, "Done")
	return((0.0 < Argument) ? (a) : (-a));
} /* end Bspline05Diff1 */

/*--------------------------------------------------------------------------*/
extern double	Bspline06Diff1
				(
					double	Argument			/* input */
				)

/* returns the value of the first derivative of a Basic spline function
	of degree 6 (order 7) evaluated at Argument */

{ /* begin Bspline06Diff1 */

	double	a, a2;

/**/DEBUG_WRITE_ENTERING(Bspline06Diff1,
/**/	"About to compute the first derivative of a B-spline of degree 6")

	a = fabs(Argument);
	if (a < 0.5) {
		a2 = a * a;
		a *= a2 * (7.0 / 12.0 - a2 * (1.0 / 6.0)) - 77.0 / 96.0;
	}
	else if (a < 1.5) {
		a = a * (a * (a * (a * (a * (1.0 / 8.0) - 35.0 / 48.0) + 21.0 / 16.0) - 35.0 / 96.0)
			- 91.0 / 128.0) - 7.0 / 768.0;
	}
	else if (a < 2.5) {
		a = a * (a * (a * (a * (7.0 / 12.0 - a * (1.0 / 20.0)) - 21.0 / 8.0) + 133.0 / 24.0)
			- 329.0 / 64.0) + 1267.0 / 960.0;
	}
	else if (a < 3.5) {
		a -= 7.0 / 2.0;
		a2 = a * a;
		a *= a2 * a2 * (1.0 / 120.0);
	}
	else {
		a = 0.0;
	}
/**/DEBUG_WRITE_LEAVING(Bspline06Diff1, "Done")
	return((0.0 < Argument) ? (a) : (-a));
} /* end Bspline06Diff1 */

/*--------------------------------------------------------------------------*/
extern double	Bspline07Diff1
				(
					double	Argument			/* input */
				)

/* returns the value of the first derivative of a Basic spline function
	of degree 7 (order 8) evaluated at Argument */

{ /* begin Bspline07Diff1 */

	double	a, a2;

/**/DEBUG_WRITE_ENTERING(Bspline07Diff1,
/**/	"About to compute the first derivative of a B-spline of degree 7")

	a = fabs(Argument);
	if (a < 1.0) {
		a2 = a * a;
		a *= a2 * (a2 * (a * (7.0 / 144.0) - 1.0 / 6.0) + 4.0 / 9.0) - 2.0 / 3.0;
	}
	else if (a < 2.0) {
		a = a * (a * (a * (a * (a * (3.0 / 10.0 - a * (7.0 / 240.0)) - 7.0 / 6.0)
			+ 2.0) - 7.0 / 6.0) - 2.0 / 10.0) - 7.0 / 90.0;
	}
	else if (a < 3.0) {
		a = a * (a * (a * (a * (a * (a * (7.0 / 720.0) - 1.0 / 6.0) + 7.0 / 6.0)
			- 38.0 / 9.0) + 49.0 / 6.0) - 23.0 / 3.0) + 217.0 / 90.0;
	}
	else if (a < 4.0) {
		a -= 4.0;
		a *= a;
		a *= a * a * (-1.0 / 720.0);
	}
	else {
		a = 0.0;
	}
/**/DEBUG_WRITE_LEAVING(Bspline07Diff1, "Done")
	return((0.0 < Argument) ? (a) : (-a));
} /* end Bspline07Diff1 */

/*--------------------------------------------------------------------------*/
extern double	Bspline08Diff1
				(
					double	Argument			/* input */
				)

/* returns the value of the first derivative of a Basic spline function
	of degree 8 (order 9) evaluated at Argument */

{ /* begin Bspline08Diff1 */

	double	a, a2;

/**/DEBUG_WRITE_ENTERING(Bspline08Diff1,
/**/	"About to compute the first derivative of a B-spline of degree 8")

	a = fabs(Argument);
	if (a < 0.5) {
		a2 = a * a;
		a *= a2 * (a2 * (a2 * (1.0 / 72.0) - 3.0 / 32) + 43.0 / 128.0) - 289.0 / 512.0;
	}
	else if (a < 1.5) {
		a = a * (a * (a * (a * (a * (a * (7.0 / 80.0 - a * (1.0 / 90.0))
			- 9.0 / 40.0) + 7.0 / 64.0) + 9.0 / 32.0) + 21.0 / 1280.0)
			- 363.0 / 640.0) + 1.0 / 5120.0;
	}
	else if (a < 2.5) {
		a = a * (a * (a * (a * (a * (a * (a * (1.0 / 180.0) - 7.0 / 80.0)
			+ 9.0 / 16.0) - 119.0 / 64.0) + 207.0 / 64.0) - 3381.0 / 1280.0)
			+ 195.0 / 256.0) - 1457.0 / 5120.0;
	}
	else if (a < 3.5) {
		a = a * (a * (a * (a * (a * (a * (3.0 / 80.0 - a * (1.0 / 630.0))
			- 3.0 / 8.0) + 131.0 / 64.0) - 209.0 / 32.0) + 15369.0 / 1280.0)
			- 1465.0 / 128.0) + 146051.0 / 35840.0;
	}
	else if (a < 4.5) {
		a -= 9.0 / 2.0;
		a2 = a * a * a;
		a2 *= a2;
		a *= a2 * (1.0 / 5040.0);
	}
	else {
		a = 0.0;
	}
/**/DEBUG_WRITE_LEAVING(Bspline08Diff1, "Done")
	return((0.0 < Argument) ? (a) : (-a));
} /* end Bspline08Diff1 */

/*--------------------------------------------------------------------------*/
extern double	Bspline09Diff1
				(
					double	Argument			/* input */
				)

/* returns the value of the first derivative of a Basic spline function
	of degree 9 (order 10) evaluated at Argument */

{ /* begin Bspline09Diff1 */

	double	a, a2;

/**/DEBUG_WRITE_ENTERING(Bspline09Diff1,
/**/	"About to compute the first derivative of a B-spline of degree 9")

	a = fabs(Argument);
	if (a < 1.0) {
		a2 = a * a;
		a *= a2 * (a2 * (a2 * (1.0 / 72.0 - a * (1.0 / 320.0)) - 5.0 / 72.0)
			+ 19.0 / 72.0) - 35.0 / 72.0;
	}
	else if (a < 2.0) {
		a = a * (a * (a * (a * (a * (a * (a * (a * (1.0 / 480.0) - 1.0 / 36.0)
			+ 7.0 / 48.0) - 13.0 / 36.0) + 35.0 / 96.0) - 1.0 / 36.0) + 7.0 / 48.0)
			- 19.0 / 36.0) + 1.0 / 192.0;
	}
	else if (a < 3.0) {
		a = a * (a * (a * (a * (a * (a * (a * (5.0 / 252.0 - a * (1.0 / 1120.0))
			- 3.0 / 16.0) +35.0 / 36.0) - 95.0 / 32.0) + 191.0 / 36.0) - 83.0 / 16.0)
			+ 635.0 / 252.0) - 339.0 / 448.0;
	}
	else if (a < 4.0) {
		a = a * (a * (a * (a * (a * (a * (a * (a * (1.0 / 4480.0) - 1.0 / 144.0)
			+ 3.0 / 32.0) - 103.0 / 144.0) + 215.0 / 64.0) - 1423.0 / 144.0)
			+ 563.0 / 32.0) - 2449.0 / 144.0) + 5883.0 / 896.0;
	}
	else if (a < 5.0) {
		a -= 5.0;
		a *= a;
		a *= a;
		a *= a * (-1.0 / 40320.0);
	}
	else {
		a = 0.0;
	}
/**/DEBUG_WRITE_LEAVING(Bspline09Diff1, "Done")
	return((0.0 < Argument) ? (a) : (-a));
} /* end Bspline09Diff1 */

/*--------------------------------------------------------------------------*/
extern double	Bspline10Diff1
				(
					double	Argument			/* input */
				)

/* returns the value of the first derivative of a Basic spline function
	of degree 10 (order 11) evaluated at Argument */

{ /* begin Bspline10Diff1 */

	double	a, a2;

/**/DEBUG_WRITE_ENTERING(Bspline10Diff1,
/**/	"About to compute the first derivative of a B-spline of degree 10")

	a = fabs(Argument);
	if (a < 0.5) {
		a2 = a * a;
		a *= a2 * (a2 * (a2 * (11.0 / 1440.0 - a2 * (1.0 / 1440.0))
			- 583.0 / 11520.0) + 14597.0 / 69120.0) - 156409.0 / 368640.0;
	}
	else if (a < 1.5) {
		a = a * (a * (a * (a * (a * (a * (a * (a * (a * (1.0 / 1728.0)
			- 11.0 / 1920.0) + 11.0 / 576.0) - 77.0 / 5760.0) - 187.0 / 4608.0)
			- 77.0 / 15360.0) + 5885.0 / 27648.0) - 11.0 / 30720.0)
			- 62557.0 / 147456.0) - 11.0 / 4423680.0;
	}
	else if (a < 2.5) {
		a = a * (a * (a * (a * (a * (a * (a * (a * (11.0 / 1680.0 - a * (1.0 / 3024.0))
			- 55.0 / 1008.0) + 11.0 / 45.0) - 715.0 / 1152.0) + 1661.0 / 1920.0)
			- 4543.0 / 6912.0) + 7513.0 / 13440.0) - 163603.0 / 258048.0)
			+ 135311.0 / 3870720.0;
	}
	else if (a < 3.5) {
		a = a * (a * (a * (a * (a * (a * (a * (a * (a * (1.0 / 8064.0) - 33.0 / 8960.0)
			+ 55.0 / 1152.0) - 451.0 / 1280.0) + 14905.0 / 9216.0) - 48433.0 / 10240.0)
			+ 479281.0 / 55296.0) - 1352153.0 / 143360.0) + 1654543.0 / 294912.0)
			- 11695211.0 / 6881280.0;
	}
	else if (a < 4.5) {
		a = a * (a * (a * (a * (a * (a * (a * (a * (11.0 / 10080.0 - a * (1.0 / 36288.0))
			- 11.0 / 576.0) + 209.0 / 1080.0) - 5753.0 / 4608.0) + 61061.0 / 11520.0)
			- 407429.0 / 27648.0) + 2070343.0 / 80640.0) - 3702215.0 / 147456.0)
			+ 237959711.0 / 23224320.0;
	}
	else if (a < 5.5) {
		a -= 11.0 / 2.0;
		a2 = a * a;
		a2 *= a2;
		a2 *= a2;
		a *= a2 * (1.0 / 362880.0);
	}
	else {
		a = 0.0;
	}
/**/DEBUG_WRITE_LEAVING(Bspline10Diff1, "Done")
	return((0.0 < Argument) ? (a) : (-a));
} /* end Bspline10Diff1 */

/*--------------------------------------------------------------------------*/
extern double	Bspline11Diff1
				(
					double	Argument			/* input */
				)

/* returns the value of the first derivative of a Basic spline function
	of degree 11 (order 12) evaluated at Argument */

{ /* begin Bspline11Diff1 */

	double	a, a2;

/**/DEBUG_WRITE_ENTERING(Bspline11Diff1,
/**/	"About to compute the first derivative of a B-spline of degree 11")

	a = fabs(Argument);
	if (a < 1.0) {
		a2 = a * a;
		a *= a2 * (a2 * (a2 * (a2 * (a * (11.0 / 86400.0) - 1.0 / 1440.0)
			+ 1.0 / 180.0) - 23.0 / 600.0) + 31.0 / 180.0) - 809.0 / 2160.0;
	}
	else if (a < 2.0) {
		a = a * (a * (a * (a * (a * (a * (a * (a * (a * (1.0 / 672.0
			- a * (11.0 / 120960.0)) - 11.0 / 1120.0) + 2.0 / 63.0) - 11.0 / 240.0)
			+ 1.0 / 60.0) - 11.0 / 240.0) + 25.0 / 126.0) - 11.0 / 1120.0)
			- 563.0 / 1512.0) - 11.0 / 50400.0;
	}
	else if (a < 3.0) {
		a = a * (a * (a * (a * (a * (a * (a * (a * (a * (a * (11.0 / 241920.0)
			- 5.0 / 4032.0) + 33.0 / 2240.0) - 25.0 / 252.0) + 33.0 / 80.0)
			- 13.0 / 12.0) + 143.0 / 80.0) - 239.0 / 126.0) + 1749.0 / 1120.0)
			- 1619.0 / 1512.0) + 781.0 / 5600.0;
	}
	else if (a < 4.0) {
		a = a * (a * (a * (a * (a * (a * (a * (a * (a * (1.0 / 1728.0
			- a * (11.0 / 725760.0)) - 11.0 / 1120.0) + 7.0 / 72.0)
			- 99.0 / 160.0) + 631.0 / 240.0) - 1199.0 / 160.0) + 1009.0 / 72.0)
			- 36597.0 / 2240.0) + 9385.0 / 864.0) - 38533.0 / 11200.0;
	}
	else if (a < 5.0) {
		a = a * (a * (a * (a * (a * (a * (a * (a * (a * (a * (11.0 / 3628800.0)
			- 1.0 / 6720.0) + 11.0 / 3360.0) - 107.0 / 2520.0) + 517.0 / 1440.0)
			- 2477.0 / 1200.0) + 11737.0 / 1440.0) - 54797.0 / 2520.0)
			+ 250657.0 / 6720.0) - 1113317.0 / 30240.0) + 4726777.0 / 302400.0;
	}
	else if (a < 6.0) {
		a -= 6.0;
		a *= a;
		a2 = a * a;
		a *= a2 * a2 * (-1.0 / 3628800.0);
	}
	else {
		a = 0.0;
	}
/**/DEBUG_WRITE_LEAVING(Bspline11Diff1, "Done")
	return((0.0 < Argument) ? (a) : (-a));
} /* end Bspline11Diff1 */

/*--------------------------------------------------------------------------*/
extern double	DodgsonDiff1
				(
					double	Argument			/* input */
				)

/* returns the value of the first derivative of the Dodgson kernel
	(order 2) evaluated at Argument */

{ /* begin DodgsonDiff1 */

	double	a;

/**/DEBUG_WRITE_ENTERING(DodgsonDiff1,
/**/	"About to compute the first derivative of a Keys")

	a = fabs(Argument);
	if (a < 0.5) {
		a *= -4.0;
	}
	else if (a == 0.5) {
		a = -7.0 / 2.0;
	}
	else if (a < 1.5) {
		a = 2.0 * a - 5.0 / 2.0;
	}
	else if (a == 1.5) {
		a = 1.0 / 4.0;
	}
	else {
		a = 0.0;
	}
/**/DEBUG_WRITE_LEAVING(DodgsonDiff1, "Done")
	return((0.0 < Argument) ? (a) : (-a));
} /* end DodgsonDiff1 */

/*--------------------------------------------------------------------------*/
extern double	German04Diff1
				(
					double	Argument			/* input */
				)

/* returns the value of the first derivative of a German function
	of degree 4 (order 5) evaluated at Argument */

{ /* begin German04Diff1 */

	double	a;

/**/DEBUG_WRITE_ENTERING(German04Diff1,
/**/	"About to compute the first derivative of a German of degree 4")

	a = fabs(Argument);
	if (a < 1.0) {
		a *= a * (a * (1.0 / 6.0) + 185.0 / 48.0) - 335.0 / 72.0;
	}
	else if (a < 2.0) {
		a = a * (a * (-29.0 / 48.0 - a * (1.0 / 4.0)) + 17.0 / 4.0) - 145.0 / 36.0;
	}
	else if (a < 3.0) {
		a = a * (a * (a * (1.0 / 12.0) - 7.0 / 16.0) + 5.0 / 9.0) + 1.0 / 36.0;
	}
	else if (a < 4.0) {
		a = a * (a * (1.0 / 48.0) - 11.0 / 72.0) + 5.0 / 18.0;
	}
	else {
		a = 0.0;
	}
/**/DEBUG_WRITE_LEAVING(German04Diff1, "Done")
	return((0.0 < Argument) ? (a) : (-a));
} /* end German04Diff1 */

/*--------------------------------------------------------------------------*/
extern double	KeysOptimalDiff1
				(
					double	Argument			/* input */
				)

/* returns the value of the first derivative of the cubic Keys optimal kernel
	(order 3) evaluated at Argument */

{ /* begin KeysOptimalDiff1 */

	double	a;

/**/DEBUG_WRITE_ENTERING(KeysOptimalDiff1,
/**/	"About to compute the first derivative of a Keys")

	a = fabs(Argument);
	if (a < 1.0) {
		a *= a * (9.0 / 2.0) - 5.0;
	}
	else if (a < 2.0) {
		a = a * (5.0 - a * (3.0 / 2.0)) - 4.0;
	}
	else {
		a = 0.0;
	}
/**/DEBUG_WRITE_LEAVING(KeysOptimalDiff1, "Done")
	return((0.0 < Argument) ? (a) : (-a));
} /* end KeysOptimalDiff1 */

/*--------------------------------------------------------------------------*/
extern int		OmomsDiff1
				(
					long	Degree,				/* degree */
					double	Argument,			/* input */
					double	*Result				/* output */
				)

/* computes the first derivative of a Blu optimum function */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin OmomsDiff1 */

	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(OmomsDiff1, Result, Status,
/**/	"No output data")
/**/DEBUG_RETURN_ON_ERROR(OmomsDiff1, Status)
/**/DEBUG_WRITE_ENTERING(OmomsDiff1,
/**/	"About to compute the first order derivative of a generic oMoms")

	switch (Degree) {
		case 0L:
			*Result = Omoms00Diff1(Argument);
			break;
		case 1L:
			*Result = Omoms01Diff1(Argument);
			break;
		case 3L:
			*Result = Omoms03Diff1(Argument);
			break;
		default:
			Status = ERROR;
			WRITE_ERROR(OmomsDiff1, "Invalid oMoms degree")
			break;
	}
/**/DEBUG_WRITE_LEAVING(OmomsDiff1, "Done")
	return(Status);
} /* end OmomsDiff1 */

/*--------------------------------------------------------------------------*/
extern double	Omoms00Diff1
				(
					double	Argument			/* input */
				)

/* returns the value of the first derivative of a Blu optimum function
	of degree 0 (order 1) evaluated at Argument */

{ /* begin Omoms00Diff1 */

/**/DEBUG_WRITE_ENTERING(Omoms00Diff1,
/**/	"About to compute the first derivative of an oMoms of degree 0")

	if (Argument == -0.5) {
		WRITE_ERROR(Omoms00Diff1, "Undefined value at Argument == -0.5")
/**/	DEBUG_WRITE_LEAVING(Omoms00Diff1, "Done")
		return(HUGE_VAL);
	}
	else if (Argument == 0.5) {
		WRITE_ERROR(Omoms00Diff1, "Undefined value at Argument == 0.5")
/**/	DEBUG_WRITE_LEAVING(Omoms00Diff1, "Done")
		return(-HUGE_VAL);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Omoms00Diff1, "Done")
		return(0.0);
	}
} /* end Omoms00Diff1 */

/*--------------------------------------------------------------------------*/
extern double	Omoms01Diff1
				(
					double	Argument			/* input */
				)

/* returns the value of the first derivative of a Blu optimum function
	of degree 1 (order 2) evaluated at Argument */

{ /* begin Omoms01Diff1 */

	double	a;

/**/DEBUG_WRITE_ENTERING(Omoms01Diff1,
/**/	"About to compute the first derivative of an oMoms of degree 1")

	a = fabs(Argument);
	if (a == 0.0) {
	}
	else if (a < 1.0) {
		a = -1.0;
	}
	else if (a == 1.0) {
		a = -1.0 / 2.0;
	}
	else {
		a = 0.0;
	}
/**/DEBUG_WRITE_LEAVING(Omoms01Diff1, "Done")
	return((0.0 < Argument) ? (a) : (-a));
} /* end Omoms01Diff1 */

/*--------------------------------------------------------------------------*/
extern double	Omoms02Diff1
				(
					double	Argument			/* input */
				)

/* returns the value of the first derivative of a Blu optimum function
	of degree 2 (order 3) evaluated at Argument */

{ /* begin Omoms02Diff1 */

	double	a;

/**/DEBUG_WRITE_ENTERING(Omoms02Diff1,
/**/	"About to compute the first derivative of an oMoms of degree 2")

	a = fabs(Argument);
	if (a < 0.5) {
		a *= -2.0;
	}
	else if (a == 0.5) {
		WRITE_ERROR(Omoms02Diff1, "Undefined value at |Argument| == 0.5")
/**/	DEBUG_WRITE_LEAVING(Omoms02Diff1, "Done")
		return(HUGE_VAL);
	}
	else if (a < 1.5) {
		a -= 3.0 / 2.0;
	}
	else if (a == 1.5) {
		WRITE_ERROR(Omoms02Diff1, "Undefined value at |Argument| == 1.5")
/**/	DEBUG_WRITE_LEAVING(Omoms02Diff1, "Done")
		return(-HUGE_VAL);
	}
	else {
		a = 0.0;
	}
/**/DEBUG_WRITE_LEAVING(Omoms02Diff1, "Done")
	return((0.0 < Argument) ? (a) : (-a));
} /* end Omoms02Diff1 */

/*--------------------------------------------------------------------------*/
extern double	Omoms03Diff1
				(
					double	Argument			/* input */
				)

/* returns the value of the first derivative of a Blu optimum function
	of degree 3 (order 4) evaluated at Argument */

{ /* begin Omoms03Diff1 */

	double	a;

/**/DEBUG_WRITE_ENTERING(Omoms03Diff1,
/**/	"About to compute the first derivative of an oMoms of degree 3")

	a = fabs(Argument);
	if (a == 0.0) {
	}
	else if (a < 1.0) {
		a = a * (a * (3.0 / 2.0) - 2.0) + 1.0 / 14.0;
	}
	else if (a == 1.0) {
		a = -10.0 / 21.0;
	}
	else if (a < 2.0) {
		a = a * (2.0 - a * (1.0 / 2.0)) - 85.0 / 42.0;
	}
	else if (a == 2.0) {
		a = -1.0 / 84.0;
	}
	else {
		a = 0.0;
	}
/**/DEBUG_WRITE_LEAVING(Omoms03Diff1, "Done")
	return((0.0 < Argument) ? (a) : (-a));
} /* end Omoms03Diff1 */

/*--------------------------------------------------------------------------*/
extern double	PositiveDiff1
				(
					double	Argument			/* input */
				)

/* returns the value of the first derivative of Positive kernel
	(order 1) evaluated at Argument */

{ /* begin PositiveDiff1 */

	double	a;

/**/DEBUG_WRITE_ENTERING(PositiveDiff1,
/**/	"About to compute the first derivative of a Positive")

	a = fabs(Argument);
	if (a < 0.5) {
		a *= a * 12.0 - 9.0;
	}
	else if (a < 1.0) {
		a = a * 3.0 - 3.0;
	}
	else if (a < 1.5) {
		a = a * (15.0 - a * 6.0) - 9.0;
	}
	else if (a < 2.0) {
		a = a * (a * 6.0 - 21.0) + 18.0;
	}
	else {
		a = 0.0;
	}
/**/DEBUG_WRITE_LEAVING(PositiveDiff1, "Done")
	return((0.0 < Argument) ? (a) : (-a));
} /* end PositiveDiff1 */

/*--------------------------------------------------------------------------*/
extern double	SincDiff1
				(
					double	Argument			/* input */
				)

/* returns the value of the first derivative of the sinc kernel evaluated at Argument */

{ /* begin SincDiff1 */

/**/DEBUG_WRITE_ENTERING(SincDiff1,
/**/	"About to compute the first derivative of a sinc")

	if (Argument == 0.0) {
/**/	DEBUG_WRITE_LEAVING(SincDiff1, "Done")
		return(0.0);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(SincDiff1, "Done")
		return((cos(PI * Argument) - Sinc(Argument)) / Argument);
	}
} /* end SincDiff1 */

