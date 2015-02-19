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
#include	"messagedisplay.h"
#include	"positivepower.h"

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
extern int		Blip
				(
					long	Degree,				/* degree */
					double	Argument,			/* input */
					double	*Result				/* output */
				)

/* computes a Blu interpolant function */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin Blip */

	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(Blip, Result, Status,
/**/	"No output data")
/**/DEBUG_RETURN_ON_ERROR(Blip, Status)
/**/DEBUG_WRITE_ENTERING(Blip,
/**/	"About to compute a generic blip")

	switch (Degree) {
		case 0L:
			*Result = Blip00(Argument);
			break;
		case 1L:
			*Result = Blip01(Argument);
			break;
		case 3L:
			*Result = Blip03(Argument);
			break;
		default:
			Status = ERROR;
			WRITE_ERROR(Blip, "Invalid blip degree")
			break;
	}
/**/DEBUG_WRITE_LEAVING(Blip, "Done")
	return(Status);
} /* end Blip */

/*--------------------------------------------------------------------------*/
extern double	Blip00
				(
					double	Argument			/* input */
				)

/* returns the value of a Blu interpolant function of degree 0 (order 1) evaluated at Argument */

{ /* begin Blip00 */

/**/DEBUG_WRITE_ENTERING(Blip00,
/**/	"About to compute a blip of degree 0")

	Argument = fabs(Argument);
	if (Argument < 0.5) {
/**/	DEBUG_WRITE_LEAVING(Blip00, "Done")
		return(1.0);
	}
	else if (Argument == 0.5) {
/**/	DEBUG_WRITE_LEAVING(Blip00, "Done")
		return(0.5);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Blip00, "Done")
		return(0.0);
	}
} /* end Blip00 */

/*--------------------------------------------------------------------------*/
extern double	Blip01
				(
					double	Argument			/* input */
				)

/* returns the value of a Blu interpolant function of degree 1 (order 2) evaluated at Argument */

{ /* begin Blip01 */

/**/DEBUG_WRITE_ENTERING(Blip01,
/**/	"About to compute a blip of degree 1")

	Argument = fabs(Argument);
	if (Argument < 1.0) {
/**/	DEBUG_WRITE_LEAVING(Blip01, "Done")
		return(1.0 - Argument);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Blip01, "Done")
		return(0.0);
	}
} /* end Blip01 */

/*--------------------------------------------------------------------------*/
extern double	Blip03
				(
					double	Argument			/* input */
				)

/* returns the value of a Blu interpolant function of degree 3 (order 4) evaluated at Argument */

{ /* begin Blip03 */

/**/DEBUG_WRITE_ENTERING(Blip03,
/**/	"About to compute a blip of degree 3")

	Argument = fabs(Argument);
	if (Argument < 1.0) {
/**/	DEBUG_WRITE_LEAVING(Blip03, "Done")
		return(Argument * (Argument * (Argument * (1.0 / 2.0) - 1.0) + 3.0 / 56.0)
			+ 53.0 / 84.0);
	}
	else if (Argument < 2.0) {
/**/	DEBUG_WRITE_LEAVING(Blip03, "Done")
		return(Argument * (Argument * (1.0 - Argument * (1.0 / 6.0)) - 113.0 / 56.0)
			+ 115.0 / 84.0);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Blip03, "Done")
		return(0.0);
	}
} /* end Blip03 */

/*--------------------------------------------------------------------------*/
extern int		Bspline
				(
					long	Degree,				/* degree */
					double	Argument,			/* input */
					double	*Result				/* output */
				)

/* computes a Basic spline function */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin Bspline */

	double	Up, Shift, Binomial;
	long	i;
	int		Status = !ERROR;

/**/DEBUG_CHECK_RANGE_LONG(Bspline, Degree, 0L, LONG_MAX, Status,
/**/	"Invalid degree (should be positive)")
/**/DEBUG_CHECK_NULL_POINTER(Bspline, Result, Status,
/**/	"No output data")
/**/DEBUG_RETURN_ON_ERROR(Bspline, Status)
/**/DEBUG_WRITE_ENTERING(Bspline,
/**/	"About to compute a generic B-spline")

	switch (Degree) {
		case 0L:
			*Result = Bspline00(Argument);
			break;
		case 1L:
			*Result = Bspline01(Argument);
			break;
		case 2L:
			*Result = Bspline02(Argument);
			break;
		case 3L:
			*Result = Bspline03(Argument);
			break;
		case 4L:
			*Result = Bspline04(Argument);
			break;
		case 5L:
			*Result = Bspline05(Argument);
			break;
		case 6L:
			*Result = Bspline06(Argument);
			break;
		case 7L:
			*Result = Bspline07(Argument);
			break;
		case 8L:
			*Result = Bspline08(Argument);
			break;
		case 9L:
			*Result = Bspline09(Argument);
			break;
		case 10L:
			*Result = Bspline10(Argument);
			break;
		case 11L:
			*Result = Bspline11(Argument);
			break;
		default:
			Argument = fabs(Argument);
			*Result = 0.0;
			Up = (double)(Degree + 1L);
			Argument += 0.5 * Up;
			Binomial = 1.0;
			for (i = 2L; (i <= Degree); i++) {
				Binomial *= (double)i;
			}
			Binomial = 1.0 / Binomial;
			Shift = 0.0;
			while ((Shift <= Up) && (Shift < Argument)) {
				*Result += Binomial * PositiveIntPower(Argument - Shift, Degree);
				Binomial *= (Shift - Up) / (Shift + 1.0);
				Shift += 1.0;
			}
			break;
	}
/**/DEBUG_WRITE_LEAVING(Bspline, "Done")
	return(Status);
} /* end Bspline */

/*--------------------------------------------------------------------------*/
extern double	Bspline00
				(
					double	Argument			/* input */
				)

/* returns the value of a Basic spline function of degree 0 (order 1) evaluated at Argument */

{ /* begin Bspline00 */

/**/DEBUG_WRITE_ENTERING(Bspline00,
/**/	"About to compute a B-spline of degree 0")

	Argument = fabs(Argument);
	if (Argument < 0.5) {
/**/	DEBUG_WRITE_LEAVING(Bspline00, "Done")
		return(1.0);
	}
	else if (Argument == 0.5) {
/**/	DEBUG_WRITE_LEAVING(Bspline00, "Done")
		return(0.5);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Bspline00, "Done")
		return(0.0);
	}
} /* end Bspline00 */

/*--------------------------------------------------------------------------*/
extern double	Bspline01
				(
					double	Argument			/* input */
				)

/* returns the value of a Basic spline function of degree 1 (order 2) evaluated at Argument */

{ /* begin Bspline01 */

/**/DEBUG_WRITE_ENTERING(Bspline01,
/**/	"About to compute a B-spline of degree 1")

	Argument = fabs(Argument);
	if (Argument < 1.0) {
/**/	DEBUG_WRITE_LEAVING(Bspline01, "Done")
		return(1.0 - Argument);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Bspline01, "Done")
		return(0.0);
	}
} /* end Bspline01 */

/*--------------------------------------------------------------------------*/
extern double	Bspline02
				(
					double	Argument			/* input */
				)

/* returns the value of a Basic spline function of degree 2 (order 3) evaluated at Argument */

{ /* begin Bspline02 */

/**/DEBUG_WRITE_ENTERING(Bspline02,
/**/	"About to compute a B-spline of degree 2")

	Argument = fabs(Argument);
	if (Argument < 0.5) {
/**/	DEBUG_WRITE_LEAVING(Bspline02, "Done")
		return(3.0 / 4.0 - Argument * Argument);
	}
	else if (Argument < 1.5) {
		Argument -= 3.0 / 2.0;
/**/	DEBUG_WRITE_LEAVING(Bspline02, "Done")
		return(Argument * Argument * (1.0 / 2.0));
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Bspline02, "Done")
		return(0.0);
	}
} /* end Bspline02 */

/*--------------------------------------------------------------------------*/
extern double	Bspline03
				(
					double	Argument			/* input */
				)

/* returns the value of a Basic spline function of degree 3 (order 4) evaluated at Argument */

{ /* begin Bspline03 */

/**/DEBUG_WRITE_ENTERING(Bspline03,
/**/	"About to compute a B-spline of degree 3")

	Argument = fabs(Argument);
	if (Argument < 1.0) {
/**/	DEBUG_WRITE_LEAVING(Bspline03, "Done")
		return(Argument * Argument * (Argument - 2.0) * (1.0 / 2.0) + 2.0 / 3.0);
	}
	else if (Argument < 2.0) {
		Argument -= 2.0;
/**/	DEBUG_WRITE_LEAVING(Bspline03, "Done")
		return(Argument * Argument * Argument * (-1.0 / 6.0));
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Bspline03, "Done")
		return(0.0);
	}
} /* end Bspline03 */

/*--------------------------------------------------------------------------*/
extern double	Bspline04
				(
					double	Argument			/* input */
				)

/* returns the value of a Basic spline function of degree 4 (order 5) evaluated at Argument */

{ /* begin Bspline04 */

/**/DEBUG_WRITE_ENTERING(Bspline04,
/**/	"About to compute a B-spline of degree 4")

	Argument = fabs(Argument);
	if (Argument < 0.5) {
		Argument *= Argument;
/**/	DEBUG_WRITE_LEAVING(Bspline04, "Done")
		return(Argument * (Argument * (1.0 / 4.0) - 5.0 / 8.0) + 115.0 / 192.0);
	}
	else if (Argument < 1.5) {
/**/	DEBUG_WRITE_LEAVING(Bspline04, "Done")
		return(Argument * (Argument * (Argument * (5.0 / 6.0 - Argument * (1.0 / 6.0))
			- 5.0 / 4.0) + 5.0 / 24.0) + 55.0 / 96.0);
	}
	else if (Argument < 2.5) {
		Argument -= 5.0 / 2.0;
		Argument *= Argument;
/**/	DEBUG_WRITE_LEAVING(Bspline04, "Done")
		return(Argument * Argument * (1.0 / 24.0));
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Bspline04, "Done")
		return(0.0);
	}
} /* end Bspline04 */

/*--------------------------------------------------------------------------*/
extern double	Bspline05
				(
					double	Argument			/* input */
				)

/* returns the value of a Basic spline function of degree 5 (order 6) evaluated at Argument */

{ /* begin Bspline05 */

	double	a;

/**/DEBUG_WRITE_ENTERING(Bspline05,
/**/	"About to compute a B-spline of degree 5")

	Argument = fabs(Argument);
	if (Argument < 1.0) {
		a = Argument * Argument;
/**/	DEBUG_WRITE_LEAVING(Bspline05, "Done")
		return(a * (a * (1.0 / 4.0 - Argument * (1.0 / 12.0)) - 1.0 / 2.0)
			+ 11.0 / 20.0);
	}
	else if (Argument < 2.0) {
/**/	DEBUG_WRITE_LEAVING(Bspline05, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument
			* (1.0 / 24.0) - 3.0 / 8.0) + 5.0 / 4.0) - 7.0 / 4.0) + 5.0 / 8.0) + 17.0 / 40.0);
	}
	else if (Argument < 3.0) {
		a = 3.0 - Argument;
		Argument = a * a;
/**/	DEBUG_WRITE_LEAVING(Bspline05, "Done")
		return(a * Argument * Argument * (1.0 / 120.0));
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Bspline05, "Done")
		return(0.0);
	}
} /* end Bspline05 */

/*--------------------------------------------------------------------------*/
extern double	Bspline06
				(
					double	Argument			/* input */
				)

/* returns the value of a Basic spline function of degree 6 (order 7) evaluated at Argument */

{ /* begin Bspline06 */

/**/DEBUG_WRITE_ENTERING(Bspline06,
/**/	"About to compute a B-spline of degree 6")

	Argument = fabs(Argument);
	if (Argument < 0.5) {
		Argument *= Argument;
/**/	DEBUG_WRITE_LEAVING(Bspline06, "Done")
		return(Argument * (Argument * (7.0 / 48.0 - Argument
			* (1.0 / 36.0)) - 77.0 / 192.0) + 5887.0 / 11520.0);
	}
	else if (Argument < 1.5) {
/**/	DEBUG_WRITE_LEAVING(Bspline06, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (1.0 / 48.0) - 7.0 / 48.0) + 21.0 / 64.0) - 35.0 / 288.0) - 91.0 / 256.0)
			- 7.0 / 768.0) + 7861.0 / 15360.0);
	}
	else if (Argument < 2.5) {
/**/	DEBUG_WRITE_LEAVING(Bspline06, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument
			* (7.0 / 60.0 - Argument * (1.0 / 120.0)) - 21.0 / 32.0) + 133.0 / 72.0)
			- 329.0 / 128.0) + 1267.0 / 960.0) + 1379.0 / 7680.0);
	}
	else if (Argument < 3.5) {
		Argument -= 7.0 / 2.0;
		Argument *= Argument * Argument;
/**/	DEBUG_WRITE_LEAVING(Bspline06, "Done")
		return(Argument * Argument * (1.0 / 720.0));
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Bspline06, "Done")
		return(0.0);
	}
} /* end Bspline06 */

/*--------------------------------------------------------------------------*/
extern double	Bspline07
				(
					double	Argument			/* input */
				)

/* returns the value of a Basic spline function of degree 7 (order 8) evaluated at Argument */

{ /* begin Bspline07 */

	double	a;

/**/DEBUG_WRITE_ENTERING(Bspline07,
/**/	"About to compute a B-spline of degree 7")

	Argument = fabs(Argument);
	if (Argument < 1.0) {
		a = Argument * Argument;
/**/	DEBUG_WRITE_LEAVING(Bspline07, "Done")
		return(a * (a * (a * (Argument * (1.0 / 144.0) - 1.0 / 36.0) + 1.0 / 9.0)
			- 1.0 / 3.0) + 151.0 / 315.0);
	}
	else if (Argument < 2.0) {
/**/	DEBUG_WRITE_LEAVING(Bspline07, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (1.0 / 20.0 - Argument * (1.0 / 240.0)) - 7.0 / 30.0) + 1.0 / 2.0)
			- 7.0 / 18.0) - 1.0 / 10.0) - 7.0 / 90.0) + 103.0 / 210.0);
	}
	else if (Argument < 3.0) {
/**/	DEBUG_WRITE_LEAVING(Bspline07, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (Argument * (1.0 / 720.0) - 1.0 / 36.0) + 7.0 / 30.0) - 19.0 / 18.0)
			+ 49.0 / 18.0) - 23.0 / 6.0) + 217.0 / 90.0) - 139.0 / 630.0);
	}
	else if (Argument < 4.0) {
		a = 4.0 - Argument;
		Argument = a * a * a;
/**/	DEBUG_WRITE_LEAVING(Bspline07, "Done")
		return(Argument * Argument * a * (1.0 / 5040.0));
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Bspline07, "Done")
		return(0.0);
	}
} /* end Bspline07 */

/*--------------------------------------------------------------------------*/
extern double	Bspline08
				(
					double	Argument			/* input */
				)

/* returns the value of a Basic spline function of degree 8 (order 9) evaluated at Argument */

{ /* begin Bspline08 */

/**/DEBUG_WRITE_ENTERING(Bspline08,
/**/	"About to compute a B-spline of degree 8")

	Argument = fabs(Argument);
	if (Argument < 0.5) {
		Argument *= Argument;
/**/	DEBUG_WRITE_LEAVING(Bspline08, "Done")
		return(Argument * (Argument * (Argument * (Argument * (1.0 / 576.0) - 1.0 / 64.0)
			+ 43.0 / 512.0) - 289.0 / 1024.0) + 259723.0 / 573440.0);
	}
	else if (Argument < 1.5) {
/**/	DEBUG_WRITE_LEAVING(Bspline08, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (Argument * (1.0 / 80.0 - Argument * (1.0 / 720.0)) - 3.0 / 80.0)
			+ 7.0 / 320.0) + 9.0 / 128.0) + 7.0 / 1280.0) - 363.0 / 1280.0)
			+ 1.0 / 5120.0) + 64929.0 / 143360.0);
	}
	else if (Argument < 2.5) {
/**/	DEBUG_WRITE_LEAVING(Bspline08, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (Argument * (Argument * (1.0 / 1440.0) - 1.0 / 80.0) + 3.0 / 32.0)
			- 119.0 / 320.0) + 207.0 / 256.0) - 1127.0 / 1280.0) + 195.0 / 512.0)
			- 1457.0 / 5120.0) + 145167.0 / 286720.0);
	}
	else if (Argument < 3.5) {
/**/	DEBUG_WRITE_LEAVING(Bspline08, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (Argument * (3.0 / 560.0 - Argument * (1.0 / 5040.0)) - 5.0 / 80.0)
			+ 131.0 / 320.0) - 209.0 / 128.0) + 5123.0 / 1280.0) - 1465.0 / 256.0)
			+ 146051.0 / 35840.0) - 122729.0 / 143360.0);
	}
	else if (Argument < 4.5) {
		Argument -= 9.0 / 2.0;
		Argument *= Argument;
		Argument *= Argument;
/**/	DEBUG_WRITE_LEAVING(Bspline08, "Done")
		return(Argument * Argument * (1.0 / 40320.0));
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Bspline08, "Done")
		return(0.0);
	}
} /* end Bspline08 */

/*--------------------------------------------------------------------------*/
extern double	Bspline09
				(
					double	Argument			/* input */
				)

/* returns the value of a Basic spline function of degree 9 (order 10) evaluated at Argument */

{ /* begin Bspline09 */

	double	a;

/**/DEBUG_WRITE_ENTERING(Bspline09,
/**/	"About to compute a B-spline of degree 9")

	Argument = fabs(Argument);
	if (Argument < 1.0) {
		a = Argument * Argument;
/**/	DEBUG_WRITE_LEAVING(Bspline09, "Done")
		return(a * (a * (a * (a * (1.0 / 576.0 - Argument * (1.0 / 2880.0))
			- 5.0 / 432.0) + 19.0 / 288.0) - 35.0 / 144.0) + 15619.0 / 36288.0);
	}
	else if (Argument < 2.0) {
/**/	DEBUG_WRITE_LEAVING(Bspline09, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (Argument * (Argument * (Argument * (1.0 / 4320.0) - 1.0 / 288.0)
			+ 1.0 / 48.0) - 13.0 / 216.0) + 7.0 / 96.0) - 1.0 / 144.0) + 7.0 / 144.0)
			- 19.0 / 72.0) + 1.0 / 192.0) + 7799.0 / 18144.0);
	}
	else if (Argument < 3.0) {
/**/	DEBUG_WRITE_LEAVING(Bspline09, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (Argument * (Argument * (5.0 / 2016.0 - Argument * (1.0 / 10080.0))
			- 3.0 / 112.0) + 35.0 / 216.0) - 19.0 / 32.0) + 191.0 / 144.0) - 83.0 / 48.0)
			+ 635.0 / 504.0) - 339.0 / 448.0) + 1553.0 / 2592.0);
	}
	else if (Argument < 4.0) {
/**/	DEBUG_WRITE_LEAVING(Bspline09, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (Argument * (Argument * (Argument * (1.0 / 40320.0) - 1.0 / 1152.0)
			+ 3.0 / 224.0) - 103.0 / 864.0) + 43.0 / 64.0) - 1423.0 / 576.0)
			+ 563.0 / 96.0) - 2449.0 / 288.0) + 5883.0 / 896.0) - 133663.0 / 72576.0);
	}
	else if (Argument < 5.0) {
		a = 5.0 - Argument;
		Argument = a * a;
		Argument *= Argument;
/**/	DEBUG_WRITE_LEAVING(Bspline09, "Done")
		return(Argument * Argument * a * (1.0 / 362880.0));
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Bspline09, "Done")
		return(0.0);
	}
} /* end Bspline09 */

/*--------------------------------------------------------------------------*/
extern double	Bspline10
				(
					double	Argument			/* input */
				)

/* returns the value of a Basic spline function of degree 10 (order 11) evaluated at Argument */

{ /* begin Bspline10 */

	double	a;

/**/DEBUG_WRITE_ENTERING(Bspline10,
/**/	"About to compute a B-spline of degree 10")

	Argument = fabs(Argument);
	if (Argument < 0.5) {
		Argument *= Argument;
/**/	DEBUG_WRITE_LEAVING(Bspline10, "Done")
		return(Argument * (Argument * (Argument * (Argument * (11.0 / 11520.0
			- Argument * (1.0 / 14400.0)) - 583.0 / 69120.0) + 14597.0 / 276480.0)
			- 156409.0 / 737280.0) + 381773117.0 / 928972800.0);
	}
	else if (Argument < 1.5) {
/**/	DEBUG_WRITE_LEAVING(Bspline10, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (Argument * (Argument * (Argument * (Argument * (1.0 / 17280.0)
			- 11.0 / 17280.0) + 11.0 / 4608.0) - 11.0 / 5760.0) - 187.0 / 27648.0)
			- 77.0 / 76800.0) + 5885.0 / 110592.0) - 11.0 / 92160.0)
			- 62557.0 / 294912.0) - 11.0 / 4423680.0) + 152709293.0 / 371589120.0);
	}
	else if (Argument < 2.5) {
/**/	DEBUG_WRITE_LEAVING(Bspline10, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (Argument * (Argument * (Argument * (11.0 / 15120.0 - Argument
			* (1.0 / 30240.0)) - 55.0 / 8064.0) + 11.0 / 315.0) - 715.0 / 6912.0)
			+ 1661.0 / 9600.0) - 4543.0 / 27648.0) + 7513.0 / 40320.0)
			- 163603.0 / 516096.0) + 135311.0 / 3870720.0) + 37690169.0 / 92897280.0);
	}
	else if (Argument < 3.5) {
/**/	DEBUG_WRITE_LEAVING(Bspline10, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (Argument * (Argument * (Argument * (Argument * (1.0 / 80640.0)
			- 11.0 / 26880.0) + 55.0 / 9216.0) - 451.0 / 8960.0) + 14905.0 / 55296.0)
			- 48433.0 / 51200.0) + 479281.0 / 221184.0) - 1352153.0 / 430080.0)
			+ 1654543.0 / 589824.0) - 11695211.0 / 6881280.0) + 623786977.0 / 743178240.0);
	}
	else if (Argument < 4.5) {
/**/	DEBUG_WRITE_LEAVING(Bspline10, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (Argument * (Argument * (Argument * (11.0 / 90720.0 - Argument
			* (1.0 / 362880.0)) - 11.0 / 4608.0) + 209.0 / 7560.0) - 5753.0 / 27648.0)
			+ 61061.0 / 57600.0) - 407429.0 / 110592.0) + 2070343.0 / 241920.0)
			- 3702215.0 / 294912.0) + 237959711.0 / 23224320.0) - 1241720381.0 / 371589120.0);
	}
	else if (Argument < 5.5) {
		Argument -= 11.0 / 2.0;
		a = Argument * Argument;
		a *= Argument * a;
/**/	DEBUG_WRITE_LEAVING(Bspline10, "Done")
		return(a * a * (1.0 / 3628800.0));
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Bspline10, "Done")
		return(0.0);
	}
} /* end Bspline10 */

/*--------------------------------------------------------------------------*/
extern double	Bspline11
				(
					double	Argument			/* input */
				)

/* returns the value of a Basic spline function of degree 11 (order 12) evaluated at Argument */

{ /* begin Bspline11 */

	double	a;

/**/DEBUG_WRITE_ENTERING(Bspline11,
/**/	"About to compute a B-spline of degree 11")

	Argument = fabs(Argument);
	if (Argument < 1.0) {
		a = Argument * Argument;
/**/	DEBUG_WRITE_LEAVING(Bspline11, "Done")
		return(a * (a * (a * (a * (a * (Argument * (1.0 / 86400.0) - 1.0 / 14400.0)
			+ 1.0 / 1440.0) - 23.0 / 3600.0) + 31.0 / 720.0) - 809.0 / 4320.0)
			+ 655177.0 / 1663200.0);
	}
	else if (Argument < 2.0) {
/**/	DEBUG_WRITE_LEAVING(Bspline11, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (Argument * (Argument * (Argument * (Argument * (1.0 / 6720.0
			- Argument * (1.0 / 120960.0)) - 11.0 / 10080.0) + 1.0 / 252.0)
			- 11.0 / 1680.0) + 1.0 / 360.0) - 11.0 / 1200.0) + 25.0 / 504.0)
			- 11.0 / 3360.0) - 563.0 / 3024.0) - 11.0 / 50400.0) + 65521.0 / 166320.0);
	}
	else if (Argument < 3.0) {
/**/	DEBUG_WRITE_LEAVING(Bspline11, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (Argument * (Argument * (Argument * (Argument * (Argument
			* (1.0 / 241920.0) - 1.0 / 8064.0) + 11.0 / 6720.0) - 25.0 / 2016.0)
			+ 33.0 / 560.0) - 13.0 / 72.0) + 143.0 / 400.0) - 239.0 / 504.0)
			+ 583.0 / 1120.0) - 1619.0 / 3024.0) + 781.0 / 5600.0) + 61297.0 / 166320.0);
	}
	else if (Argument < 4.0) {
/**/	DEBUG_WRITE_LEAVING(Bspline11, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (Argument * (Argument * (Argument * (Argument * (1.0 / 17280.0
			- Argument * (1.0 / 725760.0)) - 11.0 / 10080.0) + 7.0 / 576.0)
			- 99.0 / 1120.0) + 631.0 / 1440.0) - 1199.0 / 800.0) + 1009.0 / 288.0)
			- 12199.0 / 2240.0) + 9385.0 / 1728.0) - 38533.0 / 11200.0)
			+ 894727.0 / 665280.0);
	}
	else if (Argument < 5.0) {
/**/	DEBUG_WRITE_LEAVING(Bspline11, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (Argument * (Argument * (Argument * (Argument * (Argument
			* (1.0 / 3628800.0) - 1.0 / 67200.0) + 11.0 / 30240.0)
			- 107.0 / 20160.0) + 517.0 / 10080.0) - 2477.0 / 7200.0)
			+ 11737.0 / 7200.0) - 54797.0 / 10080.0) + 250657.0 / 20160.0)
			- 1113317.0 / 60480.0) + 4726777.0 / 302400.0) - 18595037.0 / 3326400.0);
	}
	else if (Argument < 6.0) {
		a = 6.0 - Argument;
		Argument = a * a;
		Argument *= Argument * a;
/**/	DEBUG_WRITE_LEAVING(Bspline11, "Done")
		return(Argument * Argument * a * (1.0 / 39916800.0));
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Bspline11, "Done")
		return(0.0);
	}
} /* end Bspline11 */

/*--------------------------------------------------------------------------*/
extern int		BsplineArray02
				(
					double	Argument,			/* fractional input */
					double	*b2_minus1,			/* 1st returned coefficient */
					double	*b2_plus0,			/* 2nd returned coefficient */
					double	*b2_plus1			/* 3rd returned coefficient */
				)

/* returns 3 values for a Basic spline function of degree 2 (order 3) */
/* evaluation is performed at {Argument - 1.0, Argument, Argument + 1.0} */
/* Argument must be in [-0.5, 0.5] */
/* computational load: */
/* 3 indirections */
/* 4(double)assignments, 4 (double)additions, 2 (double)multiplications */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin BsplineArray02 */

	double	a = Argument - 0.5;
	int		Status = !ERROR;

/**/DEBUG_CHECK_RANGE_DOUBLE(BsplineArray02, Argument, 0.0, 1.0, Status,
/**/	"Invalid argument (should be in [-0.5, 0.5])")
/**/DEBUG_CHECK_NULL_POINTER(BsplineArray02, b2_minus1, Status,
/**/	"No 1st output coefficient")
/**/DEBUG_CHECK_NULL_POINTER(BsplineArray02, b2_plus0, Status,
/**/	"No 2nd output coefficient")
/**/DEBUG_CHECK_NULL_POINTER(BsplineArray02, b2_plus1, Status,
/**/	"No 3rd output coefficient")
/**/DEBUG_RETURN_ON_ERROR(BsplineArray02, Status)
/**/DEBUG_WRITE_ENTERING(BsplineArray02,
/**/	"About to compute an array of B-splines of degree 2")

	*b2_plus1 = 0.5 * a * a;
	*b2_minus1 = *b2_plus1 + Argument;
	*b2_plus0 = 1.0 - *b2_plus1 - *b2_minus1;

/**/DEBUG_WRITE_LEAVING(BsplineArray02, "Done")
	return(Status);
} /* end BsplineArray02 */

/*--------------------------------------------------------------------------*/
extern int		BsplineArray03
				(
					double	Argument,			/* fractional input */
					double	*b3_minus2,			/* 1st returned coefficient */
					double	*b3_minus1,			/* 2nd returned coefficient */
					double	*b3_plus0,			/* 3rd returned coefficient */
					double	*b3_plus1			/* 4th returned coefficient */
				)

/* returns 4 values for a Basic spline function of degree 3 (order 4) */
/* evaluation is performed at {Argument - 2.0, Argument - 1.0, Argument, Argument + 1.0} */
/* Argument must be in [0.0, 1.0] */
/* computational load: */
/* 7 indirections */
/* 7(double)assignments, 6 (double)additions, 8 (double)multiplications */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin BsplineArray03 */

	double	a1 = 1.0 - Argument, a2 = Argument * Argument;
	int		Status = !ERROR;

/**/DEBUG_CHECK_RANGE_DOUBLE(BsplineArray03, Argument, 0.0, 1.0, Status,
/**/	"Invalid argument (should be in [0.0, 1.0])")
/**/DEBUG_CHECK_NULL_POINTER(BsplineArray03, b3_minus2, Status,
/**/	"No 1st output coefficient")
/**/DEBUG_CHECK_NULL_POINTER(BsplineArray03, b3_minus1, Status,
/**/	"No 2nd output coefficient")
/**/DEBUG_CHECK_NULL_POINTER(BsplineArray03, b3_plus0, Status,
/**/	"No 3rd output coefficient")
/**/DEBUG_CHECK_NULL_POINTER(BsplineArray03, b3_plus1, Status,
/**/	"No 4th output coefficient")
/**/DEBUG_RETURN_ON_ERROR(BsplineArray03, Status)
/**/DEBUG_WRITE_ENTERING(BsplineArray03,
/**/	"About to compute an array of B-splines of degree 3")

	*b3_plus1 = (1.0 / 6.0) * a1 * a1 * a1;
	*b3_plus0 = 2.0 / 3.0 - (1.0 / 2.0) * a2 * ++a1;
	*b3_minus2 = (1.0 / 6.0) * a2 * Argument;
	*b3_minus1 = 1.0 - *b3_plus1 - *b3_plus0 - *b3_minus2;

/**/DEBUG_WRITE_LEAVING(BsplineArray03, "Done")
	return(Status);
} /* end BsplineArray03 */

/*--------------------------------------------------------------------------*/
extern double	Dodgson
				(
					double	Argument			/* input */
				)

/* returns the value of a Dodgson kernel evaluated at Argument (order 2) */

{ /* begin Dodgson */

/**/DEBUG_WRITE_ENTERING(Dodgson,
/**/	"About to compute a Dodgson kernel")

	Argument = fabs(Argument);
	if (Argument < 0.5) {
/**/	DEBUG_WRITE_LEAVING(Dodgson, "Done")
		return(1.0 - 2.0 * Argument * Argument);
	}
	else if (Argument < 1.5) {
		Argument -= 5.0 / 4.0;
/**/	DEBUG_WRITE_LEAVING(Dodgson, "Done")
		return(Argument * Argument - 1.0 / 16.0);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Dodgson, "Done")
		return(0.0);
	}
} /* end Dodgson */

/*--------------------------------------------------------------------------*/
extern int		DodgsonArray
				(
					double	Argument,			/* fractional input */
					double	*d_minus1,			/* 1st returned coefficient */
					double	*d_plus0,			/* 2nd returned coefficient */
					double	*d_plus1			/* 3rd returned coefficient */
				)

/* returns 3 values for a Dodgson kernel (order 2) */
/* evaluation is performed at {Argument - 1.0, Argument, Argument + 1.0} */
/* Argument must be in [-0.5, 0.5] */
/* computational load: */
/* 5 indirections */
/* 4(double)assignments, 4 (double)additions, 3 (double)multiplications */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin DodgsonArray */

	double	a2 = Argument * Argument;
	int		Status = !ERROR;

/**/DEBUG_CHECK_RANGE_DOUBLE(DodgsonArray, Argument, -0.5, 0.5, Status,
/**/	"Invalid argument (should be in [0.0, 1.0])")
/**/DEBUG_CHECK_NULL_POINTER(DodgsonArray, d_minus1, Status,
/**/	"No 1st output coefficient")
/**/DEBUG_CHECK_NULL_POINTER(DodgsonArray, d_plus0, Status,
/**/	"No 2nd output coefficient")
/**/DEBUG_CHECK_NULL_POINTER(DodgsonArray, d_plus1, Status,
/**/	"No 3rd output coefficient")
/**/DEBUG_RETURN_ON_ERROR(DodgsonArray, Status)
/**/DEBUG_WRITE_ENTERING(DodgsonArray,
/**/	"About to compute an array of Dodgson")

	*d_plus1 = a2 - (1.0 / 2.0) * Argument;
	*d_plus0 = 1.0 - 2.0 * a2;
	*d_minus1 = 1.0 - *d_plus1 - *d_plus0;

/**/DEBUG_WRITE_LEAVING(DodgsonArray, "Done")
	return(Status);
} /* end DodgsonArray */

/*--------------------------------------------------------------------------*/
extern double	German04
				(
					double	Argument			/* input */
				)

/* returns the value of the quartic German kernel (order 5) evaluated at Argument */

{ /* begin German04 */

/**/DEBUG_WRITE_ENTERING(German04,
/**/	"About to compute a quartic German kernel")

	Argument = fabs(Argument);
	if (Argument < 1.0) {
/**/	DEBUG_WRITE_LEAVING(German04, "Done")
		return(Argument * Argument * (Argument * (Argument * (1.0 / 24.0) + 185.0 / 144.0)
			- 335.0 / 144.0) + 1.0);
	}
	else if (Argument < 2.0) {
/**/	DEBUG_WRITE_LEAVING(German04, "Done")
		return(Argument * (Argument * (Argument * (-29.0 / 144.0 - Argument * (1.0 / 16.0))
			+ 17.0 / 8.0) - 145.0 / 36.0) + 13.0 / 6.0);
	}
	else if (Argument < 3.0) {
/**/	DEBUG_WRITE_LEAVING(German04, "Done")
		return(Argument * (Argument * (Argument * (Argument * (1.0 / 48.0) - 7.0 / 48.0)
			+ 5.0 / 18.0) + 1.0 / 36.0) - 1.0 / 3.0);
	}
	else if (Argument < 4.0) {
/**/	DEBUG_WRITE_LEAVING(German04, "Done")
		return(Argument * (Argument * (Argument * (1.0 / 144.0) - 11.0 / 144.0)
			+ 5.0 / 18.0) - 1.0 / 3.0);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(German04, "Done")
		return(0.0);
	}
} /* end German04 */

/*--------------------------------------------------------------------------*/
extern double	KeysOptimal
				(
					double	Argument			/* input */
				)

/* returns the value of the cubic Keys optimal kernel (order 3) evaluated at Argument */

{ /* begin KeysOptimal */

/**/DEBUG_WRITE_ENTERING(KeysOptimal,
/**/	"About to compute an optimal cubic Keys kernel")

	Argument = fabs(Argument);
	if (Argument < 1.0) {
/**/	DEBUG_WRITE_LEAVING(KeysOptimal, "Done")
		return(Argument * Argument * (Argument * (3.0 / 2.0) - 5.0 / 2.0) + 1.0);
	}
	else if (Argument < 2.0) {
		Argument -= 2.0;
/**/	DEBUG_WRITE_LEAVING(KeysOptimal, "Done")
		return(Argument * Argument * (Argument + 1.0) * (-1.0 / 2.0));
	}
	else {
/**/	DEBUG_WRITE_LEAVING(KeysOptimal, "Done")
		return(0.0);
	}
} /* end KeysOptimal */

/*--------------------------------------------------------------------------*/
extern double	Keys
				(
					double	Argument,			/* input */
					double	a					/* tuning parameter */
				)

/* returns the value of the cubic Keys kernel evaluated at Argument */

{ /* begin Keys */

/**/DEBUG_WRITE_ENTERING(Keys,
/**/	"About to compute a cubic Keys kernel")

	Argument = fabs(Argument);
	if (Argument < 1.0) {
/**/	DEBUG_WRITE_LEAVING(Keys, "Done")
		return(Argument * Argument * (Argument * (a + 2.0) - (a + 3.0)) + 1.0);
	}
	else if (Argument < 2.0) {
		Argument -= 2.0;
/**/	DEBUG_WRITE_LEAVING(Keys, "Done")
		return(a * Argument * Argument * (Argument + 1.0));
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Keys, "Done")
		return(0.0);
	}
} /* end Keys */

/*--------------------------------------------------------------------------*/
extern int		KeysOptimalArray
				(
					double	Argument,			/* fractional input */
					double	*k_minus2,			/* 1st returned coefficient */
					double	*k_minus1,			/* 2nd returned coefficient */
					double	*k_plus0,			/* 3rd returned coefficient */
					double	*k_plus1			/* 4th returned coefficient */
				)

/* returns 4 values for a Keys kernel) */
/* evaluation is performed at {Argument - 2.0, Argument - 1.0, Argument, Argument + 1.0} */
/* Argument must be in [0.0, 1.0] */
/* computational load: */
/* 7 indirections */
/* 6(double)assignments, 6 (double)additions, 8 (double)multiplications */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin KeysOptimalArray */

	double	a1 = 1.0 - Argument, a2 = Argument * Argument;
	int		Status = !ERROR;

/**/DEBUG_CHECK_RANGE_DOUBLE(KeysOptimalArray, Argument, 0.0, 1.0, Status,
/**/	"Invalid argument (should be in [0.0, 1.0])")
/**/DEBUG_CHECK_NULL_POINTER(KeysOptimalArray, k_minus2, Status,
/**/	"No 1st output coefficient")
/**/DEBUG_CHECK_NULL_POINTER(KeysOptimalArray, k_minus1, Status,
/**/	"No 2nd output coefficient")
/**/DEBUG_CHECK_NULL_POINTER(KeysOptimalArray, k_plus0, Status,
/**/	"No 3rd output coefficient")
/**/DEBUG_CHECK_NULL_POINTER(KeysOptimalArray, k_plus1, Status,
/**/	"No 4th output coefficient")
/**/DEBUG_RETURN_ON_ERROR(KeysOptimalArray, Status)
/**/DEBUG_WRITE_ENTERING(KeysOptimalArray,
/**/	"About to compute an array of Keys")

	*k_minus2 = a1 * a2 * (-1.0 / 2.0);
	*k_plus1 = Argument * a1 * a1 * (-1.0 / 2.0);
	*k_plus0 = a2 * (Argument * (3.0 / 2.0) - 5.0 / 2.0) + 1.0;
	*k_minus1 = 1.0 - *k_minus2 - *k_plus1 - *k_plus0;

/**/DEBUG_WRITE_LEAVING(KeysOptimalArray, "Done")
	return(Status);
} /* end KeysOptimalArray */

/*--------------------------------------------------------------------------*/
extern double	Meijering05
				(
					double	Argument			/* input */
				)

/* returns the value of a Meijering function of degree 5 (order ?) evaluated at Argument */

{ /* begin Meijering05 */

	double	a;

/**/DEBUG_WRITE_ENTERING(Meijering05,
/**/	"About to compute a Meijering of degree 5")

	Argument = fabs(Argument);
	if (Argument < 1.0) {
		a = Argument * Argument;
/**/	DEBUG_WRITE_LEAVING(Meijering05, "Done")
		return(a * (a * (63.0 / 32.0 - Argument * (27.0 / 32.0)) - 17.0 / 8.0) + 1.0);
	}
	else if (Argument < 2.0) {
/**/	DEBUG_WRITE_LEAVING(Meijering05, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument
			* (13.0 / 64.0) - 21.0 / 16.0) + 85.0 / 32.0) - 7.0 / 8.0) - 165.0 / 64.0)
			+ 61.0 / 32.0);
	}
	else if (Argument < 3.0) {
		a = 3.0 - Argument;
		a *= a;
		a *= a;
/**/	DEBUG_WRITE_LEAVING(Meijering05, "Done")
		return(a * (Argument - 2.0) * (3.0 / 64.0));
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Meijering05, "Done")
		return(0.0);
	}
} /* end Meijering05 */

/*--------------------------------------------------------------------------*/
extern double	Meijering07
				(
					double	Argument			/* input */
				)

/* returns the value of a Meijering function of degree 7 (order ?) evaluated at Argument */

{ /* begin Meijering07 */

	double	a;

/**/DEBUG_WRITE_ENTERING(Meijering07,
/**/	"About to compute a Meijering of degree 7")

	Argument = fabs(Argument);
	if (Argument < 1.0) {
		a = Argument * Argument;
/**/	DEBUG_WRITE_LEAVING(Meijering07, "Done")
		return(a * (a * (a * (Argument * (22013.0 / 83232.0) - 22039.0 / 27744.0)
			+ 16775.0 / 10404.0) - 3611.0 / 1734.0) + 1.0);
	}
	else if (Argument < 2.0) {
/**/	DEBUG_WRITE_LEAVING(Meijering07, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (-477.0 / 9248.0 - Argument * (1127.0 / 83232.0)) + 9583.0 / 6936.0)
			- 5815.0 / 867.0) + 33145.0 / 2312.0) - 243.0 / 17.0) + 109739.0 / 20808.0)
			+ 259.0 / 3468.0);
	}
	else if (Argument < 3.0) {
/**/	DEBUG_WRITE_LEAVING(Meijering07, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (79.0 / 96.0 - Argument * (403.0 / 9248.0)) - 45703.0 / 6936.0)
			+ 301115.0 / 10404.0) - 521045.0 / 6936.0) + 99572.0 / 867.0)
			- 1983989.0 / 20808.0) + 114355.0 / 3468.0);
	}
	else if (Argument < 4.0) {
/**/	DEBUG_WRITE_LEAVING(Meijering07, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (213.0 / 9248.0 - Argument * (71.0 / 83232.0)) - 923.0 / 3468.0)
			+ 8875.0 / 5202.0) - 5680.0 / 867.0) + 13064.0 / 867.0) - 49984.0 / 2601.0)
			+ 9088.0 / 867.0);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Meijering07, "Done")
		return(0.0);
	}
} /* end Meijering07 */

/*--------------------------------------------------------------------------*/
extern int		Omoms
				(
					long	Degree,				/* degree */
					double	Argument,			/* input */
					double	*Result				/* output */
				)

/* computes a Blu optimum function */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin Omoms */

	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(Omoms, Result, Status,
/**/	"No output data")
/**/DEBUG_RETURN_ON_ERROR(Omoms, Status)
/**/DEBUG_WRITE_ENTERING(Omoms,
/**/	"About to compute a generic oMoms")

	switch (Degree) {
		case 0L:
			*Result = Omoms00(Argument);
			break;
		case 1L:
			*Result = Omoms01(Argument);
			break;
		case 2L:
			*Result = Omoms02(Argument);
			break;
		case 3L:
			*Result = Omoms03(Argument);
			break;
		case 4L:
			*Result = Omoms04(Argument);
			break;
		case 5L:
			*Result = Omoms05(Argument);
			break;
		case 6L:
			*Result = Omoms06(Argument);
			break;
		case 7L:
			*Result = Omoms07(Argument);
			break;
		case 8L:
			*Result = Omoms08(Argument);
			break;
		case 9L:
			*Result = Omoms09(Argument);
			break;
		case 10L:
			*Result = Omoms10(Argument);
			break;
		case 11L:
			*Result = Omoms11(Argument);
			break;
		default:
			Status = ERROR;
			WRITE_ERROR(Omoms, "Invalid oMoms degree")
			break;
	}
/**/DEBUG_WRITE_LEAVING(Omoms, "Done")
	return(Status);
} /* end Omoms */

/*--------------------------------------------------------------------------*/
extern double	Omoms00
				(
					double	Argument			/* input */
				)

/* returns the value of a Blu optimum function of degree 0 (order 1) evaluated at Argument */

{ /* begin Omoms00 */

/**/DEBUG_WRITE_ENTERING(Omoms00,
/**/	"About to compute an oMoms of degree 0")

	Argument = fabs(Argument);
	if (Argument < 0.5) {
/**/	DEBUG_WRITE_LEAVING(Omoms00, "Done")
		return(1.0);
	}
	else if (Argument == 0.5) {
/**/	DEBUG_WRITE_LEAVING(Omoms00, "Done")
		return(0.5);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Omoms00, "Done")
		return(0.0);
	}
} /* end Omoms00 */

/*--------------------------------------------------------------------------*/
extern double	Omoms01
				(
					double	Argument			/* input */
				)

/* returns the value of a Blu optimum function of degree 1 (order 2) evaluated at Argument */

{ /* begin Omoms01 */

/**/DEBUG_WRITE_ENTERING(Omoms01,
/**/	"About to compute an oMoms of degree 1")

	Argument = fabs(Argument);
	if (Argument < 1.0) {
/**/	DEBUG_WRITE_LEAVING(Omoms01, "Done")
		return(1.0 - Argument);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Omoms01, "Done")
		return(0.0);
	}
} /* end Omoms01 */

/*--------------------------------------------------------------------------*/
extern double	Omoms02
				(
					double	Argument			/* input */
				)

/* returns the value of a Blu optimum function of degree 2 (order 3) evaluated at Argument */

{ /* begin Omoms02 */

/**/DEBUG_WRITE_ENTERING(Omoms02,
/**/	"About to compute an oMoms of degree 2")

	Argument = fabs(Argument);
	if (Argument < 0.5) {
/**/	DEBUG_WRITE_LEAVING(Omoms02, "Done")
		return(43.0 / 60.0 - Argument * Argument);
	}
	else if (Argument == 0.5) {
/**/	DEBUG_WRITE_LEAVING(Omoms02, "Done")
		return(59.0 / 120.0);
	}
	else if (Argument < 1.5) {
/**/	DEBUG_WRITE_LEAVING(Omoms02, "Done")
		return(Argument * (Argument * (1.0 / 2.0) - 3.0 / 2.0) + 137.0 / 120.0);
	}
	else if (Argument == 1.5) {
/**/	DEBUG_WRITE_LEAVING(Omoms02, "Done")
		return(1.0 / 120.0);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Omoms02, "Done")
		return(0.0);
	}
} /* end Omoms02 */

/*--------------------------------------------------------------------------*/
extern double	Omoms03
				(
					double	Argument			/* input */
				)

/* returns the value of a Blu optimum function of degree 3 (order 4) evaluated at Argument */

{ /* begin Omoms03 */

/**/DEBUG_WRITE_ENTERING(Omoms03,
/**/	"About to compute an oMoms of degree 3")

	Argument = fabs(Argument);
	if (Argument < 1.0) {
/**/	DEBUG_WRITE_LEAVING(Omoms03, "Done")
		return(Argument * (Argument * (Argument * (1.0 / 2.0) - 1.0) + 1.0 / 14.0)
			+ 13.0 / 21.0);
	}
	else if (Argument < 2.0) {
/**/	DEBUG_WRITE_LEAVING(Omoms03, "Done")
		return(Argument * (Argument * (1.0 - Argument * (1.0 / 6.0)) - 85.0 / 42.0)
			+ 29.0 / 21.0);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Omoms03, "Done")
		return(0.0);
	}
} /* end Omoms03 */

/*--------------------------------------------------------------------------*/
extern double	Omoms04
				(
					double	Argument			/* input */
				)

/* returns the value of a Blu optimum function of degree 4 (order 5) evaluated at Argument */

{ /* begin Omoms04 */

/**/DEBUG_WRITE_ENTERING(Omoms04,
/**/	"About to compute an oMoms of degree 4")

	Argument = fabs(Argument);
	if (Argument < 0.5) {
		Argument *= Argument;
/**/	DEBUG_WRITE_LEAVING(Omoms04, "Done")
		return(Argument * (Argument * (1.0 / 4.0) - 13.0 / 24.0) + 11383.0 / 20160.0);
	}
	else if (Argument == 0.5) {
/**/	DEBUG_WRITE_LEAVING(Omoms04, "Done")
		return(6721.0 / 15120.0);
	}
	else if (Argument < 1.5) {
/**/	DEBUG_WRITE_LEAVING(Omoms04, "Done")
		return(Argument * (Argument * (Argument * (5.0 / 6.0 - Argument * (1.0 / 6.0))
			- 47.0 / 36.0) + 25.0 / 72.0) + 15217.0 / 30240.0);
	}
	else if (Argument == 1.5) {
/**/	DEBUG_WRITE_LEAVING(Omoms04, "Done")
		return(559.0 / 10080.0);
	}
	else if (Argument < 2.5) {
/**/	DEBUG_WRITE_LEAVING(Omoms04, "Done")
		return(Argument * (Argument * (Argument * (Argument * (1.0 / 24.0) - 5.0 / 12.0)
			+ 227.0 / 144.0) - 385.0 / 144.0) + 207383.0 / 120960.0);
	}
	else if (Argument == 2.5) {
/**/	DEBUG_WRITE_LEAVING(Omoms04, "Done")
		return(1.0 / 30240.0);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Omoms04, "Done")
		return(0.0);
	}
} /* end Omoms04 */

/*--------------------------------------------------------------------------*/
extern double	Omoms05
				(
					double	Argument			/* input */
				)

/* returns the value of a Blu optimum function of degree 5 (order 6) evaluated at Argument */

{ /* begin Omoms05 */

/**/DEBUG_WRITE_ENTERING(Omoms05,
/**/	"About to compute an oMoms of degree 5")

	Argument = fabs(Argument);
	if (Argument < 1.0) {
/**/	DEBUG_WRITE_LEAVING(Omoms05, "Done")
		return(Argument * (Argument * (Argument * (Argument * (1.0 / 4.0
			- Argument * (1.0 / 12.0)) - 5.0 / 99.0) - 9.0 / 22.0) - 1.0 / 792.0)
			+ 229.0 / 440.0);
	}
	else if (Argument < 2.0) {
/**/	DEBUG_WRITE_LEAVING(Omoms05, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (1.0 / 24.0)
			- 3.0 / 8.0) + 505.0 / 396.0) - 83.0 / 44.0) + 1351.0 / 1584.0) + 839.0 / 2640.0);
	}
	else if (Argument < 3.0) {
/**/	DEBUG_WRITE_LEAVING(Omoms05, "Done")
		return(Argument * (Argument * (Argument * (Argument * (1.0 / 8.0
			- Argument * (1.0 / 120.0)) - 299.0 / 396.0) + 101.0 / 44.0) - 27811.0 / 7920.0)
			+ 5707.0 / 2640.0);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Omoms05, "Done")
		return(0.0);
	}
} /* end Omoms05 */

/*--------------------------------------------------------------------------*/
extern double	Omoms06
				(
					double	Argument			/* input */
				)

/* returns the value of a Blu optimum function of degree 6 (order 7) evaluated at Argument */

{ /* begin Omoms06 */

/**/DEBUG_WRITE_ENTERING(Omoms06,
/**/	"About to compute an oMoms of degree 6")

	Argument = fabs(Argument);
	if (Argument < 0.5) {
		Argument *= Argument;
/**/	DEBUG_WRITE_LEAVING(Omoms06, "Done")
		return(Argument * (Argument * (223.0 / 1872.0 - Argument * (1.0 / 36.0))
			- 3173.0 / 9152.0) + 16810399.0 / 34594560.0);
	}
	else if (Argument == 0.5) {
/**/	DEBUG_WRITE_LEAVING(Omoms06, "Done")
		return(540559.0 / 1330560.0);
	}
	else if (Argument < 1.5) {
/**/	DEBUG_WRITE_LEAVING(Omoms06, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (1.0 / 48.0) - 7.0 / 48.0) + 869.0 / 2496.0) - 805.0 / 3744.0)
			- 8345.0 / 36608.0) - 11711.0 / 329472.0) + 22619117.0 / 46126080.0);
	}
	else if (Argument == 1.5) {
/**/	DEBUG_WRITE_LEAVING(Omoms06, "Done")
		return(174749.0 / 1921920.0);
	}
	else if (Argument < 2.5) {
/**/	DEBUG_WRITE_LEAVING(Omoms06, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (7.0 / 60.0
			- Argument * (1.0 / 120.0)) - 829.0 / 1248.0) + 1799.0 / 936.0)
			- 258383.0 / 91520.0) + 690851.0 / 411840.0) + 21359.0 / 1774080.0);
	}
	else if (Argument == 2.5) {
/**/	DEBUG_WRITE_LEAVING(Omoms06, "Done")
		return(4421.0 / 1572480.0);
	}
	else if (Argument < 3.5) {
/**/	DEBUG_WRITE_LEAVING(Omoms06, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (1.0 / 720.0) - 7.0 / 240.0) + 1921.0 / 7488.0) - 4529.0 / 3744.0)
			+ 590221.0 / 183040.0) - 7588511.0 / 1647360.0) + 34652699.0 / 12579840.0);
	}
	else if (Argument == 3.5) {
/**/	DEBUG_WRITE_LEAVING(Omoms06, "Done")
		return(1.0 / 17297280.0);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Omoms06, "Done")
		return(0.0);
	}
} /* end Omoms06 */

/*--------------------------------------------------------------------------*/
extern double	Omoms07
				(
					double	Argument			/* input */
				)

/* returns the value of a Blu optimum function of degree 7 (order 8) evaluated at Argument */

{ /* begin Omoms07 */

/**/DEBUG_WRITE_ENTERING(Omoms07,
/**/	"About to compute an oMoms of degree 7")

	Argument = fabs(Argument);
	if (Argument < 1.0) {
/**/	DEBUG_WRITE_LEAVING(Omoms07, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (Argument * (1.0 / 144.0) - 1.0 / 36.0) + 7.0 / 720.0) + 1.0 / 12.0)
			+ 7.0 / 5616.0) - 227.0 / 780.0) + 1.0 / 102960.0) + 247409.0 / 540540.0);
	}
	else if (Argument < 2.0) {
/**/	DEBUG_WRITE_LEAVING(Omoms07, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (1.0 / 20.0 - Argument * (1.0 / 240.0)) - 287.0 / 1200.0) + 11.0 / 20.0)
			- 567.0 / 1040.0) + 27.0 / 260.0) - 27721.0 / 171600.0) + 48687.0 / 100100.0);
	}
	else if (Argument < 3.0) {
/**/	DEBUG_WRITE_LEAVING(Omoms07, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (Argument * (1.0 / 720.0) - 1.0 / 36.0) + 847.0 / 3600.0) - 13.0 / 12.0)
			+ 16163.0 / 5616.0) - 1107.0 / 260.0) + 117277.0 / 39600.0)
			- 260329.0 / 540540.0);
	}
	else if (Argument < 4.0) {
/**/	DEBUG_WRITE_LEAVING(Omoms07, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (1.0 / 180.0 - Argument * (1.0 / 5040.0)) - 241.0 / 3600.0) + 9.0 / 20.0)
			- 51169.0 / 28080.0) + 3467.0 / 780.0) - 21787921.0 / 3603600.0)
			+ 9560851.0 / 2702700.0);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Omoms07, "Done")
		return(0.0);
	}
} /* end Omoms07 */

/*--------------------------------------------------------------------------*/
extern double	Omoms08
				(
					double	Argument			/* input */
				)

/* returns the value of a Blu optimum function of degree 8 (order 9) evaluated at Argument */

{ /* begin Omoms08 */

/**/DEBUG_WRITE_ENTERING(Omoms08,
/**/	"About to compute an oMoms of degree 8")

	Argument = fabs(Argument);
	if (Argument < 0.5) {
		Argument *= Argument;
/**/	DEBUG_WRITE_LEAVING(Omoms08, "Done")
		return(Argument * (Argument * (Argument * (Argument * (1.0 / 576.0)
			- 361.0 / 29376.0) + 5375.0 / 78336.0) - 507163.0 / 2036736.0)
			+ 49010505197.0 / 112916643840.0);
	}
	else if (Argument == 0.5) {
/**/	DEBUG_WRITE_LEAVING(Omoms08, "Done")
		return(255075379.0 / 678585600.0);
	}
	else if (Argument < 1.5) {
/**/	DEBUG_WRITE_LEAVING(Omoms08, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (Argument * (1.0 / 80.0 - Argument * (1.0 / 720.0)) - 295.0 / 7344.0)
			+ 217.0 / 5440.0) + 3049.0 / 97920.0) + 301.0 / 13056.0)
			- 656753.0 / 2545920.0) + 451.0 / 226304.0) + 3602106707.0 / 8302694400.0);
	}
	else if (Argument == 1.5) {
/**/	DEBUG_WRITE_LEAVING(Omoms08, "Done")
		return(508155113.0 / 4410806400.0);
	}
	else if (Argument < 2.5) {
/**/	DEBUG_WRITE_LEAVING(Omoms08, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (Argument * (Argument * (1.0 / 1440.0) - 1.0 / 80.0) + 6983.0 / 73440.0)
			- 2121.0 / 5440.0) + 177311.0 / 195840.0) - 14861.0 / 13056.0)
			+ 3676757.0 / 5091840.0) - 539519.0 / 1131520.0)
			+ 151654522541.0 / 282291609600.0);
	}
	else if (Argument == 2.5) {
/**/	DEBUG_WRITE_LEAVING(Omoms08, "Done")
		return(7778497.0 / 882161280.0);
	}
	else if (Argument < 3.5) {
/**/	DEBUG_WRITE_LEAVING(Omoms08, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (Argument * (3.0 / 560.0 - Argument * (1.0 / 5040.0)) - 2309.0 / 36720.0)
			+ 2269.0 / 5440.0) - 166193.0 / 97920.0) + 55937.0 / 13056.0)
			- 114066257.0 / 17821440.0) + 38899437.0 / 7920640.0)
			- 177624353243.0 / 141145804800.0);
	}
	else if (Argument == 3.5) {
/**/	DEBUG_WRITE_LEAVING(Omoms08, "Done")
		return(1462553.0 / 17643225600.0);
	}
	else if (Argument < 4.5) {
/**/	DEBUG_WRITE_LEAVING(Omoms08, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (Argument * (Argument * (1.0 / 40320.0) - 1.0 / 1120.0) + 829.0 / 58752.0)
			- 1391.0 / 10880.0) + 569033.0 / 783360.0) - 23065.0 / 8704.0)
			+ 864109733.0 / 142571520.0) - 25173235.0 / 3168256.0)
			+ 5160686670443.0 / 1129166438400.0);
	}
	else if (Argument == 4.5) {
/**/	DEBUG_WRITE_LEAVING(Omoms08, "Done")
		return(1.0 / 17643225600.0);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Omoms08, "Done")
		return(0.0);
	}
} /* end Omoms08 */

/*--------------------------------------------------------------------------*/
extern double	Omoms09
				(
					double	Argument			/* input */
				)

/* returns the value of a Blu optimum function of degree 9 (order 10) evaluated at Argument */

{ /* begin Omoms09 */

/**/DEBUG_WRITE_ENTERING(Omoms09,
/**/	"About to compute an oMoms of degree 9")

	Argument = fabs(Argument);
	if (Argument < 1.0) {
/**/	DEBUG_WRITE_LEAVING(Omoms09, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (Argument * (Argument * (1.0 / 576.0 - Argument * (1.0 / 2880.0))
			- 1.0 / 1140.0) - 67.0 / 8208.0) - 147.0 / 516800.0) + 30463.0 / 558144.0)
			- 1.0 / 77520.0) - 20129.0 / 93024.0) - 1.0 / 24186240.0)
			+ 630498007.0 / 1523733120.0);
	}
	else if (Argument < 2.0) {
/**/	DEBUG_WRITE_LEAVING(Omoms09, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (Argument * (Argument * (Argument * (1.0 / 4320.0) - 1.0 / 288.0)
			+ 293.0 / 13680.0) - 275.0 / 4104.0) + 3353.0 / 32300.0) - 20059.0 / 279072.0)
			+ 48619.0 / 465120.0) - 12685.0 / 46512.0) + 648571.0 / 36279360.0)
			+ 313316657.0 / 761866560.0);
	}
	else if (Argument < 3.0) {
/**/	DEBUG_WRITE_LEAVING(Omoms09, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (Argument * (Argument * (5.0 / 2016.0 - Argument * (1.0 / 10080.0))
			- 863.0 / 31920.0) + 685.0 / 4104.0) - 81823.0 / 129200.0) + 418073.0 / 279072.0)
			- 778483.0 / 361760.0) + 597197.0 / 325584.0) - 96512911.0 / 84651840.0)
			+ 75779819.0 / 108838080.0);
	}
	else if (Argument < 4.0) {
/**/	DEBUG_WRITE_LEAVING(Omoms09, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (Argument * (Argument * (Argument * (1.0 / 40320.0) - 1.0 / 1152.0)
			+ 859.0 / 63840.0) - 1985.0 / 16416.0) + 714871.0 / 1033600.0)
			- 2898295.0 / 1116288.0) + 13759447.0 / 2170560.0) - 1777747.0 / 186048.0)
			+ 2648715421.0 / 338607360.0) - 7480147879.0 / 3047466240.0);
	}
	else if (Argument < 5.0) {
/**/	DEBUG_WRITE_LEAVING(Omoms09, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (Argument * (Argument * (1.0 / 8064.0 - Argument * (1.0 / 362880.0))
			- 1429.0 / 574560.0) + 479.0 / 16416.0) - 684257.0 / 3100800.0)
			+ 1245313.0 / 1116288.0) - 73642277.0 / 19535040.0) + 10692427.0 / 1302336.0)
			- 31866240901.0 / 3047466240.0) + 3616276301.0 / 609493248.0);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Omoms09, "Done")
		return(0.0);
	}
} /* end Omoms09 */

/*--------------------------------------------------------------------------*/
extern double	Omoms10
				(
					double	Argument			/* input */
				)

/* returns the value of a Blu optimum function of degree 10 (order 11) evaluated at Argument */

{ /* begin Omoms10 */

/**/DEBUG_WRITE_ENTERING(Omoms10,
/**/	"About to compute an oMoms of degree 10")

	Argument = fabs(Argument);
	if (Argument < 0.5) {
		Argument *= Argument;
/**/	DEBUG_WRITE_LEAVING(Omoms10, "Done")
		return(Argument * (Argument * (Argument * (Argument * (59.0 / 80640.0
			- Argument * (1.0 / 14400.0)) - 43517.0 / 6566400.0) + 9214241.0 / 208373760.0)
			- 27203797.0 / 142884864.0) + 16999223061911.0 / 42908324659200.0);
	}
	else if (Argument == 0.5) {
/**/	DEBUG_WRITE_LEAVING(Omoms10, "Done")
		return(13852170803.0 / 39437798400.0);
	}
	else if (Argument < 1.5) {
/**/	DEBUG_WRITE_LEAVING(Omoms10, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (Argument * (Argument * (Argument * (Argument * (1.0 / 17280.0)
			- 11.0 / 17280.0) + 83.0 / 32256.0) - 143.0 / 40320.0) - 1667.0 / 875520.0)
			- 19393.0 / 4377600.0) + 3929617.0 / 83349504.0) - 119141.0 / 89303040.0)
			- 14076065.0 / 74088448.0) - 3539371.0 / 50009702400.0)
			+ 2266595190421.0 / 5721109954560.0);
	}
	else if (Argument == 1.5) {
/**/	DEBUG_WRITE_LEAVING(Omoms10, "Done")
		return(617888623367.0 / 4693098009600.0);
	}
	else if (Argument < 2.5) {
/**/	DEBUG_WRITE_LEAVING(Omoms10, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (Argument * (Argument * (Argument * (11.0 / 15120.0
			- Argument * (1.0 / 30240.0)) - 391.0 / 56448.0) + 649.0 / 17640.0)
			- 538399.0 / 4596480.0) + 865843.0 / 3830400.0) - 1934417.0 / 6945792.0)
			+ 222233.0 / 697680.0) - 1394728309.0 / 3500679168.0)
			+ 3548488471.0 / 43758489600.0) + 11469642367213.0 / 30035827261440.0);
	}
	else if (Argument == 2.5) {
/**/	DEBUG_WRITE_LEAVING(Omoms10, "Done")
		return(40687175.0 / 2453907456.0);
	}
	else if (Argument < 3.5) {
/**/	DEBUG_WRITE_LEAVING(Omoms10, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (Argument * (Argument * (Argument * (Argument * (1.0 / 80640.0)
			- 11.0 / 26880.0) + 2713.0 / 451584.0) - 3223.0 / 62720.0)
			+ 1478771.0 / 5253120.0) - 6958193.0 / 6809600.0) + 21571063.0 / 8773632.0)
			- 532328621.0 / 138915840.0) + 105364263311.0 / 28005433344.0)
			- 62407010457.0 / 25930956800.0) + 50704007329529.0 / 48057323618304.0);
	}
	else if (Argument == 3.5) {
/**/	DEBUG_WRITE_LEAVING(Omoms10, "Done")
		return(297590927.0 / 574665062400.0);
	}
	else if (Argument < 4.5) {
/**/	DEBUG_WRITE_LEAVING(Omoms10, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (Argument * (Argument * (Argument * (11.0 / 90720.0
			- Argument * (1.0 / 362880.0)) - 541.0 / 225792.0) + 2959.0 / 105840.0)
			- 1677257.0 / 7879680.0) + 25318843.0 / 22982400.0) - 977236835.0 / 250048512.0)
			+ 1092629879.0 / 117210240.0) - 594602457269.0 / 42008150016.0)
			+ 3181407449671.0 / 262550937600.0) - 1536983047451107.0 / 360429927137280.0);
	}
	else if (Argument == 4.5) {
/**/	DEBUG_WRITE_LEAVING(Omoms10, "Done")
		return(15015569.0 / 9386196019200.0);
	}
	else if (Argument < 5.5) {
/**/	DEBUG_WRITE_LEAVING(Omoms10, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (Argument * (Argument * (Argument * (Argument * (1.0 / 3628800.0)
			- 11.0 / 725760.0) + 283.0 / 752640.0) - 9383.0 / 1693440.0)
			+ 29622839.0 / 551577600.0) - 65779153.0 / 183859200.0)
			+ 4147831559.0 / 2500485120.0) - 2832070043.0 / 535818240.0)
			+ 930545778049.0 / 84016300032.0) - 28948369472011.0 / 2100407500800.0)
			+ 1639552734983279.0 / 212017604198400.0);
	}
	else if (Argument == 5.5) {
/**/	DEBUG_WRITE_LEAVING(Omoms10, "Done")
		return(1.0 / 28158588057600.0);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Omoms10, "Done")
		return(0.0);
	}
} /* end Omoms10 */

/*--------------------------------------------------------------------------*/
extern double	Omoms11
				(
					double	Argument			/* input */
				)

/* returns the value of a Blu optimum function of degree 11 (order 12) evaluated at Argument */

{ /* begin Omoms11 */

/**/DEBUG_WRITE_ENTERING(Omoms11,
/**/	"About to compute an oMoms of degree 11")

	Argument = fabs(Argument);
	if (Argument < 1.0) {
/**/	DEBUG_WRITE_LEAVING(Omoms11, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (Argument * (Argument * (Argument * (Argument * (Argument * (1.0 / 86400.0)
			- 1.0 / 14400.0) + 11.0 / 238464.0) + 31.0 / 66240.0) + 11.0 / 386400.0)
			- 79.0 / 15525.0) + 11.0 / 3146400.0) + 11473.0 / 314640.0) + 11.0 / 154047744.0)
			- 760452073.0 / 4493059200.0) + 1.0 / 8986118400.0)
			+ 56442621569.0 / 148270953600.0);
	}
	else if (Argument < 2.0) {
/**/	DEBUG_WRITE_LEAVING(Omoms11, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (Argument * (Argument * (Argument * (Argument * (1.0 / 6720.0
			- Argument * (1.0 / 120960.0)) - 9383.0 / 8346240.0) + 59.0 / 13248.0)
			- 5093.0 / 540960.0) + 3847.0 / 347760.0) - 443971.0 / 22024800.0)
			+ 3443.0 / 62928.0) - 62998331.0 / 5391671040.0) - 147576503.0 / 898611840.0)
			- 81893201.0 / 62902828800.0) + 79051603531.0 / 207579335040.0);
	}
	else if (Argument < 3.0) {
/**/	DEBUG_WRITE_LEAVING(Omoms11, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (Argument * (Argument * (Argument * (Argument * (Argument * (1.0 / 241920.0)
			- 1.0 / 8064.0) + 27599.0 / 16692480.0) - 2375.0 / 185472.0)
			+ 68387.0 / 1081920.0) - 28643.0 / 139104.0) + 19765603.0 / 44049600.0)
			- 596387.0 / 880992.0) + 8573363843.0 / 10783342080.0)
			- 273887483.0 / 359444736.0) + 33491661593.0 / 125805657600.0)
			+ 135383117087.0 / 415158670080.0);
	}
	else if (Argument < 4.0) {
/**/	DEBUG_WRITE_LEAVING(Omoms11, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (Argument * (Argument * (Argument * (Argument * (1.0 / 17280.0
			- Argument * (1.0 / 725760.0)) - 54923.0 / 50077440.0) + 109.0 / 8832.0)
			- 296153.0 / 3245760.0) + 22999.0 / 49680.0) - 215968951.0 / 132148800.0)
			+ 501677.0 / 125856.0) - 212059333211.0 / 32350026240.0)
			+ 4195676377.0 / 599074560.0) - 1766543369981.0 / 377416972800.0)
			+ 104670284681.0 / 59308381440.0);
	}
	else if (Argument < 5.0) {
/**/	DEBUG_WRITE_LEAVING(Omoms11, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (Argument * (Argument * (Argument * (Argument * (Argument * (1.0 / 3628800.0)
			- 1.0 / 67200.0) + 18271.0 / 50077440.0) - 4967.0 / 927360.0)
			+ 847781.0 / 16228800.0) - 616961.0 / 1738800.0) + 225780203.0 / 132148800.0)
			- 25605701.0 / 4404960.0) + 440867545327.0 / 32350026240.0)
			- 187002461077.0 / 8986118400.0) + 6942498798313.0 / 377416972800.0)
			- 14457478232713.0 / 2075793350400.0);
	}
	else if (Argument < 6.0) {
/**/	DEBUG_WRITE_LEAVING(Omoms11, "Done")
		return(Argument * (Argument * (Argument * (Argument * (Argument * (Argument
			* (Argument * (Argument * (Argument * (Argument * (1.0 / 604800.0
			- Argument * (1.0 / 39916800.0)) - 2489.0 / 50077440.0) + 119.0 / 132480.0)
			- 175981.0 / 16228800.0) + 35477.0 / 386400.0) - 73521109.0 / 132148800.0)
			+ 1518187.0 / 629280.0) - 237350824133.0 / 32350026240.0)
			+ 133911171673.0 / 8986118400.0) - 75527106646069.0 / 4151586700800.0)
			+ 6997050016189.0 / 691931116800.0);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Omoms11, "Done")
		return(0.0);
	}
} /* end Omoms11 */

/*--------------------------------------------------------------------------*/
extern int		OmomsArray03
				(
					double	Argument,			/* fractional input */
					double	*b3_minus2,			/* 1st returned coefficient */
					double	*b3_minus1,			/* 2nd returned coefficient */
					double	*b3_plus0,			/* 3rd returned coefficient */
					double	*b3_plus1			/* 4th returned coefficient */
				)

/* returns 4 values for an oMoms function of degree 3 (order 4) */
/* evaluation is performed at {Argument - 2.0, Argument - 1.0, Argument, Argument + 1.0} */
/* Argument must be in [0.0, 1.0] */
/* computational load: */
/* 7 indirections */
/* 6(double)assignments, 11 (double)additions, 9 (double)multiplications */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin OmomsArray03 */

	double	a1 = 1.0 - Argument, a2 = Argument * Argument;
	int		Status = !ERROR;

/**/DEBUG_CHECK_RANGE_DOUBLE(OmomsArray03, Argument, 0.0, 1.0, Status,
/**/	"Invalid argument (should be in [0.0, 1.0])")
/**/DEBUG_CHECK_NULL_POINTER(OmomsArray03, b3_minus2, Status,
/**/	"No 1st output coefficient")
/**/DEBUG_CHECK_NULL_POINTER(OmomsArray03, b3_minus1, Status,
/**/	"No 2nd output coefficient")
/**/DEBUG_CHECK_NULL_POINTER(OmomsArray03, b3_plus0, Status,
/**/	"No 3rd output coefficient")
/**/DEBUG_CHECK_NULL_POINTER(OmomsArray03, b3_plus1, Status,
/**/	"No 4th output coefficient")
/**/DEBUG_RETURN_ON_ERROR(OmomsArray03, Status)
/**/DEBUG_WRITE_ENTERING(OmomsArray03,
/**/	"About to compute an array of B-splines of degree 3")

	*b3_minus2 = Argument * (a2 * (1.0 / 6.0) + 1.0 / 42.0);
	*b3_minus1 = a1 * (a1 * (a1 * (1.0 / 2.0) - 1.0) + 1.0 / 14.0) + 13.0 / 21.0;
	*b3_plus0 = Argument * (Argument * (Argument * (1.0 / 2.0) - 1.0) + 1.0 / 14.0) + 13.0 / 21.0;
	*b3_plus1 = 1.0 - *b3_minus2 - *b3_minus1 - *b3_plus0;

/**/DEBUG_WRITE_LEAVING(OmomsArray03, "Done")
	return(Status);
} /* end OmomsArray03 */

/*--------------------------------------------------------------------------*/
extern double	Positive
				(
					double	Argument			/* input */
				)

/* returns the value of the Positive kernel of degree 3 evaluated at Argument */

{ /* begin Positive */

	double	a;

/**/DEBUG_WRITE_ENTERING(Positive,
/**/	"About to compute a Positive kernel")

	Argument = fabs(Argument);
	if (Argument < 0.5) {
		a = Argument * Argument;
/**/	DEBUG_WRITE_LEAVING(Positive, "Done")
		return(a * (Argument * 4.0 - 9.0 / 2.0) + 1.0);
	}
	if (Argument < 1.0) {
		Argument -= 1.0;
/**/	DEBUG_WRITE_LEAVING(Positive, "Done")
		return(Argument * Argument * (3.0 / 2.0));
	}
	else if (Argument < 1.5) {
		Argument -= 1.0;
/**/	DEBUG_WRITE_LEAVING(Positive, "Done")
		return(Argument * Argument * (3.0 / 2.0 - Argument * 2.0));
	}
	else if (Argument < 2.0) {
		Argument -= 2.0;
/**/	DEBUG_WRITE_LEAVING(Positive, "Done")
		return(Argument * Argument * (Argument * 2.0 + 3.0 / 2.0));
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Positive, "Done")
		return(0.0);
	}
} /* end Positive */

/*--------------------------------------------------------------------------*/
extern int		Schaum
				(
					long	Degree,				/* degree */
					double	Argument,			/* input */
					double	*Result				/* output */
				)

/* computes a Schaum interpolating function */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin Schaum */

	double	i, j, k, l;
	double	a;
	int		Status = !ERROR;

/**/DEBUG_CHECK_RANGE_LONG(Schaum, Degree, 0L, LONG_MAX, Status,
/**/	"Invalid degree (should be positive)")
/**/DEBUG_CHECK_NULL_POINTER(Schaum, Result, Status,
/**/	"No output data")
/**/DEBUG_RETURN_ON_ERROR(Schaum, Status)
/**/DEBUG_WRITE_ENTERING(Schaum,
/**/	"About to compute a generic Schaum")

	switch (Degree) {
		case 2L:
			*Result = Schaum02(Argument);
			break;
		case 3L:
			*Result = Schaum03(Argument);
			break;
		case 4L:
			*Result = Schaum04(Argument);
			break;
		default:
			Argument = fabs(Argument);
			i = -(double)((Degree + 1L) / 2L);
			if ((Degree & 1L) != 0L) {
				*Result = -1.0;
				k = floor(Argument);
				j = i + k + 1.0;
			}
			else {
				*Result = 1.0;
				k = floor(Argument + 0.5);
				j = i + k;
			}
			i -= k;
			l = (double)(Degree / 2L) - k;
			if (k <= (double)(Degree / 2L)) {
				while (i <= l) {
					if (i != 0.0) {
						*Result *= Argument + i;
					}
					if (j != 0.0) {
						*Result /= j;
					}
					i += 1.0;
					j += 1.0;
				}
				k -= 1.0;
				if (((Degree & 1L) == 0L) && ((Argument - 0.5) == k)) {
					i = -(double)((Degree + 1L) / 2L);
					a = 1.0;
					j = i + k;
					i -= k;
					l = (double)(Degree / 2L) - k;
					while (i <= l) {
						if (i != 0.0) {
							a *= Argument + i;
						}
						if (j != 0.0) {
							a /= j;
						}
						i += 1.0;
						j += 1.0;
					}
					*Result = 0.5 * (*Result + a);
				}
			}
			else {
				k -= 1.0;
				if (((Degree & 1L) == 0L) && ((Argument - 0.5) == k)
					&& (k <= (double)(Degree / 2L))) {
					i = -(double)((Degree + 1L) / 2L);
					*Result = 0.5;
					j = i + k;
					i -= k;
					l = (double)(Degree / 2L) - k;
					while (i <= l) {
						if (i != 0.0) {
							*Result *= Argument + i;
						}
						if (j != 0.0) {
							*Result /= j;
						}
						i += 1.0;
						j += 1.0;
					}
				}
				else {
					*Result *= 0.0;
				}
			}
			break;
	}
/**/DEBUG_WRITE_LEAVING(Schaum, "Done")
	return(Status);
} /* end Schaum */

/*--------------------------------------------------------------------------*/
extern double	Schaum02
				(
					double	Argument			/* input */
				)

/* returns the value of the Schaum kernel of degree 2 evaluated at Argument */

{ /* begin Schaum02 */

/**/DEBUG_WRITE_ENTERING(Schaum02,
/**/	"About to compute a Schaum of degree 2")

	Argument = fabs(Argument);
	if (Argument < 0.5) {
/**/	DEBUG_WRITE_LEAVING(Schaum02, "Done")
		return(1.0 - Argument * Argument);
	}
	else if (Argument == 0.5) {
/**/	DEBUG_WRITE_LEAVING(Schaum02, "Done")
		return(9.0 / 16.0);
	}
	else if (Argument < 1.5) {
/**/	DEBUG_WRITE_LEAVING(Schaum02, "Done")
		return(Argument * (Argument * (1.0 / 2.0) - 3.0 / 2.0) + 1.0);
	}
	else if (Argument == 1.5) {
/**/	DEBUG_WRITE_LEAVING(Schaum02, "Done")
		return(-1.0 / 16.0);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Schaum02, "Done")
		return(0.0);
	}
} /* end Schaum02 */

/*--------------------------------------------------------------------------*/
extern double	Schaum03
				(
					double	Argument			/* input */
				)

/* returns the value of the Schaum kernel of degree 3 evaluated at Argument */

{ /* begin Schaum03 */

/**/DEBUG_WRITE_ENTERING(Schaum03,
/**/	"About to compute a Schaum of degree 3")

	Argument = fabs(Argument);
	if (Argument < 1.0) {
/**/	DEBUG_WRITE_LEAVING(Schaum03, "Done")
		return(Argument * (Argument * (Argument * (1.0 / 2.0) - 1.0) - 1.0 / 2.0) + 1.0);
	}
	else if (Argument < 2.0) {
/**/	DEBUG_WRITE_LEAVING(Schaum03, "Done")
		return(Argument * (Argument * (1.0 - Argument * (1.0 / 6.0)) - 11.0 / 6.0) + 1.0);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Schaum03, "Done")
		return(0.0);
	}
} /* end Schaum03 */

/*--------------------------------------------------------------------------*/
extern double	Schaum04
				(
					double	Argument			/* input */
				)

/* returns the value of the Schaum kernel of degree 4 evaluated at Argument */

{ /* begin Schaum04 */

/**/DEBUG_WRITE_ENTERING(Schaum04,
/**/	"About to compute a Schaum of degree 4")

	Argument = fabs(Argument);
	if (Argument < 0.5) {
		Argument *= Argument;
/**/	DEBUG_WRITE_LEAVING(Schaum04, "Done")
		return(Argument * (Argument * (1.0 / 4.0) - 5.0 / 4.0) + 1.0);
	}
	else if (Argument == 0.5) {
/**/	DEBUG_WRITE_LEAVING(Schaum04, "Done")
		return(75.0 / 128.0);
	}
	else if (Argument < 1.5) {
/**/	DEBUG_WRITE_LEAVING(Schaum04, "Done")
		return(Argument * (Argument * (Argument * (5.0 / 6.0 - Argument
			* (1.0 / 6.0)) - 5.0 / 6.0) - 5.0 / 6.0) + 1.0);
	}
	else if (Argument == 1.5) {
/**/	DEBUG_WRITE_LEAVING(Schaum04, "Done")
		return(-25.0 / 256.0);
	}
	else if (Argument < 2.5) {
/**/	DEBUG_WRITE_LEAVING(Schaum04, "Done")
		return(Argument * (Argument * (Argument * (Argument * (1.0 / 24.0)
			-5.0 / 12.0) + 35.0 / 24.0) - 25.0 / 12.0) + 1.0);
	}
	else if (Argument == 2.5) {
/**/	DEBUG_WRITE_LEAVING(Schaum04, "Done")
		return(3.0 / 256.0);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Schaum04, "Done")
		return(0.0);
	}
} /* end Schaum04 */

/*--------------------------------------------------------------------------*/
extern double	Sinc
				(
					double	Argument			/* input */
				)

/* returns the value of the sinc kernel evaluated at Argument */

{ /* begin Sinc */

/**/DEBUG_WRITE_ENTERING(Sinc,
/**/	"About to compute a sinc kernel")

	if (Argument == 0.0) {
/**/	DEBUG_WRITE_LEAVING(Sinc, "Done")
		return(1.0);
	}
	else {
		Argument *= PI;
/**/	DEBUG_WRITE_LEAVING(Sinc, "Done")
		return(sin(Argument) / Argument);
	}
} /* end Sinc */

