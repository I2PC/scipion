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
#include	"kernelintegrate.h"
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
extern int		BlipIntegrate
				(
					long	Degree,				/* degree */
					double	Argument,			/* input */
					double	*Result				/* output */
				)

/* computes the integral(-Infinity, Argument) for a Blu interpolant function */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin BlipIntegrate */

	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(BlipIntegrate, Result, Status,
/**/	"No output data")
/**/DEBUG_RETURN_ON_ERROR(BlipIntegrate, Status)
/**/DEBUG_WRITE_ENTERING(BlipIntegrate,
/**/	"About to integrate a generic blip")

	switch (Degree) {
		case 0L:
			*Result = Blip00Integrate(Argument);
			break;
		case 1L:
			*Result = Blip01Integrate(Argument);
			break;
		case 3L:
			*Result = Blip03Integrate(Argument);
			break;
		default:
			Status = ERROR;
			WRITE_ERROR(BlipIntegrate, "Invalid blip degree")
			break;
	}
/**/DEBUG_WRITE_LEAVING(BlipIntegrate, "Done")
	return(Status);
} /* end BlipIntegrate */

/*--------------------------------------------------------------------------*/
extern double	Blip00Integrate
				(
					double	Argument			/* input */
				)

/* returns the integral(-Infinity, Argument) for a Blu interpolant function
	of degree 0 (order 1) */

{ /* begin Blip00Integrate */

	double	a;

/**/DEBUG_WRITE_ENTERING(Blip00Integrate,
/**/	"About to integrate a blip of degree 0")

	a = fabs(Argument);
	if (a < 0.5) {
		a = a + 1.0 / 2.0;
	}
	else {
		a = 1.0;
	}
/**/DEBUG_WRITE_LEAVING(Blip00Integrate, "Done")
	return((0.0 < Argument) ? (a) : (1.0 - a));
} /* end Blip00Integrate */

/*--------------------------------------------------------------------------*/
extern double	Blip01Integrate
				(
					double	Argument			/* input */
				)

/* returns the integral(-Infinity, Argument) for a Blu interpolant function
	of degree 1 (order 2) */

{ /* begin Blip01Integrate */

	double	a;

/**/DEBUG_WRITE_ENTERING(Blip01Integrate,
/**/	"About to integrate a blip of degree 1")

	a = fabs(Argument);
	if (a < 1.0) {
		a = a * (1.0 - a * (1.0 / 2.0)) + 1.0 / 2.0;
	}
	else {
		a = 1.0;
	}
/**/DEBUG_WRITE_LEAVING(Blip01Integrate, "Done")
	return((0.0 < Argument) ? (a) : (1.0 - a));
} /* end Blip01Integrate */

/*--------------------------------------------------------------------------*/
extern double	Blip03Integrate
				(
					double	Argument			/* input */
				)

/* returns the integral(-Infinity, Argument) for a Blu interpolant function
	of degree 3 (order 4) */

{ /* begin Blip03Integrate */

	double	a;

/**/DEBUG_WRITE_ENTERING(Blip03Integrate,
/**/	"About to integrate a blip of degree 3")

	a = fabs(Argument);
	if (a < 1.0) {
		a = a * (a * (a * (a * (1.0 / 8.0) - 1.0 / 3.0) + 3.0 / 112.0)
			+ 53.0 / 84.0) + 1.0 / 2.0;
	}
	else if (a < 2.0) {
		a = a * (a * (a * (1.0 / 3.0 - a * (1.0 / 24.0)) - 113.0 / 112.0)
			+ 115.0 / 84.0) + 25.0 / 84.0;
	}
	else {
		a = 1.0;
	}
/**/DEBUG_WRITE_LEAVING(Blip03Integrate, "Done")
	return((0.0 < Argument) ? (a) : (1.0 - a));
} /* end Blip03Integrate */

/*--------------------------------------------------------------------------*/
extern int		BsplineIntegrate
				(
					long	Degree,				/* degree */
					double	Argument,			/* input */
					double	*Result				/* output */
				)

/* computes the integral(-Infinity, Argument) for a Basic spline function of degree Degree*/
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin BsplineIntegrate */

	double	B, C;
	double	HalfWidth;
	long	i;
	int		Status = !ERROR;

/**/DEBUG_CHECK_RANGE_LONG(BsplineIntegrate, Degree, 0L, LONG_MAX, Status,
/**/	"Invalid degree (should be positive)")
/**/DEBUG_CHECK_NULL_POINTER(BsplineIntegrate, Result, Status,
/**/	"No output data")
/**/DEBUG_RETURN_ON_ERROR(BsplineIntegrate, Status)
/**/DEBUG_WRITE_ENTERING(BsplineIntegrate,
/**/	"About to integrate a generic B-spline")

	switch (Degree) {
		case 0L:
			*Result = Bspline00Integrate(Argument);
			break;
		case 1L:
			*Result = Bspline01Integrate(Argument);
			break;
		case 2L:
			*Result = Bspline02Integrate(Argument);
			break;
		case 3L:
			*Result = Bspline03Integrate(Argument);
			break;
		case 4L:
			*Result = Bspline04Integrate(Argument);
			break;
		case 5L:
			*Result = Bspline05Integrate(Argument);
			break;
		case 6L:
			*Result = Bspline06Integrate(Argument);
			break;
		case 7L:
			*Result = Bspline07Integrate(Argument);
			break;
		case 8L:
			*Result = Bspline08Integrate(Argument);
			break;
		case 9L:
			*Result = Bspline09Integrate(Argument);
			break;
		case 10L:
			*Result = Bspline10Integrate(Argument);
			break;
		case 11L:
			*Result = Bspline11Integrate(Argument);
			break;
		default:
			HalfWidth = 0.5 * (double)(Degree + 1L);
			*Result = (HalfWidth < Argument) ? (0.0)
				: ((Argument == HalfWidth) ? (0.5) : (1.0));
			C = Argument + HalfWidth;
			for (i = 0L; (i <= Degree); i++) {
				Status = Bspline(Degree - i, Argument, &B);
				if (Status == ERROR) {
/**/				DEBUG_WRITE_LEAVING(BsplineIntegrate, "Done")
					return(Status);
				}
				*Result += B * C / (double)(Degree + 1L - i);
				C -= 1.0;
				Argument -= 0.5;
			}
			break;
	}
/**/DEBUG_WRITE_LEAVING(BsplineIntegrate, "Done")
	return(Status);
} /* end BsplineIntegrate */

/*--------------------------------------------------------------------------*/
extern double	Bspline00Integrate
				(
					double	Argument			/* input */
				)

/* returns the integral(-Infinity, Argument) for a Basic spline function
	of degree 0 (order 1) */

{ /* begin Bspline00Integrate */

	double	a;

/**/DEBUG_WRITE_ENTERING(Bspline00Integrate,
/**/	"About to integrate a B-spline of degree 0")

	a = fabs(Argument);
	if (a < 0.5) {
		a = a + 1.0 / 2.0;
	}
	else {
		a = 1.0;
	}
/**/DEBUG_WRITE_LEAVING(Bspline00Integrate, "Done")
	return((0.0 < Argument) ? (a) : (1.0 - a));
} /* end Bspline00Integrate */

/*--------------------------------------------------------------------------*/
extern double	Bspline01Integrate
				(
					double	Argument			/* input */
				)

/* returns the integral(-Infinity, Argument) for a Basic spline function
	of degree 1 (order 2) */

{ /* begin Bspline01Integrate */

	double	a;

/**/DEBUG_WRITE_ENTERING(Bspline01Integrate,
/**/	"About to integrate a B-spline of degree 1")

	a = fabs(Argument);
	if (a < 1.0) {
		a = a * (1.0 - a * (1.0 / 2.0)) + 1.0 / 2.0;
	}
	else {
		a = 1.0;
	}
/**/DEBUG_WRITE_LEAVING(Bspline01Integrate, "Done")
	return((0.0 < Argument) ? (a) : (1.0 - a));
} /* end Bspline01Integrate */

/*--------------------------------------------------------------------------*/
extern double	Bspline02Integrate
				(
					double	Argument			/* input */
				)

/* returns the integral(-Infinity, Argument) for a Basic spline function
	of degree 2 (order 3) */

{ /* begin Bspline02Integrate */

	double	a;

/**/DEBUG_WRITE_ENTERING(Bspline02Integrate,
/**/	"About to integrate a B-spline of degree 2")

	a = fabs(Argument);
	if (a < 0.5) {
		a = a * (3.0 / 4.0 - a * a * (1.0 / 3.0)) + 1.0 / 2.0;
	}
	else if (a < 1.5) {
		a = a * (a * (a * (1.0 / 6.0) - 3.0 / 4.0) + 9.0 / 8.0) + 7.0 / 16.0;
	}
	else {
		a = 1.0;
	}
/**/DEBUG_WRITE_LEAVING(Bspline02Integrate, "Done")
	return((0.0 < Argument) ? (a) : (1.0 - a));
} /* end Bspline02Integrate */

/*--------------------------------------------------------------------------*/
extern double	Bspline03Integrate
				(
					double	Argument			/* input */
				)

/* returns the integral(-Infinity, Argument) for a Basic spline function
	of degree 3 (order 4) */

{ /* begin Bspline03Integrate */

	double	a;

/**/DEBUG_WRITE_ENTERING(Bspline02Integrate,
/**/	"About to integrate a B-spline of degree 3")

	a = fabs(Argument);
	if (a < 1.0) {
		a = a * (a * (a * (a * (1.0 / 8.0) - 1.0 / 3.0)) + 2.0 / 3.0) + 1.0 / 2.0;
	}
	else if (a < 2.0) {
		a = a * (a * (a * (1.0 / 3.0 - a * (1.0 / 24.0)) - 1.0) + 4.0 / 3.0) + 1.0 / 3.0;
	}
	else {
		a = 1.0;
	}
/**/DEBUG_WRITE_LEAVING(Bspline03Integrate, "Done")
	return((0.0 < Argument) ? (a) : (1.0 - a));
} /* end Bspline03Integrate */

/*--------------------------------------------------------------------------*/
extern double	Bspline04Integrate
				(
					double	Argument			/* input */
				)

/* returns the integral(-Infinity, Argument) for a Basic spline function
	of degree 4 (order 5) */

{ /* begin Bspline04Integrate */

	double	a, a2;

/**/DEBUG_WRITE_ENTERING(Bspline04Integrate,
/**/	"About to integrate a B-spline of degree 4")

	a = fabs(Argument);
	if (a < 0.5) {
		a2 = a * a;
		a = a * (a2 * (a2 * (1.0 / 20.0) - 5.0 / 24.0) + 115.0 / 192.0) + 1.0 / 2.0;
	}
	else if (a < 1.5) {
		a = a * (a * (a * (a * (5.0 / 24.0 - a * (1.0 / 30.0)) - 5.0 / 12.0)
			+ 5.0 / 48.0) + 55.0 / 96.0) + 193.0 / 384.0;
	}
	else if (a < 2.5) {
		a = a * (a * (a * (a * (a * (1.0 / 120.0) - 5.0 / 48.0) + 25.0 / 48.0)
			- 125.0 / 96.0) + 625.0 / 384.0) + 143.0 / 768.0;
	}
	else {
		a = 1.0;
	}
/**/DEBUG_WRITE_LEAVING(Bspline04Integrate, "Done")
	return((0.0 < Argument) ? (a) : (1.0 - a));
} /* end Bspline04Integrate */

/*--------------------------------------------------------------------------*/
extern double	Bspline05Integrate
				(
					double	Argument			/* input */
				)

/* returns the integral(-Infinity, Argument) for a Basic spline function
	of degree 5 (order 6) */

{ /* begin Bspline05Integrate */

	double	a, a2;

/**/DEBUG_WRITE_ENTERING(Bspline05Integrate,
/**/	"About to integrate a B-spline of degree 5")

	a = fabs(Argument);
	if (a < 1.0) {
		a2 = a * a;
		a = a * (a2 * (a2 * (1.0 / 20.0 - a * (1.0 / 72.0)) - 1.0 / 6.0)
			+ 11.0 / 20.0) + 1.0 / 2.0;
	}
	else if (a < 2.0) {
		a = a * (a * (a * (a * (a * (a * (1.0 / 144.0) - 3.0 / 40.0) + 5.0 / 16.0)
			- 7.0 / 12.0) + 5.0 / 16.0) + 17.0 / 40.0) + 25.0 / 48.0;
	}
	else if (a < 3.0) {
		a = a * (a * (a * (a * (a * (1.0 / 40.0 - a * (1.0 / 720.0)) - 3.0 / 16.0)
			+ 3.0 / 4.0) - 27.0 / 16.0) + 81.0 / 40.0) - 1.0 / 80.0;
	}
	else {
		a = 1.0;
	}
/**/DEBUG_WRITE_LEAVING(Bspline05Integrate, "Done")
	return((0.0 < Argument) ? (a) : (1.0 - a));
} /* end Bspline05Integrate */

/*--------------------------------------------------------------------------*/
extern double	Bspline06Integrate
				(
					double	Argument			/* input */
				)

/* returns the integral(-Infinity, Argument) for a Basic spline function
	of degree 6 (order 7) */

{ /* begin Bspline06Integrate */

	double	a, a2;

/**/DEBUG_WRITE_ENTERING(Bspline06Integrate,
/**/	"About to integrate a B-spline of degree 6")

	a = fabs(Argument);
	if (a < 0.5) {
		a2 = a * a;
		a = a * (a2 * (a2 * (7.0 / 240.0 - a2 * (1.0 / 252.0)) - 77.0 / 576.0)
			+ 5887.0 / 11520.0) + 1.0 / 2.0;
	}
	else if (a < 1.5) {
		a = a * (a * (a * (a * (a * (a * (a * (1.0 / 336.0) - 7.0 / 288.0)
			+ 21.0 / 320.0) - 35.0 / 1152.0) - 91.0 / 768.0) - 7.0 / 1536.0)
			+ 7861.0 / 15360.0) + 9215.0 / 18432.0;
	}
	else if (a < 2.5) {
		a = a * (a * (a * (a * (a * (a * (7.0 / 360.0 - a * (1.0 / 840.0))
			- 21.0 / 160.0) + 133.0 / 288.0) - 329.0 / 384.0) + 1267.0 / 1920.0)
			+ 1379.0 / 7680.0) + 13159.0 / 23040.0;
	}
	else if (a < 3.5) {
		a = a * (a * (a * (a * (a * (a * (a * (1.0 / 5040.0) - 7.0 / 1440.0)
			+ 49.0 / 960.0) - 343.0 / 1152.0) + 2401.0 / 2304.0) - 16807.0 / 7680.0)
			+ 117649.0 / 46080.0) - 25489.0 / 92160.0;
	}
	else {
		a = 1.0;
	}
/**/DEBUG_WRITE_LEAVING(Bspline06Integrate, "Done")
	return((0.0 < Argument) ? (a) : (1.0 - a));
} /* end Bspline06Integrate */

/*--------------------------------------------------------------------------*/
extern double	Bspline07Integrate
				(
					double	Argument			/* input */
				)

/* returns the integral(-Infinity, Argument) for a Basic spline function
	of degree 7 (order 8) */

{ /* begin Bspline07Integrate */

	double	a, a2;

/**/DEBUG_WRITE_ENTERING(Bspline07Integrate,
/**/	"About to integrate a B-spline of degree 7")

	a = fabs(Argument);
	if (a < 1.0) {
		a2 = a * a;
		a = a * (a2 * (a2 * (a2 * (a * (1.0 / 1152.0) - 1.0 / 252.0) + 1.0 / 45.0)
			- 1.0 / 9.0) + 151.0 / 315.0) + 1.0 / 2.0;
	}
	else if (a < 2.0) {
		a = a * (a * (a * (a * (a * (a * (a * (1.0 / 140.0 - a * (1.0 /1920.0))
			- 7.0 / 180.0) + 1.0 /10.0) - 7.0 / 72.0) - 1.0 / 30.0) - 7.0 / 180.0)
			+ 103.0 / 210.0) + 359.0 / 720.0;
	}
	else if (a < 3.0) {
		a = a * (a * (a * (a * (a * (a * (a * (a * (1.0 / 5760.0) - 1.0 / 252.0)
			+ 7.0 / 180.0) - 19.0 / 90.0) + 49.0 / 72.0) - 23.0 / 18.0) + 217.0 / 180.0)
			- 139.0 / 630.0) + 487.0 / 720.0;
	}
	else if (a < 4.0) {
		a = a * (a * (a * (a * (a * (a * (a * (1.0 / 1260.0 - a * (1.0 / 40320.0))
			- 1.0 / 90.0) + 4.0 / 45.0) - 4.0 / 9.0) + 64.0 / 45.0) - 128.0 / 45.0)
			+ 1024.0 / 315.0) - 197.0 / 315.0;
	}
	else {
		a = 1.0;
	}
/**/DEBUG_WRITE_LEAVING(Bspline07Integrate, "Done")
	return((0.0 < Argument) ? (a) : (1.0 - a));
} /* end Bspline07Integrate */

/*--------------------------------------------------------------------------*/
extern double	Bspline08Integrate
				(
					double	Argument			/* input */
				)

/* returns the integral(-Infinity, Argument) for a Basic spline function
	of degree 8 (order 9) */

{ /* begin Bspline08Integrate */

	double	a, a2;

/**/DEBUG_WRITE_ENTERING(Bspline08Integrate,
/**/	"About to integrate a B-spline of degree 8")

	a = fabs(Argument);
	if (a < 0.5) {
		a2 = a * a;
		a = a * (a2 * (a2 * (a2 * (a2 * (1.0 / 5184.0) - 1.0 / 448.0) + 43.0 / 2560.0)
			- 289.0 / 3072.0) + 259723.0 / 573440.0) + 1.0 / 2.0;
	}
	else if (a < 1.5) {
		a = a * (a * (a * (a * (a * (a * (a * (a * (1.0 / 640.0 - a * (1.0 / 6480.0))
			- 3.0 / 560.0) + 7.0 / 1920.0) + 9.0 / 640.0) + 7.0 / 5120.0) - 121.0 / 1280.0)
			+ 1.0 / 10240.0) + 64929.0 / 143360.0) + 737281.0 / 1474560.0;
	}
	else if (a < 2.5) {
		a = a * (a * (a * (a * (a * (a * (a * (a * (a * (1.0 / 12960.0) - 1.0 / 640.0)
			+ 3.0 / 224.0) - 119.0 / 1920.0) + 207.0 / 1280.0) - 1127.0 / 5120.0)
			+ 65.0 / 512.0) - 1457.0 / 10240.0) + 145167.0 / 286720.0) + 724159.0 / 1474560.0;
	}
	else if (a < 3.5) {
		a = a * (a * (a * (a * (a * (a * (a * (a * (3.0 / 4480.0 - a * (1.0 / 45360.0))
			- 1.0 / 112.0) + 131.0 / 1920.0) - 209.0 / 640.0) + 5123.0 / 5120.0)
			- 1465.0 / 768.0) + 146051.0 / 71680.0) - 122729.0 / 143360.0)
			+ 8975363.0 / 10321920.0;
	}
	else if (a < 4.5) {
		a = a * (a * (a * (a * (a * (a * (a * (a * (a * (1.0 / 362880.0) - 1.0 / 8960.0)
			+ 9.0 / 4480.0) - 27.0 / 1280.0) + 729.0 / 5120.0) - 6561.0 / 10240.0)
			+ 19683.0 / 10240.0) - 531441.0 / 143360.0) + 4782969.0 / 1146880.0)
			- 2489209.0 / 2293760.0;
	}
	else {
		a = 1.0;
	}
/**/DEBUG_WRITE_LEAVING(Bspline08Integrate, "Done")
	return((0.0 < Argument) ? (a) : (1.0 - a));
} /* end Bspline08Integrate */

/*--------------------------------------------------------------------------*/
extern double	Bspline09Integrate
				(
					double	Argument			/* input */
				)

/* returns the integral(-Infinity, Argument) for a Basic spline function
	of degree 9 (order 10) */

{ /* begin Bspline09Integrate */

	double	a, a2;

/**/DEBUG_WRITE_ENTERING(Bspline09Integrate,
/**/	"About to integrate a B-spline of degree 9")

	a = fabs(Argument);
	if (a < 1.0) {
		a2 = a * a;
		a = a * (a2 * (a2 * (a2 * (a2 * (1.0 / 5184.0 - a * (1.0 / 28800.0))
			- 5.0 / 3024.0) + 19.0 / 1440.0) - 35.0 / 432.0) + 15619.0 / 36288.0) + 1.0 / 2.0;
	}
	else if (a < 2.0) {
		a = a * (a * (a * (a * (a * (a * (a * (a * (a * (a * (1.0 / 43200.0)
			- 1.0 / 2592.0) + 1.0 / 384.0) - 13.0 / 1512.0) + 7.0 / 576.0)
			- 1.0 / 720.0) + 7.0 / 576.0) - 19.0 / 216.0) + 1.0 / 384.0)
			+ 7799.0 / 18144.0) + 8641.0 / 17280.0;
	}
	else if (a < 3.0) {
		a = a * (a * (a * (a * (a * (a * (a * (a * (a * (5.0 / 18144.0
			- a * (1.0 / 100800.0)) - 3.0 / 896.0) + 5.0 / 216.0) - 19.0 / 192.0)
			+ 191.0 / 720.0) - 83.0 / 192.0) + 635.0 / 1512.0) - 339.0 / 896.0)
			+ 1553.0 / 2592.0) + 18797.0 / 40320.0;
	}
	else if (a < 4.0) {
		a = a * (a * (a * (a * (a * (a * (a * (a * (a * (a * (1.0 / 403200.0)
			- 1.0 / 10368.0) + 3.0 / 1792.0) - 103.0 / 6048.0) + 43.0 / 384.0)
			- 1423.0 / 2880.0) + 563.0 / 384.0) - 2449.0 / 864.0) + 5883.0 / 1792.0)
			- 133663.0 / 72576.0) + 96643.0 / 80640.0;
	}
	else if (a < 5.0) {
		a = a * (a * (a * (a * (a * (a * (a * (a * (a * (1.0 / 72576.0
			- a * (1.0 / 3628800.0)) - 5.0 / 16128.0) + 25.0 / 6048.0)
			- 125.0 / 3456.0) + 125.0 / 576.0) - 3125.0 / 3456.0) + 15625.0 / 6048.0)
			- 78125.0 / 16128.0) + 390625.0 / 72576.0) - 245473.0 / 145152.0;
	}
	else {
		a = 1.0;
	}
/**/DEBUG_WRITE_LEAVING(Bspline09Integrate, "Done")
	return((0.0 < Argument) ? (a) : (1.0 - a));
} /* end Bspline09Integrate */

/*--------------------------------------------------------------------------*/
extern double	Bspline10Integrate
				(
					double	Argument			/* input */
				)

/* returns the integral(-Infinity, Argument) for a Basic spline function
	of degree 10 (order 11) */

{ /* begin Bspline10Integrate */

	double	a, a2;

/**/DEBUG_WRITE_ENTERING(Bspline10Integrate,
/**/	"About to integrate a B-spline of degree 10")

	a = fabs(Argument);
	if (a < 0.5) {
		a2 = a * a;
		a = a * (a2 * (a2 * (a2 * (a2 * (11.0 / 103680.0 - a2 * (1.0 / 158400.0))
			- 583.0 / 483840.0) + 14597.0 / 1382400.0) - 156409.0 / 2211840.0)
			+ 381773117.0 / 928972800.0) + 1.0 / 2.0;
	}
	else if (a < 1.5) {
		a = a * (a * (a * (a * (a * (a * (a * (a * (a * (a * (a * (1.0 / 190080.0)
			- 11.0 / 172800.0) + 11.0 / 41472.0) - 11.0 / 46080.0) - 187.0 / 193536.0)
			- 77.0 / 460800.0) + 1177.0 / 110592.0) - 11.0 / 368640.0) - 62557.0 / 884736.0)
			- 11.0 / 8847360.0) + 152709293.0 / 371589120.0) + 88473599.0 / 176947200.0;
	}
	else if (a < 2.5) {
		a = a * (a * (a * (a * (a * (a * (a * (a * (a * (a * (11.0 / 151200.0
			- a * (1.0 / 332640.0)) - 55.0 / 72576.0) + 11.0 / 2520.0)
			- 715.0 / 48384.0) + 1661.0 / 57600.0) - 4543.0 / 138240.0)
			+ 7513.0 / 161280.0) - 163603.0 / 1548288.0) + 135311.0 / 7741440.0)
			+ 37690169.0 / 92897280.0) + 19381279.0 / 38707200.0;
	}
	else if (a < 3.5) {
		a = a * (a * (a * (a * (a * (a * (a * (a * (a * (a * (a * (1.0 / 887040.0)
			- 11.0 / 268800.0) + 55.0 / 82944.0) - 451.0 / 71680.0) + 14905.0 / 387072.0)
			- 48433.0 / 307200.0) + 479281.0 / 1105920.0) - 1352153.0 / 1720320.0)
			+ 1654543.0 / 1769472.0) - 11695211.0 / 13762560.0) + 623786977.0 / 743178240.0)
			+ 47441011.0 / 117964800.0;
	}
	else if (a < 4.5) {
		a = a * (a * (a * (a * (a * (a * (a * (a * (a * (a * (11.0 / 907200.0
			- a * (1.0 / 3991680.0)) - 11.0 / 41472.0) + 209.0 / 60480.0)
			- 5753.0 / 193536.0) + 61061.0 / 345600.0) - 407429.0 / 552960.0)
			+ 2070343.0 / 967680.0) - 3702215.0 / 884736.0) + 237959711.0 / 46448640.0)
			- 1241720381.0 / 371589120.0) + 28739771.0 / 16588800.0;
	}
	else if (a < 5.5) {
		a = a * (a * (a * (a * (a * (a * (a * (a * (a * (a * (a * (1.0 / 39916800.0)
			- 11.0 / 7257600.0) + 121.0 / 2903040.0) - 1331.0 / 1935360.0)
			+ 14641.0 / 1935360.0) - 161051.0 / 2764800.0) + 1771561.0 / 5529600.0)
			- 19487171.0 / 15482880.0) + 214358881.0 / 61931520.0)
			- 2357947691.0 / 371589120.0) + 25937424601.0 / 3715891200.0)
			- 18505642201.0 / 7431782400.0;
	}
	else {
		a = 1.0;
	}
/**/DEBUG_WRITE_LEAVING(Bspline10Integrate, "Done")
	return((0.0 < Argument) ? (a) : (1.0 - a));
} /* end Bspline10Integrate */

/*--------------------------------------------------------------------------*/
extern double	Bspline11Integrate
				(
					double	Argument			/* input */
				)

/* returns the integral(-Infinity, Argument) for a Basic spline function
	of degree 11 (order 12) */

{ /* begin Bspline11Integrate */

	double	a, a2;

/**/DEBUG_WRITE_ENTERING(Bspline11Integrate,
/**/	"About to integrate a B-spline of degree 11")

	a = fabs(Argument);
	if (a < 1.0) {
		a2 = a * a;
		a = a * (a2 * (a2 * (a2 * (a2 * (a2 * (a * (1.0 / 1036800.0)
			- 1.0 / 158400.0) + 1.0 / 12960.0) - 23.0 / 25200.0) + 31.0 / 3600.0)
			- 809.0 / 12960.0) + 655177.0 / 1663200.0) + 1.0 / 2.0;
	}
	else if (a < 2.0) {
		a = a * (a * (a * (a * (a * (a * (a * (a * (a * (a * (a * (1.0 / 73920.0
			- a * (1.0 / 1451520.0)) - 11.0 / 100800.0) + 1.0 / 2268.0)
			- 11.0 / 13440.0) + 1.0 / 2520.0) - 11.0 / 7200.0) + 5.0 / 504.0)
			- 11.0 / 13440.0) - 563.0 / 9072.0) - 11.0 / 100800.0) + 65521.0 / 166320.0)
			+ 302399.0 / 604800.0;
	}
	else if (a < 3.0) {
		a = a * (a * (a * (a * (a * (a * (a * (a * (a * (a * (a * (a * (1.0 / 2903040.0)
			- 1.0 / 88704.0) + 11.0 / 67200.0) - 25.0 / 18144.0) + 33.0 / 4480.0)
			- 13.0 / 504.0) + 143.0 / 2400.0) - 239.0 / 2520.0) + 583.0 / 4480.0)
			- 1619.0 / 9072.0) + 781.0 / 11200.0) + 61297.0 / 166320.0)
			+ 101653.0 / 201600.0;
	}
	else if (a < 4.0) {
		a = a * (a * (a * (a * (a * (a * (a * (a * (a * (a * (a * (1.0 / 190080.0
			- a * (1.0 / 8709120.0)) - 11.0 / 100800.0) + 7.0 / 5184.0) - 99.0 / 8960.0)
			+ 631.0 / 10080.0) - 1199.0 /4800.0) + 1009.0 / 1440.0) - 12199.0 / 8960.0)
			+ 9385.0 / 5184.0) - 38533.0 / 22400.0) + 894727.0 / 665280.0)
			+ 104891.0 / 403200.0;
	}
	else if (a < 5.0) {
		a = a * (a * (a * (a * (a * (a * (a * (a * (a * (a * (a * (a * (1.0 / 43545600.0)
			- 1.0 / 739200.0) + 11.0 / 302400.0) - 107.0 / 181440.0) + 517.0 / 80640.0)
			- 2477.0 / 50400.0) + 11737.0 / 43200.0) - 54797.0 / 50400.0)
			+ 250657.0 / 80640.0) - 1113317.0 / 181440.0) + 4726777.0 / 604800.0)
			- 18595037.0 / 3326400.0) + 9332627.0 / 3628800.0;
	}
	else if (a < 6.0) {
		a = a * (a * (a * (a * (a * (a * (a * (a * (a * (a * (a * (1.0 / 6652800.0
			- a * (1.0 / 479001600.0)) - 1.0 / 201600.0) + 1.0 / 10080.0) - 3.0 / 2240.0)
			+ 9.0 / 700.0) - 9.0 / 100.0) + 81.0 / 175.0) - 243.0 / 140.0)
			+ 162.0 / 35.0) - 1458.0 / 175.0) + 17496.0 / 1925.0) - 6823.0 / 1925.0;
	}
	else {
		a = 1.0;
	}
/**/DEBUG_WRITE_LEAVING(Bspline11Integrate, "Done")
	return((0.0 < Argument) ? (a) : (1.0 - a));
} /* end Bspline11Integrate */

/*--------------------------------------------------------------------------*/
extern double	DodgsonIntegrate
				(
					double	Argument			/* input */
				)

/* returns the integral(-Infinity, Argument) for a Dodgson kernel
	(order 2) evaluated at Argument */

{ /* begin DodgsonIntegrate */

	double	a;

/**/DEBUG_WRITE_ENTERING(DodgsonIntegrate,
/**/	"About to integrate a Dodgson kernel")

	a = fabs(Argument);
	if (a < 0.5) {
		a = a * (1.0 - a * a * (2.0 / 3.0)) + 1.0 / 2.0;
	}
	else if (a < 1.5) {
		a = a * (a * (a * (1.0 / 3.0) - 5.0 / 4.0) + 3.0 / 2) + 7.0 / 16;
	}
	else {
		a = 1.0;
	}
/**/DEBUG_WRITE_LEAVING(DodgsonIntegrate, "Done")
	return((0.0 < Argument) ? (a) : (1.0 - a));
} /* end DodgsonIntegrate */

/*--------------------------------------------------------------------------*/
extern double	KeysOptimalIntegrate
				(
					double	Argument			/* input */
				)

/* returns the integral(-Infinity, Argument) for a cubic Keys optimal kernel
	(order 3) evaluated at Argument */

{ /* begin KeysOptimalIntegrate */

	double	a;

/**/DEBUG_WRITE_ENTERING(KeysOptimalIntegrate,
/**/	"About to integrate a cubic Keys kernel")

	a = fabs(Argument);
	if (a < 1.0) {
		a = a * (a * (a * (a * (3.0 / 8.0) - 5.0 / 6.0)) + 1.0) + 1.0 / 2.0;
	}
	else if (a < 2.0) {
		a = a * (a * (a * (5.0 / 6.0 - a * (1.0 / 8.0)) - 2.0) + 2.0) + 1.0 / 3.0;
	}
	else {
		a = 1.0;
	}
/**/DEBUG_WRITE_LEAVING(KeysOptimalIntegrate, "Done")
	return((0.0 < Argument) ? (a) : (1.0 - a));
} /* end KeysOptimalIntegrate */

/*--------------------------------------------------------------------------*/
extern int		OmomsIntegrate
				(
					long	Degree,				/* degree */
					double	Argument,			/* input */
					double	*Result				/* output */
				)

/* computes the integral(-Infinity, Argument) for a Blu optimum function */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin OmomsIntegrate */

	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(OmomsIntegrate, Result, Status,
/**/	"No output data")
/**/DEBUG_RETURN_ON_ERROR(OmomsIntegrate, Status)
/**/DEBUG_WRITE_ENTERING(OmomsIntegrate,
/**/	"About to integrate a generic oMoms")

	switch (Degree) {
		case 0L:
			*Result = Omoms00Integrate(Argument);
			break;
		case 1L:
			*Result = Omoms01Integrate(Argument);
			break;
		case 2L:
			*Result = Omoms02Integrate(Argument);
			break;
		case 3L:
			*Result = Omoms03Integrate(Argument);
			break;
		default:
			Status = ERROR;
			WRITE_ERROR(OmomsIntegrate, "Invalid oMoms degree")
			break;
	}
/**/DEBUG_WRITE_LEAVING(OmomsIntegrate, "Done")
	return(Status);
} /* end OmomsIntegrate */

/*--------------------------------------------------------------------------*/
extern double	Omoms00Integrate
				(
					double	Argument			/* input */
				)

/* returns the integral(-Infinity, Argument) for a Blu interpolant function
	of degree 0 (order 1) */

{ /* begin Omoms00Integrate */

	double	a;

/**/DEBUG_WRITE_ENTERING(Omoms00Integrate,
/**/	"About to integrate an oMoms of degree 0")

	a = fabs(Argument);
	if (a < 0.5) {
		a = a + 1.0 / 2.0;
	}
	else {
		a = 1.0;
	}
/**/DEBUG_WRITE_LEAVING(Omoms00Integrate, "Done")
	return((0.0 < Argument) ? (a) : (1.0 - a));
} /* end Omoms00Integrate */

/*--------------------------------------------------------------------------*/
extern double	Omoms01Integrate
				(
					double	Argument			/* input */
				)

/* returns the integral(-Infinity, Argument) for a Blu interpolant function
	of degree 1 (order 2) */

{ /* begin Omoms01Integrate */

	double	a;

/**/DEBUG_WRITE_ENTERING(Omoms01Integrate,
/**/	"About to integrate an oMoms of degree 1")

	a = fabs(Argument);
	if (a < 1.0) {
		a = a * (1.0 - a * (1.0 / 2.0)) + 1.0 / 2.0;
	}
	else {
		a = 1.0;
	}
/**/DEBUG_WRITE_LEAVING(Omoms01Integrate, "Done")
	return((0.0 < Argument) ? (a) : (1.0 - a));
} /* end Omoms01Integrate */

/*--------------------------------------------------------------------------*/
extern double	Omoms02Integrate
				(
					double	Argument			/* input */
				)

/* returns the integral(-Infinity, Argument) for a Blu optimum function
	of degree 2 (order 3) */

{ /* begin Omoms02Integrate */

	double	a;

/**/DEBUG_WRITE_ENTERING(Omoms02Integrate,
/**/	"About to integrate an oMoms of degree 2")

	a = fabs(Argument);
	if (a < 0.5) {
		a = a * (43.0 / 60.0 - a * a * (1.0 / 3.0)) + 1.0 / 2.0;
	}
	else if (a < 1.5) {
		a = a * (a * (a * (1.0 / 6.0) - 3.0 / 4.0) + 137.0 / 120.0) + 33.0 / 80.0;
	}
	else {
		a = 1.0;
	}
/**/DEBUG_WRITE_LEAVING(Omoms02Integrate, "Done")
	return((0.0 < Argument) ? (a) : (1.0 - a));
} /* end Omoms02Integrate */

/*--------------------------------------------------------------------------*/
extern double	Omoms03Integrate
				(
					double	Argument			/* input */
				)

/* returns the integral(-Infinity, Argument) for a Blu optimum function
	of degree 3 (order 4) */

{ /* begin Omoms03Integrate */

	double	a;

/**/DEBUG_WRITE_ENTERING(Omoms03Integrate,
/**/	"About to integrate an oMoms of degree 3")

	a = fabs(Argument);
	if (a < 1.0) {
		a = a * (a * (a * (a * (1.0 / 8.0) - 1.0 / 3.0) + 1.0 / 28.0)
			+ 13.0 / 21.0) + 1.0 / 2.0;
	}
	else if (a < 2.0) {
		a = a * (a * (a * (1.0 / 3.0 - a * (1.0 / 24.0)) - 85.0 / 84.0)
			+ 29.0 / 21.0) + 2.0 / 7.0;
	}
	else {
		a = 1.0;
	}
/**/DEBUG_WRITE_LEAVING(Omoms03Integrate, "Done")
	return((0.0 < Argument) ? (a) : (-a));
} /* end Omoms03Integrate */

/*--------------------------------------------------------------------------*/
extern double	SincIntegrate
				(
					double	Argument			/* input */
				)

/* returns the approximate integral(-Infinity, Argument) for a sinc kernel
	evaluated at Argument */
/* maximum error is about 2.0E-7 */

{ /* begin SincIntegrate */

	double	a, a2;

/**/DEBUG_WRITE_ENTERING(SincIntegrate,
/**/	"About to integrate a sinc kernel")

	a = fabs(Argument);
	if (a < 0.5) {
		a2 = (PI * PI) * a * a;
		a *= a2 * (a2 * (a2 * (a2 * (1.0 /3265920.0) - 1.0 / 35280.0)
			+ 1.0 / 600.0) - 1.0 /18.0) + 1.0;
	}
	else {
		a *= PI;
		a2 = a * a;
		a = 0.5 - (1.0 / PI) * ((1.0 / a) * (a2 * (a2 * (a2 * (a2 + 38.027264)
			+ 265.187033) + 335.677320) + 38.102495) / (a2 * (a2 * (a2 * (a2
			+ 40.021433) + 322.624911) + 570.236280) + 157.105423) * cos(a)
			+ (1.0 / a2) * (a2 * (a2 * (a2 * (a2 + 42.242855) + 302.757865)
			+ 352.018498) + 21.821899) / (a2 * (a2 * (a2 * (a2 + 48.196927)
			+ 482.485984) + 1114.978885) + 449.690326) * sin(a));
	}
/**/DEBUG_WRITE_LEAVING(SincIntegrate, "Done")
	return((0.0 < Argument) ? (a) : (1.0 - a));
} /* end SincIntegrate */

