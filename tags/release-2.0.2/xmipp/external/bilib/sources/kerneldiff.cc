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
#include	"kerneldiff.h"
#include	"kerneldiff1.h"
#include	"kerneldiff2.h"
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
extern int		BlipDiff
				(
					long	Degree,				/* degree */
					long	Derivative,			/* order of the derivative */
					double	Argument,			/* input */
					double	*Result				/* output */
				)

/* computes a derivative of order Derivative for a Blu interpolant function */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin BlipDiff */

	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(BlipDiff, Result, Status,
/**/	"No output data")
/**/DEBUG_RETURN_ON_ERROR(BlipDiff, Status)
/**/DEBUG_WRITE_ENTERING(BlipDiff,
/**/	"About to compute the derivative of a blip")

	switch (Degree) {
		case 0L:
			switch (Derivative) {
				case 0L:
					*Result = Blip00(Argument);
					break;
				case 1L:
					*Result = Blip00Diff1(Argument);
					break;
				case 2L:
					*Result = Blip00Diff2(Argument);
					break;
				default:
					Status = ERROR;
					WRITE_ERROR(BlipDiff,
						"Invalid derivative order for a blip of degree 0")
					break;
			}
			break;
		case 1L:
			switch (Derivative) {
				case 0L:
					*Result = Blip01(Argument);
					break;
				case 1L:
					*Result = Blip01Diff1(Argument);
					break;
				case 2L:
					*Result = Blip01Diff2(Argument);
					break;
				default:
					Status = ERROR;
					WRITE_ERROR(BlipDiff,
						"Invalid derivative order for a blip of degree 1")
					break;
			}
			break;
		case 3L:
			switch (Derivative) {
				case 0L:
					*Result = Blip03(Argument);
					break;
				case 1L:
					*Result = Blip03Diff1(Argument);
					break;
				case 2L:
					*Result = Blip03Diff2(Argument);
					break;
				default:
					Status = ERROR;
					WRITE_ERROR(BlipDiff,
						"Invalid derivative order for a blip of degree 3")
					break;
			}
			break;
		default:
			Status = ERROR;
			WRITE_ERROR(BlipDiff, "Invalid blip degree")
			break;
	}
/**/DEBUG_WRITE_LEAVING(BlipDiff, "Done")
	return(Status);
} /* end BlipDiff */

/*--------------------------------------------------------------------------*/
extern int		BsplineDiff
				(
					long	Degree,				/* degree */
					long	Derivative,			/* order of the derivative */
					double	Argument,			/* input */
					double	*Result				/* output */
				)

/* computes a derivative of order Derivative for a Basic spline function of degree Degree*/
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin BsplineDiff */

	double	Left, Right;
	int		Status = !ERROR;
	int		Done = TRUE;

/**/DEBUG_CHECK_RANGE_LONG(BsplineDiff, Degree, 0L, LONG_MAX, Status,
/**/	"Invalid degree (should be positive)")
/**/DEBUG_CHECK_NULL_POINTER(BsplineDiff, Result, Status,
/**/	"No output data")
/**/DEBUG_RETURN_ON_ERROR(BsplineDiff, Status)
/**/DEBUG_WRITE_ENTERING(BsplineDiff,
/**/	"About to compute the derivative of a B-spline")

	switch (Degree) {
		case 0L:
			switch (Derivative) {
				case 0L:
					*Result = Bspline00(Argument);
					break;
				case 1L:
					*Result = Bspline00Diff1(Argument);
					break;
				case 2L:
					*Result = Bspline00Diff2(Argument);
					break;
				default:
					Status = ERROR;
					WRITE_ERROR(BsplineDiff,
						"Invalid derivative order for a B-spline of degree 0")
					break;
			}
			break;
		case 1L:
			switch (Derivative) {
				case 0L:
					*Result = Bspline01(Argument);
					break;
				case 1L:
					*Result = Bspline01Diff1(Argument);
					break;
				case 2L:
					*Result = Bspline01Diff2(Argument);
					break;
				default:
					Status = ERROR;
					WRITE_ERROR(BsplineDiff,
						"Invalid derivative order for a B-spline of degree 1")
					break;
			}
			break;
		case 2L:
			switch (Derivative) {
				case 0L:
					*Result = Bspline02(Argument);
					break;
				case 1L:
					*Result = Bspline02Diff1(Argument);
					break;
				case 2L:
					*Result = Bspline02Diff2(Argument);
					break;
				default:
					Status = ERROR;
					WRITE_ERROR(BsplineDiff,
						"Invalid derivative order for a B-spline of degree 2")
					break;
			}
			break;
		case 3L:
			switch (Derivative) {
				case 0L:
					*Result = Bspline03(Argument);
					break;
				case 1L:
					*Result = Bspline03Diff1(Argument);
					break;
				case 2L:
					*Result = Bspline03Diff2(Argument);
					break;
				default:
					Done = FALSE;
					break;
			}
			break;
		case 4L:
			switch (Derivative) {
				case 0L:
					*Result = Bspline04(Argument);
					break;
				case 1L:
					*Result = Bspline04Diff1(Argument);
					break;
				case 2L:
					*Result = Bspline04Diff2(Argument);
					break;
				default:
					Done = FALSE;
					break;
			}
			break;
		case 5L:
			switch (Derivative) {
				case 0L:
					*Result = Bspline05(Argument);
					break;
				case 1L:
					*Result = Bspline05Diff1(Argument);
					break;
				case 2L:
					*Result = Bspline05Diff2(Argument);
					break;
				default:
					Done = FALSE;
					break;
			}
			break;
		case 6L:
			switch (Derivative) {
				case 0L:
					*Result = Bspline06(Argument);
					break;
				case 1L:
					*Result = Bspline06Diff1(Argument);
					break;
				case 2L:
					*Result = Bspline06Diff2(Argument);
					break;
				default:
					Done = FALSE;
					break;
			}
			break;
		case 7L:
			switch (Derivative) {
				case 0L:
					*Result = Bspline07(Argument);
					break;
				case 1L:
					*Result = Bspline07Diff1(Argument);
					break;
				case 2L:
					*Result = Bspline07Diff2(Argument);
					break;
				default:
					Done = FALSE;
					break;
			}
			break;
		case 8L:
			switch (Derivative) {
				case 0L:
					*Result = Bspline08(Argument);
					break;
				case 1L:
					*Result = Bspline08Diff1(Argument);
					break;
				case 2L:
					*Result = Bspline08Diff2(Argument);
					break;
				default:
					Done = FALSE;
					break;
			}
			break;
		case 9L:
			switch (Derivative) {
				case 0L:
					*Result = Bspline09(Argument);
					break;
				case 1L:
					*Result = Bspline09Diff1(Argument);
					break;
				case 2L:
					*Result = Bspline09Diff2(Argument);
					break;
				default:
					Done = FALSE;
					break;
			}
			break;
		case 10L:
			switch (Derivative) {
				case 0L:
					*Result = Bspline10(Argument);
					break;
				case 1L:
					*Result = Bspline10Diff1(Argument);
					break;
				case 2L:
					*Result = Bspline10Diff2(Argument);
					break;
				default:
					Done = FALSE;
					break;
			}
			break;
		case 11L:
			switch (Derivative) {
				case 0L:
					*Result = Bspline11(Argument);
					break;
				case 1L:
					*Result = Bspline11Diff1(Argument);
					break;
				case 2L:
					*Result = Bspline11Diff2(Argument);
					break;
				default:
					Done = FALSE;
					break;
			}
			break;
		default:
			Done = FALSE;
			break;
	}
	if (!Done) {
/**/	DEBUG_CHECK_RANGE_LONG(BsplineDiff, Derivative, 0L, Degree, Status,
/**/		"Invalid derivative order (should be positive and should not exceed Degree)")
/**/	DEBUG_RETURN_ON_ERROR(BsplineDiff, Status)
		switch (Derivative) {
			case 0L:
				Status = Bspline(Degree, Argument, Result);
				break;
			case 1L:
				Status = BsplineDiff1(Degree, Argument, Result);
				break;
			case 2L:
				Status = BsplineDiff2(Degree, Argument, Result);
				break;
			default:
				Status = BsplineDiff(Degree - 1L, Derivative - 1L, Argument - 0.5, &Left);
				if (Status == ERROR) {
/**/				DEBUG_WRITE_LEAVING(BsplineDiff, "Done")
					return(Status);
				}
				Status = BsplineDiff(Degree - 1L, Derivative - 1L, Argument + 0.5, &Right);
				if (Status == ERROR) {
/**/				DEBUG_WRITE_LEAVING(BsplineDiff, "Done")
					return(Status);
				}
				*Result = Right - Left;
				break;
		}
	}
/**/DEBUG_WRITE_LEAVING(BsplineDiff, "Done")
	return(Status);
} /* end BsplineDiff */

/*--------------------------------------------------------------------------*/
extern int		DodgsonDiff
				(
					long	Derivative,			/* order of the derivative */
					double	Argument,			/* input */
					double	*Result				/* output */
				)

/* computes a derivative of order Derivative for a Dodgson kernel
	(order 2) evaluated at Argument */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin DodgsonDiff */

	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(DodgsonDiff, Result, Status,
/**/	"No output data")
/**/DEBUG_RETURN_ON_ERROR(DodgsonDiff, Status)
/**/DEBUG_WRITE_ENTERING(DodgsonDiff,
/**/	"About to compute the derivative of a Dodgson kernel")

	switch (Derivative) {
		case 0L:
			*Result = Dodgson(Argument);
			break;
		case 1L:
			*Result = DodgsonDiff1(Argument);
			break;
		case 2L:
			*Result = DodgsonDiff2(Argument);
			break;
		default:
			Status = ERROR;
			WRITE_ERROR(KeysOptimalDiff,
				"Invalid derivative order for a Dodgson kernel")
			break;
	}
/**/DEBUG_WRITE_LEAVING(KeysOptimalDiff, "Done")
	return(Status);
} /* end DodgsonDiff */

/*--------------------------------------------------------------------------*/
extern int		KeysOptimalDiff
				(
					long	Derivative,			/* order of the derivative */
					double	Argument,			/* input */
					double	*Result				/* output */
				)

/* computes a derivative of order Derivative for a cubic Keys optimal kernel
	(order 3) evaluated at Argument */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin KeysOptimalDiff */

	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(KeysOptimalDiff, Result, Status,
/**/	"No output data")
/**/DEBUG_RETURN_ON_ERROR(KeysOptimalDiff, Status)
/**/DEBUG_WRITE_ENTERING(KeysOptimalDiff,
/**/	"About to compute the derivative of a cubic Keys optimal kernel")

	switch (Derivative) {
		case 0L:
			*Result = KeysOptimal(Argument);
			break;
		case 1L:
			*Result = KeysOptimalDiff1(Argument);
			break;
		case 2L:
			*Result = KeysOptimalDiff2(Argument);
			break;
		default:
			Status = ERROR;
			WRITE_ERROR(KeysOptimalDiff,
				"Invalid derivative order for a cubic Keys optimal kernel")
			break;
	}
/**/DEBUG_WRITE_LEAVING(KeysOptimalDiff, "Done")
	return(Status);
} /* end KeysOptimalDiff */

/*--------------------------------------------------------------------------*/
extern int		OmomsDiff
				(
					long	Degree,				/* degree */
					long	Derivative,			/* order of the derivative */
					double	Argument,			/* input */
					double	*Result				/* output */
				)

/* computes a derivative of order Derivative for a Blu optimum function */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin OmomsDiff */

	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(OmomsDiff, Result, Status,
/**/	"No output data")
/**/DEBUG_RETURN_ON_ERROR(OmomsDiff, Status)
/**/DEBUG_WRITE_ENTERING(OmomsDiff,
/**/	"About to compute the derivative of an oMoms")

	switch (Degree) {
		case 0L:
			switch (Derivative) {
				case 0L:
					*Result = Omoms00(Argument);
					break;
				case 1L:
					*Result = Omoms00Diff1(Argument);
					break;
				case 2L:
					*Result = Omoms00Diff2(Argument);
					break;
				default:
					Status = ERROR;
					WRITE_ERROR(OmomsDiff,
						"Invalid derivative order for an oMoms of degree 0")
					break;
			}
			break;
		case 1L:
			switch (Derivative) {
				case 0L:
					*Result = Omoms01(Argument);
					break;
				case 1L:
					*Result = Omoms01Diff1(Argument);
					break;
				case 2L:
					*Result = Omoms01Diff2(Argument);
					break;
				default:
					Status = ERROR;
					WRITE_ERROR(OmomsDiff,
						"Invalid derivative order for an oMoms of degree 1")
					break;
			}
			break;
		case 2L:
			switch (Derivative) {
				case 0L:
					*Result = Omoms02(Argument);
					break;
				case 1L:
					*Result = Omoms02Diff1(Argument);
					break;
				case 2L:
					*Result = Omoms02Diff2(Argument);
					break;
				default:
					Status = ERROR;
					WRITE_ERROR(OmomsDiff,
						"Invalid derivative order for an oMoms of degree 2")
					break;
			}
			break;
		case 3L:
			switch (Derivative) {
				case 0L:
					*Result = Omoms03(Argument);
					break;
				case 1L:
					*Result = Omoms03Diff1(Argument);
					break;
				case 2L:
					*Result = Omoms03Diff2(Argument);
					break;
				default:
					Status = ERROR;
					WRITE_ERROR(OmomsDiff,
						"Invalid derivative order for an oMoms of degree 3")
					break;
			}
			break;
		default:
			Status = ERROR;
			WRITE_ERROR(OmomsDiff, "Invalid oMoms degree")
			break;
	}
/**/DEBUG_WRITE_LEAVING(OmomsDiff, "Done")
	return(Status);
} /* end OmomsDiff */

/*--------------------------------------------------------------------------*/
extern int		PositiveDiff
				(
					long	Derivative,			/* order of the derivative */
					double	Argument,			/* input */
					double	*Result				/* output */
				)

/* computes a derivative of order Derivative for a Positive kernel
	(order 1) evaluated at Argument */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin PositiveDiff */

	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(PositiveDiff, Result, Status,
/**/	"No output data")
/**/DEBUG_RETURN_ON_ERROR(PositiveDiff, Status)
/**/DEBUG_WRITE_ENTERING(PositiveDiff,
/**/	"About to compute the derivative of a Positive kernel")

	switch (Derivative) {
		case 0L:
			*Result = Positive(Argument);
			break;
		case 1L:
			*Result = PositiveDiff1(Argument);
			break;
		case 2L:
			*Result = PositiveDiff2(Argument);
			break;
		default:
			Status = ERROR;
			WRITE_ERROR(PositiveDiff,
				"Invalid derivative order for a Positive kernel")
			break;
	}
/**/DEBUG_WRITE_LEAVING(PositiveDiff, "Done")
	return(Status);
} /* end PositiveDiff */

/*--------------------------------------------------------------------------*/
extern int		SincDiff
				(
					long	Derivative,			/* order of the derivative */
					double	Argument,			/* input */
					double	*Result				/* output */
				)

/* computes a derivative of order Derivative for a sinc kernel evaluated at Argument */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin SincDiff */

	double	SumSin, SumCos, Sign, Fact;
	long	k;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(SincDiff, Result, Status,
/**/	"No output data")
/**/DEBUG_RETURN_ON_ERROR(SincDiff, Status)
/**/DEBUG_WRITE_ENTERING(SincDiff,
/**/	"About to compute the derivative of a sinc kernel")

	switch (Derivative) {
		case 0L:
			*Result = Sinc(Argument);
			break;
		case 1L:
			*Result = SincDiff1(Argument);
			break;
		case 2L:
			*Result = SincDiff2(Argument);
			break;
		default:
			if (Argument == 0.0) {
				*Result = cos((double)Derivative * (PI / 2.0))
					* PositiveIntPower(PI, Derivative) / (double)(Derivative + 1L);
			}
			else {
				SumSin = 0.0;
				Sign = 1.0 - (double)(2L * (Derivative & 1L));
				Fact = 1.0;
				k = 0L;
				while (k <= (2L * (Derivative / 2L))) {
					SumSin += Sign * PositiveIntPower(PI * Argument, k) / Fact;
					Sign *= -1.0;
					Fact *= (double)++k;
					Fact *= (double)++k;
				}
				SumSin *= Sinc(Argument);
				SumCos = 0.0;
				Sign = 1.0 - (double)(2L * (Derivative & 1L));
				Fact = 1.0;
				k = 1L;
				while (k <= (2L * ((Derivative + 1L) / 2L))) {
					SumCos += Sign * PositiveIntPower(PI * Argument, k - 1L) / Fact;
					Sign *= -1.0;
					Fact *= (double)++k;
					Fact *= (double)++k;
				}
				SumCos *= cos(PI * Argument);
				Fact = 1.0;
				for (k = 2L; (k <= Derivative); k++) {
					Fact *= (double)k;
				}
				*Result = Fact * (SumSin - SumCos) / PositiveIntPower(Argument, Derivative);
			}
			break;
	}
/**/DEBUG_WRITE_LEAVING(SincDiff, "Done")
	return(Status);
} /* end SincDiff */

