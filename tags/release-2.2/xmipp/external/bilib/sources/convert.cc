/*****************************************************************************
 *	System includes
 ****************************************************************************/
#include	<float.h>
#include	<limits.h>
#include	<math.h>
#include	<stdio.h>
#include	<string.h>

/*****************************************************************************
 *	Toolbox defines
 ****************************************************************************/
#include	"configs.h"
#include	"debug.h"
#include	"error.h"

/*****************************************************************************
 *	Toolbox includes
 ****************************************************************************/
#include	"convert.h"
#include	"messagedisplay.h"

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
extern float	ConvertDoubleToFloat
				(
					double	Argument			/* value to convert */
				)

/* converts a double to a float */
/* saturation occurs when Argument generates overflow */

{ /* begin ConvertDoubleToFloat */

/**/DEBUG_WRITE_ENTERING(ConvertDoubleToFloat,
/**/	"About to convert a (double)argument to a (float)")

	if ((double)FLT_MAX < Argument) {
		WRITE_WARNING(ConvertDoubleToFloat, "Positive saturation occured")
/**/	DEBUG_WRITE_LEAVING(ConvertDoubleToFloat, "Done")
		return(FLT_MAX);
	}
	else if (Argument < (double)(-FLT_MAX)) {
		WRITE_WARNING(ConvertDoubleToFloat, "Negative saturation occured")
/**/	DEBUG_WRITE_LEAVING(ConvertDoubleToFloat, "Done")
		return(-FLT_MAX);
	}
	else {
#ifdef DEBUG
/**/	if ((fabs(Argument) < (double)FLT_MIN) && (Argument != 0.0)) {
/**/		WRITE_WARNING(ConvertDoubleToFloat, "Underflow occured")
/**/	}
#endif
/**/	DEBUG_WRITE_LEAVING(ConvertDoubleToFloat, "Done")
		return((float)Argument);
	}
} /* end ConvertDoubleToFloat */

/*--------------------------------------------------------------------------*/
extern int		ConvertDoubleToInt
				(
					double	Argument			/* value to round */
				)

/* rounds a double to an int */
/* saturation occurs when Argument generates overflow */

{ /* begin ConvertDoubleToInt */

	double	Rounded;

/**/DEBUG_WRITE_ENTERING(ConvertDoubleToInt,
/**/	"About to round a (double)argument to an (int)")

	if ((double)INT_MAX < Argument) {
		WRITE_WARNING(ConvertDoubleToInt, "Positive saturation occured")
/**/	DEBUG_WRITE_LEAVING(ConvertDoubleToInt, "Done")
		return(INT_MAX);
	}
	else if (Argument < (double)INT_MIN) {
		WRITE_WARNING(ConvertDoubleToInt, "Negative saturation occured")
/**/	DEBUG_WRITE_LEAVING(ConvertDoubleToInt, "Done")
		return(INT_MIN);
	}
	Argument += 0.5;
	Rounded = (double)((long)Argument);
	if ((Argument < 0.0) && (Rounded != Argument)) {
		Rounded -= 1.0;
	}
/**/DEBUG_WRITE_LEAVING(ConvertDoubleToInt, "Done")
	return((int)Rounded);
} /* end ConvertDoubleToInt */

/*--------------------------------------------------------------------------*/
extern long		ConvertDoubleToLong
				(
					double	Argument			/* value to round */
				)

/* rounds a double to a long */
/* saturation occurs when Argument generates overflow */

{ /* begin ConvertDoubleToLong */

	double	Rounded;

/**/DEBUG_WRITE_ENTERING(ConvertDoubleToLong,
/**/	"About to round a (double)argument to a (long)")

	if ((double)LONG_MAX < Argument) {
		WRITE_WARNING(ConvertDoubleToLong, "Positive saturation occured")
/**/	DEBUG_WRITE_LEAVING(ConvertDoubleToLong, "Done")
		return(LONG_MAX);
	}
	else if (Argument < (double)LONG_MIN) {
		WRITE_WARNING(ConvertDoubleToLong, "Negative saturation occured")
/**/	DEBUG_WRITE_LEAVING(ConvertDoubleToLong, "Done")
		return(LONG_MIN);
	}
	Argument += 0.5;
	Rounded = (double)((long)Argument);
	if ((Argument < 0.0) && (Rounded != Argument)) {
		Rounded -= 1.0;
	}
/**/DEBUG_WRITE_LEAVING(ConvertDoubleToLong, "Done")
	return((long)Rounded);
} /* end ConvertDoubleToLong */

/*--------------------------------------------------------------------------*/
extern short	ConvertDoubleToShort
				(
					double	Argument			/* value to round */
				)

/* rounds a double to a short */
/* saturation occurs when Argument generates overflow */

{ /* begin ConvertDoubleToShort */

	double	Rounded;

/**/DEBUG_WRITE_ENTERING(ConvertDoubleToShort,
/**/	"About to round a (double)argument to a (short)")

	if ((double)SHRT_MAX < Argument) {
		WRITE_WARNING(ConvertDoubleToShort, "Positive saturation occured")
/**/	DEBUG_WRITE_LEAVING(ConvertDoubleToShort, "Done")
		return((short)SHRT_MAX);
	}
	else if (Argument < (double)SHRT_MIN) {
		WRITE_WARNING(ConvertDoubleToShort, "Negative saturation occured")
/**/	DEBUG_WRITE_LEAVING(ConvertDoubleToShort, "Done")
		return((short)SHRT_MIN);
	}
	Argument += 0.5;
	Rounded = (double)((long)Argument);
	if ((Argument < 0.0) && (Rounded != Argument)) {
		Rounded -= 1.0;
	}
/**/DEBUG_WRITE_LEAVING(ConvertDoubleToShort, "Done")
	return((short)Rounded);
} /* end ConvertDoubleToShort */

/*--------------------------------------------------------------------------*/
extern int		ConvertFloatToInt
				(
					float	Argument			/* value to round */
				)

/* rounds a float to an int */
/* saturation occurs when Argument generates overflow */

{ /* begin ConvertFloatToInt */

	float	Rounded;

/**/DEBUG_WRITE_ENTERING(ConvertFloatToInt,
/**/	"About to round a (float)argument to an (int)")

	if ((double)INT_MAX < (double)Argument) {
		WRITE_WARNING(ConvertFloatToInt, "Positive saturation occured")
/**/	DEBUG_WRITE_LEAVING(ConvertFloatToInt, "Done")
		return(INT_MAX);
	}
	else if ((double)Argument < (double)INT_MIN) {
		WRITE_WARNING(ConvertFloatToInt, "Negative saturation occured")
/**/	DEBUG_WRITE_LEAVING(ConvertFloatToInt, "Done")
		return(INT_MIN);
	}
	Argument += 0.5F;
	Rounded = (float)((long)Argument);
	if ((Argument < 0.0F) && (Rounded != Argument)) {
		Rounded -= 1.0F;
	}
/**/DEBUG_WRITE_LEAVING(ConvertFloatToInt, "Done")
	return((int)Rounded);
} /* end ConvertFloatToInt */

/*--------------------------------------------------------------------------*/
extern long		ConvertFloatToLong
				(
					float	Argument			/* value to round */
				)

/* rounds a float to a long */
/* saturation occurs when Argument generates overflow */

{ /* begin ConvertFloatToLong */

	float	Rounded;

/**/DEBUG_WRITE_ENTERING(ConvertFloatToLong,
/**/	"About to round a (float)argument to a (long)")

	if ((double)LONG_MAX < (double)Argument) {
		WRITE_WARNING(ConvertFloatToLong, "Positive saturation occured")
/**/	DEBUG_WRITE_LEAVING(ConvertFloatToLong, "Done")
		return(LONG_MAX);
	}
	else if ((double)Argument < (double)LONG_MIN) {
		WRITE_WARNING(ConvertFloatToLong, "Negative saturation occured")
/**/	DEBUG_WRITE_LEAVING(ConvertFloatToLong, "Done")
		return(LONG_MIN);
	}
	Argument += 0.5F;
	Rounded = (float)((long)Argument);
	if ((Argument < 0.0) && (Rounded != Argument)) {
		Rounded -= 1.0F;
	}
/**/DEBUG_WRITE_LEAVING(ConvertFloatToLong, "Done")
	return((long)Rounded);
} /* end ConvertFloatToLong */

/*--------------------------------------------------------------------------*/
extern short	ConvertFloatToShort
				(
					float	Argument			/* value to round */
				)

/* rounds a float to a short */
/* saturation occurs when Argument generates overflow */

{ /* begin ConvertFloatToShort */

	float	Rounded;

/**/DEBUG_WRITE_ENTERING(ConvertFloatToShort,
/**/	"About to round a (float)argument to a (short)")

	if ((float)SHRT_MAX < Argument) {
		WRITE_WARNING(ConvertFloatToShort, "Positive saturation occured")
/**/	DEBUG_WRITE_LEAVING(ConvertFloatToShort, "Done")
		return((short)SHRT_MAX);
	}
	else if (Argument < (float)SHRT_MIN) {
		WRITE_WARNING(ConvertFloatToShort, "Negative saturation occured")
/**/	DEBUG_WRITE_LEAVING(ConvertFloatToShort, "Done")
		return((short)SHRT_MIN);
	}
	Argument += 0.5F;
	Rounded = (float)((long)Argument);
	if ((Argument < 0.0) && (Rounded != Argument)) {
		Rounded -= 1.0F;
	}
/**/DEBUG_WRITE_LEAVING(ConvertFloatToShort, "Done")
	return((short)Rounded);
} /* end ConvertFloatToShort */

/*--------------------------------------------------------------------------*/
extern short	ConvertIntToShort
				(
					int		Argument			/* value to convert */
				)

/* converts an int to a short */
/* saturation occurs when Argument generates overflow */

{ /* begin ConvertIntToShort */

/**/DEBUG_WRITE_ENTERING(ConvertIntToShort,
/**/	"About to convert an (int)argument to a (short)")

	if (SHRT_MAX < Argument) {
		WRITE_WARNING(ConvertIntToShort, "Positive saturation occured")
/**/	DEBUG_WRITE_LEAVING(ConvertIntToShort, "Done")
		return((short)SHRT_MAX);
	}
	else if (Argument < SHRT_MIN) {
		WRITE_WARNING(ConvertIntToShort, "Negative saturation occured")
/**/	DEBUG_WRITE_LEAVING(ConvertIntToShort, "Done")
		return((short)SHRT_MIN);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(ConvertIntToShort, "Done")
		return((short)Argument);
	}
} /* end ConvertIntToShort */

/*--------------------------------------------------------------------------*/
extern int		ConvertLongToInt
				(
					long	Argument			/* value to convert */
				)

/* converts a long to an int */
/* saturation occurs when Argument generates overflow */

{ /* begin ConvertLongToInt */

/**/DEBUG_WRITE_ENTERING(ConvertLongToInt,
/**/	"About to convert a (long)argument to an (int)")

	if ((long)INT_MAX < Argument) {
		WRITE_WARNING(ConvertLongToInt, "Positive saturation occured")
/**/	DEBUG_WRITE_LEAVING(ConvertLongToInt, "Done")
		return(INT_MAX);
	}
	else if (Argument < (long)INT_MIN) {
		WRITE_WARNING(ConvertLongToInt, "Negative saturation occured")
/**/	DEBUG_WRITE_LEAVING(ConvertLongToInt, "Done")
		return(INT_MIN);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(ConvertLongToInt, "Done")
		return((int)Argument);
	}
} /* end ConvertLongToInt */

/*--------------------------------------------------------------------------*/
extern short	ConvertLongToShort
				(
					long	Argument			/* value to convert */
				)

/* converts a long to a short */
/* saturation occurs when Argument generates overflow */

{ /* begin ConvertLongToShort */

/**/DEBUG_WRITE_ENTERING(ConvertLongToShort,
/**/	"About to convert a (long)argument to a (short)")

	if ((long)SHRT_MAX < Argument) {
		WRITE_WARNING(ConvertLongToShort, "Positive saturation occured")
/**/	DEBUG_WRITE_LEAVING(ConvertLongToShort, "Done")
		return((short)SHRT_MAX);
	}
	else if (Argument < (long)SHRT_MIN) {
		WRITE_WARNING(ConvertLongToShort, "Negative saturation occured")
/**/	DEBUG_WRITE_LEAVING(ConvertLongToShort, "Done")
		return((short)SHRT_MIN);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(ConvertLongToShort, "Done")
		return((short)Argument);
	}
} /* end ConvertLongToShort */

