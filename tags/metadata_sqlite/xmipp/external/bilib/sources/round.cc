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
#include	"../headers/convert.h"
#include	"round.h"

/*****************************************************************************
 *	Conditional includes
 ****************************************************************************/
#ifdef		DEBUG
#include	<stdio.h>
#include	<string.h>
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
extern double	RoundDoubleToDouble
				(
					double	Argument			/* value to round */
				)

/* rounds a double to a double */
/* returns the nearest integer value */

{ /* begin RoundDoubleToDouble */

/**/DEBUG_WRITE_ENTERING(RoundDoubleToDouble,
/**/	"About to round a (double)argument")

	Argument += 0.5;

/**/DEBUG_WRITE_LEAVING(RoundDoubleToDouble, "Done")
	return(floor(Argument));
} /* end RoundDoubleToDouble */

/*--------------------------------------------------------------------------*/
extern float	RoundDoubleToFloat
				(
					double	Argument			/* value to round */
				)

/* rounds a double to a float */
/* returns the nearest integer value */

{ /* begin RoundDoubleToFloat */

/**/DEBUG_WRITE_ENTERING(RoundDoubleToFloat,
/**/	"About to round a (double)argument to a (float)")

	Argument += 0.5;

/**/DEBUG_WRITE_LEAVING(RoundDoubleToFloat, "Done")
	return(ConvertDoubleToFloat(floor(Argument)));
} /* end RoundDoubleToFloat */

/*--------------------------------------------------------------------------*/
extern float	RoundFloatToFloat
				(
					float	Argument			/* value to round */
				)

/* rounds a foat to a foat */
/* returns the nearest integer value */

{ /* begin RoundFloatToFloat */

/**/DEBUG_WRITE_ENTERING(RoundFloatToFloat,
/**/	"About to round a (float)argument")

	Argument += 0.5F;

/**/DEBUG_WRITE_LEAVING(RoundFloatToFloat, "Done")
	return(ConvertDoubleToFloat(floor((double)Argument)));
} /* end RoundFloatToFloat */

