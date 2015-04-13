/*****************************************************************************
 *	System includes
 ****************************************************************************/
/* None */

/*****************************************************************************
 *	Toolbox defines
 ****************************************************************************/
#include	"configs.h"
#include	"debug.h"

/*****************************************************************************
 *	Toolbox includes
 ****************************************************************************/
#include	"minmax.h"

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
extern double	MaxDouble
				(
					double	Argument1,			/* first argument */
					double	Argument2			/* second argument */
				)

/* returns the most positive of a pair of double arguments */

{ /* begin MaxDouble */

/**/DEBUG_WRITE_ENTERING(MaxDouble,
/**/	"About to return the most positive of two (double)arguments")
/**/DEBUG_WRITE_LEAVING(MaxDouble, "Done")
	return((Argument1 < Argument2) ? (Argument2) : (Argument1));
} /* end MaxDouble */

/*--------------------------------------------------------------------------*/
extern float	MaxFloat
				(
					float	Argument1,			/* first argument */
					float	Argument2			/* second argument */
				)

/* returns the most positive of a pair of float arguments */

{ /* begin MaxFloat */

/**/DEBUG_WRITE_ENTERING(MaxFloat,
/**/	"About to return the most positive of two (float)arguments")
/**/DEBUG_WRITE_LEAVING(MaxFloat, "Done")
	return((Argument1 < Argument2) ? (Argument2) : (Argument1));
} /* end MaxFloat */

/*--------------------------------------------------------------------------*/
extern int		MaxInt
				(
					int		Argument1,			/* first argument */
					int		Argument2			/* second argument */
				)

/* returns the most positive of a pair of int arguments */

{ /* begin MaxInt */

/**/DEBUG_WRITE_ENTERING(MaxInt,
/**/	"About to return the most positive of two (int)arguments")
/**/DEBUG_WRITE_LEAVING(MaxInt, "Done")
	return((Argument1 < Argument2) ? (Argument2) : (Argument1));
} /* end MaxInt */

/*--------------------------------------------------------------------------*/
extern long		MaxLong
				(
					long	Argument1,			/* first argument */
					long	Argument2			/* second argument */
				)

/* returns the most positive of a pair of long arguments */

{ /* begin MaxLong */

/**/DEBUG_WRITE_ENTERING(MaxLong,
/**/	"About to return the most positive of two (long)arguments")
/**/DEBUG_WRITE_LEAVING(MaxLong, "Done")
	return((Argument1 < Argument2) ? (Argument2) : (Argument1));
} /* end MaxLong */

/*--------------------------------------------------------------------------*/
extern short	MaxShort
				(
					short	Argument1,			/* first argument */
					short	Argument2)			/* second argument */

/* returns the most positive of a pair of short arguments */

{ /* begin MaxShort */

/**/DEBUG_WRITE_ENTERING(MaxShort,
/**/	"About to return the most positive of two (short)arguments")
/**/DEBUG_WRITE_LEAVING(MaxShort, "Done")
	return((Argument1 < Argument2) ? (Argument2) : (Argument1));
} /* end MaxShort */

/*--------------------------------------------------------------------------*/
extern double	MinDouble
				(
					double	Argument1,			/* first argument */
					double	Argument2			/* second argument */
				)

/* returns the most negative of a pair of double arguments */

{ /* begin MinDouble */

/**/DEBUG_WRITE_ENTERING(MinDouble,
/**/	"About to return the most negative of two (double)arguments")
/**/DEBUG_WRITE_LEAVING(MinDouble, "Done")
	return((Argument1 < Argument2) ? (Argument1) : (Argument2));
} /* end MinDouble */

/*--------------------------------------------------------------------------*/
extern float	MinFloat
				(
					float	Argument1,			/* first argument */
					float	Argument2			/* second argument */
				)

/* returns the most negative of a pair of float arguments */

{ /* begin MinFloat */

/**/DEBUG_WRITE_ENTERING(MinFloat,
/**/	"About to return the most negative of two (float)arguments")
/**/DEBUG_WRITE_LEAVING(MinFloat, "Done")
	return((Argument1 < Argument2) ? (Argument1) : (Argument2));
} /* end MinFloat */

/*--------------------------------------------------------------------------*/
extern int		MinInt
				(
					int		Argument1,			/* first argument */
					int		Argument2			/* second argument */
				)

/* returns the most negative of a pair of int arguments */

{ /* begin MinInt */

/**/DEBUG_WRITE_ENTERING(MinInt,
/**/	"About to return the most negative of two (int)arguments")
/**/DEBUG_WRITE_LEAVING(MinInt, "Done")
	return((Argument1 < Argument2) ? (Argument1) : (Argument2));
} /* end MinInt */

/*--------------------------------------------------------------------------*/
extern long		MinLong
				(
					long	Argument1,			/* first argument */
					long	Argument2			/* second argument */
				)

/* returns the most negative of a pair of long arguments */

{ /* begin MinLong */

/**/DEBUG_WRITE_ENTERING(MinLong,
/**/	"About to return the most negative of two (long)arguments")
/**/DEBUG_WRITE_LEAVING(MinLong, "Done")
	return((Argument1 < Argument2) ? (Argument1) : (Argument2));
} /* end MinLong */

/*--------------------------------------------------------------------------*/
extern short	MinShort
				(
					short	Argument1,			/* first argument */
					short	Argument2			/* second argument */
				)

/* returns the most negative of a pair of short arguments */

{ /* begin MinShort */

/**/DEBUG_WRITE_ENTERING(MinShort,
/**/	"About to return the most negative of two (short)arguments")
/**/DEBUG_WRITE_LEAVING(MinShort, "Done")
	return((Argument1 < Argument2) ? (Argument1) : (Argument2));
} /* end MinShort */

