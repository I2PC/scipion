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
#include	"messagedisplay.h"
#include	"findroot.h"

/*****************************************************************************
 *	Conditional includes
 ****************************************************************************/
#ifdef		DEBUG
#include	<float.h>
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
extern int		RootBracket
				(
					int		(*Function)(double, void *, double *),
												/* function to bracket */
					void	*AuxilliaryData,	/* parameters used by the function */
					double	*LowerBound,		/* lower interval bound to be updated */
					double	*UpperBound,		/* upper interval bound to be updated */
					double	*LowerSample,		/* value of Function for argument LowerBound */
					double	*UpperSample,		/* value of Function for argument UpperBound */
					double	Tolerance,			/* admissible relative error */
					int		*ValidBracket,		/* whether or not a root could be bracketed */
					int		*Status				/* error management */
				)

/* bracketing of a root of (Function) inside the interval [LowerBound, UpperBound] */
/* the purpose is to search for a pair of arguments for which 'Function' differs in sign */
/* the search is conducted within [LowerBound, UpperBound] */
/* at most [(UpperBound - LowerBound) / Tolerance] function evaluations are performed */
/* even-order roots cannot be bracketed */
/* the evaluation of Function at LowerBound (UpperBound) */
/* is returned in LowerSample (UpperSample) */
/* success: return(!ERROR); failure: return(ERROR) */
/* the function 'Function' must be declared as follows: */
/*
		extern int myFunction(double myArgument, void *AuxilliaryData, double *myResult);
*/
/*		it must return ERROR upon failure, and !ERROR upon success */
/*		it is evaluated for the value of the variable 'myArgument' */
/*		the result of the function evaluation must be returned in 'myResult' */
/*		the generic pointer 'AuxilliaryData' can be used to pass additional parameters */
/* what follows is a developed example of the function f(x) = a * x^2 + b * x + c */
/* a, b, and c are free parameters */
/*
..............................................................................

struct myStruct
{
	double	a, b, c;
};

extern int		myFunction
				(
					double	myArgument,
					void	*AuxilliaryData,
					double	*myResult
				);
extern int		myFunction
				(
					double	myArgument,
					void	*AuxilliaryData,
					double	*myResult
				)

{
	struct myStruct
			myData;

	myData = *((struct myStruct *)AuxilliaryData);
	*myResult = myArgument * (myArgument * myData.a + myData.b) + myData.c;
	return(!ERROR);
}

int				main
				(
					void
				)

{
	struct myStruct
			myData;
	double	LowerBound = -100.0, UpperBound = 100.0;
	double	LowerSample, UpperSample;
	double	Tolerance = FLT_EPSILON;
	int		ValidBracket;
	int		Status;

	myData.a = 1.0;
	myData.b = 5.0;
	myData.c = 4.0;
	RootBracket(*myFunction, (void *)&myData, &LowerBound, &UpperBound,
		&LowerSample, &UpperSample, Tolerance, &ValidBracket, &Status);
	return(0);
}

..............................................................................
*/

{ /* begin RootBracket */

	double	Attempt, Sample, Sign;
	long	i, n;

	*Status = !ERROR;
	*ValidBracket = FALSE;

/**/DEBUG_CHECK_NULL_POINTER(RootBracket, Function, *Status,
/**/	"Missing function specification")
/**/DEBUG_CHECK_NULL_POINTER(RootBracket, AuxilliaryData, *Status,
/**/	"Missing AuxilliaryData")
/**/DEBUG_CHECK_NULL_POINTER(RootBracket, LowerBound, *Status,
/**/	"Missing LowerBound")
/**/DEBUG_CHECK_NULL_POINTER(RootBracket, UpperBound, *Status,
/**/	"Missing UpperBound")
/**/DEBUG_CHECK_NULL_POINTER(RootBracket, LowerSample, *Status,
/**/	"Missing LowerSample")
/**/DEBUG_CHECK_NULL_POINTER(RootBracket, UpperSample, *Status,
/**/	"Missing UpperSample")
/**/DEBUG_CHECK_RANGE_DOUBLE(RootBracket, Tolerance, 0.0, DBL_MAX, *Status,
/**/	"Invalid Tolerance (should be positive)")
/**/DEBUG_RETURN_ON_ERROR(RootBracket, *Status)
/**/DEBUG_WRITE_ENTERING(RootBracket,
/**/	"About to execute RootBracket")

	if (*UpperBound <= *LowerBound) {
		*Status = ERROR;
		WRITE_ERROR(RootBracket, "LowerBound must be strictly smaller than UpperBound")
/**/	DEBUG_WRITE_LEAVING(RootBracket, "Done")
		return(*Status);
	}
	*Status = (*Function)(*LowerBound, AuxilliaryData, LowerSample);
	if (*Status == ERROR) {
		WRITE_ERROR(RootBracket, "Function call failed")
/**/	DEBUG_WRITE_LEAVING(RootBracket, "Done")
		return(*Status);
	}
	Sign = *LowerSample;
	*Status = (*Function)(*UpperBound, AuxilliaryData, UpperSample);
	if (*Status == ERROR) {
		WRITE_ERROR(RootBracket, "Function call failed")
/**/	DEBUG_WRITE_LEAVING(RootBracket, "Done")
		return(*Status);
	}
	if ((*LowerSample * *UpperSample) < 0.0) {
		*ValidBracket = TRUE;
/**/	DEBUG_WRITE_LEAVING(RootBracket, "Done")
		return(*Status);
	}
	n = 1L;
	while (((double)n * Tolerance) <= 1.0) {
		n *= 2L;
		for (i = 1L; (i < n); i += 2L) {
			Attempt = *LowerBound + (double)i * (*UpperBound - *LowerBound) / (double)n;
			*Status = (*Function)(Attempt, AuxilliaryData, &Sample);
			if (*Status == ERROR) {
				WRITE_ERROR(RootBracket, "Function call failed")
/**/			DEBUG_WRITE_LEAVING(RootBracket, "Done")
				return(*Status);
			}
			if ((Sign * Sample) < 0.0) {
				*ValidBracket = TRUE;
				*LowerBound = Attempt - (*UpperBound - *LowerBound) / (double)n;
				*UpperBound = Attempt;
				*Status = (*Function)(*LowerBound, AuxilliaryData, LowerSample);
				if (*Status == ERROR) {
					WRITE_ERROR(RootBracket, "Function call failed")
/**/				DEBUG_WRITE_LEAVING(RootBracket, "Done")
					return(*Status);
				}
				*UpperSample = Sample;
/**/			DEBUG_WRITE_LEAVING(RootBracket, "Done")
				return(*Status);
			}
		}
	}
/**/DEBUG_WRITE_LEAVING(RootBracket, "Done")
	return(*Status);
} /* end RootBracket */

/*--------------------------------------------------------------------------*/
extern int		RootFindBisection
				(
					int		(*Function)(double, void *, double *),
												/* function, of which a root is sought */
					void	*AuxilliaryData,	/* parameters used by the function */
					double	*Root,				/* returned root */
					double	LowerBound,			/* lower bound of an interval containing a root */
					double	UpperBound,			/* upper bound of an interval containing a root */
					double	Tolerance,			/* admissible relative error */
					int		*Status				/* error management */
				)

/* search for a root of (Function) inside the bracketing interval [LowerBound, UpperBound] */
/* the strategy proceeds by iteratively cutting the interval in two equal parts */
/* even-order roots generally cannot be found */
/* only one root is returned, even if there are several ones */
/* Tolerance is relative to the size of the bracketing interval */
/* success: return(!ERROR); failure: return(ERROR) */
/* the function 'Function' must be declared as follows: */
/*
		extern int myFunction(double myArgument, void *AuxilliaryData, double *myResult);
*/
/*		it must return ERROR upon failure, and !ERROR upon success */
/*		it is evaluated for the value of the variable 'myArgument' */
/*		the result of the function evaluation must be returned in 'myResult' */
/*		the generic pointer 'AuxilliaryData' can be used to pass additional parameters */
/* what follows is a developed example of the function f(x) = a * x^2 + b * x + c */
/* a, b, and c are free parameters */
/*
..............................................................................

struct myStruct
{
	double	a, b, c;
};

extern int		myFunction
				(
					double	myArgument,
					void	*AuxilliaryData,
					double	*myResult
				);
extern int		myFunction
				(
					double	myArgument,
					void	*AuxilliaryData,
					double	*myResult
				)

{
	struct myStruct
			myData;

	myData = *((struct myStruct *)AuxilliaryData);
	*myResult = myArgument * (myArgument * myData.a + myData.b) + myData.c;
	return(!ERROR);
}

int				main
				(
					void
				)

{
	struct myStruct
			myData;
	double	LowerBound = -100.0, UpperBound = 100.0;
	double	Root;
	double	Tolerance = FLT_EPSILON;
	int		Status;

	myData.a = 1.0;
	myData.b = 5.0;
	myData.c = 4.0;
	RootBracket(*myFunction, (void *)&myData, &LowerBound, &UpperBound, Tolerance, &Status);
	RootFindBisection(*myFunction, (void *)&myData, &Root, LowerBound, UpperBound, Tolerance, &Status);
	return(0);
}

..............................................................................
*/

{ /* begin RootFindBisection */

	const double
			dx = Tolerance * (UpperBound - LowerBound);
	double	x0 = LowerBound, x1, x2 = UpperBound, old;
	double	fx0, fx1, fx2;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(RootFindBisection, Function, *Status,
/**/	"Missing function specification")
/**/DEBUG_CHECK_NULL_POINTER(RootFindBisection, AuxilliaryData, *Status,
/**/	"Missing AuxilliaryData")
/**/DEBUG_CHECK_NULL_POINTER(RootFindBisection, Root, *Status,
/**/	"Missing UpperBound")
/**/DEBUG_RETURN_ON_ERROR(RootFindBisection, *Status)
/**/DEBUG_WRITE_ENTERING(RootFindBisection,
/**/	"About to execute RootFindBisection")

	if (UpperBound <= LowerBound) {
		*Status = ERROR;
		WRITE_ERROR(RootFindBisection, "LowerBound must be strictly smaller than UpperBound")
/**/	DEBUG_WRITE_LEAVING(RootFindBisection, "Done")
		return(*Status);
	}
	if (Tolerance <= 0.0) {
		*Status = ERROR;
		WRITE_ERROR(RootFindBisection, "Invalid Tolerance (should be strictly positive)")
/**/	DEBUG_WRITE_LEAVING(RootFindBisection, "Done")
		return(*Status);
	}
	if (1.0 <= Tolerance) {
		*Status = ERROR;
		WRITE_ERROR(RootFindBisection, "Invalid Tolerance (should be strictly less than unity)")
/**/	DEBUG_WRITE_LEAVING(RootFindBisection, "Done")
		return(*Status);
	}
	*Status = (*Function)(x0, AuxilliaryData, &fx0);
	if (*Status == ERROR) {
		WRITE_ERROR(RootFindBisection, "Function call failed")
/**/	DEBUG_WRITE_LEAVING(RootFindBisection, "Done")
		return(*Status);
	}
	if (fx0 == 0.0) {
		*Root = x0;
/**/	DEBUG_WRITE_LEAVING(RootFindBisection, "Done")
		return(*Status);
	}
	*Status = (*Function)(x2, AuxilliaryData, &fx2);
	if (*Status == ERROR) {
		WRITE_ERROR(RootFindBisection, "Function call failed")
/**/	DEBUG_WRITE_LEAVING(RootFindBisection, "Done")
		return(*Status);
	}
	if (fx2 == 0.0) {
		*Root = x2;
/**/	DEBUG_WRITE_LEAVING(RootFindBisection, "Done")
		return(*Status);
	}
	if (0.0 < (fx0 * fx2)) {
		*Status = ERROR;
		WRITE_ERROR(RootFindBisection, "Non-bracketing interval")
/**/	DEBUG_WRITE_LEAVING(RootFindBisection, "Done")
		return(*Status);
	}
	old = x0;
	while (dx < (x2 - x0)) {
		x1 = 0.5 * (x0 + x2);
		if (x1 == old) {
			*Root = x1;
			WRITE_WARNING(RootFindBisection, "Unable to reach the desired precision")
/**/		DEBUG_WRITE_LEAVING(RootFindBisection, "Done")
			return(*Status);
		}
		old = x1;
		*Status = (*Function)(x1, AuxilliaryData, &fx1);
		if (*Status == ERROR) {
			WRITE_ERROR(RootFindBisection, "Function call failed")
/**/		DEBUG_WRITE_LEAVING(RootFindBisection, "Done")
			return(*Status);
		}
		if (fx1 == 0.0) {
			*Root = x1;
/**/		DEBUG_WRITE_LEAVING(RootFindBisection, "Done")
			return(*Status);
		}
		if ((fx0 * fx1) < 0.0) {
			x2 = x1;
			fx2 = fx1;
		}
		else {
			x0 = x1;
			fx0 = fx1;
		}
	}
	x1 = 0.5 * (x0 + x2);
	*Status = (*Function)(x1, AuxilliaryData, &fx1);
	if (*Status == ERROR) {
		WRITE_ERROR(RootFindBisection, "Function call failed")
/**/	DEBUG_WRITE_LEAVING(RootFindBisection, "Done")
		return(*Status);
	}
	*Root = (fabs(fx0) < fabs(fx1)) ? ((fabs(fx0) < fabs(fx2)) ? (x0) : (x2))
		: (fabs(fx1) < fabs(fx2) ? (x1) : (x2));
/**/DEBUG_WRITE_LEAVING(RootFindBisection, "Done")
	return(*Status);
} /* end RootFindBisection */

/*--------------------------------------------------------------------------*/
extern int		RootFindBrent
				(
					int		(*Function)(double, void *, double *),
												/* function, of which a root is sought */
					void	*AuxilliaryData,	/* parameters used by the function */
					double	*Root,				/* returned root */
					double	LowerBound,			/* lower bound of an interval containing a root */
					double	UpperBound,			/* upper bound of an interval containing a root */
					double	Tolerance,			/* admissible relative error */
					int		*Status				/* error management */
				)

/* search for a root of (Function) inside the bracketing interval [LowerBound, UpperBound] */
/* the strategy proceeds by fitting an inverse quadratic function and by using the root that */
/* belong to the interval, if any */
/* even-order roots generally cannot be found */
/* only one root is returned, even if there are several ones */
/* Tolerance is relative to the size of the bracketing interval */
/* success: return(!ERROR); failure: return(ERROR) */
/* the function 'Function' must be declared as follows: */
/*
		extern int myFunction(double myArgument, void *AuxilliaryData, double *myResult);
*/
/*		it must return ERROR upon failure, and !ERROR upon success */
/*		it is evaluated for the value of the variable 'myArgument' */
/*		the result of the function evaluation must be returned in 'myResult' */
/*		the generic pointer 'AuxilliaryData' can be used to pass additional parameters */
/* what follows is a developed example of the function f(x) = a * x^2 + b * x + c */
/* a, b, and c are free parameters */
/*
..............................................................................

struct myStruct
{
	double	a, b, c;
};

extern int		myFunction
				(
					double	myArgument,
					void	*AuxilliaryData,
					double	*myResult
				);
extern int		myFunction
				(
					double	myArgument,
					void	*AuxilliaryData,
					double	*myResult
				)

{
	struct myStruct
			myData;

	myData = *((struct myStruct *)AuxilliaryData);
	*myResult = myArgument * (myArgument * myData.a + myData.b) + myData.c;
	return(!ERROR);
}

int				main
				(
					void
				)

{
	struct myStruct
			myData;
	double	LowerBound = -100.0, UpperBound = 100.0;
	double	Root;
	double	Tolerance = FLT_EPSILON;
	int		Status;

	myData.a = 1.0;
	myData.b = 5.0;
	myData.c = 4.0;
	RootBracket(*myFunction, (void *)&myData, &LowerBound, &UpperBound, Tolerance, &Status);
	RootFindBrent(*myFunction, (void *)&myData, &Root, LowerBound, UpperBound, Tolerance, &Status);
	return(0);
}

..............................................................................
*/

{ /* begin RootFindBrent */

	const double
			dx = Tolerance * (UpperBound - LowerBound);
	double	x, x0 = LowerBound, x1, x2 = UpperBound, old;
	double	fx0, fx1, fx2;
	double	p, q, r, s;
	double	d = 0.0, e= 0.0;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(RootFindBrent, Function, *Status,
/**/	"Missing function specification")
/**/DEBUG_CHECK_NULL_POINTER(RootFindBrent, AuxilliaryData, *Status,
/**/	"Missing AuxilliaryData")
/**/DEBUG_CHECK_NULL_POINTER(RootFindBrent, Root, *Status,
/**/	"Missing UpperBound")
/**/DEBUG_RETURN_ON_ERROR(RootFindBrent, *Status)
/**/DEBUG_WRITE_ENTERING(RootFindBrent,
/**/	"About to execute RootFindBrent")

	if (UpperBound <= LowerBound) {
		*Status = ERROR;
		WRITE_ERROR(RootFindBrent, "LowerBound must be strictly smaller than UpperBound")
/**/	DEBUG_WRITE_LEAVING(RootFindBrent, "Done")
		return(*Status);
	}
	if (Tolerance <= 0.0) {
		*Status = ERROR;
		WRITE_ERROR(RootFindBrent, "Invalid Tolerance (should be strictly positive)")
/**/	DEBUG_WRITE_LEAVING(RootFindBrent, "Done")
		return(*Status);
	}
	if (1.0 <= Tolerance) {
		*Status = ERROR;
		WRITE_ERROR(RootFindBrent, "Invalid Tolerance (should be strictly less than unity)")
/**/	DEBUG_WRITE_LEAVING(RootFindBrent, "Done")
		return(*Status);
	}
	*Status = (*Function)(x0, AuxilliaryData, &fx0);
	if (*Status == ERROR) {
		WRITE_ERROR(RootFindBrent, "Function call failed")
/**/	DEBUG_WRITE_LEAVING(RootFindBrent, "Done")
		return(*Status);
	}
	if (fx0 == 0.0) {
		*Root = x0;
/**/	DEBUG_WRITE_LEAVING(RootFindBrent, "Done")
		return(*Status);
	}
	*Status = (*Function)(x2, AuxilliaryData, &fx2);
	if (*Status == ERROR) {
		WRITE_ERROR(RootFindBrent, "Function call failed")
/**/	DEBUG_WRITE_LEAVING(RootFindBrent, "Done")
		return(*Status);
	}
	if (fx2 == 0.0) {
		*Root = x2;
/**/	DEBUG_WRITE_LEAVING(RootFindBrent, "Done")
		return(*Status);
	}
	if (0.0 < (fx0 * fx2)) {
		*Status = ERROR;
		WRITE_ERROR(RootFindBrent, "Non-bracketing interval")
/**/	DEBUG_WRITE_LEAVING(RootFindBrent, "Done")
		return(*Status);
	}
	x1 = x2;
	fx1 = fx2;
	old = x0;
	while (TRUE) {
		if (0.0 < (fx1 * fx2)) {
			x2 = x0;
			fx2 = fx0;
			d = x1 - x0;
			e = d;
		}
		if (fabs(fx2) < fabs(fx1)) {
			x0 = x1;
			fx0 = fx1;
			x1 = x2;
			fx1 = fx2;
			x2 = x0;
			fx2 = fx0;
		}
		x = 0.5 * (x2 - x1);
		if (fabs(x) <= dx) {
			*Status = (*Function)(x1 + x, AuxilliaryData, &fx0);
			if (*Status == ERROR) {
				WRITE_ERROR(RootFindBisection, "Function call failed")
/**/			DEBUG_WRITE_LEAVING(RootFindBisection, "Done")
				return(*Status);
			}
			*Root = (fabs(fx0) < fabs(fx1)) ? ((fabs(fx0) < fabs(fx2)) ? (x0) : (x2))
				: (fabs(fx1) < fabs(fx2) ? (x1) : (x2));
/**/		DEBUG_WRITE_LEAVING(RootFindBrent, "Done")
			return(*Status);
		}
		if ((dx <= fabs(e)) && (fabs(fx1) < fabs(fx0))) {
			s = fx1 / fx0;
			if (x0 == x2) {
				p = 2.0 * x * s;
				q = 1.0 - s;
			}
			else {
				q = fx0 / fx2;
				r = fx1 / fx2;
				p = s * (2.0 * x * q * (q - r) - (x1 - x0) * (r - 1.0));
				q = (q - 1.0) * (r - 1.0) * (s - 1.0);
			}
			if (0.0 < p) {
				q = -q;
			}
			p = fabs(p);
			if (((2.0 * p) < (3.0 * x * q - dx * fabs(q))) && ((2.0 * p) < fabs(e * q))) {
				e = d;
				d = p / q;
			}
			else {
				d = x;
				e = d;
			}
		}
		else {
			d = x;
			e = d;
		}
		x0 = x1;
		fx0 = fx1;
		if (dx < fabs(d)) {
			x1 += d;
		}
		else {
			x1 += (0.0 <= x) ? (dx) : (-dx);
		}
		if (x1 == old) {
			*Root = x1;
			WRITE_WARNING(RootFindBrent, "Unable to reach the desired precision")
/**/		DEBUG_WRITE_LEAVING(RootFindBrent, "Done")
			return(*Status);
		}
		old = x1;
		*Status = (*Function)(x1, AuxilliaryData, &fx1);
		if (*Status == ERROR) {
			WRITE_ERROR(RootFindBrent, "Function call failed")
/**/		DEBUG_WRITE_LEAVING(RootFindBrent, "Done")
			return(*Status);
		}
		if (fx1 == 0.0) {
			*Root = x1;
/**/		DEBUG_WRITE_LEAVING(RootFindBrent, "Done")
			return(*Status);
		}
	}
} /* end RootFindBrent */

