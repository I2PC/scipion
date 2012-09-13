/*****************************************************************************
 *	System includes
 ****************************************************************************/
#include	<stdio.h>
#include	<stdlib.h>

/*****************************************************************************
 *	Toolbox defines
 ****************************************************************************/
#include	"configs.h"
#include	"debug.h"
#include	"error.h"
#include	"tboundaryconvention.h"

/*****************************************************************************
 *	Toolbox includes
 ****************************************************************************/
#include	"../headers/convert.h"
#include	"fold.h"
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
extern int		GetFoldedIndex
				(
					long	InputIndex,			/* index to fold back */
					long	*OutputIndex,		/* folded index */
					long	SignalLength,		/* length of the signal */
					enum TBoundaryConvention
							Convention			/* boundary convention */
				)

/* for (Periodic, MirrorOnBounds, MirrorOffBounds) boundary conditions:
	set the output index in the range [0L..SignalLength-1L] */
/* see (GetFoldedValueDouble, GetFoldedValueShort) functions
	for (FiniteDataSupport, AntiMirrorOnBounds) */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin GetFoldedIndex */

	long	N;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(GetFoldedIndex, OutputIndex, Status,
/**/	"No output index")
/**/DEBUG_CHECK_RANGE_LONG(GetFoldedIndex, SignalLength, 1L, (LONG_MAX - 1L) / 2L, Status,
/**/	"Invalid length (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(GetFoldedIndex, Status)
/**/DEBUG_WRITE_ENTERING(GetFoldedIndex,
/**/	"About to fold back an index into its legitimate range")

	switch (Convention) {
		case MirrorOffBounds:
			N = SignalLength << 1L;
			if (InputIndex < 0L) {
				InputIndex += N * ((N - InputIndex - 1L) / N);
			}
			else {
				InputIndex -= N * (InputIndex / N);
			}
			*OutputIndex = (InputIndex < SignalLength) ? (InputIndex)
				: (N - InputIndex - 1L);
			break;
		case MirrorOnBounds:
			if (SignalLength == 1L) {
				*OutputIndex = 0L;
			}
			else {
				N = SignalLength;
				N = --N << 1L;
				InputIndex = labs(InputIndex);
				InputIndex -= N * (InputIndex / N);
				*OutputIndex = (InputIndex < SignalLength) ? (InputIndex)
					: (N - InputIndex);
			}
			break;
		case Periodic:
			*OutputIndex = InputIndex;
			if (InputIndex < 0L) {
				*OutputIndex += SignalLength * ((SignalLength - InputIndex - 1L)
					/ SignalLength);
			}
			else {
				*OutputIndex -= SignalLength * (InputIndex / SignalLength);
			}
			break;
		default:
			Status = ERROR;
			WRITE_ERROR(GetFoldedIndex, "Invalid boundary convention")
			break;
	}
/**/DEBUG_WRITE_LEAVING(GetFoldedIndex, "Done")
	return(Status);
} /* end GetFoldedIndex */

/*--------------------------------------------------------------------------*/
extern int		GetFoldedValueDouble
				(
					double	Signal[],			/* double input data */
					long	InputIndex,			/* index to fold back */
					double	*OutputValue,		/* double output value */
					long	SignalLength,		/* length of the input data */
					enum TBoundaryConvention
							Convention			/* boundary convention */
				)

/* for (FiniteDataSupport, AntiMirrorOnBounds) boundary conditions:
	set a double output value */
/* see (GetFoldedValueShort) function for short */
/* see (GetFoldedIndex) function for (Periodic, MirrorOnBounds, MirrorOffBounds) */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin GetFoldedValueDouble */

	double	F0, Fn;
	long	K, N;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(GetFoldedValueDouble, Signal, Status,
/**/	"No input data")
/**/DEBUG_CHECK_NULL_POINTER(GetFoldedValueDouble, OutputValue, Status,
/**/	"No output value")
/**/DEBUG_CHECK_RANGE_LONG(GetFoldedValueDouble, SignalLength, 1L, (LONG_MAX - 1L) / 2L, Status,
/**/	"Invalid length (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(GetFoldedValueDouble, Status)
/**/DEBUG_WRITE_ENTERING(GetFoldedValueDouble,
/**/	"About to fold back an index into its legitimate range and to get the corresponding value")

	switch (Convention) {
		case AntiMirrorOnBounds:
			if (SignalLength == 1L) {
				*OutputValue = Signal[0L];
			}
			else {
				F0 = Signal[0L];
				N = SignalLength;
				Fn = Signal[--N];
				N = N << 1L;
				K = InputIndex - ((InputIndex < 0L) ? (N * ((1L + InputIndex - N) / N))
					: ((InputIndex < N) ? (0L) : (N * (InputIndex / N))));
				*OutputValue = ((K < SignalLength) ? (Signal[K]) : (2.0 * Fn - Signal[N - K]))
					+ (Fn - F0) * (double)((InputIndex - K) / (SignalLength - 1L));
			}
			break;
		case FiniteDataSupport:
			*OutputValue = (InputIndex < 0L) ? (0.0)
				: ((InputIndex < SignalLength) ? (Signal[InputIndex]) : (0.0));
			break;
		case MirrorOffBounds:
			Status = GetFoldedIndex(InputIndex, &K, SignalLength, Convention);
			if (Status == ERROR) {
/**/			DEBUG_WRITE_LEAVING(GetFoldedValueDouble, "Done")
				return(Status);
			}
			*OutputValue = Signal[K];
			break;
		case MirrorOnBounds:
			Status = GetFoldedIndex(InputIndex, &K, SignalLength, Convention);
			if (Status == ERROR) {
/**/			DEBUG_WRITE_LEAVING(GetFoldedValueDouble, "Done")
				return(Status);
			}
			*OutputValue = Signal[K];
			break;
		case Periodic:
			Status = GetFoldedIndex(InputIndex, &K, SignalLength, Convention);
			if (Status == ERROR) {
/**/			DEBUG_WRITE_LEAVING(GetFoldedValueDouble, "Done")
				return(Status);
			}
			*OutputValue = Signal[K];
			break;
		default:
			Status = ERROR;
			WRITE_ERROR(GetFoldedValueDouble, "Invalid boundary convention")
			break;
	}
/**/DEBUG_WRITE_LEAVING(GetFoldedValueDouble, "Done")
	return(Status);
} /* end GetFoldedValueDouble */

/*--------------------------------------------------------------------------*/
extern int		GetFoldedValueShort
				(
					short	Signal[],			/* short input data */
					long	InputIndex,			/* index to fold back */
					short	*OutputValue,		/* short output value */
					long	SignalLength,		/* length of the input data */
					enum TBoundaryConvention
							Convention			/* boundary convention */
				)

/* for (FiniteDataSupport, AntiMirrorOnBounds) boundary conditions:
	set a float output value */
/* see (GetFoldedValueShort) function for short */
/* see (GetFoldedIndex) function for (Periodic, MirrorOnBounds, MirrorOffBounds) */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin GetFoldedValueShort */

	double	F0, Fn;
	long	K, N;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(GetFoldedValueShort, Signal, Status,
/**/	"No input data")
/**/DEBUG_CHECK_NULL_POINTER(GetFoldedValueShort, OutputValue, Status,
/**/	"No output value")
/**/DEBUG_CHECK_RANGE_LONG(GetFoldedValueShort, SignalLength, 1L, (LONG_MAX - 1L) / 2L, Status,
/**/	"Invalid length (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(GetFoldedValueShort, Status)
/**/DEBUG_WRITE_ENTERING(GetFoldedValueShort,
/**/	"About to fold back an index into its legitimate range and to get the corresponding value")

	switch (Convention) {
		case AntiMirrorOnBounds:
			if (SignalLength == 1L) {
				*OutputValue = Signal[0L];
			}
			else {
				F0 = (double)Signal[0L];
				N = SignalLength;
				Fn = (double)Signal[--N];
				N = N << 1L;
				K = InputIndex - ((InputIndex < 0L) ? (N * ((1L + InputIndex - N) / N))
					: ((InputIndex < N) ? (0L) : (N * (InputIndex / N))));
				*OutputValue = ConvertDoubleToShort(((K < SignalLength) ? ((double)Signal[K])
					: (2.0 * Fn - (double)Signal[N - K])) + (Fn - F0)
					* (double)((InputIndex - K) / (SignalLength - 1L)));
			}
			break;
		case FiniteDataSupport:
			*OutputValue = (InputIndex < 0L) ? ((short)0)
				: ((SignalLength <= InputIndex) ? ((short)0) : (Signal[InputIndex]));
			break;
		case MirrorOffBounds:
			Status = GetFoldedIndex(InputIndex, &K, SignalLength, Convention);
			if (Status == ERROR) {
/**/			DEBUG_WRITE_LEAVING(GetFoldedValueShort, "Done")
				return(Status);
			}
			*OutputValue = Signal[K];
			break;
		case MirrorOnBounds:
			Status = GetFoldedIndex(InputIndex, &K, SignalLength, Convention);
			if (Status == ERROR) {
/**/			DEBUG_WRITE_LEAVING(GetFoldedValueShort, "Done")
				return(Status);
			}
			*OutputValue = Signal[K];
			break;
		case Periodic:
			Status = GetFoldedIndex(InputIndex, &K, SignalLength, Convention);
			if (Status == ERROR) {
/**/			DEBUG_WRITE_LEAVING(GetFoldedValueShort, "Done")
				return(Status);
			}
			*OutputValue = Signal[K];
			break;
		default:
			Status = ERROR;
			WRITE_ERROR(GetFoldedValueShort, "Invalid boundary convention")
			break;
	}
/**/DEBUG_WRITE_LEAVING(GetFoldedValueShort, "Done")
	return(Status);
} /* end GetFoldedValueShort */

