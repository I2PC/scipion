/*****************************************************************************
 *	System includes
 ****************************************************************************/
#include	<math.h>
#include	<stddef.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>

/*****************************************************************************
 *	Toolbox defines
 ****************************************************************************/
#include	"configs.h"
#include	"debug.h"
#include        "error.h"

/*****************************************************************************
 *	Toolbox includes
 ****************************************************************************/
#include	"dft.h"
#include	"dht.h"
#include        "getputd.h"

/*****************************************************************************
 *	Conditional includes
 ****************************************************************************/
#include	<limits.h>
#include	<stdio.h>
#include	"messagedisplay.h"

/*****************************************************************************
 *	Local defines
 ****************************************************************************/
/*--------------------------------------------------------------------------*/
#undef		FORWARD
#define		FORWARD				(TRUE)
#undef		BACKWARD
#define		BACKWARD			(!FORWARD)

/*****************************************************************************
 *	Declaration of static procedures
 ****************************************************************************/
/*--------------------------------------------------------------------------*/
static int		bufferAllocR
				(
					double	*(H[]),
					double	*(tmp[]),
					double	*(CaS[]),
					long	SignalLength
				);

/*--------------------------------------------------------------------------*/
static int		bufferAllocRI
				(
					double	*(Re[]),
					double	*(Im[]),
					double	*(TmpRe[]),
					double	*(TmpIm[]),
					double	*(CaS[]),
					long	SignalLength
				);

/*--------------------------------------------------------------------------*/
static int		bufferFreeR
				(
					double	*(H[]),
					double	*(tmp[]),
					double	*(CaS[])
				);

/*--------------------------------------------------------------------------*/
static int		bufferFreeRI
				(
					double	*(Re[]),
					double	*(Im[]),
					double	*(TmpRe[]),
					double	*(TmpIm[]),
					double	*(CaS[])
				);

/*****************************************************************************
 *	Definition of static procedures
 ****************************************************************************/
/*--------------------------------------------------------------------------*/
static int		bufferAllocR
				(
					double	*(H[]),
					double	*(tmp[]),
					double	*(CaS[]),
					long	SignalLength
				)

{ /* Begin bufferAllocR */

	int		Status = !ERROR;
/**/DEBUG_WRITE_ENTERING(bufferAllocR,
/**/	"About to execute bufferAllocR")
	if (*H != (double *)NULL) {
		Status = ERROR;
		WRITE_ERROR(bufferAllocR, "Unexpected nonempty H")
/**/	DEBUG_WRITE_LEAVING(bufferAllocR, "Done")
		return(Status);
	}
	*H = (double *)malloc((size_t)SignalLength * sizeof(double));
	if (*H == (double *)NULL) {
		Status = ERROR;
		WRITE_ERROR(bufferAllocR, "Not enough memory for H")
/**/	DEBUG_WRITE_LEAVING(bufferAllocR, "Done")
		return(Status);
	}
	if (*tmp != (double *)NULL) {
		Status = ERROR;
		WRITE_ERROR(bufferAllocR, "Unexpected nonempty tmp")
/**/	DEBUG_WRITE_LEAVING(bufferAllocR, "Done")
		return(Status);
	}
	*tmp = (double *)malloc((size_t)SignalLength * sizeof(double));
	if (*tmp == (double *)NULL) {
		Status = ERROR;
		WRITE_ERROR(bufferAllocR, "Not enough memory for tmp")
/**/	DEBUG_WRITE_LEAVING(bufferAllocR, "Done")
		return(Status);
	}
	if (*CaS != (double *)NULL) {
		Status = ERROR;
		WRITE_ERROR(bufferAllocR, "Unexpected nonempty CaS")
/**/	DEBUG_WRITE_LEAVING(bufferAllocR, "Done")
		return(Status);
	}
	*CaS = (double *)malloc((size_t)SignalLength * sizeof(double));
	if (*CaS == (double *)NULL) {
		Status = ERROR;
		WRITE_ERROR(bufferAllocR, "Not enough memory for CaS")
/**/	DEBUG_WRITE_LEAVING(bufferAllocR, "Done")
		return(Status);
	}
/**/DEBUG_WRITE_LEAVING(bufferAllocR, "Done")
	return(Status);
} /* End bufferAllocR */

/*--------------------------------------------------------------------------*/
static int		bufferAllocRI
				(
					double	*(Re[]),
					double	*(Im[]),
					double	*(TmpRe[]),
					double	*(TmpIm[]),
					double	*(CaS[]),
					long	SignalLength
				)

{ /* Begin bufferAllocRI */

	int		Status = !ERROR;

/**/DEBUG_WRITE_ENTERING(bufferAllocRI,
/**/	"About to execute bufferAllocRI")
	if (*Re != (double *)NULL) {
		Status = ERROR;
		WRITE_ERROR(bufferAllocRI, "Unexpected nonempty Re")
/**/	DEBUG_WRITE_LEAVING(bufferAllocRI, "Done")
		return(Status);
	}
	*Re = (double *)malloc((size_t)SignalLength * sizeof(double));
	if (*Re == (double *)NULL) {
		Status = ERROR;
		WRITE_ERROR(bufferAllocRI, "Not enough memory for Re")
/**/	DEBUG_WRITE_LEAVING(bufferAllocRI, "Done")
		return(Status);
	}
	if (*Im != (double *)NULL) {
		Status = ERROR;
		WRITE_ERROR(bufferAllocRI, "Unexpected nonempty Im")
/**/	DEBUG_WRITE_LEAVING(bufferAllocRI, "Done")
		return(Status);
	}
	*Im = (double *)malloc((size_t)SignalLength * sizeof(double));
	if (*Im == (double *)NULL) {
		Status = ERROR;
		WRITE_ERROR(bufferAllocRI, "Not enough memory for Im")
/**/	DEBUG_WRITE_LEAVING(bufferAllocRI, "Done")
		return(Status);
	}
	if (*TmpRe != (double *)NULL) {
		Status = ERROR;
		WRITE_ERROR(bufferAllocRI, "Unexpected nonempty TmpRe")
/**/	DEBUG_WRITE_LEAVING(bufferAllocRI, "Done")
		return(Status);
	}
	*TmpRe = (double *)malloc((size_t)SignalLength * sizeof(double));
	if (*TmpRe == (double *)NULL) {
		Status = ERROR;
		WRITE_ERROR(bufferAllocRI, "Not enough memory for TmpRe")
/**/	DEBUG_WRITE_LEAVING(bufferAllocRI, "Done")
		return(Status);
	}
	if (*TmpIm != (double *)NULL) {
		Status = ERROR;
		WRITE_ERROR(bufferAllocRI, "Unexpected nonempty TmpIm")
/**/	DEBUG_WRITE_LEAVING(bufferAllocRI, "Done")
		return(Status);
	}
	*TmpIm = (double *)malloc((size_t)SignalLength * sizeof(double));
	if (*TmpIm == (double *)NULL) {
		Status = ERROR;
		WRITE_ERROR(bufferAllocRI, "Not enough memory for TmpIm")
/**/	DEBUG_WRITE_LEAVING(bufferAllocRI, "Done")
		return(Status);
	}
	if (*CaS != (double *)NULL) {
		Status = ERROR;
		WRITE_ERROR(bufferAllocRI, "Unexpected nonempty CaS")
/**/	DEBUG_WRITE_LEAVING(bufferAllocRI, "Done")
		return(Status);
	}
	*CaS = (double *)malloc((size_t)SignalLength * sizeof(double));
	if (*CaS == (double *)NULL) {
		Status = ERROR;
		WRITE_ERROR(bufferAllocRI, "Not enough memory for CaS")
/**/	DEBUG_WRITE_LEAVING(bufferAllocRI, "Done")
		return(Status);
	}
/**/DEBUG_WRITE_LEAVING(bufferAllocRI, "Done")
	return(Status);
} /* End bufferAllocRI */

/*--------------------------------------------------------------------------*/
static int		bufferFreeR
				(
					double	*(H[]),
					double	*(tmp[]),
					double	*(CaS[])
				)

{ /* Begin bufferFreeR */

	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(bufferFreeR, H, Status,
/**/	"Missing H")
/**/DEBUG_CHECK_NULL_POINTER(bufferFreeR, tmp, Status,
/**/	"Missing tmp")
/**/DEBUG_CHECK_NULL_POINTER(bufferFreeR, CaS, Status,
/**/	"Missing CaS")
/**/DEBUG_RETURN_ON_ERROR(bufferFreeR, Status)
/**/DEBUG_WRITE_ENTERING(bufferFreeR,
/**/	"About to execute bufferFreeR")

	free(*H);
	*H = (double *)NULL;
	free(*tmp);
	*tmp = (double *)NULL;
	free(*CaS);
	*CaS = (double *)NULL;
/**/DEBUG_WRITE_LEAVING(bufferFreeR, "Done")
	return(Status);
} /* End bufferFreeR */

/*--------------------------------------------------------------------------*/
static int		bufferFreeRI
				(
					double	*(Re[]),
					double	*(Im[]),
					double	*(TmpRe[]),
					double	*(TmpIm[]),
					double	*(CaS[])
				)

{ /* Begin bufferFreeRI */

	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(bufferFreeRI, Re, Status,
/**/	"Missing Re")
/**/DEBUG_CHECK_NULL_POINTER(bufferFreeRI, Im, Status,
/**/	"Missing Im")
/**/DEBUG_CHECK_NULL_POINTER(bufferFreeRI, TmpRe, Status,
/**/	"Missing TmpRe")
/**/DEBUG_CHECK_NULL_POINTER(bufferFreeRI, TmpIm, Status,
/**/	"Missing TmpIm")
/**/DEBUG_CHECK_NULL_POINTER(bufferFreeRI, CaS, Status,
/**/	"Missing CaS")
/**/DEBUG_RETURN_ON_ERROR(bufferFreeRI, Status)
/**/DEBUG_WRITE_ENTERING(bufferFreeRI,
/**/	"About to execute bufferFreeRI")

	free(*Re);
	*Re = (double *)NULL;
	free(*Im);
	*Im = (double *)NULL;
	free(*TmpRe);
	*TmpRe = (double *)NULL;
	free(*TmpIm);
	*TmpIm = (double *)NULL;
	free(*CaS);
	*CaS = (double *)NULL;
/**/DEBUG_WRITE_LEAVING(bufferFreeRI, "Done")
	return(Status);
} /* End bufferFreeRI */

/*****************************************************************************
 *	Definition of extern procedures
 ****************************************************************************/
/*--------------------------------------------------------------------------*/
extern int		AmplitudePhaseToRealImaginary
				(
					double	Am2Re[],			/* (amplitude -> real) */
					double	Ph2Im[],			/* (phase -> imaginary) */
					long	SignalLength		/* signal length */
				)

/* converts an (amplitude, phase) representation of a complex signal
	into a (real, imaginary) representation */
/* the input phase is in [rad] */
/* in-place processing */
/* the input signal (amplitude Am2Re and phase Ph2Im) is
	replaced by the output signal (real Am2Re and imaginary Ph2Im) */
/* SignalLength is the signal length */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin AmplitudePhaseToRealImaginary */

	double	Phase;
	long	i;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(AmplitudePhaseToRealImaginary, Am2Re, Status,
/**/	"Missing Am2Re ")
/**/DEBUG_CHECK_NULL_POINTER(AmplitudePhaseToRealImaginary, Ph2Im, Status,
/**/	"Missing Ph2Im ")
/**/DEBUG_CHECK_RANGE_LONG(AmplitudePhaseToRealImaginary, SignalLength, 1L, LONG_MAX, Status,
/**/	"Invalid length (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(AmplitudePhaseToRealImaginary, Status)
/**/DEBUG_WRITE_ENTERING(AmplitudePhaseToRealImaginary,
/**/	"About to convert an amplitude/phase representation to a real/imaginary representation")

	for (i = -SignalLength; (i < 0L); i++) {
		Phase = *Ph2Im;
		*Ph2Im++ = *Am2Re * sin(Phase);
		*Am2Re++ *= cos(Phase);
	}
/**/DEBUG_WRITE_LEAVING(AmplitudePhaseToRealImaginary, "Done")
	return(Status);
} /* end AmplitudePhaseToRealImaginary */

/*--------------------------------------------------------------------------*/
extern int		DftAmplitudePhaseToAmplitudePhase
				(
					double	Am2Am[],			/* amplitude -> amplitude */
					double	Ph2Ph[],			/* phase -> phase */
					double	*TmpRe,				/* first scratch workspace */
					double	*TmpIm,				/* second scratch workspace */
					double	CaS[],				/* Hartley transform coefficients */
					long	SignalLength		/* signal length */
				)

/* computes the direct DFT of a complex signal given in (amplitude, phase) representation
	and returns an (amplitude, phase) representation */
/* the input phase is in [rad] */
/* the output phase is in [rad]; its domain is (-PI, PI) */
/* the origin is at index [0] */
/* the highest coordinate (SignalLength-2)/2 is at index [(SignalLength-2)/2] for SignalLength even */
/* the highest coordinate (SignalLength-1)/2 is at index [(SignalLength-1)/2] for SignalLength odd */
/* the lowest coordinate -SignalLength/2 is at index [SignalLength/2] for SignalLength even */
/* the lowest coordinate -(SignalLength-1)/2 is at index [(SignalLength+1)/2] for SignalLength odd */
/* the coordinate -1 is at index [SignalLength-1] for SignalLength even or odd */
/* in-place processing */
/* the input signal (amplitude Am2Am and phase Ph2Ph) is
	replaced by the output signal (amplitude Am2Am and phase Ph2Ph) */
/* (TmpRe, TmpIm) are pre-allocated workspaces of size SignalLength each */
/* the values returned in (TmpRe, TmpIm) are meaningless */
/* CaS is an input array of coefficients of size SignalLength (see GetCaS function) */
/* CaS is modified internally, but is restored when DftAmplitudePhaseToAmplitudePhase returns */
/* SignalLength is the length of the signal */
/* no restriction on SignalLength, but best efficiency for radix 2, 3, 4, 5 */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin DftAmplitudePhaseToAmplitudePhase */

	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(DftAmplitudePhaseToAmplitudePhase, Am2Am, Status,
/**/	"Missing Am2Am ")
/**/DEBUG_CHECK_NULL_POINTER(DftAmplitudePhaseToAmplitudePhase, Ph2Ph, Status,
/**/	"Missing Ph2Ph ")
/**/DEBUG_CHECK_NULL_POINTER(DftAmplitudePhaseToAmplitudePhase, TmpRe, Status,
/**/	"Missing TmpRe ")
/**/DEBUG_CHECK_NULL_POINTER(DftAmplitudePhaseToAmplitudePhase, TmpIm, Status,
/**/	"Missing TmpIm ")
/**/DEBUG_CHECK_NULL_POINTER(DftAmplitudePhaseToAmplitudePhase, CaS, Status,
/**/	"Missing CaS ")
/**/DEBUG_CHECK_RANGE_LONG(DftAmplitudePhaseToAmplitudePhase, SignalLength, 1L, LONG_MAX, Status,
/**/	"Invalid length (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(DftAmplitudePhaseToAmplitudePhase, Status)
/**/DEBUG_WRITE_ENTERING(DftAmplitudePhaseToAmplitudePhase,
/**/	"About to compute an amplitude/phase to amplitude/phase direct Fourier transform")

	Status = AmplitudePhaseToRealImaginary(Am2Am, Ph2Ph, SignalLength);
	if (Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(DftAmplitudePhaseToAmplitudePhase, "Done")
		return(Status);
	}
	Status = DftRealImaginaryToAmplitudePhase(Am2Am, Ph2Ph, TmpRe, TmpIm, CaS,
		SignalLength);
/**/DEBUG_WRITE_LEAVING(DftAmplitudePhaseToAmplitudePhase, "Done")
	return(Status);
} /* end DftAmplitudePhaseToAmplitudePhase */

/*--------------------------------------------------------------------------*/
extern int		DftAmplitudePhaseToRealImaginary
				(
					double	Am2Re[],			/* amplitude -> real */
					double	Ph2Im[],			/* phase -> imaginary */
					double	*TmpRe,				/* first scratch workspace */
					double	*TmpIm,				/* second scratch workspace */
					double	CaS[],				/* Hartley transform coefficients */
					long	SignalLength		/* signal length */
				)

/* computes the direct DFT of a complex signal given in (amplitude, phase) representation
	and returns a (real, imaginary) representation */
/* the input phase is in [rad] */
/* the origin is at index [0] */
/* the highest coordinate (SignalLength-2)/2 is at index [(SignalLength-2)/2] for SignalLength even */
/* the highest coordinate (SignalLength-1)/2 is at index [(SignalLength-1)/2] for SignalLength odd */
/* the lowest coordinate -SignalLength/2 is at index [SignalLength/2] for SignalLength even */
/* the lowest coordinate -(SignalLength-1)/2 is at index [(SignalLength+1)/2] for SignalLength odd */
/* the coordinate -1 is at index [SignalLength-1] for SignalLength even or odd */
/* in-place processing */
/* the input signal (amplitude Am2Re and phase Ph2Im) is
	replaced by the output signal (real Am2Re and imaginary Ph2Im) */
/* (TmpRe, TmpIm) are pre-allocated workspaces of size SignalLength each */
/* the values returned in (TmpRe, TmpIm) are meaningless */
/* CaS is an input array of coefficients of size SignalLength (see GetCaS function) */
/* CaS is modified internally, but is restored when DftAmplitudePhaseToRealImaginary returns */
/* SignalLength is the length of the signal */
/* no restriction on SignalLength, but best efficiency for radix 2, 3, 4, 5 */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin DftAmplitudePhaseToRealImaginary */

	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(DftAmplitudePhaseToRealImaginary, Am2Re, Status,
/**/	"Missing Am2Re ")
/**/DEBUG_CHECK_NULL_POINTER(DftAmplitudePhaseToRealImaginary, Ph2Im, Status,
/**/	"Missing Ph2Im ")
/**/DEBUG_CHECK_NULL_POINTER(DftAmplitudePhaseToRealImaginary, TmpRe, Status,
/**/	"Missing TmpRe ")
/**/DEBUG_CHECK_NULL_POINTER(DftAmplitudePhaseToRealImaginary, TmpIm, Status,
/**/	"Missing TmpIm ")
/**/DEBUG_CHECK_NULL_POINTER(DftAmplitudePhaseToRealImaginary, CaS, Status,
/**/	"Missing CaS ")
/**/DEBUG_CHECK_RANGE_LONG(DftAmplitudePhaseToRealImaginary, SignalLength, 1L, LONG_MAX, Status,
/**/	"Invalid length (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(DftAmplitudePhaseToRealImaginary, Status)
/**/DEBUG_WRITE_ENTERING(DftAmplitudePhaseToRealImaginary,
/**/	"About to compute an amplitude/phase to real/imaginary direct Fourier transform")

	Status = AmplitudePhaseToRealImaginary(Am2Re, Ph2Im, SignalLength);
	if (Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(DftAmplitudePhaseToRealImaginary, "Done")
		return(Status);
	}
	Status = DftRealImaginaryToRealImaginary(Am2Re, Ph2Im, TmpRe, TmpIm, CaS,
		SignalLength);
/**/DEBUG_WRITE_LEAVING(DftAmplitudePhaseToRealImaginary, "Done")
	return(Status);
} /* end DftAmplitudePhaseToRealImaginary */

/*--------------------------------------------------------------------------*/
extern int		DftRealImaginaryToAmplitudePhase
				(
					double	Re2Am[],			/* real -> amplitude */
					double	Im2Ph[],			/* imaginary -> phase */
					double	*TmpRe,				/* first scratch workspace */
					double	*TmpIm,				/* second scratch workspace */
					double	CaS[],				/* Hartley transform coefficients */
					long	SignalLength		/* signal length */
				)

/* computes the direct DFT of a complex signal given in (real, imaginary) representation
	and returns an (amplitude, phase) representation */
/* the output phase is in [rad]; its domain is (-PI, PI) */
/* the origin is at index [0] */
/* the highest coordinate (SignalLength-2)/2 is at index [(SignalLength-2)/2] for SignalLength even */
/* the highest coordinate (SignalLength-1)/2 is at index [(SignalLength-1)/2] for SignalLength odd */
/* the lowest coordinate -SignalLength/2 is at index [SignalLength/2] for SignalLength even */
/* the lowest coordinate -(SignalLength-1)/2 is at index [(SignalLength+1)/2] for SignalLength odd */
/* the coordinate -1 is at index [SignalLength-1] for SignalLength even or odd */
/* in-place processing */
/* the input signal (real Re2Am and imaginary Im2Ph) is
	replaced by the output signal (amplitude Re2Am and phase Im2Ph) */
/* (TmpRe, TmpIm) are pre-allocated workspaces of size SignalLength each */
/* the values returned in (TmpRe, TmpIm) are meaningless */
/* CaS is an input array of coefficients of size SignalLength (see GetCaS function) */
/* CaS is modified internally, but is restored when DftRealImaginaryToAmplitudePhase returns */
/* SignalLength is the length of the signal */
/* no restriction on SignalLength, but best efficiency for radix 2, 3, 4, 5 */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin DftRealImaginaryToAmplitudePhase */

	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(DftRealImaginaryToAmplitudePhase, Re2Am, Status,
/**/	"Missing Re2Am ")
/**/DEBUG_CHECK_NULL_POINTER(DftRealImaginaryToAmplitudePhase, Im2Ph, Status,
/**/	"Missing Im2Ph ")
/**/DEBUG_CHECK_NULL_POINTER(DftRealImaginaryToAmplitudePhase, TmpRe, Status,
/**/	"Missing TmpRe ")
/**/DEBUG_CHECK_NULL_POINTER(DftRealImaginaryToAmplitudePhase, TmpIm, Status,
/**/	"Missing TmpIm ")
/**/DEBUG_CHECK_NULL_POINTER(DftRealImaginaryToAmplitudePhase, CaS, Status,
/**/	"Missing CaS ")
/**/DEBUG_CHECK_RANGE_LONG(DftRealImaginaryToAmplitudePhase, SignalLength, 1L, LONG_MAX, Status,
/**/	"Invalid length (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(DftRealImaginaryToAmplitudePhase, Status)
/**/DEBUG_WRITE_ENTERING(DftRealImaginaryToAmplitudePhase,
/**/	"About to compute a real/imaginary to amplitude/phase direct Fourier transform")

	Status = DftRealImaginaryToRealImaginary(Re2Am, Im2Ph, TmpRe, TmpIm, CaS,
		SignalLength);
	if (Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(DftRealImaginaryToAmplitudePhase, "Done")
		return(Status);
	}
	Status = RealImaginaryToAmplitudePhase(Re2Am, Im2Ph, SignalLength);
/**/DEBUG_WRITE_LEAVING(DftRealImaginaryToAmplitudePhase, "Done")
	return(Status);
} /* end DftRealImaginaryToAmplitudePhase */

/*--------------------------------------------------------------------------*/
extern int		DftRealImaginaryToRealImaginary
				(
					double	Re2Re[],			/* real -> real */
					double	Im2Im[],			/* imaginary -> imaginary */
					double	*TmpRe,				/* first scratch workspace */
					double	*TmpIm,				/* second scratch workspace */
					double	CaS[],				/* Hartley transform coefficients */
					long	SignalLength		/* signal length */
				)

/* computes the direct DFT of a complex signal given in (real, imaginary) representation
	and returns a (real, imaginary) representation */
/* the origin is at index [0] */
/* the highest coordinate (SignalLength-2)/2 is at index [(SignalLength-2)/2] for SignalLength even */
/* the highest coordinate (SignalLength-1)/2 is at index [(SignalLength-1)/2] for SignalLength odd */
/* the lowest coordinate -SignalLength/2 is at index [SignalLength/2] for SignalLength even */
/* the lowest coordinate -(SignalLength-1)/2 is at index [(SignalLength+1)/2] for SignalLength odd */
/* the coordinate -1 is at index [SignalLength-1] for SignalLength even or odd */
/* in-place processing */
/* the input signal (real Re2Re and imaginary Im2Im) is
	replaced by the output signal (real Re2Re and imaginary Im2Im) */
/* (TmpRe, TmpIm) are pre-allocated workspaces of size SignalLength each */
/* the values returned in (TmpRe, TmpIm) are meaningless */
/* CaS is an input array of coefficients of size SignalLength (see GetCaS function) */
/* CaS is modified internally, but is restored when DftRealImaginaryToRealImaginary returns */
/* SignalLength is the length of the signal */
/* no restriction on SignalLength, but best efficiency for radix 2, 3, 4, 5 */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin DftRealImaginaryToRealImaginary */

	double	*ReBack, *ImBack;
	long	i;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(DftRealImaginaryToRealImaginary, Re2Re, Status,
/**/	"Missing Re2Re ")
/**/DEBUG_CHECK_NULL_POINTER(DftRealImaginaryToRealImaginary, Im2Im, Status,
/**/	"Missing Im2Im ")
/**/DEBUG_CHECK_NULL_POINTER(DftRealImaginaryToRealImaginary, TmpRe, Status,
/**/	"Missing TmpRe ")
/**/DEBUG_CHECK_NULL_POINTER(DftRealImaginaryToRealImaginary, TmpIm, Status,
/**/	"Missing TmpIm ")
/**/DEBUG_CHECK_NULL_POINTER(DftRealImaginaryToRealImaginary, CaS, Status,
/**/	"Missing CaS ")
/**/DEBUG_CHECK_RANGE_LONG(DftRealImaginaryToRealImaginary, SignalLength, 1L, LONG_MAX, Status,
/**/	"Invalid length (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(DftRealImaginaryToRealImaginary, Status)
/**/DEBUG_WRITE_ENTERING(DftRealImaginaryToRealImaginary,
/**/	"About to compute a real/imaginary to real/imaginary direct Fourier transform")

	TmpRe = (double *)memcpy(TmpRe, Re2Re, (size_t)(SignalLength
		* (long)sizeof(double)));
	TmpIm = (double *)memcpy(TmpIm, Im2Im, (size_t)(SignalLength
		* (long)sizeof(double)));
	Status = DiscreteHartleyTransform(TmpRe, Re2Re, CaS, SignalLength, FORWARD);
	if (Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(DftRealImaginaryToRealImaginary, "Done")
		return(Status);
	}
	Status = DiscreteHartleyTransform(TmpIm, Im2Im, CaS, SignalLength, FORWARD);
	if (Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(DftRealImaginaryToRealImaginary, "Done")
		return(Status);
	}
	ReBack = TmpRe + (ptrdiff_t)(SignalLength - 1L);
	ImBack = TmpIm + (ptrdiff_t)(SignalLength - 1L);
	*Re2Re = *TmpRe;
	*Im2Im = *TmpIm;
	Re2Re += (ptrdiff_t)(SignalLength - 1L);
	Im2Im += (ptrdiff_t)(SignalLength - 1L);
	TmpRe++;
	TmpIm++;
	for (i = 1L - SignalLength; (i < 0L); i++) {
		*Re2Re-- = (*TmpRe + *ReBack - *TmpIm + *ImBack) * 0.5;
		*Im2Im-- = (*TmpRe++ - *ReBack-- + *TmpIm++ + *ImBack--) * 0.5;
	}
/**/DEBUG_WRITE_LEAVING(DftRealImaginaryToRealImaginary, "Done")
	return(Status);
} /* end DftRealImaginaryToRealImaginary */

/*--------------------------------------------------------------------------*/
extern int		DftRealToAmplitudePhase
				(
					double	Re2Am[],			/* real -> amplitude */
					double	PhOut[],			/* output phase */
					double	*Tmp,				/* scratch workspace */
					double	CaS[],				/* Hartley transform coefficients */
					long	SignalLength		/* signal length */
				)

/* computes the direct DFT of a real signal
	and returns an (amplitude, phase) representation */
/* the output phase is in [rad]; its domain is (-PI, PI) */
/* the origin is at index [0] */
/* the highest coordinate (SignalLength-2)/2 is at index [(SignalLength-2)/2] for SignalLength even */
/* the highest coordinate (SignalLength-1)/2 is at index [(SignalLength-1)/2] for SignalLength odd */
/* the lowest coordinate -SignalLength/2 is at index [SignalLength/2] for SignalLength even */
/* the lowest coordinate -(SignalLength-1)/2 is at index [(SignalLength+1)/2] for SignalLength odd */
/* the coordinate -1 is at index [SignalLength-1] for SignalLength even or odd */
/* in-place processing */
/* the input signal (real Re2Am) is
	replaced by the output signal (amplitude Re2Am and phase PhOut) */
/* Tmp is a pre-allocated workspace of size SignalLength */
/* the values returned in Tmp are meaningless */
/* CaS is an input array of coefficients of size SignalLength (see GetCaS function) */
/* CaS is modified internally, but is restored when DftRealToAmplitudePhase returns */
/* SignalLength is the length of the signal */
/* no restriction on SignalLength, but best efficiency for radix 2, 3, 4, 5 */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin DftRealToAmplitudePhase */

	double	*AmForw, *PhForw;
	double	*AmBack, *PhBack;
	double	Re, Im;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(DftRealToAmplitudePhase, Re2Am, Status,
/**/	"Missing Re2Am ")
/**/DEBUG_CHECK_NULL_POINTER(DftRealToAmplitudePhase, PhOut, Status,
/**/	"Missing PhOut ")
/**/DEBUG_CHECK_NULL_POINTER(DftRealToAmplitudePhase, Tmp, Status,
/**/	"Missing Tmp ")
/**/DEBUG_CHECK_NULL_POINTER(DftRealToAmplitudePhase, CaS, Status,
/**/	"Missing CaS ")
/**/DEBUG_CHECK_RANGE_LONG(DftRealToAmplitudePhase, SignalLength, 1L, LONG_MAX, Status,
/**/	"Invalid length (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(DftRealToAmplitudePhase, Status)
/**/DEBUG_WRITE_ENTERING(DftRealToAmplitudePhase,
/**/	"About to compute a real to amplitude/phase direct Fourier transform")

	Status = DiscreteHartleyTransform(Re2Am, Tmp, CaS, SignalLength, FORWARD);
	if (Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(DftRealToAmplitudePhase, "Done")
		return(Status);
	}
	AmForw = Re2Am;
	AmBack = Re2Am + (ptrdiff_t)(SignalLength - 1L);
	PhForw = PhOut;
	PhBack = PhOut + (ptrdiff_t)(SignalLength - 1L);
	if (0.0 <= *AmForw) {
		*PhForw = 0.0;
	}
	else {
		*AmForw = -*AmForw;
		*PhForw = PI;
	}
	PhForw++;
	AmForw++;
	while (AmForw < AmBack) {
		Im = *AmBack - *AmForw;
		Re = *AmBack + *AmForw;
		*PhForw = atan2(Im, Re);
		*PhBack-- = -(*PhForw++);
		*AmForw = sqrt((Re * Re + Im * Im) * 0.25);
		*AmBack-- = *AmForw++;
	}
	if (AmForw == AmBack) {
		if (0.0 <= *AmForw) {
			*PhForw = 0.0;
		}
		else {
			*AmForw = -*AmForw;
			*PhForw = PI;
		}
	}
/**/DEBUG_WRITE_LEAVING(DftRealToAmplitudePhase, "Done")
	return(Status);
} /* end DftRealToAmplitudePhase */

/*--------------------------------------------------------------------------*/
extern int		DftRealToRealImaginary
				(
					double	Re2Re[],			/* real -> real */
					double	ImOut[],			/* output imaginary */
					double	*Tmp,				/* scratch workspace */
					double	CaS[],				/* Hartley transform coefficients */
					long	SignalLength		/* signal length */
				)

/* computes the direct DFT of a real signal
	and returns a (real, imaginary) representation */
/* the origin is at index [0] */
/* the highest coordinate (SignalLength-2)/2 is at index [(SignalLength-2)/2] for SignalLength even */
/* the highest coordinate (SignalLength-1)/2 is at index [(SignalLength-1)/2] for SignalLength odd */
/* the lowest coordinate -SignalLength/2 is at index [SignalLength/2] for SignalLength even */
/* the lowest coordinate -(SignalLength-1)/2 is at index [(SignalLength+1)/2] for SignalLength odd */
/* the coordinate -1 is at index [SignalLength-1] for SignalLength even or odd */
/* in-place processing */
/* the input signal (real Re2Re) is
	replaced by the output signal (real Re2Re and imaginary ImOut) */
/* Tmp is a pre-allocated workspace of size SignalLength */
/* the values returned in Tmp are meaningless */
/* CaS is an input array of coefficients of size SignalLength (see GetCaS function) */
/* CaS is modified internally, but is restored when DftRealToRealImaginary returns */
/* SignalLength is the length of the signal */
/* no restriction on SignalLength, but best efficiency for radix 2, 3, 4, 5 */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin DftRealToRealImaginary */

	double	*ReForw, *ImForw;
	double	*ReBack, *ImBack;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(DftRealToRealImaginary, Re2Re, Status,
/**/	"Missing Re2Re ")
/**/DEBUG_CHECK_NULL_POINTER(DftRealToRealImaginary, ImOut, Status,
/**/	"Missing ImOut ")
/**/DEBUG_CHECK_NULL_POINTER(DftRealToRealImaginary, Tmp, Status,
/**/	"Missing Tmp ")
/**/DEBUG_CHECK_NULL_POINTER(DftRealToRealImaginary, CaS, Status,
/**/	"Missing CaS ")
/**/DEBUG_CHECK_RANGE_LONG(DftRealToRealImaginary, SignalLength, 1L, LONG_MAX, Status,
/**/	"Invalid length (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(DftRealToRealImaginary, Status)
/**/DEBUG_WRITE_ENTERING(DftRealToRealImaginary,
/**/	"About to compute a real to real/imaginary direct Fourier transform")

	Status = DiscreteHartleyTransform(Re2Re, Tmp, CaS, SignalLength, FORWARD);
	if (Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(DftRealToRealImaginary, "Done")
		return(Status);
	}
	ReForw = Re2Re;
	ReBack = Re2Re + (ptrdiff_t)(SignalLength - 1L);
	ImForw = ImOut;
	ImBack = ImOut + (ptrdiff_t)(SignalLength - 1L);
	*ImForw++ = 0.0;
	ReForw++;
	while (ReForw < ReBack) {
		*ImForw = (*ReBack - *ReForw) * 0.5;
		*ImBack-- = -(*ImForw++);
		*ReForw = (*ReBack + *ReForw) * 0.5;
		*ReBack-- = *ReForw++;
	}
	if (ReForw == ReBack) {
		*ImForw = 0.0;
	}
/**/DEBUG_WRITE_LEAVING(DftRealToRealImaginary, "Done")
	return(Status);
} /* end DftRealToRealImaginary */

/*--------------------------------------------------------------------------*/
extern int		GetCaS
				(
					double	CaS[],				/* Hartley transform coefficients */
					long	SignalLength		/* signal length */
				)

/* computes an array of coefficients of size SignalLength */
/* these coefficients are necessary for performing a Hartley transform */
/* the Hartley transform is an alternate representation of Fourier for real signals */
/* Hartley computations are more accurate (less roundoff errors) than Fourier */
/* the same coefficients are used for direct and inverse transforms */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin GetCaS */

	const double
			s2 = sqrt(2.0), pi2 = 2.0 * PI, pi4 = PI * 0.25;
	long	i;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(GetCaS, CaS, Status,
/**/	"Missing CaS ")
/**/DEBUG_CHECK_RANGE_LONG(GetCaS, SignalLength, 1L, LONG_MAX, Status,
/**/	"Invalid length (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(GetCaS, Status)
/**/DEBUG_WRITE_ENTERING(GetCaS,
/**/	"About to compute Hartley coefficients")

	for (i = 0L; (i < SignalLength); i++) {
		*CaS++ = s2 * sin(pi4 + pi2 * (double)i / (double)SignalLength);
	}
/**/DEBUG_WRITE_LEAVING(GetCaS, "Done")
	return(Status);
} /* end GetCaS */

/*--------------------------------------------------------------------------*/
extern int		InvDftAmplitudePhaseToAmplitudePhase
				(
					double	Am2Am[],			/* amplitude -> amplitude */
					double	Ph2Ph[],			/* phase -> phase */
					double	*TmpRe,				/* first scratch workspace */
					double	*TmpIm,				/* second scratch workspace */
					double	CaS[],				/* Hartley transform coefficients */
					long	SignalLength		/* signal length */
				)

/* computes the inverse DFT of a complex signal given in (amplitude, phase) representation
	and returns an (amplitude, phase) representation */
/* the input phase is in [rad] */
/* the output phase is in [rad]; its domain is (-PI, PI) */
/* the origin is at index [0] */
/* the highest coordinate (SignalLength-2)/2 is at index [(SignalLength-2)/2] for SignalLength even */
/* the highest coordinate (SignalLength-1)/2 is at index [(SignalLength-1)/2] for SignalLength odd */
/* the lowest coordinate -SignalLength/2 is at index [SignalLength/2] for SignalLength even */
/* the lowest coordinate -(SignalLength-1)/2 is at index [(SignalLength+1)/2] for SignalLength odd */
/* the coordinate -1 is at index [SignalLength-1] for SignalLength even or odd */
/* in-place processing */
/* the input signal (amplitude Am2Am and phase Ph2Ph) is
	replaced by the output signal (amplitude Am2Am and phase Ph2Ph) */
/* (TmpRe, TmpIm) are pre-allocated workspaces of size SignalLength each */
/* the values returned in (TmpRe, TmpIm) are meaningless */
/* CaS is an input array of coefficients of size SignalLength (see GetCaS function) */
/* CaS is modified internally, but is restored when DftAmplitudePhaseToAmplitudePhase returns */
/* SignalLength is the length of the signal */
/* no restriction on SignalLength, but best efficiency for radix 2, 3, 4, 5 */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin InvDftAmplitudePhaseToAmplitudePhase */

	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(InvDftAmplitudePhaseToAmplitudePhase, Am2Am, Status,
/**/	"Missing Am2Am ")
/**/DEBUG_CHECK_NULL_POINTER(InvDftAmplitudePhaseToAmplitudePhase, Ph2Ph, Status,
/**/	"Missing Ph2Ph ")
/**/DEBUG_CHECK_NULL_POINTER(InvDftAmplitudePhaseToAmplitudePhase, TmpRe, Status,
/**/	"Missing TmpRe ")
/**/DEBUG_CHECK_NULL_POINTER(InvDftAmplitudePhaseToAmplitudePhase, TmpIm, Status,
/**/	"Missing TmpIm ")
/**/DEBUG_CHECK_NULL_POINTER(InvDftAmplitudePhaseToAmplitudePhase, CaS, Status,
/**/	"Missing CaS ")
/**/DEBUG_CHECK_RANGE_LONG(InvDftAmplitudePhaseToAmplitudePhase, SignalLength, 1L, LONG_MAX, Status,
/**/	"Invalid length (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(InvDftAmplitudePhaseToAmplitudePhase, Status)
/**/DEBUG_WRITE_ENTERING(InvDftAmplitudePhaseToAmplitudePhase,
/**/	"About to compute an amplitude/phase to amplitude/phase inverse Fourier transform")

	Status = AmplitudePhaseToRealImaginary(Am2Am, Ph2Ph, SignalLength);
	if (Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(InvDftAmplitudePhaseToAmplitudePhase, "Done")
		return(Status);
	}
	Status = InvDftRealImaginaryToAmplitudePhase(Am2Am, Ph2Ph, TmpRe, TmpIm, CaS,
		SignalLength);
/**/DEBUG_WRITE_LEAVING(InvDftAmplitudePhaseToAmplitudePhase, "Done")
	return(Status);
} /* end InvDftAmplitudePhaseToAmplitudePhase */

/*--------------------------------------------------------------------------*/
extern int		InvDftAmplitudePhaseToReal
				(
					double	Am2Re[],			/* amplitude -> real */
					double	PhIn[],				/* input phase */
					double	*Tmp,				/* scratch workspace */
					double	CaS[],				/* Hartley transform coefficients */
					long	SignalLength		/* signal length */
				)

/* computes the inverse DFT of a complex signal given in (amplitude, phase) representation
	and returns a real signal */
/* the complex Fourier signal is symmetrized before the inverse transformation is applied */
/* the input phase is in [rad] */
/* the origin is at index [0] */
/* the highest coordinate (SignalLength-2)/2 is at index [(SignalLength-2)/2] for SignalLength even */
/* the highest coordinate (SignalLength-1)/2 is at index [(SignalLength-1)/2] for SignalLength odd */
/* the lowest coordinate -SignalLength/2 is at index [SignalLength/2] for SignalLength even */
/* the lowest coordinate -(SignalLength-1)/2 is at index [(SignalLength+1)/2] for SignalLength odd */
/* the coordinate -1 is at index [SignalLength-1] for SignalLength even or odd */
/* in-place processing */
/* the input signal (amplitude Am2Re and phase PhIn) is
	replaced by the output signal (real Am2Re) */
/* Tmp is a pre-allocated workspace of size SignalLength */
/* the values returned in Tmp are meaningless */
/* CaS is an input array of coefficients of size SignalLength (see GetCaS function) */
/* CaS is modified internally, but is restored when InvDftAmplitudePhaseToReal returns */
/* SignalLength is the length of the signal */
/* no restriction on SignalLength, but best efficiency for radix 2, 3, 4, 5 */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin InvDftAmplitudePhaseToReal */

	double	*AmForw, *PhForw;
	double	*AmBack, *PhBack;
	double	Re, Im;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(InvDftAmplitudePhaseToReal, Am2Re, Status,
/**/	"Missing Am2Re ")
/**/DEBUG_CHECK_NULL_POINTER(InvDftAmplitudePhaseToReal, PhIn, Status,
/**/	"Missing PhIn ")
/**/DEBUG_CHECK_NULL_POINTER(InvDftAmplitudePhaseToReal, Tmp, Status,
/**/	"Missing Tmp ")
/**/DEBUG_CHECK_NULL_POINTER(InvDftAmplitudePhaseToReal, CaS, Status,
/**/	"Missing CaS ")
/**/DEBUG_CHECK_RANGE_LONG(InvDftAmplitudePhaseToReal, SignalLength, 1L, LONG_MAX, Status,
/**/	"Invalid length (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(InvDftAmplitudePhaseToReal, Status)
/**/DEBUG_WRITE_ENTERING(InvDftAmplitudePhaseToReal,
/**/	"About to compute an amplitude/phase to real inverse Fourier transform")

	AmForw = Am2Re;
	AmBack = Am2Re + (ptrdiff_t)(SignalLength - 1L);
	PhForw = PhIn;
	PhBack = PhIn + (ptrdiff_t)(SignalLength - 1L);
	*Am2Re = *Am2Re * cos(*PhIn);
	PhForw++;
	AmForw++;
	while (AmForw < AmBack) {
		Re = *AmForw * cos(*PhForw) + *AmBack * cos(*PhBack);
		Im = *AmForw * sin(*PhForw++) - *AmBack * sin(*PhBack--);
		*AmForw++ = (Re - Im) * 0.5;
		*AmBack-- = (Re + Im) * 0.5;
	}
	if (AmForw == AmBack) {
		*AmForw *= cos(*PhForw);
	}
	Status = DiscreteHartleyTransform(Am2Re, Tmp, CaS, SignalLength, FALSE);
/**/DEBUG_WRITE_LEAVING(InvDftAmplitudePhaseToReal, "Done")
	return(Status);
} /* end InvDftAmplitudePhaseToReal */

/*--------------------------------------------------------------------------*/
extern int		InvDftAmplitudePhaseToRealImaginary
				(
					double	Am2Re[],			/* amplitude -> real */
					double	Ph2Im[],			/* phase -> imaginary */
					double	*TmpRe,				/* first scratch workspace */
					double	*TmpIm,				/* second scratch workspace */
					double	CaS[],				/* Hartley transform coefficients */
					long	SignalLength		/* signal length */
				)

/* computes the inverse DFT of a complex signal given in (amplitude, phase) representation
	and returns a (real, imaginary) representation */
/* the input phase is in [rad] */
/* the origin is at index [0] */
/* the highest coordinate (SignalLength-2)/2 is at index [(SignalLength-2)/2] for SignalLength even */
/* the highest coordinate (SignalLength-1)/2 is at index [(SignalLength-1)/2] for SignalLength odd */
/* the lowest coordinate -SignalLength/2 is at index [SignalLength/2] for SignalLength even */
/* the lowest coordinate -(SignalLength-1)/2 is at index [(SignalLength+1)/2] for SignalLength odd */
/* the coordinate -1 is at index [SignalLength-1] for SignalLength even or odd */
/* in-place processing */
/* the input signal (amplitude Am2Re and phase Ph2Im) is
	replaced by the output signal (real Am2Re and imaginary Ph2Im) */
/* (TmpRe, TmpIm) are pre-allocated workspaces of size SignalLength each */
/* the values returned in (TmpRe, TmpIm) are meaningless */
/* CaS is an input array of coefficients of size SignalLength (see GetCaS function) */
/* CaS is modified internally, but is restored when InvDftAmplitudePhaseToRealImaginary returns */
/* SignalLength is the length of the signal */
/* no restriction on SignalLength, but best efficiency for radix 2, 3, 4, 5 */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin InvDftAmplitudePhaseToRealImaginary */

	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(InvDftAmplitudePhaseToRealImaginary, Am2Re, Status,
/**/	"Missing Am2Re ")
/**/DEBUG_CHECK_NULL_POINTER(InvDftAmplitudePhaseToRealImaginary, Ph2Im, Status,
/**/	"Missing Ph2Im ")
/**/DEBUG_CHECK_NULL_POINTER(InvDftAmplitudePhaseToRealImaginary, TmpRe, Status,
/**/	"Missing TmpRe ")
/**/DEBUG_CHECK_NULL_POINTER(InvDftAmplitudePhaseToRealImaginary, TmpIm, Status,
/**/	"Missing TmpIm ")
/**/DEBUG_CHECK_NULL_POINTER(InvDftAmplitudePhaseToRealImaginary, CaS, Status,
/**/	"Missing CaS ")
/**/DEBUG_CHECK_RANGE_LONG(InvDftAmplitudePhaseToRealImaginary, SignalLength, 1L, LONG_MAX, Status,
/**/	"Invalid length (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(InvDftAmplitudePhaseToRealImaginary, Status)
/**/DEBUG_WRITE_ENTERING(InvDftAmplitudePhaseToRealImaginary,
/**/	"About to compute an amplitude/phase to real/imaginary inverse Fourier transform")

	Status = AmplitudePhaseToRealImaginary(Am2Re, Ph2Im, SignalLength);
	if (Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(InvDftAmplitudePhaseToRealImaginary, "Done")
		return(Status);
	}
	Status = InvDftRealImaginaryToRealImaginary(Am2Re, Ph2Im, TmpRe, TmpIm, CaS,
		SignalLength);
/**/DEBUG_WRITE_LEAVING(InvDftAmplitudePhaseToRealImaginary, "Done")
	return(Status);
} /* end InvDftAmplitudePhaseToRealImaginary */

/*--------------------------------------------------------------------------*/
extern int		InvDftRealImaginaryToAmplitudePhase
				(
					double	Re2Am[],			/* real -> amplitude */
					double	Im2Ph[],			/* imaginary -> phase */
					double	*TmpRe,				/* first scratch workspace */
					double	*TmpIm,				/* second scratch workspace */
					double	CaS[],				/* Hartley transform coefficients */
					long	SignalLength		/* signal length */
				)

/* computes the inverse DFT of a complex signal given in (real, imaginary) representation
	and returns an (amplitude, phase) representation */
/* the output phase is in [rad]; its domain is (-PI, PI) */
/* the origin is at index [0] */
/* the highest coordinate (SignalLength-2)/2 is at index [(SignalLength-2)/2] for SignalLength even */
/* the highest coordinate (SignalLength-1)/2 is at index [(SignalLength-1)/2] for SignalLength odd */
/* the lowest coordinate -SignalLength/2 is at index [SignalLength/2] for SignalLength even */
/* the lowest coordinate -(SignalLength-1)/2 is at index [(SignalLength+1)/2] for SignalLength odd */
/* the coordinate -1 is at index [SignalLength-1] for SignalLength even or odd */
/* in-place processing */
/* the input signal (real Re2Am and imaginary Im2Ph) is
	replaced by the output signal (amplitude Re2Am and phase Im2Ph) */
/* (TmpRe, TmpIm) are pre-allocated workspaces of size SignalLength each */
/* the values returned in (TmpRe, TmpIm) are meaningless */
/* CaS is an input array of coefficients of size SignalLength (see GetCaS function) */
/* CaS is modified internally, but is restored when InvDftRealImaginaryToAmplitudePhase returns */
/* SignalLength is the length of the signal */
/* no restriction on SignalLength, but best efficiency for radix 2, 3, 4, 5 */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin InvDftRealImaginaryToAmplitudePhase */

	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(InvDftRealImaginaryToAmplitudePhase, Re2Am, Status,
/**/	"Missing Re2Am ")
/**/DEBUG_CHECK_NULL_POINTER(InvDftRealImaginaryToAmplitudePhase, Im2Ph, Status,
/**/	"Missing Im2Ph ")
/**/DEBUG_CHECK_NULL_POINTER(InvDftRealImaginaryToAmplitudePhase, TmpRe, Status,
/**/	"Missing TmpRe ")
/**/DEBUG_CHECK_NULL_POINTER(InvDftRealImaginaryToAmplitudePhase, TmpIm, Status,
/**/	"Missing TmpIm ")
/**/DEBUG_CHECK_NULL_POINTER(InvDftRealImaginaryToAmplitudePhase, CaS, Status,
/**/	"Missing CaS ")
/**/DEBUG_CHECK_RANGE_LONG(InvDftRealImaginaryToAmplitudePhase, SignalLength, 1L, LONG_MAX, Status,
/**/	"Invalid length (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(InvDftRealImaginaryToAmplitudePhase, Status)
/**/DEBUG_WRITE_ENTERING(InvDftRealImaginaryToAmplitudePhase,
/**/	"About to compute a real/imaginary to amplitude/phase inverse Fourier transform")

	Status = InvDftRealImaginaryToRealImaginary(Re2Am, Im2Ph, TmpRe, TmpIm, CaS,
		SignalLength);
	if (Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(InvDftRealImaginaryToAmplitudePhase, "Done")
		return(Status);
	}
	Status = RealImaginaryToAmplitudePhase(Re2Am, Im2Ph, SignalLength);
/**/DEBUG_WRITE_LEAVING(InvDftRealImaginaryToAmplitudePhase, "Done")
	return(Status);
} /* end InvDftRealImaginaryToAmplitudePhase */

/*--------------------------------------------------------------------------*/
extern int		InvDftRealImaginaryToReal
				(
					double	Re2Re[],			/* real -> real */
					double	ImIn[],				/* input imaginary */
					double	*Tmp,				/* scratch workspace */
					double	CaS[],				/* Hartley transform coefficients */
					long	SignalLength		/* signal length */
				)

/* computes the inverse DFT of a complex signal given in (real, imaginary) representation
	and returns a real signal */
/* the complex Fourier signal is symmetrized before the inverse transformation is applied */
/* the origin is at index [0] */
/* the highest coordinate (SignalLength-2)/2 is at index [(SignalLength-2)/2] for SignalLength even */
/* the highest coordinate (SignalLength-1)/2 is at index [(SignalLength-1)/2] for SignalLength odd */
/* the lowest coordinate -SignalLength/2 is at index [SignalLength/2] for SignalLength even */
/* the lowest coordinate -(SignalLength-1)/2 is at index [(SignalLength+1)/2] for SignalLength odd */
/* the coordinate -1 is at index [SignalLength-1] for SignalLength even or odd */
/* in-place processing */
/* the input signal (real Re2Re and imaginary ImIn) is
	replaced by the output signal (real Re2Re) */
/* Tmp is a pre-allocated workspace of size SignalLength */
/* the values returned in Tmp are meaningless */
/* CaS is an input array of coefficients of size SignalLength (see GetCaS function) */
/* CaS is modified internally, but is restored when InvDftRealImaginaryToReal returns */
/* SignalLength is the length of the signal */
/* no restriction on SignalLength, but best efficiency for radix 2, 3, 4, 5 */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin InvDftRealImaginaryToReal */

	double	*ReForw, *ImForw;
	double	*ReBack, *ImBack;
	double	Re, Im;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(InvDftRealImaginaryToReal, Re2Re, Status,
/**/	"Missing Re2Re ")
/**/DEBUG_CHECK_NULL_POINTER(InvDftRealImaginaryToReal, ImIn, Status,
/**/	"Missing ImIn ")
/**/DEBUG_CHECK_NULL_POINTER(InvDftRealImaginaryToReal, Tmp, Status,
/**/	"Missing Tmp ")
/**/DEBUG_CHECK_NULL_POINTER(InvDftRealImaginaryToReal, CaS, Status,
/**/	"Missing CaS ")
/**/DEBUG_CHECK_RANGE_LONG(InvDftRealImaginaryToReal, SignalLength, 1L, LONG_MAX, Status,
/**/	"Invalid length (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(InvDftRealImaginaryToReal, Status)
/**/DEBUG_WRITE_ENTERING(InvDftRealImaginaryToReal,
/**/	"About to compute a real/imaginary to real inverse Fourier transform")

	ReForw = Re2Re;
	ReBack = Re2Re + (ptrdiff_t)(SignalLength - 1L);
	ImForw = ImIn;
	ImBack = ImIn + (ptrdiff_t)(SignalLength - 1L);
	ImForw++;
	ReForw++;
	while (ReForw < ReBack) {
		Re = *ReForw + *ReBack;
		Im = *ImForw++ - *ImBack--;
		*ReForw++ = (Re - Im) * 0.5;
		*ReBack-- = (Re + Im) * 0.5;
	}
	Status = DiscreteHartleyTransform(Re2Re, Tmp, CaS, SignalLength, BACKWARD);
/**/DEBUG_WRITE_LEAVING(InvDftRealImaginaryToReal, "Done")
	return(Status);
} /* end InvDftRealImaginaryToReal */

/*--------------------------------------------------------------------------*/
extern int		InvDftRealImaginaryToRealImaginary
				(
					double	Re2Re[],			/* real -> real */
					double	Im2Im[],			/* imaginary -> imaginary */
					double	*TmpRe,				/* first scratch workspace */
					double	*TmpIm,				/* second scratch workspace */
					double	CaS[],				/* Hartley transform coefficients */
					long	SignalLength		/* signal length */
				)

/* computes the inverse DFT of a complex signal given in (real, imaginary) representation
	and returns a (real, imaginary) representation */
/* the origin is at index [0] */
/* the highest coordinate (SignalLength-2)/2 is at index [(SignalLength-2)/2] for SignalLength even */
/* the highest coordinate (SignalLength-1)/2 is at index [(SignalLength-1)/2] for SignalLength odd */
/* the lowest coordinate -SignalLength/2 is at index [SignalLength/2] for SignalLength even */
/* the lowest coordinate -(SignalLength-1)/2 is at index [(SignalLength+1)/2] for SignalLength odd */
/* the coordinate -1 is at index [SignalLength-1] for SignalLength even or odd */
/* in-place processing */
/* the input signal (real Re2Re and imaginary Im2Im) is
	replaced by the output signal (real Re2Re and imaginary Im2Im) */
/* (TmpRe, TmpIm) are pre-allocated workspaces of size SignalLength each */
/* the values returned in (TmpRe, TmpIm) are meaningless */
/* CaS is an input array of coefficients of size SignalLength (see GetCaS function) */
/* CaS is modified internally, but is restored when InvDftRealImaginaryToRealImaginary returns */
/* SignalLength is the length of the signal */
/* no restriction on SignalLength, but best efficiency for radix 2, 3, 4, 5 */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin InvDftRealImaginaryToRealImaginary */

	double	*ReBack, *ImBack;
	long	i;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(InvDftRealImaginaryToRealImaginary, Re2Re, Status,
/**/	"Missing Re2Re ")
/**/DEBUG_CHECK_NULL_POINTER(InvDftRealImaginaryToRealImaginary, Im2Im, Status,
/**/	"Missing Im2Im ")
/**/DEBUG_CHECK_NULL_POINTER(InvDftRealImaginaryToRealImaginary, TmpRe, Status,
/**/	"Missing TmpRe ")
/**/DEBUG_CHECK_NULL_POINTER(InvDftRealImaginaryToRealImaginary, TmpIm, Status,
/**/	"Missing TmpIm ")
/**/DEBUG_CHECK_NULL_POINTER(InvDftRealImaginaryToRealImaginary, CaS, Status,
/**/	"Missing CaS ")
/**/DEBUG_CHECK_RANGE_LONG(InvDftRealImaginaryToRealImaginary, SignalLength, 1L, LONG_MAX, Status,
/**/	"Invalid length (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(InvDftRealImaginaryToRealImaginary, Status)
/**/DEBUG_WRITE_ENTERING(InvDftRealImaginaryToRealImaginary,
/**/	"About to compute a real/imaginary to real/imaginary inverse Fourier transform")

	TmpRe = (double *)memcpy(TmpRe, Re2Re, (size_t)(SignalLength
		* (long)sizeof(double)));
	TmpIm = (double *)memcpy(TmpIm, Im2Im, (size_t)(SignalLength
		* (long)sizeof(double)));
	Status = DiscreteHartleyTransform(TmpRe, Re2Re, CaS, SignalLength, BACKWARD);
	if (Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(InvDftRealImaginaryToRealImaginary, "Done")
		return(Status);
	}
	Status = DiscreteHartleyTransform(TmpIm, Im2Im, CaS, SignalLength, BACKWARD);
	if (Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(InvDftRealImaginaryToRealImaginary, "Done")
		return(Status);
	}
	ReBack = TmpRe + (ptrdiff_t)(SignalLength - 1L);
	ImBack = TmpIm + (ptrdiff_t)(SignalLength - 1L);
	*Re2Re = *TmpRe;
	*Im2Im = *TmpIm;
	Re2Re += (ptrdiff_t)(SignalLength - 1L);
	Im2Im += (ptrdiff_t)(SignalLength - 1L);
	TmpRe++;
	TmpIm++;
	for (i = 1L - SignalLength; (i < 0L); i++) {
		*Re2Re-- = (*ReBack + *TmpRe - *ImBack + *TmpIm) * 0.5;
		*Im2Im-- = (*ReBack-- - *TmpRe++ + *ImBack-- + *TmpIm++) * 0.5;
	}
/**/DEBUG_WRITE_LEAVING(InvDftRealImaginaryToRealImaginary, "Done")
	return(Status);
} /* end InvDftRealImaginaryToRealImaginary */

/*--------------------------------------------------------------------------*/
extern int		RealImaginaryToAmplitudePhase
				(
					double	Re2Am[],			/* real -> amplitude */
					double	Im2Ph[],			/* imaginary -> phase */
					long	SignalLength		/* signal length */
				)

/* converts a (real, imaginary) representation of a complex signal
	into an (amplitude, phase) representation */
/* the output phase is in [rad]; its domain is (-PI, PI) */
/* in-place processing */
/* the input signal (real Re2Am and imaginary Im2Ph) is
	replaced by the output signal (amplitude Re2Am and phase Im2Ph) */
/* SignalLength is the signal length */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin RealImaginaryToAmplitudePhase */

	double	Phase;
	long	i;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(RealImaginaryToAmplitudePhase, Re2Am, Status,
/**/	"Missing Re2Am ")
/**/DEBUG_CHECK_NULL_POINTER(RealImaginaryToAmplitudePhase, Im2Ph, Status,
/**/	"Missing Im2Ph ")
/**/DEBUG_CHECK_RANGE_LONG(RealImaginaryToAmplitudePhase, SignalLength, 1L, LONG_MAX, Status,
/**/	"Invalid length (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(RealImaginaryToAmplitudePhase, Status)
/**/DEBUG_WRITE_ENTERING(RealImaginaryToAmplitudePhase,
/**/	"About to convert a real/imaginary representation to an amplitude/phase representation")

	for (i = -SignalLength; (i < 0L); Re2Am++, i++) {
		Phase = atan2(*Im2Ph, *Re2Am);
		*Re2Am = sqrt(*Re2Am * *Re2Am + *Im2Ph * *Im2Ph);
		*Im2Ph++ = Phase;
	}
/**/DEBUG_WRITE_LEAVING(RealImaginaryToAmplitudePhase, "Done")
	return(Status);
} /* end RealImaginaryToAmplitudePhase */

/*--------------------------------------------------------------------------*/
extern int		VolumeAmplitudePhaseToRealImaginary
				(
					double	Am2Re[],			/* (amplitude -> real) */
					double	Ph2Im[],			/* (phase -> imaginary) */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz					/* depth of the volume */
				)

/* converts an (amplitude, phase) representation of a complex signal
	into a (real, imaginary) representation */
/* the input phase is in [rad] */
/* in-place processing */
/* the input signal (amplitude Am2Re and phase Ph2Im) is
	replaced by the output signal (real Am2Re and imaginary Ph2Im) */
/* Nx is the width of the volume */
/* Ny is the height of the volume */
/* Nz is the depth of the volume */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin VolumeAmplitudePhaseToRealImaginary */

	double	Phase;
	long	i;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(VolumeAmplitudePhaseToRealImaginary, Am2Re, Status,
/**/	"Missing Am2Re ")
/**/DEBUG_CHECK_NULL_POINTER(VolumeAmplitudePhaseToRealImaginary, Ph2Im, Status,
/**/	"Missing Ph2Im ")
/**/DEBUG_CHECK_RANGE_LONG(VolumeAmplitudePhaseToRealImaginary, Nx, 1L, LONG_MAX, Status,
/**/	"Invalid length (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(VolumeAmplitudePhaseToRealImaginary, Ny, 1L, LONG_MAX, Status,
/**/	"Invalid length (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(VolumeAmplitudePhaseToRealImaginary, Nz, 1L, LONG_MAX, Status,
/**/	"Invalid length (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(VolumeAmplitudePhaseToRealImaginary, Status)
/**/DEBUG_WRITE_ENTERING(VolumeAmplitudePhaseToRealImaginary,
/**/	"About to convert an amplitude/phase representation to a real/imaginary representation")

	for (i = -Nx * Ny * Nz; (i < 0L); i++) {
		Phase = *Ph2Im;
		*Ph2Im++ = (double)(*Am2Re * sin(Phase));
		*Am2Re++ *= (double)cos(Phase);
	}
/**/DEBUG_WRITE_LEAVING(VolumeAmplitudePhaseToRealImaginary, "Done")
	return(Status);
} /* end VolumeAmplitudePhaseToRealImaginary */

/*--------------------------------------------------------------------------*/
extern int		VolumeDftAmplitudePhaseToAmplitudePhase
				(
					double	*Am2Am,				/* amplitude -> amplitude */
					double	*Ph2Ph,				/* phase -> phase */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					int		*Status				/* error management */
				)

/* computes the direct DFT of a complex signal given in (amplitude, phase) representation
	and returns an (amplitude, phase) representation */
/* the input phase is in [rad] */
/* the output phase is in [rad]; its domain is (-PI, PI) */
/* the origin is at index [0] */
/* the highest coordinate (SignalLength-2)/2 is at index [(SignalLength-2)/2] for SignalLength even */
/* the highest coordinate (SignalLength-1)/2 is at index [(SignalLength-1)/2] for SignalLength odd */
/* the lowest coordinate -SignalLength/2 is at index [SignalLength/2] for SignalLength even */
/* the lowest coordinate -(SignalLength-1)/2 is at index [(SignalLength+1)/2] for SignalLength odd */
/* the coordinate -1 is at index [SignalLength-1] for SignalLength even or odd */
/* no restriction on SignalLength, but best efficiency for radix 2, 3, 4, 5 */
/* in the explanations above, SignalLength has to be replaced by Nx, Ny, Nz */
/* Nx is the width of the volume */
/* Ny is the height of the volume */
/* Nz is the depth of the volume */
/* in-place processing */
/* the input signal (amplitude Am2Am and phase Ph2Ph) is
	replaced by the output signal (amplitude Am2Am and phase Ph2Ph) */
/* success: return(!ERROR); failure: return(ERROR); */
/* the returned value is duplicated in Status */

{ /* begin VolumeDftAmplitudePhaseToAmplitudePhase */

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(VolumeDftAmplitudePhaseToAmplitudePhase, Am2Am, *Status,
/**/	"Missing Am2Am ")
/**/DEBUG_CHECK_NULL_POINTER(VolumeDftAmplitudePhaseToAmplitudePhase, Ph2Ph, *Status,
/**/	"Missing Ph2Ph ")
/**/DEBUG_CHECK_RANGE_LONG(VolumeDftAmplitudePhaseToAmplitudePhase, Nx, 1L, LONG_MAX, *Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(VolumeDftAmplitudePhaseToAmplitudePhase, Ny, 1L, LONG_MAX, *Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(VolumeDftAmplitudePhaseToAmplitudePhase, Nz, 1L, LONG_MAX, *Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(VolumeDftAmplitudePhaseToAmplitudePhase, *Status)
/**/DEBUG_WRITE_ENTERING(VolumeDftAmplitudePhaseToAmplitudePhase,
/**/	"About to compute a 3D amplitude/phase to amplitude/phase direct Fourier transform")

	*Status = VolumeAmplitudePhaseToRealImaginary(Am2Am, Ph2Ph, Nx, Ny, Nz);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(VolumeDftAmplitudePhaseToAmplitudePhase, "Done")
		return(*Status);
	}
	VolumeDftRealImaginaryToAmplitudePhase(Am2Am, Ph2Ph, Nx, Ny, Nz, Status);
/**/DEBUG_WRITE_LEAVING(VolumeDftAmplitudePhaseToAmplitudePhase, "Done")
	return(*Status);
} /* end VolumeDftAmplitudePhaseToAmplitudePhase */

/*--------------------------------------------------------------------------*/
extern int		VolumeDftAmplitudePhaseToRealImaginary
				(
					double	*Am2Re,				/* amplitude -> real */
					double	*Ph2Im,				/* phase -> imaginary */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					int		*Status				/* error management */
				)

/* computes the direct DFT of a complex signal given in (amplitude, phase) representation
	and returns a (real, imaginary) representation */
/* the input phase is in [rad] */
/* the origin is at index [0] */
/* the highest coordinate (SignalLength-2)/2 is at index [(SignalLength-2)/2] for SignalLength even */
/* the highest coordinate (SignalLength-1)/2 is at index [(SignalLength-1)/2] for SignalLength odd */
/* the lowest coordinate -SignalLength/2 is at index [SignalLength/2] for SignalLength even */
/* the lowest coordinate -(SignalLength-1)/2 is at index [(SignalLength+1)/2] for SignalLength odd */
/* the coordinate -1 is at index [SignalLength-1] for SignalLength even or odd */
/* no restriction on SignalLength, but best efficiency for radix 2, 3, 4, 5 */
/* in the explanations above, SignalLength has to be replaced by Nx, Ny, Nz */
/* Nx is the width of the volume */
/* Ny is the height of the volume */
/* Nz is the depth of the volume */
/* in-place processing */
/* the input signal (amplitude Am2Re and phase Ph2Im) is
	replaced by the output signal (real Am2Re and imaginary Ph2Im) */
/* success: return(!ERROR); failure: return(ERROR); */
/* the returned value is duplicated in Status */

{ /* begin VolumeDftAmplitudePhaseToRealImaginary */

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(VolumeDftAmplitudePhaseToRealImaginary, Am2Re, *Status,
/**/	"Missing Am2Re ")
/**/DEBUG_CHECK_NULL_POINTER(VolumeDftAmplitudePhaseToRealImaginary, Ph2Im, *Status,
/**/	"Missing Ph2Im ")
/**/DEBUG_CHECK_RANGE_LONG(VolumeDftAmplitudePhaseToRealImaginary, Nx, 1L, LONG_MAX, *Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(VolumeDftAmplitudePhaseToRealImaginary, Ny, 1L, LONG_MAX, *Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(VolumeDftAmplitudePhaseToRealImaginary, Nz, 1L, LONG_MAX, *Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(VolumeDftAmplitudePhaseToRealImaginary, *Status)
/**/DEBUG_WRITE_ENTERING(VolumeDftAmplitudePhaseToRealImaginary,
/**/	"About to compute a 3D amplitude/phase to real/imaginary direct Fourier transform")

	*Status = VolumeAmplitudePhaseToRealImaginary(Am2Re, Ph2Im, Nx, Ny, Nz);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(VolumeDftAmplitudePhaseToRealImaginary, "Done")
		return(*Status);
	}
	VolumeDftRealImaginaryToRealImaginary(Am2Re, Ph2Im, Nx, Ny, Nz, Status);
/**/DEBUG_WRITE_LEAVING(VolumeDftAmplitudePhaseToRealImaginary, "Done")
	return(*Status);
} /* end VolumeDftAmplitudePhaseToRealImaginary */

/*--------------------------------------------------------------------------*/
extern int		VolumeDftRealImaginaryToAmplitudePhase
				(
					double	*Re2Am,				/* real -> amplitude */
					double	*Im2Ph,				/* imaginary -> phase */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					int		*Status				/* error management */
				)

/* computes the direct DFT of a complex signal given in (real, imaginary) representation
	and returns an (amplitude, phase) representation */
/* the output phase is in [rad]; its domain is (-PI, PI) */
/* the origin is at index [0] */
/* the highest coordinate (SignalLength-2)/2 is at index [(SignalLength-2)/2] for SignalLength even */
/* the highest coordinate (SignalLength-1)/2 is at index [(SignalLength-1)/2] for SignalLength odd */
/* the lowest coordinate -SignalLength/2 is at index [SignalLength/2] for SignalLength even */
/* the lowest coordinate -(SignalLength-1)/2 is at index [(SignalLength+1)/2] for SignalLength odd */
/* the coordinate -1 is at index [SignalLength-1] for SignalLength even or odd */
/* no restriction on SignalLength, but best efficiency for radix 2, 3, 4, 5 */
/* in the explanations above, SignalLength has to be replaced by Nx, Ny, Nz */
/* Nx is the width of the volume */
/* Ny is the height of the volume */
/* Nz is the depth of the volume */
/* in-place processing */
/* the input signal (real Re2Am and imaginary Im2Ph) is
	replaced by the output signal (amplitude Re2Am and phase Im2Ph) */
/* success: return(!ERROR); failure: return(ERROR); */
/* the returned value is duplicated in Status */

{ /* begin VolumeDftRealImaginaryToAmplitudePhase */

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(VolumeDftRealImaginaryToAmplitudePhase, Re2Am, *Status,
/**/	"Missing Re2Am ")
/**/DEBUG_CHECK_NULL_POINTER(VolumeDftRealImaginaryToAmplitudePhase, Im2Ph, *Status,
/**/	"Missing Im2Ph ")
/**/DEBUG_CHECK_RANGE_LONG(VolumeDftRealImaginaryToAmplitudePhase, Nx, 1L, LONG_MAX, *Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(VolumeDftRealImaginaryToAmplitudePhase, Ny, 1L, LONG_MAX, *Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(VolumeDftRealImaginaryToAmplitudePhase, Nz, 1L, LONG_MAX, *Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(VolumeDftRealImaginaryToAmplitudePhase, *Status)
/**/DEBUG_WRITE_ENTERING(VolumeDftRealImaginaryToAmplitudePhase,
/**/	"About to compute a 3D real/imaginary to amplitude/phase direct Fourier transform")

	VolumeDftRealImaginaryToRealImaginary(Re2Am, Im2Ph, Nx, Ny, Nz, Status);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(VolumeDftRealImaginaryToAmplitudePhase, "Done")
		return(*Status);
	}
	*Status = VolumeRealImaginaryToAmplitudePhase(Re2Am, Im2Ph, Nx, Ny, Nz);
/**/DEBUG_WRITE_LEAVING(VolumeDftRealImaginaryToAmplitudePhase, "Done")
	return(*Status);
} /* end VolumeDftRealImaginaryToAmplitudePhase */

/*--------------------------------------------------------------------------*/
extern int		VolumeDftRealImaginaryToRealImaginary
				(
					double	*Re2Re,				/* real -> real */
					double	*Im2Im,				/* imaginary -> imaginary */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					int		*Status				/* error management */
				)

/* computes the direct DFT of a complex signal given in (real, imaginary) representation
	and returns a (real, imaginary) representation */
/* the origin is at index [0] */
/* the highest coordinate (SignalLength-2)/2 is at index [(SignalLength-2)/2] for SignalLength even */
/* the highest coordinate (SignalLength-1)/2 is at index [(SignalLength-1)/2] for SignalLength odd */
/* the lowest coordinate -SignalLength/2 is at index [SignalLength/2] for SignalLength even */
/* the lowest coordinate -(SignalLength-1)/2 is at index [(SignalLength+1)/2] for SignalLength odd */
/* the coordinate -1 is at index [SignalLength-1] for SignalLength even or odd */
/* no restriction on SignalLength, but best efficiency for radix 2, 3, 4, 5 */
/* in the explanations above, SignalLength has to be replaced by Nx, Ny, Nz */
/* Nx is the width of the volume */
/* Ny is the height of the volume */
/* Nz is the depth of the volume */
/* in-place processing */
/* the input signal (real Re2Re and imaginary Im2Im) is
	replaced by the output signal (real Re2Re and imaginary Im2Ph) */
/* success: return(!ERROR); failure: return(ERROR); */
/* the returned value is duplicated in Status */

{ /* begin VolumeDftRealImaginaryToRealImaginary */

	double	*Re = (double *)NULL, *Im = (double *)NULL;
	double	*TmpRe = (double *)NULL, *TmpIm = (double *)NULL;
	double	*CaS = (double *)NULL;
	long	x, y, z;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(VolumeDftRealImaginaryToRealImaginary, Re2Re, *Status,
/**/	"Missing Re2Re ")
/**/DEBUG_CHECK_NULL_POINTER(VolumeDftRealImaginaryToRealImaginary, Im2Im, *Status,
/**/	"Missing Im2Im ")
/**/DEBUG_CHECK_RANGE_LONG(VolumeDftRealImaginaryToRealImaginary, Nx, 1L, LONG_MAX, *Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(VolumeDftRealImaginaryToRealImaginary, Ny, 1L, LONG_MAX, *Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(VolumeDftRealImaginaryToRealImaginary, Nz, 1L, LONG_MAX, *Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(VolumeDftRealImaginaryToRealImaginary, *Status)
/**/DEBUG_WRITE_ENTERING(VolumeDftRealImaginaryToRealImaginary,
/**/	"About to compute a 3D real/imaginary to real/imaginary direct Fourier transform")

	if (Nz > 1L) {
		*Status = bufferAllocRI(&Re, &Im, &TmpRe, &TmpIm, &CaS, Nz);
		if (*Status == ERROR) {
			bufferFreeRI(&Re, &Im, &TmpRe, &TmpIm, &CaS);
/**/		DEBUG_WRITE_LEAVING(VolumeDftRealImaginaryToRealImaginary, "Done")
			return(*Status);
		}
		*Status = GetCaS(CaS, Nz);
		if (*Status == ERROR) {
			bufferFreeRI(&Re, &Im, &TmpRe, &TmpIm, &CaS);
/**/		DEBUG_WRITE_LEAVING(VolumeDftRealImaginaryToRealImaginary, "Done")
			return(*Status);
		}
		for (y = 0L; (y < Ny); y++) {
			for (x = 0L; (x < Nx); x++) {
				*Status = GetzDoubleToDouble(Re2Re, Nx, Ny, Nz, x, y, 0L, Re, Nz);
				if (*Status == ERROR) {
					bufferFreeRI(&Re, &Im, &TmpRe, &TmpIm, &CaS);
/**/				DEBUG_WRITE_LEAVING(VolumeDftRealImaginaryToRealImaginary, "Done")
					return(*Status);
				}
				*Status = GetzDoubleToDouble(Im2Im, Nx, Ny, Nz, x, y, 0L, Im, Nz);
				if (*Status == ERROR) {
					bufferFreeRI(&Re, &Im, &TmpRe, &TmpIm, &CaS);
/**/				DEBUG_WRITE_LEAVING(VolumeDftRealImaginaryToRealImaginary, "Done")
					return(*Status);
				}
				*Status = DftRealImaginaryToRealImaginary(Re, Im, TmpRe, TmpIm, CaS, Nz);
				if (*Status == ERROR) {
					bufferFreeRI(&Re, &Im, &TmpRe, &TmpIm, &CaS);
/**/				DEBUG_WRITE_LEAVING(VolumeDftRealImaginaryToRealImaginary, "Done")
					return(*Status);
				}
				*Status = PutzDoubleToDouble(Re2Re, Nx, Ny, Nz, x, y, 0L, Re, Nz);
				if (*Status == ERROR) {
					bufferFreeRI(&Re, &Im, &TmpRe, &TmpIm, &CaS);
/**/				DEBUG_WRITE_LEAVING(VolumeDftRealImaginaryToRealImaginary, "Done")
					return(*Status);
				}
				*Status = PutzDoubleToDouble(Im2Im, Nx, Ny, Nz, x, y, 0L, Im, Nz);
				if (*Status == ERROR) {
					bufferFreeRI(&Re, &Im, &TmpRe, &TmpIm, &CaS);
/**/				DEBUG_WRITE_LEAVING(VolumeDftRealImaginaryToRealImaginary, "Done")
					return(*Status);
				}
			}
		}
		*Status = bufferFreeRI(&Re, &Im, &TmpRe, &TmpIm, &CaS);
		if (*Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(VolumeDftRealImaginaryToRealImaginary, "Done")
			return(*Status);
		}
	}
	if (Ny > 1L) {
		*Status = bufferAllocRI(&Re, &Im, &TmpRe, &TmpIm, &CaS, Ny);
		if (*Status == ERROR) {
			bufferFreeRI(&Re, &Im, &TmpRe, &TmpIm, &CaS);
/**/		DEBUG_WRITE_LEAVING(VolumeDftRealImaginaryToRealImaginary, "Done")
			return(*Status);
		}
		*Status = GetCaS(CaS, Ny);
		if (*Status == ERROR) {
			bufferFreeRI(&Re, &Im, &TmpRe, &TmpIm, &CaS);
/**/		DEBUG_WRITE_LEAVING(VolumeDftRealImaginaryToRealImaginary, "Done")
			return(*Status);
		}
		for (z = 0L; (z < Nz); z++) {
			for (x = 0L; (x < Nx); x++) {
				*Status = GetyDoubleToDouble(Re2Re, Nx, Ny, Nz, x, 0L, z, Re, Ny);
				if (*Status == ERROR) {
					bufferFreeRI(&Re, &Im, &TmpRe, &TmpIm, &CaS);
/**/				DEBUG_WRITE_LEAVING(VolumeDftRealImaginaryToRealImaginary, "Done")
					return(*Status);
				}
				*Status = GetyDoubleToDouble(Im2Im, Nx, Ny, Nz, x, 0L, z, Im, Ny);
				if (*Status == ERROR) {
					bufferFreeRI(&Re, &Im, &TmpRe, &TmpIm, &CaS);
/**/				DEBUG_WRITE_LEAVING(VolumeDftRealImaginaryToRealImaginary, "Done")
					return(*Status);
				}
				*Status = DftRealImaginaryToRealImaginary(Re, Im, TmpRe, TmpIm, CaS, Ny);
				if (*Status == ERROR) {
					bufferFreeRI(&Re, &Im, &TmpRe, &TmpIm, &CaS);
/**/				DEBUG_WRITE_LEAVING(VolumeDftRealImaginaryToRealImaginary, "Done")
					return(*Status);
				}
				*Status = PutyDoubleToDouble(Re2Re, Nx, Ny, Nz, x, 0L, z, Re, Ny);
				if (*Status == ERROR) {
					bufferFreeRI(&Re, &Im, &TmpRe, &TmpIm, &CaS);
/**/				DEBUG_WRITE_LEAVING(VolumeDftRealImaginaryToRealImaginary, "Done")
					return(*Status);
				}
				*Status = PutyDoubleToDouble(Im2Im, Nx, Ny, Nz, x, 0L, z, Im, Ny);
				if (*Status == ERROR) {
					bufferFreeRI(&Re, &Im, &TmpRe, &TmpIm, &CaS);
/**/				DEBUG_WRITE_LEAVING(VolumeDftRealImaginaryToRealImaginary, "Done")
					return(*Status);
				}
			}
		}
		*Status = bufferFreeRI(&Re, &Im, &TmpRe, &TmpIm, &CaS);
		if (*Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(VolumeDftRealImaginaryToRealImaginary, "Done")
			return(*Status);
		}
	}
	if (Nx > 1L) {
		*Status = bufferAllocRI(&Re, &Im, &TmpRe, &TmpIm, &CaS, Nx);
		if (*Status == ERROR) {
			bufferFreeRI(&Re, &Im, &TmpRe, &TmpIm, &CaS);
/**/		DEBUG_WRITE_LEAVING(VolumeDftRealImaginaryToRealImaginary, "Done")
			return(*Status);
		}
		*Status = GetCaS(CaS, Nx);
		if (*Status == ERROR) {
			bufferFreeRI(&Re, &Im, &TmpRe, &TmpIm, &CaS);
/**/		DEBUG_WRITE_LEAVING(VolumeDftRealImaginaryToRealImaginary, "Done")
			return(*Status);
		}
		for (z = 0L; (z < Nz); z++) {
			for (y = 0L; (y < Ny); y++) {
				*Status = GetxDoubleToDouble(Re2Re, Nx, Ny, Nz, 0L, y, z, Re, Nx);
				if (*Status == ERROR) {
					bufferFreeRI(&Re, &Im, &TmpRe, &TmpIm, &CaS);
/**/				DEBUG_WRITE_LEAVING(VolumeDftRealImaginaryToRealImaginary, "Done")
					return(*Status);
				}
				*Status = GetxDoubleToDouble(Im2Im, Nx, Ny, Nz, 0L, y, z, Im, Nx);
				if (*Status == ERROR) {
					bufferFreeRI(&Re, &Im, &TmpRe, &TmpIm, &CaS);
/**/				DEBUG_WRITE_LEAVING(VolumeDftRealImaginaryToRealImaginary, "Done")
					return(*Status);
				}
				*Status = DftRealImaginaryToRealImaginary(Re, Im, TmpRe, TmpIm, CaS, Nx);
				if (*Status == ERROR) {
					bufferFreeRI(&Re, &Im, &TmpRe, &TmpIm, &CaS);
/**/				DEBUG_WRITE_LEAVING(VolumeDftRealImaginaryToRealImaginary, "Done")
					return(*Status);
				}
				*Status = PutxDoubleToDouble(Re2Re, Nx, Ny, Nz, 0L, y, z, Re, Nx);
				if (*Status == ERROR) {
					bufferFreeRI(&Re, &Im, &TmpRe, &TmpIm, &CaS);
/**/				DEBUG_WRITE_LEAVING(VolumeDftRealImaginaryToRealImaginary, "Done")
					return(*Status);
				}
				*Status = PutxDoubleToDouble(Im2Im, Nx, Ny, Nz, 0L, y, z, Im, Nx);
				if (*Status == ERROR) {
					bufferFreeRI(&Re, &Im, &TmpRe, &TmpIm, &CaS);
/**/				DEBUG_WRITE_LEAVING(VolumeDftRealImaginaryToRealImaginary, "Done")
					return(*Status);
				}
			}
		}
		*Status = bufferFreeRI(&Re, &Im, &TmpRe, &TmpIm, &CaS);
		if (*Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(VolumeDftRealImaginaryToRealImaginary, "Done")
			return(*Status);
		}
	}
/**/DEBUG_WRITE_LEAVING(VolumeDftRealImaginaryToRealImaginary, "Done")
	return(*Status);
} /* end VolumeDftRealImaginaryToRealImaginary */

/*--------------------------------------------------------------------------*/
extern int		VolumeDftRealToAmplitudePhase
				(
					double	*Re2Am,				/* real -> amplitude */
					double	*PhOut,				/* output phase */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					int		*Status				/* error management */
				)

/* computes the direct DFT of a real signal
	and returns an (amplitude, phase) representation */
/* the output phase is in [rad]; its domain is (-PI, PI) */
/* the origin is at index [0] */
/* the highest coordinate (SignalLength-2)/2 is at index [(SignalLength-2)/2] for SignalLength even */
/* the highest coordinate (SignalLength-1)/2 is at index [(SignalLength-1)/2] for SignalLength odd */
/* the lowest coordinate -SignalLength/2 is at index [SignalLength/2] for SignalLength even */
/* the lowest coordinate -(SignalLength-1)/2 is at index [(SignalLength+1)/2] for SignalLength odd */
/* the coordinate -1 is at index [SignalLength-1] for SignalLength even or odd */
/* no restriction on SignalLength, but best efficiency for radix 2, 3, 4, 5 */
/* in the explanations above, SignalLength has to be replaced by Nx, Ny, Nz */
/* Nx is the width of the volume */
/* Ny is the height of the volume */
/* Nz is the depth of the volume */
/* in-place processing */
/* the input signal Re2Am is
	replaced by the output signal (amplitude Re2Am and phase PhOut) */
/* success: return(!ERROR); failure: return(ERROR); */
/* the returned value is duplicated in Status */

{ /* begin VolumeDftRealToAmplitudePhase */

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(VolumeDftRealToAmplitudePhase, Re2Am, *Status,
/**/	"Missing Re2Am ")
/**/DEBUG_CHECK_NULL_POINTER(VolumeDftRealToAmplitudePhase, PhOut, *Status,
/**/	"Missing PhOut ")
/**/DEBUG_CHECK_RANGE_LONG(VolumeDftRealToAmplitudePhase, Nx, 1L, LONG_MAX, *Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(VolumeDftRealToAmplitudePhase, Ny, 1L, LONG_MAX, *Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(VolumeDftRealToAmplitudePhase, Nz, 1L, LONG_MAX, *Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(VolumeDftRealToAmplitudePhase, *Status)
/**/DEBUG_WRITE_ENTERING(VolumeDftRealToAmplitudePhase,
/**/	"About to compute a 3D real to amplitude/phase direct Fourier transform")

	VolumeDftRealToRealImaginary(Re2Am, PhOut, Nx, Ny, Nz, Status);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(VolumeDftRealToAmplitudePhase, "Done")
		return(*Status);
	}
	*Status = VolumeRealImaginaryToAmplitudePhase(Re2Am, PhOut, Nx, Ny, Nz);
/**/DEBUG_WRITE_LEAVING(VolumeDftRealToAmplitudePhase, "Done")
	return(*Status);
} /* end VolumeDftRealToAmplitudePhase */

/*--------------------------------------------------------------------------*/
extern int		VolumeDftRealToRealImaginary
				(
					double	*Re2Re,				/* real -> real */
					double	*ImOut,				/* output imaginary */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					int		*Status				/* error management */
				)

/* computes the direct DFT of a real signal
	and returns a (real, imaginary) representation */
/* the origin is at index [0] */
/* the highest coordinate (SignalLength-2)/2 is at index [(SignalLength-2)/2] for SignalLength even */
/* the highest coordinate (SignalLength-1)/2 is at index [(SignalLength-1)/2] for SignalLength odd */
/* the lowest coordinate -SignalLength/2 is at index [SignalLength/2] for SignalLength even */
/* the lowest coordinate -(SignalLength-1)/2 is at index [(SignalLength+1)/2] for SignalLength odd */
/* the coordinate -1 is at index [SignalLength-1] for SignalLength even or odd */
/* no restriction on SignalLength, but best efficiency for radix 2, 3, 4, 5 */
/* in the explanations above, SignalLength has to be replaced by Nx, Ny, Nz */
/* Nx is the width of the volume */
/* Ny is the height of the volume */
/* Nz is the depth of the volume */
/* in-place processing */
/* the input signal Re2Re is
	replaced by the output signal (real Re2Re and imaginary ImOut) */
/* success: return(!ERROR); failure: return(ERROR); */
/* the returned value is duplicated in Status */

{ /* begin VolumeDftRealToRealImaginary */

	double	*H = (double *)NULL, *tmp = (double *)NULL;
	double	*CaS = (double *)NULL;
	double	*V = (double *)NULL;
	double	*Rppp, *Rppm, *Rpmp, *Rpmm, *Rmpp, *Rmpm, *Rmmp, *Rmmm;
	long	x, y, z;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(VolumeDftRealToRealImaginary, Re2Re, *Status,
/**/	"Missing Re2Re ")
/**/DEBUG_CHECK_NULL_POINTER(VolumeDftRealToRealImaginary, ImOut, *Status,
/**/	"Missing ImOut ")
/**/DEBUG_CHECK_RANGE_LONG(VolumeDftRealToRealImaginary, Nx, 1L, LONG_MAX, *Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(VolumeDftRealToRealImaginary, Ny, 1L, LONG_MAX, *Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(VolumeDftRealToRealImaginary, Nz, 1L, LONG_MAX, *Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(VolumeDftRealToRealImaginary, *Status)
/**/DEBUG_WRITE_ENTERING(VolumeDftRealToRealImaginary,
/**/	"About to compute a 3D real to amplitude/phase direct Fourier transform")

	AllocateVolumeDouble(&V, Nx, Ny, Nz, Status);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(VolumeDftRealToRealImaginary, "Done")
		return(*Status);
	}
	if (Nz > 1L) {
		*Status = bufferAllocR(&H, &tmp, &CaS, Nz);
		if (*Status == ERROR) {
			bufferFreeR(&H, &tmp, &CaS);
			FreeVolumeDouble(&V);
/**/		DEBUG_WRITE_LEAVING(VolumeDftRealToRealImaginary, "Done")
			return(*Status);
		}
		*Status = GetCaS(CaS, Nz);
		if (*Status == ERROR) {
			bufferFreeR(&H, &tmp, &CaS);
			FreeVolumeDouble(&V);
/**/		DEBUG_WRITE_LEAVING(VolumeDftRealToRealImaginary, "Done")
			return(*Status);
		}
		for (y = 0L; (y < Ny); y++) {
			for (x = 0L; (x < Nx); x++) {
				*Status = GetzDoubleToDouble(Re2Re, Nx, Ny, Nz, x, y, 0L, H, Nz);
				if (*Status == ERROR) {
					bufferFreeR(&H, &tmp, &CaS);
					FreeVolumeDouble(&V);
/**/				DEBUG_WRITE_LEAVING(VolumeDftRealToRealImaginary, "Done")
					return(*Status);
				}
				*Status = DiscreteHartleyTransform(H, tmp, CaS, Nz, FORWARD);
				if (*Status == ERROR) {
					bufferFreeR(&H, &tmp, &CaS);
					FreeVolumeDouble(&V);
/**/				DEBUG_WRITE_LEAVING(VolumeDftRealToRealImaginary, "Done")
					return(*Status);
				}
				*Status = PutzDoubleToDouble(Re2Re, Nx, Ny, Nz, x, y, 0L, H, Nz);
				if (*Status == ERROR) {
					bufferFreeR(&H, &tmp, &CaS);
					FreeVolumeDouble(&V);
/**/				DEBUG_WRITE_LEAVING(VolumeDftRealToRealImaginary, "Done")
					return(*Status);
				}
			}
		}
		*Status = bufferFreeR(&H, &tmp, &CaS);
		if (*Status == ERROR) {
			FreeVolumeDouble(&V);
/**/		DEBUG_WRITE_LEAVING(VolumeDftRealToRealImaginary, "Done")
			return(*Status);
		}
	}
	if (Ny > 1L) {
		*Status = bufferAllocR(&H, &tmp, &CaS, Ny);
		if (*Status == ERROR) {
			bufferFreeR(&H, &tmp, &CaS);
			FreeVolumeDouble(&V);
/**/		DEBUG_WRITE_LEAVING(VolumeDftRealToRealImaginary, "Done")
			return(*Status);
		}
		*Status = GetCaS(CaS, Ny);
		if (*Status == ERROR) {
			bufferFreeR(&H, &tmp, &CaS);
			FreeVolumeDouble(&V);
/**/		DEBUG_WRITE_LEAVING(VolumeDftRealToRealImaginary, "Done")
			return(*Status);
		}
		for (z = 0L; (z < Nz); z++) {
			for (x = 0L; (x < Nx); x++) {
				*Status = GetyDoubleToDouble(Re2Re, Nx, Ny, Nz, x, 0L, z, H, Ny);
				if (*Status == ERROR) {
					bufferFreeR(&H, &tmp, &CaS);
					FreeVolumeDouble(&V);
/**/				DEBUG_WRITE_LEAVING(VolumeDftRealToRealImaginary, "Done")
					return(*Status);
				}
				*Status = DiscreteHartleyTransform(H, tmp, CaS, Ny, FORWARD);
				if (*Status == ERROR) {
					bufferFreeR(&H, &tmp, &CaS);
					FreeVolumeDouble(&V);
/**/				DEBUG_WRITE_LEAVING(VolumeDftRealToRealImaginary, "Done")
					return(*Status);
				}
				*Status = PutyDoubleToDouble(Re2Re, Nx, Ny, Nz, x, 0L, z, H, Ny);
				if (*Status == ERROR) {
					bufferFreeR(&H, &tmp, &CaS);
					FreeVolumeDouble(&V);
/**/				DEBUG_WRITE_LEAVING(VolumeDftRealToRealImaginary, "Done")
					return(*Status);
				}
			}
		}
		*Status = bufferFreeR(&H, &tmp, &CaS);
		if (*Status == ERROR) {
			FreeVolumeDouble(&V);
/**/		DEBUG_WRITE_LEAVING(VolumeDftRealToRealImaginary, "Done")
			return(*Status);
		}
	}
	if (Nx > 1L) {
		*Status = bufferAllocR(&H, &tmp, &CaS, Nx);
		if (*Status == ERROR) {
			bufferFreeR(&H, &tmp, &CaS);
			FreeVolumeDouble(&V);
/**/		DEBUG_WRITE_LEAVING(VolumeDftRealToRealImaginary, "Done")
			return(*Status);
		}
		*Status = GetCaS(CaS, Nx);
		if (*Status == ERROR) {
			bufferFreeR(&H, &tmp, &CaS);
			FreeVolumeDouble(&V);
/**/		DEBUG_WRITE_LEAVING(VolumeDftRealToRealImaginary, "Done")
			return(*Status);
		}
		for (z = 0L; (z < Nz); z++) {
			for (y = 0L; (y < Ny); y++) {
				*Status = GetxDoubleToDouble(Re2Re, Nx, Ny, Nz, 0L, y, z, H, Nx);
				if (*Status == ERROR) {
					bufferFreeR(&H, &tmp, &CaS);
					FreeVolumeDouble(&V);
/**/				DEBUG_WRITE_LEAVING(VolumeDftRealToRealImaginary, "Done")
					return(*Status);
				}
				*Status = DiscreteHartleyTransform(H, tmp, CaS, Nx, FORWARD);
				if (*Status == ERROR) {
					bufferFreeR(&H, &tmp, &CaS);
					FreeVolumeDouble(&V);
/**/				DEBUG_WRITE_LEAVING(VolumeDftRealToRealImaginary, "Done")
					return(*Status);
				}
				*Status = PutxDoubleToDouble(V, Nx, Ny, Nz, 0L, y, z, H, Nx);
				if (*Status == ERROR) {
					bufferFreeR(&H, &tmp, &CaS);
					FreeVolumeDouble(&V);
/**/				DEBUG_WRITE_LEAVING(VolumeDftRealToRealImaginary, "Done")
					return(*Status);
				}
			}
		}
		*Status = bufferFreeR(&H, &tmp, &CaS);
		if (*Status == ERROR) {
			FreeVolumeDouble(&V);
/**/		DEBUG_WRITE_LEAVING(VolumeDftRealToRealImaginary, "Done")
			return(*Status);
		}
	}
	else {
		V = (double *)memcpy(V, Re2Re, (size_t)(Nx * Ny * Nz * (long)sizeof(double)));
	}
	*ImOut++ = 0.0F;
	*Re2Re++ = *V;
	Rppp = V + (ptrdiff_t)1;
	Rmpp = V + (ptrdiff_t)(Nx - 1L);
	x = Nx - 1L;
	while (x-- > 0L) {
		*Re2Re++ = (double)(((double)*Rppp + (double)*Rmpp) * 0.5);
		*ImOut++ = (double)((-(double)*Rppp++ + (double)*Rmpp--) * 0.5);
	}
	Rpmp = Rmmp = V + (ptrdiff_t)(Nx * (Ny - 1L));
	Rmpp = Rppp;
	y = Ny - 1L;
	while (y-- > 0L) {
		*Re2Re++ = (double)(((double)*Rpmp++ + (double)*Rmpp) * 0.5);
		*ImOut++ = (double)((-(double)*Rppp++ + (double)*Rmmp) * 0.5);
		Rmpp += (ptrdiff_t)(Nx - 1L);
		Rmmp += (ptrdiff_t)(Nx - 1L);
		x = Nx - 1L;
		while (x-- > 0L) {
			*Re2Re++ = (double)(((double)*Rpmp++ + (double)*Rmpp--) * 0.5);
			*ImOut++ = (double)((-(double)*Rppp++ + (double)*Rmmp--) * 0.5);
		}
		Rmmp -= (ptrdiff_t)Nx;
		Rmpp = Rppp;
		Rpmp = Rmmp;
	}
	Rppm = Rpmm = Rmpm = Rmmm = V + (ptrdiff_t)(Nx * Ny * (Nz - 1L));
	Rpmp = Rmmp = Rppp;
	z = Nz - 1L;
	while (z-- > 0L) {
		*Re2Re++ = (double)((-(double)*Rppp + (double)*Rppm
			+ (double)*Rpmp + (double)*Rpmm
			+ (double)*Rmpp + (double)*Rmpm
			+ (double)*Rmmp - (double)*Rmmm) * 0.25);
		*ImOut++ = (double)((-(double)*Rppp++ - (double)*Rppm++
			- (double)*Rpmp++ + (double)*Rpmm++
			- (double)*Rmpp + (double)*Rmpm
			+ (double)*Rmmp + (double)*Rmmm) * 0.25);
		Rmpp += (ptrdiff_t)(Nx - 1L);
		Rmpm += (ptrdiff_t)(Nx - 1L);
		Rmmp += (ptrdiff_t)(Nx - 1L);
		Rmmm += (ptrdiff_t)(Nx - 1L);
		x = Nx - 1L;
		while (x-- > 0L) {
			*Re2Re++ = (double)((-(double)*Rppp + (double)*Rppm
				+ (double)*Rpmp + (double)*Rpmm
				+ (double)*Rmpp + (double)*Rmpm
				+ (double)*Rmmp - (double)*Rmmm) * 0.25);
			*ImOut++ = (double)((-(double)*Rppp++ - (double)*Rppm++
				- (double)*Rpmp++ + (double)*Rpmm++
				- (double)*Rmpp-- + (double)*Rmpm--
				+ (double)*Rmmp-- + (double)*Rmmm--) * 0.25);
		}
		Rpmp += (ptrdiff_t)(Nx * (Ny - 2L));
		Rpmm += (ptrdiff_t)(Nx * (Ny - 2L));
		Rmpp = Rppp;
		Rmpm = Rppm;
		Rmmp = Rpmp;
		Rmmm = Rpmm;
		y = Ny - 1L;
		while (y-- > 0L) {
			*Re2Re++ = (double)((-(double)*Rppp + (double)*Rppm
				+ (double)*Rpmp + (double)*Rpmm
				+ (double)*Rmpp + (double)*Rmpm
				+ (double)*Rmmp - (double)*Rmmm) * 0.25);
			*ImOut++ = (double)((-(double)*Rppp++ - (double)*Rppm++
				- (double)*Rpmp++ + (double)*Rpmm++
				- (double)*Rmpp + (double)*Rmpm
				+ (double)*Rmmp + (double)*Rmmm) * 0.25);
			Rmpp += (ptrdiff_t)(Nx - 1L);
			Rmpm += (ptrdiff_t)(Nx - 1L);
			Rmmp += (ptrdiff_t)(Nx - 1L);
			Rmmm += (ptrdiff_t)(Nx - 1L);
			x = Nx - 1L;
			while (x-- > 0L) {
				*Re2Re++ = (double)((-(double)*Rppp + (double)*Rppm
					+ (double)*Rpmp + (double)*Rpmm
					+ (double)*Rmpp + (double)*Rmpm
					+ (double)*Rmmp - (double)*Rmmm) * 0.25);
				*ImOut++ = (double)((-(double)*Rppp++ - (double)*Rppm++
					- (double)*Rpmp++ + (double)*Rpmm++
					- (double)*Rmpp-- + (double)*Rmpm--
					+ (double)*Rmmp-- + (double)*Rmmm--) * 0.25);
			}
			Rmmp -= (ptrdiff_t)Nx;
			Rmmm -= (ptrdiff_t)Nx;
			Rmpp = Rppp;
			Rmpm = Rppm;
			Rpmp = Rmmp;
			Rpmm = Rmmm;
		}
		Rmmm -= (ptrdiff_t)(Nx * Ny);
		Rpmp = Rmmp = Rppp;
		Rppm = Rpmm = Rmpm = Rmmm;
	}
	*Status = FreeVolumeDouble(&V);
/**/DEBUG_WRITE_LEAVING(VolumeDftRealToRealImaginary, "Done")
	return(*Status);
} /* end VolumeDftRealToRealImaginary */

/*--------------------------------------------------------------------------*/
extern int		VolumeInvDftAmplitudePhaseToAmplitudePhase
				(
					double	*Am2Am,				/* amplitude -> amplitude */
					double	*Ph2Ph,				/* phase -> phase */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					int		*Status				/* error management */
				)

/* computes the inverse DFT of a complex signal given in (amplitude, phase) representation
	and returns an (amplitude, phase) representation */
/* the input phase is in [rad] */
/* the output phase is in [rad]; its domain is (-PI, PI) */
/* the origin is at index [0] */
/* the highest coordinate (SignalLength-2)/2 is at index [(SignalLength-2)/2] for SignalLength even */
/* the highest coordinate (SignalLength-1)/2 is at index [(SignalLength-1)/2] for SignalLength odd */
/* the lowest coordinate -SignalLength/2 is at index [SignalLength/2] for SignalLength even */
/* the lowest coordinate -(SignalLength-1)/2 is at index [(SignalLength+1)/2] for SignalLength odd */
/* the coordinate -1 is at index [SignalLength-1] for SignalLength even or odd */
/* no restriction on SignalLength, but best efficiency for radix 2, 3, 4, 5 */
/* in the explanations above, SignalLength has to be replaced by Nx, Ny, Nz */
/* Nx is the width of the volume */
/* Ny is the height of the volume */
/* Nz is the depth of the volume */
/* in-place processing */
/* the input signal (amplitude Am2Am and phase Ph2Ph) is
	replaced by the output signal (amplitude Am2Am and phase Ph2Ph) */
/* success: return(!ERROR); failure: return(ERROR); */
/* the returned value is duplicated in Status */

{ /* begin VolumeInvDftAmplitudePhaseToAmplitudePhase */

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(VolumeInvDftAmplitudePhaseToAmplitudePhase, Am2Am, *Status,
/**/	"Missing Am2Am ")
/**/DEBUG_CHECK_NULL_POINTER(VolumeInvDftAmplitudePhaseToAmplitudePhase, Ph2Ph, *Status,
/**/	"Missing Ph2Ph ")
/**/DEBUG_CHECK_RANGE_LONG(VolumeInvDftAmplitudePhaseToAmplitudePhase, Nx, 1L, LONG_MAX, *Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(VolumeInvDftAmplitudePhaseToAmplitudePhase, Ny, 1L, LONG_MAX, *Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(VolumeInvDftAmplitudePhaseToAmplitudePhase, Nz, 1L, LONG_MAX, *Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(VolumeInvDftAmplitudePhaseToAmplitudePhase, *Status)
/**/DEBUG_WRITE_ENTERING(VolumeInvDftAmplitudePhaseToAmplitudePhase,
/**/	"About to compute a 3D amplitude/phase to amplitude/phase inverse Fourier transform")

	*Status = VolumeAmplitudePhaseToRealImaginary(Am2Am, Ph2Ph, Nx, Ny, Nz);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(VolumeInvDftAmplitudePhaseToAmplitudePhase, "Done")
		return(*Status);
	}
	VolumeInvDftRealImaginaryToAmplitudePhase(Am2Am, Ph2Ph, Nx, Ny, Nz, Status);
/**/DEBUG_WRITE_LEAVING(VolumeInvDftAmplitudePhaseToAmplitudePhase, "Done")
	return(*Status);
} /* end VolumeInvDftAmplitudePhaseToAmplitudePhase */

/*--------------------------------------------------------------------------*/
extern int		VolumeInvDftAmplitudePhaseToReal
				(
					double	*Am2Re,				/* amplitude -> real */
					double	*PhIn,				/* input phase */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					int		*Status				/* error management */
				)

/* computes the inverse DFT of a complex signal given in (amplitude, phase) representation
	and returns a real signal */
/* the complex Fourier signal is symmetrized before the inverse transformation is applied */
/* the input phase is in [rad] */
/* the origin is at index [0] */
/* the highest coordinate (SignalLength-2)/2 is at index [(SignalLength-2)/2] for SignalLength even */
/* the highest coordinate (SignalLength-1)/2 is at index [(SignalLength-1)/2] for SignalLength odd */
/* the lowest coordinate -SignalLength/2 is at index [SignalLength/2] for SignalLength even */
/* the lowest coordinate -(SignalLength-1)/2 is at index [(SignalLength+1)/2] for SignalLength odd */
/* the coordinate -1 is at index [SignalLength-1] for SignalLength even or odd */
/* no restriction on SignalLength, but best efficiency for radix 2, 3, 4, 5 */
/* in the explanations above, SignalLength has to be replaced by Nx, Ny, Nz */
/* Nx is the width of the volume */
/* Ny is the height of the volume */
/* Nz is the depth of the volume */
/* in-place processing */
/* the input signal (amplitude Am2Re and phase PhIn) is
	replaced by the output signal (real Am2Re) */
/* PhIn is destroyed */
/* success: return(!ERROR); failure: return(ERROR); */
/* the returned value is duplicated in Status */

{ /* begin VolumeInvDftAmplitudePhaseToReal */

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(VolumeInvDftAmplitudePhaseToReal, Am2Re, *Status,
/**/	"Missing Am2Re ")
/**/DEBUG_CHECK_NULL_POINTER(VolumeInvDftAmplitudePhaseToReal, PhIn, *Status,
/**/	"Missing PhIn ")
/**/DEBUG_CHECK_RANGE_LONG(VolumeInvDftAmplitudePhaseToReal, Nx, 1L, LONG_MAX, *Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(VolumeInvDftAmplitudePhaseToReal, Ny, 1L, LONG_MAX, *Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(VolumeInvDftAmplitudePhaseToReal, Nz, 1L, LONG_MAX, *Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(VolumeInvDftAmplitudePhaseToReal, *Status)
/**/DEBUG_WRITE_ENTERING(VolumeInvDftAmplitudePhaseToReal,
/**/	"About to compute a 3D amplitude/phase to real inverse Fourier transform")

	*Status = VolumeAmplitudePhaseToRealImaginary(Am2Re, PhIn, Nx, Ny, Nz);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(VolumeInvDftAmplitudePhaseToReal, "Done")
		return(*Status);
	}
	VolumeInvDftRealImaginaryToReal(Am2Re, PhIn, Nx, Ny, Nz, Status);
/**/DEBUG_WRITE_LEAVING(VolumeInvDftAmplitudePhaseToReal, "Done")
	return(*Status);
} /* end VolumeInvDftAmplitudePhaseToReal */

/*--------------------------------------------------------------------------*/
extern int		VolumeInvDftAmplitudePhaseToRealImaginary
				(
					double	*Am2Re,				/* amplitude -> real */
					double	*Ph2Im,				/* phase -> imaginary */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					int		*Status				/* error management */
				)

/* computes the inverse DFT of a complex signal given in (amplitude, phase) representation
	and returns a (real, imaginary) representation */
/* the input phase is in [rad] */
/* the origin is at index [0] */
/* the highest coordinate (SignalLength-2)/2 is at index [(SignalLength-2)/2] for SignalLength even */
/* the highest coordinate (SignalLength-1)/2 is at index [(SignalLength-1)/2] for SignalLength odd */
/* the lowest coordinate -SignalLength/2 is at index [SignalLength/2] for SignalLength even */
/* the lowest coordinate -(SignalLength-1)/2 is at index [(SignalLength+1)/2] for SignalLength odd */
/* the coordinate -1 is at index [SignalLength-1] for SignalLength even or odd */
/* no restriction on SignalLength, but best efficiency for radix 2, 3, 4, 5 */
/* in the explanations above, SignalLength has to be replaced by Nx, Ny, Nz */
/* Nx is the width of the volume */
/* Ny is the height of the volume */
/* Nz is the depth of the volume */
/* in-place processing */
/* the input signal (amplitude Am2Re and phase Ph2Im) is
	replaced by the output signal (real Am2Re and imaginary Ph2Im) */
/* success: return(!ERROR); failure: return(ERROR); */
/* the returned value is duplicated in Status */

{ /* begin VolumeInvDftAmplitudePhaseToRealImaginary */

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(VolumeInvDftAmplitudePhaseToRealImaginary, Am2Re, *Status,
/**/	"Missing Am2Re ")
/**/DEBUG_CHECK_NULL_POINTER(VolumeInvDftAmplitudePhaseToRealImaginary, Ph2Im, *Status,
/**/	"Missing Ph2Im ")
/**/DEBUG_CHECK_RANGE_LONG(VolumeInvDftAmplitudePhaseToRealImaginary, Nx, 1L, LONG_MAX, *Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(VolumeInvDftAmplitudePhaseToRealImaginary, Ny, 1L, LONG_MAX, *Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(VolumeInvDftAmplitudePhaseToRealImaginary, Nz, 1L, LONG_MAX, *Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(VolumeInvDftAmplitudePhaseToRealImaginary, *Status)
/**/DEBUG_WRITE_ENTERING(VolumeInvDftAmplitudePhaseToRealImaginary,
/**/	"About to compute a 3D amplitude/phase to real/imaginary inverse Fourier transform")

	*Status = VolumeAmplitudePhaseToRealImaginary(Am2Re, Ph2Im, Nx, Ny, Nz);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(VolumeInvDftAmplitudePhaseToRealImaginary, "Done")
		return(*Status);
	}
	VolumeInvDftRealImaginaryToRealImaginary(Am2Re, Ph2Im, Nx, Ny, Nz, Status);
/**/DEBUG_WRITE_LEAVING(VolumeInvDftAmplitudePhaseToRealImaginary, "Done")
	return(*Status);
} /* end VolumeInvDftAmplitudePhaseToRealImaginary */

/*--------------------------------------------------------------------------*/
extern int		VolumeInvDftRealImaginaryToAmplitudePhase
				(
					double	*Re2Am,				/* real -> amplitude */
					double	*Im2Ph,				/* imaginary -> phase */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					int		*Status				/* error management */
				)

/* computes the inverse DFT of a complex signal given in (real, imaginary) representation
	and returns an (amplitude, phase) representation */
/* the output phase is in [rad]; its domain is (-PI, PI) */
/* the origin is at index [0] */
/* the highest coordinate (SignalLength-2)/2 is at index [(SignalLength-2)/2] for SignalLength even */
/* the highest coordinate (SignalLength-1)/2 is at index [(SignalLength-1)/2] for SignalLength odd */
/* the lowest coordinate -SignalLength/2 is at index [SignalLength/2] for SignalLength even */
/* the lowest coordinate -(SignalLength-1)/2 is at index [(SignalLength+1)/2] for SignalLength odd */
/* the coordinate -1 is at index [SignalLength-1] for SignalLength even or odd */
/* no restriction on SignalLength, but best efficiency for radix 2, 3, 4, 5 */
/* in the explanations above, SignalLength has to be replaced by Nx, Ny, Nz */
/* Nx is the width of the volume */
/* Ny is the height of the volume */
/* Nz is the depth of the volume */
/* in-place processing */
/* the input signal (real Re2Am and imaginary Im2Ph) is
	replaced by the output signal (amplitude Re2Am and phase Im2Ph) */
/* success: return(!ERROR); failure: return(ERROR); */
/* the returned value is duplicated in Status */

{ /* begin VolumeInvDftRealImaginaryToAmplitudePhase */

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(VolumeInvDftRealImaginaryToAmplitudePhase, Re2Am, *Status,
/**/	"Missing Re2Am ")
/**/DEBUG_CHECK_NULL_POINTER(VolumeInvDftRealImaginaryToAmplitudePhase, Im2Ph, *Status,
/**/	"Missing Im2Ph ")
/**/DEBUG_CHECK_RANGE_LONG(VolumeInvDftRealImaginaryToAmplitudePhase, Nx, 1L, LONG_MAX, *Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(VolumeInvDftRealImaginaryToAmplitudePhase, Ny, 1L, LONG_MAX, *Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(VolumeInvDftRealImaginaryToAmplitudePhase, Nz, 1L, LONG_MAX, *Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(VolumeInvDftRealImaginaryToAmplitudePhase, *Status)
/**/DEBUG_WRITE_ENTERING(VolumeInvDftRealImaginaryToAmplitudePhase,
/**/	"About to compute a 3D real/imaginary to amplitude/phase inverse Fourier transform")

	VolumeInvDftRealImaginaryToRealImaginary(Re2Am, Im2Ph, Nx, Ny, Nz, Status);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(VolumeInvDftRealImaginaryToAmplitudePhase, "Done")
		return(*Status);
	}
	*Status = VolumeRealImaginaryToAmplitudePhase(Re2Am, Im2Ph, Nx, Ny, Nz);
/**/DEBUG_WRITE_LEAVING(VolumeInvDftRealImaginaryToAmplitudePhase, "Done")
	return(*Status);
} /* end VolumeInvDftRealImaginaryToAmplitudePhase */

/*--------------------------------------------------------------------------*/
extern int		VolumeInvDftRealImaginaryToReal
				(
					double	*Re2Re,				/* real -> real */
					double	*ImIn,				/* input imaginary */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					int		*Status				/* error management */
				)

/* computes the inverse DFT of a complex signal given in (real, imaginary) representation
	and returns a real signal */
/* the complex Fourier signal is symmetrized before the inverse transformation is applied */
/* the origin is at index [0] */
/* the highest coordinate (SignalLength-2)/2 is at index [(SignalLength-2)/2] for SignalLength even */
/* the highest coordinate (SignalLength-1)/2 is at index [(SignalLength-1)/2] for SignalLength odd */
/* the lowest coordinate -SignalLength/2 is at index [SignalLength/2] for SignalLength even */
/* the lowest coordinate -(SignalLength-1)/2 is at index [(SignalLength+1)/2] for SignalLength odd */
/* the coordinate -1 is at index [SignalLength-1] for SignalLength even or odd */
/* no restriction on SignalLength, but best efficiency for radix 2, 3, 4, 5 */
/* in the explanations above, SignalLength has to be replaced by Nx, Ny, Nz */
/* Nx is the width of the volume */
/* Ny is the height of the volume */
/* Nz is the depth of the volume */
/* in-place processing */
/* the input signal (real Re2Re and imaginary ImIn) is
	replaced by the output signal (real Re2Re) */
/* success: return(!ERROR); failure: return(ERROR); */
/* the returned value is duplicated in Status */

{ /* begin VolumeInvDftRealImaginaryToReal */

	double	*H = (double *)NULL, *tmp = (double *)NULL;
	double	*CaS = (double *)NULL;
	double	*V = (double *)NULL;
	double	*Rppp, *Rppm, *Rpmp, *Rpmm, *Rmpp, *Rmpm, *Rmmp, *Rmmm;
	double	*Ippp, *Ippm, *Ipmp, *Ipmm, *Impp, *Impm, *Immp, *Immm;
	double	re, im;
	long	x, y, z;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(VolumeInvDftRealImaginaryToReal, Re2Re, *Status,
/**/	"Missing Re2Re ")
/**/DEBUG_CHECK_NULL_POINTER(VolumeInvDftRealImaginaryToReal, ImIn, *Status,
/**/	"Missing ImIn ")
/**/DEBUG_CHECK_RANGE_LONG(VolumeInvDftRealImaginaryToReal, Nx, 1L, LONG_MAX, *Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(VolumeInvDftRealImaginaryToReal, Ny, 1L, LONG_MAX, *Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(VolumeInvDftRealImaginaryToReal, Nz, 1L, LONG_MAX, *Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(VolumeInvDftRealImaginaryToReal, *Status)
/**/DEBUG_WRITE_ENTERING(VolumeInvDftRealImaginaryToReal,
/**/	"About to compute a 3D real/imaginary to real inverse Fourier transform")

	AllocateVolumeDouble(&V, Nx, Ny, Nz, Status);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(VolumeInvDftRealImaginaryToReal, "Done")
		return(*Status);
	}
	*V = (double)sqrt((double)*Re2Re * (double)*Re2Re + (double)*ImIn * (double)*ImIn);
	*V = (*Re2Re < 0.0) ? (-*V) : (*V);
	V++;
	Rppp = Re2Re + (ptrdiff_t)1;
	Rmpp = Re2Re + (ptrdiff_t)(Nx - 1L);
	Ippp = ImIn + (ptrdiff_t)1;
	Impp = ImIn + (ptrdiff_t)(Nx - 1L);
	x = Nx - 1L;
	while (x-- > 0L) {
		re = (double)*Rmpp-- + (double)*Rppp++;
		im = (double)*Impp-- - (double)*Ippp++;
		*V++ = (double)((re + im) * 0.5);
	}
	Rpmp = Rmmp = Re2Re + (ptrdiff_t)(Nx * (Ny - 1L));
	Rmpp = Rppp;
	Immp = ImIn + (ptrdiff_t)(Nx * (Ny - 1L));
	y = Ny - 1L;
	while (y-- > 0L) {
		re = (double)*Rpmp++ + (double)*Rmpp;
		im = (double)*Immp - (double)*Ippp++;
		*V++ = (double)((re + im) * 0.5);
		Rmpp += (ptrdiff_t)(Nx - 1L);
		Immp += (ptrdiff_t)(Nx - 1L);
		x = Nx - 1L;
		while (x-- > 0L) {
			re = (double)*Rpmp++ + (double)*Rmpp--;
			im = (double)*Immp-- - (double)*Ippp++;
			*V++ = (double)((re + im) * 0.5);
		}
		Rmmp -= (ptrdiff_t)Nx;
		Rmpp += (ptrdiff_t)Nx;
		Rpmp = Rmmp;
		Immp -= (ptrdiff_t)Nx;
	}
	Rppm = Rpmm = Rmpm = Rmmm = Re2Re + (ptrdiff_t)(Nx * Ny * (Nz - 1L));
	Rppp += (ptrdiff_t)(Nx * (Ny - 1L));
	Rpmp = Rmmp = Rppp;
	Ippm = Ipmm = Impm = Immm = ImIn + (ptrdiff_t)(Nx * Ny * (Nz - 1L));
	Ipmp = Impp = Immp = Ippp;
	z = Nz - 1L;
	while (z-- > 0L) {
		re = -(double)*Rppp++ + (double)*Rppm++ + (double)*Rpmp++ + (double)*Rpmm++
			+ (double)*Rmpp + (double)*Rmpm + (double)*Rmmp - (double)*Rmmm;
		im = -(double)*Ippp++ - (double)*Ippm++ - (double)*Ipmp++ + (double)*Ipmm++
			- (double)*Impp + (double)*Impm + (double)*Immp + (double)*Immm;
		*V++ = (double)((re + im) * 0.25);
		Rmpp += (ptrdiff_t)(Nx - 1L);
		Rmpm += (ptrdiff_t)(Nx - 1L);
		Rmmp += (ptrdiff_t)(Nx - 1L);
		Rmmm += (ptrdiff_t)(Nx - 1L);
		Impp += (ptrdiff_t)(Nx - 1L);
		Impm += (ptrdiff_t)(Nx - 1L);
		Immp += (ptrdiff_t)(Nx - 1L);
		Immm += (ptrdiff_t)(Nx - 1L);
		x = Nx - 1L;
		while (x-- > 0L) {
			re = -(double)*Rppp++ + (double)*Rppm++ + (double)*Rpmp++ + (double)*Rpmm++
				+ (double)*Rmpp-- + (double)*Rmpm-- + (double)*Rmmp-- - (double)*Rmmm--;
			im = -(double)*Ippp++ - (double)*Ippm++ - (double)*Ipmp++ + (double)*Ipmm++
				- (double)*Impp-- + (double)*Impm-- + (double)*Immp-- + (double)*Immm--;
			*V++ = (double)((re + im) * 0.25);
		}
		Rpmp += (ptrdiff_t)(Nx * (Ny - 2L));
		Rpmm += (ptrdiff_t)(Nx * (Ny - 2L));
		Rmpp = Rppp;
		Rmpm = Rppm;
		Rmmp = Rpmp;
		Rmmm = Rpmm;
		Ipmp += (ptrdiff_t)(Nx * (Ny - 2L));
		Ipmm += (ptrdiff_t)(Nx * (Ny - 2L));
		Impp = Ippp;
		Impm = Ippm;
		Immp = Ipmp;
		Immm = Ipmm;
		y = Ny - 1L;
		while (y-- > 0L) {
			re = -(double)*Rppp++ + (double)*Rppm++ + (double)*Rpmp++ + (double)*Rpmm++
				+ (double)*Rmpp + (double)*Rmpm + (double)*Rmmp - (double)*Rmmm;
			im = -(double)*Ippp++ - (double)*Ippm++ - (double)*Ipmp++ + (double)*Ipmm++
				- (double)*Impp + (double)*Impm + (double)*Immp + (double)*Immm;
			*V++ = (double)((re + im) * 0.25);
			Rmpp += (ptrdiff_t)(Nx - 1L);
			Rmpm += (ptrdiff_t)(Nx - 1L);
			Rmmp += (ptrdiff_t)(Nx - 1L);
			Rmmm += (ptrdiff_t)(Nx - 1L);
			Impp += (ptrdiff_t)(Nx - 1L);
			Impm += (ptrdiff_t)(Nx - 1L);
			Immp += (ptrdiff_t)(Nx - 1L);
			Immm += (ptrdiff_t)(Nx - 1L);
			x = Nx - 1L;
			while (x-- > 0L) {
				re = -(double)*Rppp++ + (double)*Rppm++ + (double)*Rpmp++ + (double)*Rpmm++
					+ (double)*Rmpp-- + (double)*Rmpm-- + (double)*Rmmp-- - (double)*Rmmm--;
				im = -(double)*Ippp++ - (double)*Ippm++ - (double)*Ipmp++ + (double)*Ipmm++
					- (double)*Impp-- + (double)*Impm-- + (double)*Immp-- + (double)*Immm--;
				*V++ = (double)((re + im) * 0.25);
			}
			Rmmp -= (ptrdiff_t)Nx;
			Rmmm -= (ptrdiff_t)Nx;
			Rmpp = Rppp;
			Rmpm = Rppm;
			Rpmp = Rmmp;
			Rpmm = Rmmm;
			Immp -= (ptrdiff_t)Nx;
			Immm -= (ptrdiff_t)Nx;
			Impp = Ippp;
			Impm = Ippm;
			Ipmp = Immp;
			Ipmm = Immm;
			}
		Rmmm -= (ptrdiff_t)(Nx * Ny);
		Rpmp = Rmmp = Rppp;
		Rppm = Rpmm = Rmpm = Rmmm;
		Immm -= (ptrdiff_t)(Nx * Ny);
		Ipmp = Immp = Ippp;
		Ippm = Ipmm = Impm = Immm;
	}
	V -= (ptrdiff_t)(Nx * Ny * Nz);
	if (Nz > 1L) {
		*Status = bufferAllocR(&H, &tmp, &CaS, Nz);
		if (*Status == ERROR) {
			bufferFreeR(&H, &tmp, &CaS);
			FreeVolumeDouble(&V);
/**/		DEBUG_WRITE_LEAVING(VolumeInvDftRealImaginaryToReal, "Done")
			return(*Status);
		}
		*Status = GetCaS(CaS, Nz);
		if (*Status == ERROR) {
			bufferFreeR(&H, &tmp, &CaS);
			FreeVolumeDouble(&V);
/**/		DEBUG_WRITE_LEAVING(VolumeInvDftRealImaginaryToReal, "Done")
			return(*Status);
		}
		for (y = 0L; (y < Ny); y++) {
			for (x = 0L; (x < Nx); x++) {
				*Status = GetzDoubleToDouble(V, Nx, Ny, Nz, x, y, 0L, H, Nz);
				if (*Status == ERROR) {
					bufferFreeR(&H, &tmp, &CaS);
					FreeVolumeDouble(&V);
/**/				DEBUG_WRITE_LEAVING(VolumeInvDftRealImaginaryToReal, "Done")
					return(*Status);
				}
				*Status = DiscreteHartleyTransform(H, tmp, CaS, Nz, BACKWARD);
				if (*Status == ERROR) {
					bufferFreeR(&H, &tmp, &CaS);
					FreeVolumeDouble(&V);
/**/				DEBUG_WRITE_LEAVING(VolumeInvDftRealImaginaryToReal, "Done")
					return(*Status);
				}
				*Status = PutzDoubleToDouble(V, Nx, Ny, Nz, x, y, 0L, H, Nz);
				if (*Status == ERROR) {
					bufferFreeR(&H, &tmp, &CaS);
					FreeVolumeDouble(&V);
/**/				DEBUG_WRITE_LEAVING(VolumeInvDftRealImaginaryToReal, "Done")
					return(*Status);
				}
			}
		}
		*Status = bufferFreeR(&H, &tmp, &CaS);
		if (*Status == ERROR) {
			FreeVolumeDouble(&V);
/**/		DEBUG_WRITE_LEAVING(VolumeInvDftRealImaginaryToReal, "Done")
			return(*Status);
		}
	}
	if (Ny > 1L) {
		*Status = bufferAllocR(&H, &tmp, &CaS, Ny);
		if (*Status == ERROR) {
			bufferFreeR(&H, &tmp, &CaS);
			FreeVolumeDouble(&V);
/**/		DEBUG_WRITE_LEAVING(VolumeInvDftRealImaginaryToReal, "Done")
			return(*Status);
		}
		*Status = GetCaS(CaS, Ny);
		if (*Status == ERROR) {
			bufferFreeR(&H, &tmp, &CaS);
			FreeVolumeDouble(&V);
/**/		DEBUG_WRITE_LEAVING(VolumeInvDftRealImaginaryToReal, "Done")
			return(*Status);
		}
		for (z = 0L; (z < Nz); z++) {
			for (x = 0L; (x < Nx); x++) {
				*Status = GetyDoubleToDouble(V, Nx, Ny, Nz, x, 0L, z, H, Ny);
				if (*Status == ERROR) {
					bufferFreeR(&H, &tmp, &CaS);
					FreeVolumeDouble(&V);
/**/				DEBUG_WRITE_LEAVING(VolumeInvDftRealImaginaryToReal, "Done")
					return(*Status);
				}
				*Status = DiscreteHartleyTransform(H, tmp, CaS, Ny, BACKWARD);
				if (*Status == ERROR) {
					bufferFreeR(&H, &tmp, &CaS);
					FreeVolumeDouble(&V);
/**/				DEBUG_WRITE_LEAVING(VolumeInvDftRealImaginaryToReal, "Done")
					return(*Status);
				}
				*Status = PutyDoubleToDouble(V, Nx, Ny, Nz, x, 0L, z, H, Ny);
				if (*Status == ERROR) {
					bufferFreeR(&H, &tmp, &CaS);
					FreeVolumeDouble(&V);
/**/				DEBUG_WRITE_LEAVING(VolumeInvDftRealImaginaryToReal, "Done")
					return(*Status);
				}
			}
		}
	        *Status = bufferFreeR(&H, &tmp, &CaS);
		if (*Status == ERROR) {
			FreeVolumeDouble(&V);
/**/		DEBUG_WRITE_LEAVING(VolumeInvDftRealImaginaryToReal, "Done")
			return(*Status);
		}
	}
	if (Nx > 1L) {
		*Status = bufferAllocR(&H, &tmp, &CaS, Nx);
		if (*Status == ERROR) {
			bufferFreeR(&H, &tmp, &CaS);
			FreeVolumeDouble(&V);
/**/		DEBUG_WRITE_LEAVING(VolumeInvDftRealImaginaryToReal, "Done")
			return(*Status);
		}
		*Status = GetCaS(CaS, Nx);
		if (*Status == ERROR) {
			bufferFreeR(&H, &tmp, &CaS);
			FreeVolumeDouble(&V);
/**/		DEBUG_WRITE_LEAVING(VolumeInvDftRealImaginaryToReal, "Done")
			return(*Status);
		}
		for (z = 0L; (z < Nz); z++) {
			for (y = 0L; (y < Ny); y++) {
				*Status = GetxDoubleToDouble(V, Nx, Ny, Nz, 0L, y, z, H, Nx);
				if (*Status == ERROR) {
					bufferFreeR(&H, &tmp, &CaS);
					FreeVolumeDouble(&V);
/**/				DEBUG_WRITE_LEAVING(VolumeInvDftRealImaginaryToReal, "Done")
					return(*Status);
				}
				*Status = DiscreteHartleyTransform(H, tmp, CaS, Nx, BACKWARD);
				if (*Status == ERROR) {
					bufferFreeR(&H, &tmp, &CaS);
					FreeVolumeDouble(&V);
/**/				DEBUG_WRITE_LEAVING(VolumeInvDftRealImaginaryToReal, "Done")
					return(*Status);
				}
				*Status = PutxDoubleToDouble(Re2Re, Nx, Ny, Nz, 0L, y, z, H, Nx);
				if (*Status == ERROR) {
					bufferFreeR(&H, &tmp, &CaS);
					FreeVolumeDouble(&V);
/**/				DEBUG_WRITE_LEAVING(VolumeInvDftRealImaginaryToReal, "Done")
					return(*Status);
				}
			}
		}
		*Status = bufferFreeR(&H, &tmp, &CaS);
		if (*Status == ERROR) {
			FreeVolumeDouble(&V);
/**/		DEBUG_WRITE_LEAVING(VolumeInvDftRealImaginaryToReal, "Done")
			return(*Status);
		}
	}
	else {
		Re2Re = (double *)memcpy(Re2Re, V, (size_t)(Nx * Ny * Nz * (long)sizeof(double)));
	}
	*Status = FreeVolumeDouble(&V);
/**/DEBUG_WRITE_LEAVING(VolumeInvDftRealImaginaryToReal, "Done")
	return(*Status);
} /* end VolumeInvDftRealImaginaryToReal */

/*--------------------------------------------------------------------------*/
extern int		VolumeInvDftRealImaginaryToRealImaginary
				(
					double	*Re2Re,				/* real -> real */
					double	*Im2Im,				/* imaginary -> imaginary */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					int		*Status				/* error management */
				)

/* computes the inverse DFT of a complex signal given in (real, imaginary) representation
	and returns a (real, imaginary) representation */
/* the origin is at index [0] */
/* the highest coordinate (SignalLength-2)/2 is at index [(SignalLength-2)/2] for SignalLength even */
/* the highest coordinate (SignalLength-1)/2 is at index [(SignalLength-1)/2] for SignalLength odd */
/* the lowest coordinate -SignalLength/2 is at index [SignalLength/2] for SignalLength even */
/* the lowest coordinate -(SignalLength-1)/2 is at index [(SignalLength+1)/2] for SignalLength odd */
/* the coordinate -1 is at index [SignalLength-1] for SignalLength even or odd */
/* no restriction on SignalLength, but best efficiency for radix 2, 3, 4, 5 */
/* in the explanations above, SignalLength has to be replaced by Nx, Ny, Nz */
/* Nx is the width of the volume */
/* Ny is the height of the volume */
/* Nz is the depth of the volume */
/* in-place processing */
/* the input signal (real Re2Re and imaginary Im2Im) is
	replaced by the output signal (real Re2Re and imaginary Im2Ph) */
/* success: return(!ERROR); failure: return(ERROR); */
/* the returned value is duplicated in Status */

{ /* begin VolumeInvDftRealImaginaryToRealImaginary */

	double	*Re = (double *)NULL, *Im = (double *)NULL;
	double	*TmpRe = (double *)NULL, *TmpIm = (double *)NULL;
	double	*CaS = (double *)NULL;
	long	x, y, z;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(VolumeInvDftRealImaginaryToRealImaginary, Re2Re, *Status,
/**/	"Missing Re2Re ")
/**/DEBUG_CHECK_NULL_POINTER(VolumeInvDftRealImaginaryToRealImaginary, Im2Im, *Status,
/**/	"Missing Im2Im ")
/**/DEBUG_CHECK_RANGE_LONG(VolumeInvDftRealImaginaryToRealImaginary, Nx, 1L, LONG_MAX, *Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(VolumeInvDftRealImaginaryToRealImaginary, Ny, 1L, LONG_MAX, *Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(VolumeInvDftRealImaginaryToRealImaginary, Nz, 1L, LONG_MAX, *Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(VolumeInvDftRealImaginaryToRealImaginary, *Status)
/**/DEBUG_WRITE_ENTERING(VolumeInvDftRealImaginaryToRealImaginary,
/**/	"About to compute a 3D real/imaginary to real/imaginary direct Fourier transform")

	if (Nz > 1L) {
		*Status = bufferAllocRI(&Re, &Im, &TmpRe, &TmpIm, &CaS, Nz);
		if (*Status == ERROR) {
			bufferFreeRI(&Re, &Im, &TmpRe, &TmpIm, &CaS);
/**/		DEBUG_WRITE_LEAVING(VolumeInvDftRealImaginaryToRealImaginary, "Done")
			return(*Status);
		}
		*Status = GetCaS(CaS, Nz);
		if (*Status == ERROR) {
			bufferFreeRI(&Re, &Im, &TmpRe, &TmpIm, &CaS);
/**/		DEBUG_WRITE_LEAVING(VolumeInvDftRealImaginaryToRealImaginary, "Done")
			return(*Status);
		}
		for (y = 0L; (y < Ny); y++) {
			for (x = 0L; (x < Nx); x++) {
				*Status = GetzDoubleToDouble(Re2Re, Nx, Ny, Nz, x, y, 0L, Re, Nz);
				if (*Status == ERROR) {
					bufferFreeRI(&Re, &Im, &TmpRe, &TmpIm, &CaS);
/**/				DEBUG_WRITE_LEAVING(VolumeInvDftRealImaginaryToRealImaginary, "Done")
					return(*Status);
				}
				*Status = GetzDoubleToDouble(Im2Im, Nx, Ny, Nz, x, y, 0L, Im, Nz);
				if (*Status == ERROR) {
					bufferFreeRI(&Re, &Im, &TmpRe, &TmpIm, &CaS);
/**/				DEBUG_WRITE_LEAVING(VolumeInvDftRealImaginaryToRealImaginary, "Done")
					return(*Status);
				}
				*Status = InvDftRealImaginaryToRealImaginary(Re, Im, TmpRe, TmpIm, CaS, Nz);
				if (*Status == ERROR) {
					bufferFreeRI(&Re, &Im, &TmpRe, &TmpIm, &CaS);
/**/				DEBUG_WRITE_LEAVING(VolumeInvDftRealImaginaryToRealImaginary, "Done")
					return(*Status);
				}
				*Status = PutzDoubleToDouble(Re2Re, Nx, Ny, Nz, x, y, 0L, Re, Nz);
				if (*Status == ERROR) {
					bufferFreeRI(&Re, &Im, &TmpRe, &TmpIm, &CaS);
/**/				DEBUG_WRITE_LEAVING(VolumeInvDftRealImaginaryToRealImaginary, "Done")
					return(*Status);
				}
				*Status = PutzDoubleToDouble(Im2Im, Nx, Ny, Nz, x, y, 0L, Im, Nz);
				if (*Status == ERROR) {
					bufferFreeRI(&Re, &Im, &TmpRe, &TmpIm, &CaS);
/**/				DEBUG_WRITE_LEAVING(VolumeInvDftRealImaginaryToRealImaginary, "Done")
					return(*Status);
				}
			}
		}
		*Status = bufferFreeRI(&Re, &Im, &TmpRe, &TmpIm, &CaS);
		if (*Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(VolumeInvDftRealImaginaryToRealImaginary, "Done")
			return(*Status);
		}
	}
	if (Ny > 1L) {
		*Status = bufferAllocRI(&Re, &Im, &TmpRe, &TmpIm, &CaS, Ny);
		if (*Status == ERROR) {
			bufferFreeRI(&Re, &Im, &TmpRe, &TmpIm, &CaS);
/**/		DEBUG_WRITE_LEAVING(VolumeInvDftRealImaginaryToRealImaginary, "Done")
			return(*Status);
		}
		*Status = GetCaS(CaS, Ny);
		if (*Status == ERROR) {
			bufferFreeRI(&Re, &Im, &TmpRe, &TmpIm, &CaS);
/**/		DEBUG_WRITE_LEAVING(VolumeInvDftRealImaginaryToRealImaginary, "Done")
			return(*Status);
		}
		for (z = 0L; (z < Nz); z++) {
			for (x = 0L; (x < Nx); x++) {
				*Status = GetyDoubleToDouble(Re2Re, Nx, Ny, Nz, x, 0L, z, Re, Ny);
				if (*Status == ERROR) {
					bufferFreeRI(&Re, &Im, &TmpRe, &TmpIm, &CaS);
/**/				DEBUG_WRITE_LEAVING(VolumeInvDftRealImaginaryToRealImaginary, "Done")
					return(*Status);
				}
				*Status = GetyDoubleToDouble(Im2Im, Nx, Ny, Nz, x, 0L, z, Im, Ny);
				if (*Status == ERROR) {
					bufferFreeRI(&Re, &Im, &TmpRe, &TmpIm, &CaS);
/**/				DEBUG_WRITE_LEAVING(VolumeInvDftRealImaginaryToRealImaginary, "Done")
					return(*Status);
				}
				*Status = InvDftRealImaginaryToRealImaginary(Re, Im, TmpRe, TmpIm, CaS, Ny);
				if (*Status == ERROR) {
					bufferFreeRI(&Re, &Im, &TmpRe, &TmpIm, &CaS);
/**/				DEBUG_WRITE_LEAVING(VolumeInvDftRealImaginaryToRealImaginary, "Done")
					return(*Status);
				}
				*Status = PutyDoubleToDouble(Re2Re, Nx, Ny, Nz, x, 0L, z, Re, Ny);
				if (*Status == ERROR) {
					bufferFreeRI(&Re, &Im, &TmpRe, &TmpIm, &CaS);
/**/				DEBUG_WRITE_LEAVING(VolumeInvDftRealImaginaryToRealImaginary, "Done")
					return(*Status);
				}
				*Status = PutyDoubleToDouble(Im2Im, Nx, Ny, Nz, x, 0L, z, Im, Ny);
				if (*Status == ERROR) {
					bufferFreeRI(&Re, &Im, &TmpRe, &TmpIm, &CaS);
/**/				DEBUG_WRITE_LEAVING(VolumeInvDftRealImaginaryToRealImaginary, "Done")
					return(*Status);
				}
			}
		}
		*Status = bufferFreeRI(&Re, &Im, &TmpRe, &TmpIm, &CaS);
		if (*Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(VolumeInvDftRealImaginaryToRealImaginary, "Done")
			return(*Status);
		}
	}
	if (Nx > 1L) {
		*Status = bufferAllocRI(&Re, &Im, &TmpRe, &TmpIm, &CaS, Nx);
		if (*Status == ERROR) {
			bufferFreeRI(&Re, &Im, &TmpRe, &TmpIm, &CaS);
/**/		DEBUG_WRITE_LEAVING(VolumeInvDftRealImaginaryToRealImaginary, "Done")
			return(*Status);
		}
		*Status = GetCaS(CaS, Nx);
		if (*Status == ERROR) {
			bufferFreeRI(&Re, &Im, &TmpRe, &TmpIm, &CaS);
/**/		DEBUG_WRITE_LEAVING(VolumeInvDftRealImaginaryToRealImaginary, "Done")
			return(*Status);
		}
		for (z = 0L; (z < Nz); z++) {
			for (y = 0L; (y < Ny); y++) {
				*Status = GetxDoubleToDouble(Re2Re, Nx, Ny, Nz, 0L, y, z, Re, Nx);
				if (*Status == ERROR) {
					bufferFreeRI(&Re, &Im, &TmpRe, &TmpIm, &CaS);
/**/				DEBUG_WRITE_LEAVING(VolumeInvDftRealImaginaryToRealImaginary, "Done")
					return(*Status);
				}
				*Status = GetxDoubleToDouble(Im2Im, Nx, Ny, Nz, 0L, y, z, Im, Nx);
				if (*Status == ERROR) {
					bufferFreeRI(&Re, &Im, &TmpRe, &TmpIm, &CaS);
/**/				DEBUG_WRITE_LEAVING(VolumeInvDftRealImaginaryToRealImaginary, "Done")
					return(*Status);
				}
				*Status = InvDftRealImaginaryToRealImaginary(Re, Im, TmpRe, TmpIm, CaS, Nx);
				if (*Status == ERROR) {
					bufferFreeRI(&Re, &Im, &TmpRe, &TmpIm, &CaS);
/**/				DEBUG_WRITE_LEAVING(VolumeInvDftRealImaginaryToRealImaginary, "Done")
					return(*Status);
				}
				*Status = PutxDoubleToDouble(Re2Re, Nx, Ny, Nz, 0L, y, z, Re, Nx);
				if (*Status == ERROR) {
					bufferFreeRI(&Re, &Im, &TmpRe, &TmpIm, &CaS);
/**/				DEBUG_WRITE_LEAVING(VolumeInvDftRealImaginaryToRealImaginary, "Done")
					return(*Status);
				}
				*Status = PutxDoubleToDouble(Im2Im, Nx, Ny, Nz, 0L, y, z, Im, Nx);
				if (*Status == ERROR) {
					bufferFreeRI(&Re, &Im, &TmpRe, &TmpIm, &CaS);
/**/				DEBUG_WRITE_LEAVING(VolumeInvDftRealImaginaryToRealImaginary, "Done")
					return(*Status);
				}
			}
		}
		*Status = bufferFreeRI(&Re, &Im, &TmpRe, &TmpIm, &CaS);
		if (*Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(VolumeInvDftRealImaginaryToRealImaginary, "Done")
			return(*Status);
		}
	}
/**/DEBUG_WRITE_LEAVING(VolumeInvDftRealImaginaryToRealImaginary, "Done")
	return(*Status);
} /* end VolumeInvInvDftRealImaginaryToRealImaginary */

/*--------------------------------------------------------------------------*/
extern int		VolumeRealImaginaryToAmplitudePhase
				(
					double	Re2Am[],			/* real -> amplitude */
					double	Im2Ph[],			/* imaginary -> phase */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz					/* depth of the volume */
				)

/* converts a (real, imaginary) representation of a complex signal
	into an (amplitude, phase) representation */
/* the output phase is in [rad]; its domain is (-PI, PI) */
/* in-place processing */
/* the input signal (real Re2Am and imaginary Im2Ph) is
	replaced by the output signal (amplitude Re2Am and phase Im2Ph) */
/* Nx is the width of the volume */
/* Ny is the height of the volume */
/* Nz is the depth of the volume */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin VolumeRealImaginaryToAmplitudePhase */

	double	Phase;
	long	i;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(VolumeRealImaginaryToAmplitudePhase, Re2Am, Status,
/**/	"Missing Re2Am ")
/**/DEBUG_CHECK_NULL_POINTER(VolumeRealImaginaryToAmplitudePhase, Im2Ph, Status,
/**/	"Missing Im2Ph ")
/**/DEBUG_CHECK_RANGE_LONG(VolumeRealImaginaryToAmplitudePhase, Nx, 1L, LONG_MAX, Status,
/**/	"Invalid length (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(VolumeRealImaginaryToAmplitudePhase, Ny, 1L, LONG_MAX, Status,
/**/	"Invalid length (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(VolumeRealImaginaryToAmplitudePhase, Nz, 1L, LONG_MAX, Status,
/**/	"Invalid length (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(VolumeRealImaginaryToAmplitudePhase, Status)
/**/DEBUG_WRITE_ENTERING(VolumeRealImaginaryToAmplitudePhase,
/**/	"About to convert a real/imaginary representation to an amplitude/phase representation")

	for (i = -Nx * Ny * Nz; (i < 0L); Re2Am++, i++) {
		Phase = (double)atan2(*Im2Ph, *Re2Am);
		*Re2Am = (double)sqrt(*Re2Am * *Re2Am + *Im2Ph * *Im2Ph);
		*Im2Ph++ = Phase;
	}
/**/DEBUG_WRITE_LEAVING(VolumeRealImaginaryToAmplitudePhase, "Done")
	return(Status);
} /* end VolumeRealImaginaryToAmplitudePhase */

