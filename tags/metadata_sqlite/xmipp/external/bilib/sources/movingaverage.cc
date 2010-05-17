/*****************************************************************************
 *	System includes
 ****************************************************************************/
#include	<stddef.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>

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
#include	"fold.h"
#include	"getput.h"
#include	"messagedisplay.h"
#include	"movingaverage.h"

/*****************************************************************************
 *	Conditional includes
 ****************************************************************************/
#ifdef		DEBUG
#include	<limits.h>
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
extern int		MovingAverage
				(
					double	InputData[],		/* data to process */
					double	OutputData[],		/* result */
					long	SignalLength,		/* length of the 1D data array */
					long	KernelOrigin,		/* center of the kernel */
					long	KernelLength,		/* length of the 1D kernel */
					enum TBoundaryConvention
							Convention			/* boundary convention */
				)

/* the specified boundary convention applies to the input data only, not to the kernel */
/* the boundary convention applied to the kernel is FiniteDataSupport */
/* the kernel has a constant value and sums up to 1.0 */
/* the input and the output have the same length */
/* the origin for the kernel is given with respect to the leftmost sample [0] */
/* success: return(!ERROR); failure: return(ERROR) */
/* general structure is as follows:
	for (i = 0L; (i < SignalLength); i++) {
		Sum = 0.0;
		for (j = -Infinity; (j <= Infinity); j++) {
			Sum += InputData[j] * Kernel[KernelOrigin + i - j];
		}
		OutputData[i] = Sum;
	}
*/

{ /* begin MovingAverage */

	double	*p, *q;
	double	Sum, Norm;
	long	i, j;
	long	km, kp;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(MovingAverage, InputData, Status,
/**/	"No input data")
/**/DEBUG_CHECK_NULL_POINTER(MovingAverage, OutputData, Status,
/**/	"No output data")
/**/DEBUG_CHECK_RANGE_LONG(MovingAverage, SignalLength, 1L, LONG_MAX, Status,
/**/	"Invalid signal length (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MovingAverage, KernelLength, 1L, LONG_MAX, Status,
/**/	"Invalid kernel length (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(MovingAverage, Status)
/**/DEBUG_WRITE_ENTERING(MovingAverage,
/**/	"About to perform a convolution by moving average")

	Norm = 1.0 / (double)KernelLength;
	switch (Convention) {
		case AntiMirrorOnBounds:
			if (SignalLength == 1L) {
				*OutputData = *InputData;
				break;
			}
			Status = ERROR;
			WRITE_ERROR(MovingAverage, "Not yet implemented")
			break;
		case FiniteDataSupport:
			Sum = 0.0;
			KernelOrigin = KernelLength - KernelOrigin - 1L;
			km = (KernelOrigin < 0L) ? (-KernelOrigin) : (0L);
			kp = ((KernelLength - KernelOrigin) < SignalLength) ? (KernelLength - KernelOrigin)
				: (SignalLength);
			p = InputData + (ptrdiff_t)km;
			for (j = km - kp; (j < 0L); j++) {
				Sum += *p++;
			}
			*OutputData++ = Sum * Norm;
			km = -KernelOrigin;
			kp = km + KernelLength;
			p = InputData + (ptrdiff_t)km;
			q = InputData + (ptrdiff_t)kp;
			i = 1L - SignalLength;
			if (i < kp) {
				if (kp < 0L) {
					p -= (ptrdiff_t)kp;
					q -= (ptrdiff_t)kp;
					km -= kp;
					i -= kp;
					while (kp++ < 0L) {
						*OutputData++ = Sum * Norm;
					}
				}
			}
			else {
				Sum *= Norm;
				while (i++ < 0L) {
					*OutputData++ = Sum;
				}
				break;
			}
			if (i < km) {
				if (km < (kp - SignalLength)) {
					if (kp < SignalLength) {
						p += (ptrdiff_t)(SignalLength - kp);
						km += SignalLength - kp;
						i += SignalLength - kp;
						while (kp++ < SignalLength) {
							Sum += *q++;
							*OutputData++ = Sum * Norm;
						}
					}
				}
				else {
					if (km < 0L) {
						p -= (ptrdiff_t)km;
						kp -= km;
						i -= km;
						while (km++ < 0L) {
							Sum += *q++;
							*OutputData++ = Sum * Norm;
						}
					}
				}
			}
			else {
				if (i < (kp - SignalLength)) {
					if (kp < SignalLength) {
						p += (ptrdiff_t)(SignalLength - kp);
						km += SignalLength - kp;
						i += SignalLength - kp;
						while (kp++ < SignalLength) {
							Sum += *q++;
							*OutputData++ = Sum * Norm;
						}
					}
				}
				else {
					while (i++ < 0L) {
						Sum += *q++;
						*OutputData++ = Sum * Norm;
					}
					break;
				}
			}
			if (i < km) {
				if (km < 0L) {
					p -= (ptrdiff_t)km;
					q -= (ptrdiff_t)km;
					kp -= km;
					i -= km;
					while (km++ < 0L) {
						*OutputData++ = Sum * Norm;
					}
				}
			}
			else {
				Sum *= Norm;
				while (i++ < 0L) {
					*OutputData++ = Sum;
				}
				break;
			}
			if (i < (kp - SignalLength)) {
				if (kp < SignalLength) {
					km += SignalLength - kp;
					i += SignalLength - kp;
					while (kp++ < SignalLength) {
						Sum -= *p++;
						Sum += *q++;
						*OutputData++ = Sum * Norm;
					}
				}
			}
			else {
				while (i++ < 0L) {
					Sum -= *p++;
					Sum += *q++;
					*OutputData++ = Sum * Norm;
				}
				break;
			}
			if (i < (km - SignalLength)) {
				i += SignalLength - km;
				while (km++ < SignalLength) {
					Sum -= *p++;
					*OutputData++ = Sum * Norm;
				}
			}
			else {
				while (i++ < 0L) {
					Sum -= *p++;
					*OutputData++ = Sum * Norm;
				}
				break;
			}
			Sum *= Norm;
			while (i++ < 0L) {
				*OutputData++ = Sum;
			}
			break;
		case MirrorOffBounds:
			if (SignalLength == 1L) {
				*OutputData = *InputData;
				break;
			}
			Status = ERROR;
			WRITE_ERROR(MovingAverage, "Not yet implemented")
			break;
		case MirrorOnBounds:
			if (SignalLength == 1L) {
				*OutputData = *InputData;
				break;
			}
			Status = ERROR;
			WRITE_ERROR(MovingAverage, "Not yet implemented")
			break;
		case Periodic:
			if (SignalLength == 1L) {
				*OutputData = *InputData;
				break;
			}
			Sum = 0.0;
			if (KernelLength < SignalLength) {
				km = KernelOrigin - KernelLength + 1L;
				km += (km < 0L)
					? (SignalLength * ((SignalLength - km - 1L) / SignalLength))
					: (-SignalLength * (km / SignalLength));
				kp = KernelOrigin;
				kp += (kp < 0L)
					? (SignalLength * ((SignalLength - kp - 1L) / SignalLength))
					: (-SignalLength * (kp / SignalLength));
				if (km < kp) {
					p = InputData + (ptrdiff_t)km;
					q = p;
					for (j = -KernelLength; (j < 0L); j++) {
						Sum += *q++;
					}
					*OutputData++ = Sum * Norm;
					i = -kp;
					j = -km;
					km += SignalLength - kp - 1L;
					while (++kp < SignalLength) {
						Sum -= *p++;
						Sum += *q++;
						*OutputData++ = Sum * Norm;
					}
					q = InputData;
					kp = 0L;
					if ((km - SignalLength) <= i) {
						while (i++ < 0L) {
							Sum -= *p++;
							Sum += *q++;
							*OutputData++ = Sum * Norm;
						}
						break;
					}
					else {
						while (km++ < SignalLength) {
							Sum -= *p++;
							Sum += *q++;
							*OutputData++ = Sum * Norm;
						}
						p = InputData;
						while (++j < 0L) {
							Sum -= *p++;
							Sum += *q++;
							*OutputData++ = Sum * Norm;
						}
						break;
					}
				}
				else if (km == kp) {
					kp = SignalLength - km;
					OutputData = (double *)memcpy(OutputData, InputData + (ptrdiff_t)km,
						(size_t)(kp * (long)sizeof(double)));
					OutputData += (ptrdiff_t)kp;
					OutputData = (double *)memcpy(OutputData, InputData,
						(size_t)((SignalLength - kp) * (long)sizeof(double)));
				}
				else {
					q = InputData + (ptrdiff_t)kp;
					p = q;
					do {
						Sum += *p--;
					} while (InputData <= p);
					p = InputData + (ptrdiff_t)km;
					for (j = km - SignalLength; (j < 0L); j++) {
						Sum += *p++;
					}
					*OutputData++ = Sum * Norm;
					i = -km;
					j = -kp;
					p = InputData + (ptrdiff_t)km;
					kp += SignalLength - km + 1L;
					while (km++ < SignalLength) {
						Sum -= *p++;
						Sum += *++q;
						*OutputData++ = Sum * Norm;
					}
					p = InputData;
					km = 0L;
					while (kp++ < SignalLength) {
						Sum -= *p++;
						Sum += *++q;
						*OutputData++ = Sum * Norm;
					}
					q = InputData;
					while (j++ < 0L) {
						Sum -= *p++;
						Sum += *q++;
						*OutputData++ = Sum * Norm;
					}
					break;
				}
			}
			if (KernelLength == SignalLength) {
				for (j = -SignalLength; (j < 0L); j++) {
					Sum += *InputData++;
				}
				Sum *= Norm;
				for (j = -SignalLength; (j < 0L); j++) {
					*OutputData++ = Sum;
				}
				break;
			}
			if (SignalLength < KernelLength) {
/*
@
*/
			}
			Status = ERROR;
			WRITE_ERROR(MovingAverage, "Not yet implemented")
			break;
		default:
			Status = ERROR;
			WRITE_ERROR(MovingAverage, "Invalid boundary convention")
			break;
	}
/**/DEBUG_WRITE_LEAVING(MovingAverage, "Done")
	return(Status);
} /* end MovingAverage */

/*--------------------------------------------------------------------------*/
extern int		MovingAverageVolume
				(
					float	*VolumeSource,		/* data to process */
					float	*VolumeDestination,	/* result */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					long	KernelOrigin,		/* center of the kernel */
					long	KernelLength,		/* length of the 1D kernel */
					enum TBoundaryConvention
							Convention			/* boundary convention */
				)

/* the specified boundary convention applies to the input data only, not to the kernel */
/* the boundary convention applied to the kernel is FiniteDataSupport */
/* VolumeSource is a (float)volume of size (Nx x Ny x Nz) */
/* OutputData is a (float)volume of size (Nx x Ny x Nz) */
/* the origin for the kernel is given with respect to the leftmost sample [0] */
/* the 1D kernel is applied successively to each principal direction in a separable fashion */
/* success: return(!ERROR); failure: return(ERROR) */
/* general structure is as follows:
	for (i = 0L; (i < SignalLength); i++) {
		Sum = 0.0;
		for (j = i - KernelOrigin; (j < (KernelLength + i - KernelOrigin)); j++) {
			Sum += InputData[j];
		}
		OutputData[i] = Sum / KernelLength;
	}
*/

{ /* begin MovingAverageVolume */

	double	*InBuffer = (double *)NULL, *OutBuffer = (double *)NULL;
	long	i, j, k;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(MovingAverageVolume, VolumeSource, Status,
/**/	"No VolumeSource")
/**/DEBUG_CHECK_NULL_POINTER(MovingAverageVolume, VolumeDestination, Status,
/**/	"No OutputData")
/**/DEBUG_CHECK_RANGE_LONG(MovingAverageVolume, Nx, 1L, LONG_MAX, Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MovingAverageVolume, Ny, 1L, LONG_MAX, Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MovingAverageVolume, Nz, 1L, LONG_MAX, Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MovingAverageVolume, KernelLength, 1L, LONG_MAX, Status,
/**/	"Invalid kernel length (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(MovingAverageVolume, Status)
/**/DEBUG_WRITE_ENTERING(MovingAverageVolume,
/**/	"About to perform IIR recursive convolution for a volume")

	switch (Convention) {
		case AntiMirrorOnBounds:
		case FiniteCoefficientSupport:
		case FiniteDataSupport:
		case MirrorOffBounds:
		case MirrorOnBounds:
		case Periodic:
			break;
		default:
			Status = ERROR;
			WRITE_ERROR(MovingAverageVolume, "Invalid boundary convention")
/**/		DEBUG_WRITE_LEAVING(MovingAverageVolume, "Done")
			return(Status);
	}
	if (1L < Nx) {
		AllocateLineDouble(&InBuffer, Nx, &Status);
		if (Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(MovingAverageVolume, "Done")
			return(Status);
		}
		AllocateLineDouble(&OutBuffer, Nx, &Status);
		if (Status == ERROR) {
			FreeLineDouble(&InBuffer);
/**/		DEBUG_WRITE_LEAVING(MovingAverageVolume, "Done")
			return(Status);
		}
		for (k = 0L; (k < Nz); k++) {
			for (j = 0L; (j < Ny); j++) {
				Status = GetxFloatToDouble(VolumeSource, Nx, Ny, Nz, 0L, j, k, InBuffer, Nx);
				if (Status == ERROR) {
					FreeLineDouble(&OutBuffer);
					FreeLineDouble(&InBuffer);
/**/				DEBUG_WRITE_LEAVING(MovingAverageVolume, "Done")
					return(Status);
				}
				Status = MovingAverage(InBuffer, OutBuffer, Nx,
					KernelOrigin, KernelLength, Convention);
				if (Status == ERROR) {
					FreeLineDouble(&OutBuffer);
					FreeLineDouble(&InBuffer);
/**/				DEBUG_WRITE_LEAVING(MovingAverageVolume, "Done")
					return(Status);
				}
				Status = PutxDoubleToFloat(VolumeDestination, Nx, Ny, Nz, 0L, j, k, OutBuffer, Nx);
				if (Status == ERROR) {
					FreeLineDouble(&OutBuffer);
					FreeLineDouble(&InBuffer);
/**/				DEBUG_WRITE_LEAVING(MovingAverageVolume, "Done")
					return(Status);
				}
			}
		}
		Status = FreeLineDouble(&OutBuffer);
		if (Status == ERROR) {
			FreeLineDouble(&InBuffer);
/**/		DEBUG_WRITE_LEAVING(MovingAverageVolume, "Done")
			return(Status);
		}
		Status = FreeLineDouble(&InBuffer);
		if (Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(MovingAverageVolume, "Done")
			return(Status);
		}
	}
	else {
		VolumeDestination = (float *)memcpy(VolumeDestination, VolumeSource,
			(size_t)(Nx * Ny * Nz * (long)sizeof(float)));
	}
	if (1L < Ny) {
		AllocateLineDouble(&InBuffer, Ny, &Status);
		if (Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(MovingAverageVolume, "Done")
			return(Status);
		}
		AllocateLineDouble(&OutBuffer, Ny, &Status);
		if (Status == ERROR) {
			FreeLineDouble(&InBuffer);
/**/		DEBUG_WRITE_LEAVING(MovingAverageVolume, "Done")
			return(Status);
		}
		for (k = 0L; (k < Nz); k++) {
			for (i = 0L; (i < Nx); i++) {
				Status = GetyFloatToDouble(VolumeDestination, Nx, Ny, Nz, i, 0L, k, InBuffer, Ny);
				if (Status == ERROR) {
					FreeLineDouble(&OutBuffer);
					FreeLineDouble(&InBuffer);
/**/				DEBUG_WRITE_LEAVING(MovingAverageVolume, "Done")
					return(Status);
				}
				Status = MovingAverage(InBuffer, OutBuffer, Ny,
					KernelOrigin, KernelLength, Convention);
				if (Status == ERROR) {
					FreeLineDouble(&OutBuffer);
					FreeLineDouble(&InBuffer);
/**/				DEBUG_WRITE_LEAVING(MovingAverageVolume, "Done")
					return(Status);
				}
				Status = PutyDoubleToFloat(VolumeDestination, Nx, Ny, Nz, i, 0L, k, OutBuffer, Ny);
				if (Status == ERROR) {
					FreeLineDouble(&OutBuffer);
					FreeLineDouble(&InBuffer);
/**/				DEBUG_WRITE_LEAVING(MovingAverageVolume, "Done")
					return(Status);
				}
			}
		}
		Status = FreeLineDouble(&OutBuffer);
		if (Status == ERROR) {
			FreeLineDouble(&InBuffer);
/**/		DEBUG_WRITE_LEAVING(MovingAverageVolume, "Done")
			return(Status);
		}
		Status = FreeLineDouble(&InBuffer);
		if (Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(MovingAverageVolume, "Done")
			return(Status);
		}
	}
	if (1L < Nz) {
		AllocateLineDouble(&InBuffer, Nz, &Status);
		if (Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(MovingAverageVolume, "Done")
			return(Status);
		}
		AllocateLineDouble(&OutBuffer, Nz, &Status);
		if (Status == ERROR) {
			FreeLineDouble(&InBuffer);
/**/		DEBUG_WRITE_LEAVING(MovingAverageVolume, "Done")
			return(Status);
		}
		for (j = 0L; (j < Ny); j++) {
			for (i = 0L; (i < Nx); i++) {
				Status = GetzFloatToDouble(VolumeDestination, Nx, Ny, Nz, i, j, 0L, InBuffer, Nz);
				if (Status == ERROR) {
					FreeLineDouble(&OutBuffer);
					FreeLineDouble(&InBuffer);
/**/				DEBUG_WRITE_LEAVING(MovingAverageVolume, "Done")
					return(Status);
				}
				Status = MovingAverage(InBuffer, OutBuffer, Nz,
					KernelOrigin, KernelLength, Convention);
				if (Status == ERROR) {
					FreeLineDouble(&OutBuffer);
					FreeLineDouble(&InBuffer);
/**/				DEBUG_WRITE_LEAVING(MovingAverageVolume, "Done")
					return(Status);
				}
				Status = PutzDoubleToFloat(VolumeDestination, Nx, Ny, Nz, i, j, 0L, OutBuffer, Nz);
				if (Status == ERROR) {
					FreeLineDouble(&OutBuffer);
					FreeLineDouble(&InBuffer);
/**/				DEBUG_WRITE_LEAVING(MovingAverageVolume, "Done")
					return(Status);
				}
			}
		}
		Status = FreeLineDouble(&OutBuffer);
		if (Status == ERROR) {
			FreeLineDouble(&InBuffer);
/**/		DEBUG_WRITE_LEAVING(MovingAverageVolume, "Done")
			return(Status);
		}
		Status = FreeLineDouble(&InBuffer);
		if (Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(MovingAverageVolume, "Done")
			return(Status);
		}
	}
/**/DEBUG_WRITE_LEAVING(MovingAverageVolume, "Done")
	return(Status);
} /* end MovingAverageVolume */

