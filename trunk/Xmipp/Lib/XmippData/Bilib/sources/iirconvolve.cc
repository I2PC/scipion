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
#include	"error.h"
#include	"tboundaryconvention.h"

/*****************************************************************************
 *	Toolbox includes
 ****************************************************************************/
#include	"../headers/convert.h"
#include	"getput.h"
#include	"getputd.h"
#include	"iirconvolve.h"
#include	"messagedisplay.h"
#include	"positivepower.h"

/*****************************************************************************
 *	Conditional includes
 ****************************************************************************/
#ifdef		DEBUG
#include	<float.h>
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
extern int		IirConvolveCanonicProgressive
				(
					double	InputData[],		/* data to process */
					double	OutputData[],		/* result */
					long	SignalLength,		/* length of the 1D data array */
					double	Kernel[],			/* kernel */
					double	RightInit[],		/* progressive recursion initialization */
					long	KernelLength		/* length of the 1D kernel */
				)

/* computes OutputData[i] = InputData[i] + SUM(k): OutputData[i + k + 1] * Kernel[k] */
/* the input and the output have the same length */
/* the kernel and the initial output values have the same length */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin IirConvolveCanonicProgressive */

	double	*p, *q;
	double	Sum;
	long	i, j;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(IirConvolveCanonicProgressive, InputData, Status,
/**/	"No input data")
/**/DEBUG_CHECK_NULL_POINTER(IirConvolveCanonicProgressive, OutputData, Status,
/**/	"No output data")
/**/DEBUG_CHECK_RANGE_LONG(IirConvolveCanonicProgressive, SignalLength, 1L, LONG_MAX, Status,
/**/	"Invalid length (should be strictly positive)")
/**/DEBUG_CHECK_NULL_POINTER(IirConvolveCanonicProgressive, Kernel, Status,
/**/	"No kernel")
/**/DEBUG_CHECK_NULL_POINTER(IirConvolveCanonicProgressive, RightInit, Status,
/**/	"No initialization")
/**/DEBUG_CHECK_RANGE_LONG(IirConvolveCanonicProgressive, KernelLength, 1L, LONG_MAX, Status,
/**/	"Invalid length (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(IirConvolveCanonicProgressive, Status)
/**/DEBUG_WRITE_ENTERING(IirConvolveCanonicProgressive,
/**/	"About to perform IIR progressive canonic convolution")

/*
	if (KernelLength <= SignalLength) {
		for (i = SignalLength - 1L; ((SignalLength - KernelLength) <= i); i--) {
			Sum = InputData[i];
			for (j = KernelLength - 1L; ((SignalLength - i - 1L) <= j); j--) {
				Sum += RightInit[i - SignalLength + j + 1L] * Kernel[j];
			}
			while (0L <= j) {
				Sum += OutputData[i + j + 1L] * Kernel[j];
				j--;
			}
			OutputData[i] = Sum;
		}
		while (0L <= i) {
			Sum = InputData[i];
			for (j = KernelLength - 1L; (0L <= j); j--) {
				Sum += OutputData[i + j + 1L] * Kernel[j];
			}
			OutputData[i] = Sum;
			i--;
		}
	}
	else {
		for (i = SignalLength - 1L; (0L <= i); i--) {
			Sum = InputData[i];
			for (j = KernelLength - 1L; ((SignalLength - i - 1L) <= j); j--) {
				Sum += RightInit[i - SignalLength + j + 1L] * Kernel[j];
			}
			while (0L <= j) {
				Sum += OutputData[i + j + 1L] * Kernel[j];
				j--;
			}
			OutputData[i] = Sum;
		}
	}
*/
	InputData += (ptrdiff_t)(SignalLength - 1L);
	OutputData += (ptrdiff_t)(SignalLength - 1L);
	RightInit += (ptrdiff_t)(KernelLength - 1L);
	Kernel += (ptrdiff_t)(KernelLength - 1L);
	if (KernelLength <= SignalLength) {
		for (i = SignalLength - 1L; ((SignalLength - KernelLength) <= i); i--) {
			Sum = *InputData--;
			p = RightInit--;
			q = Kernel;
			for (j = KernelLength - 1L; ((SignalLength - i - 1L) <= j); j--) {
				Sum += *p-- * *q--;
			}
			p = OutputData + (ptrdiff_t)(j + 1L);
			while (0L <= j--) {
				Sum += *p-- * *q--;
			}
			*OutputData-- = Sum;
		}
		while (0L <= i--) {
			Sum = *InputData--;
			p = OutputData + (ptrdiff_t)KernelLength;
			q = Kernel;
			for (j = KernelLength - 1L; (0L <= j); j--) {
				Sum += *p-- * *q--;
			}
			*OutputData-- = Sum;
		}
	}
	else {
		for (i = SignalLength - 1L; (0L <= i); i--) {
			Sum = *InputData--;
			p = RightInit--;
			q = Kernel;
			for (j = KernelLength - 1L; ((SignalLength - i - 1L) <= j); j--) {
				Sum += *p-- * *q--;
			}
			p = OutputData + (ptrdiff_t)(j + 1L);
			while (0L <= j--) {
				Sum += *p-- * *q--;
			}
			*OutputData-- = Sum;
		}
	}
/**/DEBUG_WRITE_LEAVING(IirConvolveCanonicProgressive, "Done")
	return(Status);
} /* end IirConvolveCanonicProgressive */

/*--------------------------------------------------------------------------*/
extern int		IirConvolveCanonicRegressive
				(
					double	InputData[],		/* data to process */
					double	OutputData[],		/* result */
					long	SignalLength,		/* length of the 1D data array */
					double	Kernel[],			/* kernel */
					double	LeftInit[],			/* regressive recursion initialization */
					long	KernelLength		/* length of the 1D kernel */
				)

/* computes OutputData[i] = InputData[i] + SUM(k): OutputData[i - k - 1] * Kernel[k] */
/* the input and the output have the same length */
/* the kernel and the initial output values have the same length */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin IirConvolveCanonicRegressive */

	double	*p, *q;
	double	Sum;
	long	i, j;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(IirConvolveCanonicRegressive, InputData, Status,
/**/	"No input data")
/**/DEBUG_CHECK_NULL_POINTER(IirConvolveCanonicRegressive, OutputData, Status,
/**/	"No output data")
/**/DEBUG_CHECK_RANGE_LONG(IirConvolveCanonicRegressive, SignalLength, 1L, LONG_MAX, Status,
/**/	"Invalid length (should be strictly positive)")
/**/DEBUG_CHECK_NULL_POINTER(IirConvolveCanonicRegressive, Kernel, Status,
/**/	"No kernel")
/**/DEBUG_CHECK_NULL_POINTER(IirConvolveCanonicRegressive, LeftInit, Status,
/**/	"No initialization")
/**/DEBUG_CHECK_RANGE_LONG(IirConvolveCanonicRegressive, KernelLength, 1L, LONG_MAX, Status,
/**/	"Invalid length (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(IirConvolveCanonicRegressive, Status)
/**/DEBUG_WRITE_ENTERING(IirConvolveCanonicRegressive,
/**/	"About to perform IIR regressive canonic convolution")

/*
	if (KernelLength < SignalLength) {
		for (i = 0L; (i <= KernelLength); i++) {
			Sum = InputData[i];
			for (j = 0L; (j < i); j++) {
				Sum += OutputData[i - j - 1L] * Kernel[j];
			}
			while (j < KernelLength) {
				Sum += LeftInit[KernelLength + i - j - 1L] * Kernel[j];
				j++;
			}
			OutputData[i] = Sum;
		}
		while (i < SignalLength) {
			Sum = InputData[i];
			for (j = 0L; (j < KernelLength); j++) {
				Sum += OutputData[i - j - 1L] * Kernel[j];
			}
			OutputData[i] = Sum;
			i++;
		}
	}
	else {
		for (i = 0L; (i < SignalLength); i++) {
			Sum = InputData[i];
			for (j = 0L; (j < i); j++) {
				Sum += OutputData[i - j - 1L] * Kernel[j];
			}
			while (j < KernelLength) {
				Sum += LeftInit[KernelLength + i - j - 1L] * Kernel[j];
				j++;
			}
			OutputData[i] = Sum;
		}
	}
*/
	LeftInit += (ptrdiff_t)(KernelLength - 1L);
	if (KernelLength < SignalLength) {
		for (i = 0L; (i <= KernelLength); i++) {
			Sum = *InputData++;
			p = OutputData - (ptrdiff_t)1L;
			q = Kernel;
			for (j = 0L; (j < i); j++) {
				Sum += *p-- * *q++;
			}
			p = LeftInit;
			while (j++ < KernelLength) {
				Sum += *p-- * *q++;
			}
			*OutputData++ = Sum;
		}
		while (i++ < SignalLength) {
			Sum = *InputData++;
			p = OutputData - (ptrdiff_t)1L;
			q = Kernel;
			for (j = 0L; (j < KernelLength); j++) {
				Sum += *p-- * *q++;
			}
			*OutputData++ = Sum;
		}
	}
	else {
		for (i = 0L; (i < SignalLength); i++) {
			Sum = *InputData++;
			p = OutputData - (ptrdiff_t)1L;
			q = Kernel;
			for (j = 0L; (j < i); j++) {
				Sum += *p-- * *q++;
			}
			p = LeftInit;
			while (j++ < KernelLength) {
				Sum += *p-- * *q++;
			}
			*OutputData++ = Sum;
		}
	}
/**/DEBUG_WRITE_LEAVING(IirConvolveCanonicRegressive, "Done")
	return(Status);
} /* end IirConvolveCanonicRegressive */

/*--------------------------------------------------------------------------*/
extern int		IirConvolvePoles
				(
					double	InputData[],		/* data to process */
					double	OutputData[],		/* result */
					long	SignalLength,		/* length of the 1D data array */
					double	RealPoles[],		/* array of real poles */
					long	PoleNumber,			/* number of poles */
					enum TBoundaryConvention
							Convention,			/* boundary convention */
					double	Tolerance			/* admissible relative error */
				)

/* the input and the output have the same length */
/* in-place processing is allowed */
/* no more than two poles are allowed for a finite support boundary */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin IirConvolvePoles */

	double	*p, *q;
	double	Sum, Gain;
	double	z, z0, z1, z2, iz;
	long	Horizon;
	long	i, j;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(IirConvolvePoles, InputData, Status,
/**/	"No input data")
/**/DEBUG_CHECK_NULL_POINTER(IirConvolvePoles, OutputData, Status,
/**/	"No output data")
/**/DEBUG_CHECK_RANGE_LONG(IirConvolvePoles, SignalLength, 1L, LONG_MAX, Status,
/**/	"Invalid length (should be strictly positive)")
/**/DEBUG_CHECK_NULL_POINTER(IirConvolvePoles, RealPoles, Status,
/**/	"No poles")
/**/DEBUG_CHECK_RANGE_LONG(IirConvolvePoles, PoleNumber, 1L, LONG_MAX, Status,
/**/	"Invalid number of poles (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_DOUBLE(IirConvolvePoles, Tolerance, 0.0, DBL_MAX, Status,
/**/	"Invalid Tolerance (should be positive)")
/**/DEBUG_RETURN_ON_ERROR(IirConvolvePoles, Status)
/**/DEBUG_WRITE_ENTERING(IirConvolvePoles,
/**/	"About to perform IIR recursive convolution")

	Gain = 1.0;
	for (i = 0L; (i < PoleNumber); i++) {
		z = RealPoles[i];
		if (fabs(z) < Tolerance) {
			Status = ERROR;
			WRITE_ERROR(IirConvolvePoles, "Invalid pole (too small)")
/**/		DEBUG_WRITE_LEAVING(IirConvolvePoles, "Done")
			return(Status);
		}
		z1 = 1.0 - z;
		Gain *= -(z1 * z1) / z;
	}
	p = InputData;
	q = OutputData;
	for (i = -SignalLength; (i < 0L); i++) {
		*q++ = *p++ * Gain;
	}
	switch (Convention) {
		case AntiMirrorOnBounds:
			if (SignalLength == 1L) {
				*OutputData = *InputData;
/**/			DEBUG_WRITE_LEAVING(IirConvolvePoles, "Done")
				return(Status);
			}
			for (i = 0L; (i < PoleNumber); i++) {
				z = RealPoles[i];
				Sum = ((1.0 + z) / (1.0 - z)) * (OutputData[0]
					- PositiveIntPower(z, SignalLength - 1L) * OutputData[SignalLength - 1L]);
				z1 = z;
				z2 = PositiveIntPower(z, 2L * SignalLength - 3L);
				iz = 1.0 / z;
				for (j = 1L; (j < (SignalLength - 1L)); j++) {
					Sum += (z2 - z1) * OutputData[j];
					z1 *= z;
					z2 *= iz;
				}
				OutputData[0] = Sum / (1.0 - PositiveIntPower(z, 2L * SignalLength - 2L));
				for (j = 1L; (j < SignalLength); j++) {
					OutputData[j] += z * OutputData[j - 1L];
				}
				OutputData[SignalLength - 1L] = (-z / ((1.0 - z) * (1.0 - z)))
					* (OutputData[SignalLength - 1L] - z * OutputData[SignalLength - 2L]);
				for (j = SignalLength - 2L; (0L <= j); j--) {
					OutputData[j] = z * (OutputData[j + 1L] - OutputData[j]);
				}
			}
			break;
		case FiniteCoefficientSupport:
			switch (PoleNumber) {
				case 1L:
					z = RealPoles[0];
					Sum = 0.0;
					z1 = z;
					z2 = PositiveIntPower(z, 2L * SignalLength - 1L);
					iz = 1.0 / z;
					for (j = 1L; (j < SignalLength); j++) {
						Sum += (z1 - z2) * OutputData[j];
						z1 *= z;
						z2 *= iz;
					}
					z2 = z * z;
					OutputData[0] -= Sum * z2 / (1.0 - z2);
					OutputData[0] *= (1.0 - z2)
						/ (1.0 - PositiveIntPower(z2, SignalLength + 1L));
					for (j = 1L; (j < SignalLength); j++) {
						OutputData[j] += z * OutputData[j - 1L];
					}
					OutputData[SignalLength - 1L] *= -z;
					for (j = SignalLength - 2L; (0L <= j); j--) {
						OutputData[j] = z * (OutputData[j + 1L] - OutputData[j]);
					}
					break;
				default:
					Status = ERROR;
					WRITE_ERROR(IirConvolvePoles,
						"Invalid number of poles (should be 1 for FiniteCoefficientSupport)")
					break;
			}
			break;
		case FiniteDataSupport:
			switch (PoleNumber) {
				case 1L:
					z = RealPoles[0];
					for (j = 1L; (j < SignalLength); j++) {
						OutputData[j] += z * OutputData[j - 1L];
					}
					OutputData[SignalLength - 1L] *= z / (z * z - 1.0);
					for (j = SignalLength - 2L; (0L <= j); j--) {
						OutputData[j] = z * (OutputData[j + 1L] - OutputData[j]);
					}
					break;
				case 2L:
					z0 = RealPoles[0];
					for (j = 1L; (j < SignalLength); j++) {
						OutputData[j] += z0 * OutputData[j - 1L];
					}
					OutputData[SignalLength - 1L] *= z0 / (z0 * z0 - 1.0);
					z = OutputData[SignalLength - 1L];
					for (j = SignalLength - 2L; (0L <= j); j--) {
						OutputData[j] = z0 * (OutputData[j + 1L] - OutputData[j]);
					}
					z1 = RealPoles[1];
					OutputData[0] /= 1.0 - z0 * z1;
					for (j = 1L; (j < SignalLength); j++) {
						OutputData[j] += z1 * OutputData[j - 1L];
					}
					OutputData[SignalLength - 1L] = (z1 / (z1 * z1 - 1.0))
						* (OutputData[SignalLength - 1L] - (z0 * z1 / (z0 * z1 - 1.0)) * z);
					for (j = SignalLength - 2L; (0L <= j); j--) {
						OutputData[j] = z1 * (OutputData[j + 1L] - OutputData[j]);
					}
					break;
				default:
					Status = ERROR;
					WRITE_ERROR(IirConvolvePoles,
						"Invalid number of poles (should be 1 or 2 for FiniteDataSupport)")
					break;
			}
			break;
		case MirrorOffBounds:
			if (SignalLength == 1L) {
				*OutputData = *InputData;
/**/			DEBUG_WRITE_LEAVING(IirConvolvePoles, "Done")
				return(Status);
			}
			for (i = 0L; (i < PoleNumber); i++) {
				z = RealPoles[i];
				Sum = (OutputData[0] + PositiveIntPower(z, SignalLength)
					* OutputData[SignalLength - 1L]) * (1.0 + z) / z;
				z1 = z;
				z2 = PositiveIntPower(z, 2L * SignalLength - 2L);
				iz = 1.0 / z;
				for (j = 1L; (j < (SignalLength - 1L)); j++) {
					Sum += (z2 + z1) * OutputData[j];
					z1 *= z;
					z2 *= iz;
				}
				OutputData[0] = Sum * z / (1.0 - PositiveIntPower(z, 2L * SignalLength));
				for (j = 1L; (j < SignalLength); j++) {
					OutputData[j] += z * OutputData[j - 1L];
				}
				OutputData[SignalLength - 1L] *= z / (z - 1.0);
				for (j = SignalLength - 2L; (0L <= j); j--) {
					OutputData[j] = z * (OutputData[j + 1L] - OutputData[j]);
				}
			}
			break;
		case MirrorOnBounds:
			if (SignalLength == 1L) {
				*OutputData = *InputData;
/**/			DEBUG_WRITE_LEAVING(IirConvolvePoles, "Done")
				return(Status);
			}
/*
			for (i = 0L; (i < PoleNumber); i++) {
				z = RealPoles[i];
				Sum = OutputData[0] + PositiveIntPower(z, SignalLength - 1L)
					* OutputData[SignalLength - 1L];
				for (j = 1L; (j < (SignalLength - 1L)); j++) {
					Sum += (PositiveIntPower(z, 2L * SignalLength - 2L - j)
						+ PositiveIntPower(z, j)) * OutputData[j];
				}
				OutputData[0] = Sum / (1.0 - PositiveIntPower(z, 2L * SignalLength - 2L));
				for (j = 1L; (j < SignalLength); j++) {
					OutputData[j] += z * OutputData[j - 1L];
				}
				OutputData[SignalLength - 1L] = (z / (z * z - 1.0))
					* (z * OutputData[SignalLength - 2L] + OutputData[SignalLength - 1L]);
				for (j = SignalLength - 2L; (0L <= j); j--) {
					OutputData[j] = z * (OutputData[j + 1L] - OutputData[j]);
				}
			}
*/
			for (i = 0L; (i < PoleNumber); i++) {
				z = *RealPoles++;
				if (Tolerance == 0.0) {
					Horizon = SignalLength;
				}
				else {
					Horizon = ConvertDoubleToLong(ceil(log(Tolerance) / log(fabs(z))));
				}
				if (Horizon < SignalLength) {
					z1 = z;
					Sum = *OutputData;
					p = OutputData;
					for (j = -Horizon; (j < 0L); j++) {
						Sum += z1 * *++p;
						z1 *= z;
					}
					p = OutputData;
					q = OutputData;
					*p++ = Sum;
					for (j = 1L - SignalLength; (j < 0L); j++) {
						*p++ += z * *q++;
					}
					p = q--;
					*p = (z * *q + *p) * z / (z * z - 1.0);
					for (j = SignalLength - 2L; (0L <= j); j--) {
						*q-- = z * (*p-- - *q);
					}
				}
				else {
					z1 = z;
					z2 = PositiveIntPower(z, SignalLength - 1L);
					iz = 1.0 / z;
					Sum = *OutputData + z2 * OutputData[SignalLength - 1L];
					z2 *= z2 * iz;
					p = OutputData;
					for (j = 2L - SignalLength; (j < 0L); j++) {
						Sum += (z1 + z2) * *++p;
						z1 *= z;
						z2 *= iz;
					}
					p = OutputData;
					q = OutputData;
					*p++ = Sum / (1.0 - z2 * z2);
					for (j = 1L - SignalLength; (j < 0L); j++) {
						*p++ += z * *q++;
					}
					p = q--;
					*p = (z * *q + *p) * z / (z * z - 1.0);
					for (j = SignalLength - 2L; (0L <= j); j--) {
						*q-- = z * (*p-- - *q);
					}
				}
			}
			break;
		case Periodic:
			if (SignalLength == 1L) {
				*OutputData = *InputData;
/**/			DEBUG_WRITE_LEAVING(IirConvolvePoles, "Done")
				return(Status);
			}
			for (i = 0L; (i < PoleNumber); i++) {
				z = RealPoles[i];
				if (Tolerance == 0.0) {
					Horizon = SignalLength;
				}
				else {
					Horizon = ConvertDoubleToLong(ceil(log(Tolerance) / log(fabs(z))));
				}
				if (Horizon < SignalLength) {
					z1 = z;
					Sum = OutputData[0];
					for (j = 1L; (j < Horizon); j++) {
						Sum += z1 * OutputData[SignalLength - j];
						z1 *= z;
					}
					OutputData[0] = Sum / (1.0 - PositiveIntPower(z, SignalLength));
					for (j = 1L; (j < SignalLength); j++) {
						OutputData[j] += z * OutputData[j - 1L];
					}
					z1 = z;
					Sum = OutputData[0] + OutputData[SignalLength - 1L] / z;
					for (j = 1L; (j < Horizon); j++) {
						Sum += z1 * OutputData[j];
						z1 *= z;
					}
					OutputData[SignalLength - 1L] = Sum * z * z
						/ (PositiveIntPower(z, SignalLength) - 1.0);
					for (j = SignalLength - 2L; (0L <= j); j--) {
						OutputData[j] = z * (OutputData[j + 1L] - OutputData[j]);
					}
				}
				else {
					z1 = PositiveIntPower(z, SignalLength - 1L);
					iz = 1.0 / z;
					Sum = OutputData[0];
					for (j = 1L; (j < SignalLength); j++) {
						Sum += z1 * OutputData[j];
						z1 *= iz;
					}
					OutputData[0] = Sum / (1.0 - PositiveIntPower(z, SignalLength));
					for (j = 1L; (j < SignalLength); j++) {
						OutputData[j] += z * OutputData[j - 1L];
					}
					z1 = z;
					Sum = OutputData[0] + OutputData[SignalLength - 1L] / z;
					for (j = 1L; (j < (SignalLength - 1L)); j++) {
						Sum += z1 * OutputData[j];
						z1 *= z;
					}
					OutputData[SignalLength - 1L] = Sum * z * z
						/ (PositiveIntPower(z, SignalLength) - 1.0);
					for (j = SignalLength - 2L; (0L <= j); j--) {
						OutputData[j] = z * (OutputData[j + 1L] - OutputData[j]);
					}
				}
			}
			break;
		default:
			Status = ERROR;
			WRITE_ERROR(IirConvolvePoles, "Invalid boundary convention")
			break;
	}
/**/DEBUG_WRITE_LEAVING(IirConvolvePoles, "Done")
	return(Status);
} /* end IirConvolvePoles */

/*--------------------------------------------------------------------------*/
extern int		IirConvolvePolesVolume
				(
					double	*VolumeSource,		/* data to process */
					double	*VolumeDestination,	/* result */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					double	RealPoles[],		/* array of real poles */
					long	PoleNumber,			/* number of poles */
					enum TBoundaryConvention
							Convention,			/* boundary convention */
					double	Tolerance,			/* admissible relative error */
					int		*Status				/* error management */
				)

/* no more than two poles are allowed for a finite support boundary */
/* VolumeSource is a (double)volume of size (Nx x Ny x Nz) */
/* OutputData is a (double)volume of size (Nx x Ny x Nz) */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin IirConvolvePolesVolume */

	double	*Buffer = (double *)NULL, *cBuffer = (double *)NULL;
	long	i, j, k;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(IirConvolvePolesVolume, VolumeSource, *Status,
/**/	"No VolumeSource")
/**/DEBUG_CHECK_NULL_POINTER(IirConvolvePolesVolume, VolumeDestination, *Status,
/**/	"No OutputData")
/**/DEBUG_CHECK_RANGE_LONG(IirConvolvePolesVolume, Nx, 1L, LONG_MAX, *Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(IirConvolvePolesVolume, Ny, 1L, LONG_MAX, *Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(IirConvolvePolesVolume, Nz, 1L, LONG_MAX, *Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_CHECK_NULL_POINTER(IirConvolvePolesVolume, RealPoles, *Status,
/**/	"No poles")
/**/DEBUG_CHECK_RANGE_LONG(IirConvolvePolesVolume, PoleNumber, 1L, LONG_MAX, *Status,
/**/	"Invalid number of poles (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_DOUBLE(IirConvolvePolesVolume, Tolerance, 0.0, DBL_MAX, *Status,
/**/	"Invalid Tolerance (should be positive)")
/**/DEBUG_RETURN_ON_ERROR(IirConvolvePolesVolume, *Status)
/**/DEBUG_WRITE_ENTERING(IirConvolvePolesVolume,
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
			*Status = ERROR;
			WRITE_ERROR(IirConvolvePolesVolume, "Invalid boundary convention")
/**/		DEBUG_WRITE_LEAVING(IirConvolvePolesVolume, "Done")
			return(*Status);
	}
	if (1L < Nx) {
		AllocateLineDouble(&Buffer, Nx, Status);
		if (*Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(IirConvolvePolesVolume, "Done")
			return(*Status);
		}
		AllocateLineDouble(&cBuffer, Nx, Status);
		if (*Status == ERROR) {
			FreeLineDouble(&Buffer);
/**/		DEBUG_WRITE_LEAVING(IirConvolvePolesVolume, "Done")
			return(*Status);
		}
		for (k = 0L; (k < Nz); k++) {
			for (j = 0L; (j < Ny); j++) {
				*Status = GetxDoubleToDouble(VolumeSource, Nx, Ny, Nz, 0L, j, k, Buffer, Nx);
				if (*Status == ERROR) {
					FreeLineDouble(&cBuffer);
					FreeLineDouble(&Buffer);
/**/				DEBUG_WRITE_LEAVING(IirConvolvePolesVolume, "Done")
					return(*Status);
				}
				*Status = IirConvolvePoles(Buffer, cBuffer, Nx,
					RealPoles, PoleNumber, Convention, Tolerance);
				if (*Status == ERROR) {
					FreeLineDouble(&cBuffer);
					FreeLineDouble(&Buffer);
/**/				DEBUG_WRITE_LEAVING(IirConvolvePolesVolume, "Done")
					return(*Status);
				}
				*Status = PutxDoubleToDouble(VolumeDestination, Nx, Ny, Nz, 0L, j, k, cBuffer, Nx);
				if (*Status == ERROR) {
					FreeLineDouble(&cBuffer);
					FreeLineDouble(&Buffer);
/**/				DEBUG_WRITE_LEAVING(IirConvolvePolesVolume, "Done")
					return(*Status);
				}
			}
		}
		*Status = FreeLineDouble(&cBuffer);
		if (*Status == ERROR) {
			FreeLineDouble(&Buffer);
/**/		DEBUG_WRITE_LEAVING(IirConvolvePolesVolume, "Done")
			return(*Status);
		}
		*Status = FreeLineDouble(&Buffer);
		if (*Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(IirConvolvePolesVolume, "Done")
			return(*Status);
		}
	}
	else {
		VolumeDestination = (double *)memcpy(VolumeDestination, VolumeSource,
			(size_t)(Nx * Ny * Nz * (long)sizeof(double)));
	}
	if (1L < Ny) {
		AllocateLineDouble(&Buffer, Ny, Status);
		if (*Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(IirConvolvePolesVolume, "Done")
			return(*Status);
		}
		AllocateLineDouble(&cBuffer, Ny, Status);
		if (*Status == ERROR) {
			FreeLineDouble(&Buffer);
/**/		DEBUG_WRITE_LEAVING(IirConvolvePolesVolume, "Done")
			return(*Status);
		}
		for (k = 0L; (k < Nz); k++) {
			for (i = 0L; (i < Nx); i++) {
				*Status = GetyDoubleToDouble(VolumeDestination, Nx, Ny, Nz, i, 0L, k, Buffer, Ny);
				if (*Status == ERROR) {
					FreeLineDouble(&cBuffer);
					FreeLineDouble(&Buffer);
/**/				DEBUG_WRITE_LEAVING(IirConvolvePolesVolume, "Done")
					return(*Status);
				}
				*Status = IirConvolvePoles(Buffer, cBuffer, Ny,
					RealPoles, PoleNumber, Convention, Tolerance);
				if (*Status == ERROR) {
					FreeLineDouble(&cBuffer);
					FreeLineDouble(&Buffer);
/**/				DEBUG_WRITE_LEAVING(IirConvolvePolesVolume, "Done")
					return(*Status);
				}
				*Status = PutyDoubleToDouble(VolumeDestination, Nx, Ny, Nz, i, 0L, k, cBuffer, Ny);
				if (*Status == ERROR) {
					FreeLineDouble(&cBuffer);
					FreeLineDouble(&Buffer);
/**/				DEBUG_WRITE_LEAVING(IirConvolvePolesVolume, "Done")
					return(*Status);
				}
			}
		}
		*Status = FreeLineDouble(&cBuffer);
		if (*Status == ERROR) {
			FreeLineDouble(&Buffer);
/**/		DEBUG_WRITE_LEAVING(IirConvolvePolesVolume, "Done")
			return(*Status);
		}
		*Status = FreeLineDouble(&Buffer);
		if (*Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(IirConvolvePolesVolume, "Done")
			return(*Status);
		}
	}
	if (1L < Nz) {
		AllocateLineDouble(&Buffer, Nz, Status);
		if (*Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(IirConvolvePolesVolume, "Done")
			return(*Status);
		}
		AllocateLineDouble(&cBuffer, Nz, Status);
		if (*Status == ERROR) {
			FreeLineDouble(&Buffer);
/**/		DEBUG_WRITE_LEAVING(IirConvolvePolesVolume, "Done")
			return(*Status);
		}
		for (j = 0L; (j < Ny); j++) {
			for (i = 0L; (i < Nx); i++) {
				*Status = GetzDoubleToDouble(VolumeDestination, Nx, Ny, Nz, i, j, 0L, Buffer, Nz);
				if (*Status == ERROR) {
					FreeLineDouble(&cBuffer);
					FreeLineDouble(&Buffer);
/**/				DEBUG_WRITE_LEAVING(IirConvolvePolesVolume, "Done")
					return(*Status);
				}
				*Status = IirConvolvePoles(Buffer, cBuffer, Nz,
					RealPoles, PoleNumber, Convention, Tolerance);
				if (*Status == ERROR) {
					FreeLineDouble(&cBuffer);
					FreeLineDouble(&Buffer);
/**/				DEBUG_WRITE_LEAVING(IirConvolvePolesVolume, "Done")
					return(*Status);
				}
				*Status = PutzDoubleToDouble(VolumeDestination, Nx, Ny, Nz, i, j, 0L, cBuffer, Nz);
				if (*Status == ERROR) {
					FreeLineDouble(&cBuffer);
					FreeLineDouble(&Buffer);
/**/				DEBUG_WRITE_LEAVING(IirConvolvePolesVolume, "Done")
					return(*Status);
				}
			}
		}
		*Status = FreeLineDouble(&cBuffer);
		if (*Status == ERROR) {
			FreeLineDouble(&Buffer);
/**/		DEBUG_WRITE_LEAVING(IirConvolvePolesVolume, "Done")
			return(*Status);
		}
		*Status = FreeLineDouble(&Buffer);
		if (*Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(IirConvolvePolesVolume, "Done")
			return(*Status);
		}
	}
/**/DEBUG_WRITE_LEAVING(IirConvolvePolesVolume, "Done")
	return(*Status);
} /* end IirConvolvePolesVolume */

