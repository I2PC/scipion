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
#include	"firconvolve.h"
#include	"getput.h"
#include	"messagedisplay.h"

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
extern int		FirConvolve
				(
					double	InputData[],		/* data to process */
					double	OutputData[],		/* result */
					long	SignalLength,		/* length of the 1D data array */
					double	Kernel[],			/* kernel */
					long	KernelOrigin,		/* center of the kernel */
					long	KernelLength,		/* length of the 1D kernel */
					enum TBoundaryConvention
							Convention			/* boundary convention */
				)

/* the specified boundary convention applies to the input data only, not to the kernel */
/* the boundary convention applied to the kernel is FiniteDataSupport */
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

{ /* begin FirConvolve */

	double	*p, *q;
	double	f0, fk, fn, df;
	double	Sum;
	long	i, j, k;
	long	m, n;
	long	n2;
	long	kp, km;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(FirConvolve, InputData, Status,
/**/	"No input data")
/**/DEBUG_CHECK_NULL_POINTER(FirConvolve, OutputData, Status,
/**/	"No output data")
/**/DEBUG_CHECK_RANGE_LONG(FirConvolve, SignalLength, 1L, LONG_MAX, Status,
/**/	"Invalid signal length (should be strictly positive)")
/**/DEBUG_CHECK_NULL_POINTER(FirConvolve, Kernel, Status,
/**/	"No kernel")
/**/DEBUG_CHECK_RANGE_LONG(FirConvolve, KernelLength, 1L, LONG_MAX, Status,
/**/	"Invalid kernel length (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(FirConvolve, Status)
/**/DEBUG_WRITE_ENTERING(FirConvolve,
/**/	"About to perform FIR convolution")

	switch (Convention) {
		case AntiMirrorOnBounds:
/* More optimization needed for better efficiency in AntiMirrorOnBounds */
			Kernel += (ptrdiff_t)(KernelLength - 1L);
			if (SignalLength == 1L) {
				Sum = 0.0;
				for (j = -KernelLength; (j < 0L); j++) {
					Sum += *Kernel--;
				}
				*OutputData = *InputData * Sum;
			}
			else {
				n2 = 2L * (SignalLength - 1L);
				f0 = *InputData;
				fn = InputData[SignalLength - 1L];
				df = fn - f0;
				m = 1L + KernelOrigin - KernelLength;
				km = m;
				m -= (m < 0L) ? (n2 * ((m + 1L - n2) / n2)) : (n2 * (m / n2));
				for (i = 0L; (i < SignalLength); km++, i++) {
					q = Kernel;
					k = m;
					kp = km;
					Sum = 0.0;
					for (j = -KernelLength; (j < 0L); kp++, j++) {
						fk = (k < SignalLength) ? (InputData[k]) : (2.0 * fn - InputData[n2 - k]);
						Sum += ((double)((kp - k) / (SignalLength - 1L)) * df + fk) * *q--;
						if (++k == n2) {
							k = 0L;
						}
					}
					if (++m == n2) {
						m = 0L;
					}
					*OutputData++ = Sum;
				}
			}
			break;
		case FiniteDataSupport:
			Kernel += (ptrdiff_t)(KernelLength - 1L);
			k = KernelOrigin - (KernelLength - 1L);
			for (i = -SignalLength; (i < 0L); k++, i++) {
				kp = (0L < k) ? (k) : (0L);
				km = k - kp;
				p = InputData + (ptrdiff_t)kp;
				q = Kernel + (ptrdiff_t)km;
				Sum = 0.0;
				for (j = kp - ((SignalLength < (k + KernelLength)) ? (SignalLength)
					: (k + KernelLength)); (j < 0L); j++) {
					Sum += *p++ * *q--;
				}
				*OutputData++ = Sum;
			}
			break;
		case MirrorOffBounds:
			Kernel += (ptrdiff_t)(KernelLength - 1L);
			if (SignalLength == 1L) {
				Sum = 0.0;
				for (j = -KernelLength; (j < 0L); j++) {
					Sum += *Kernel--;
				}
				*OutputData = *InputData * Sum;
			}
			else {
				n2 = 2L * SignalLength;
				m = 1L + KernelOrigin - KernelLength;
				m -= (m < 0L) ? (n2 * ((m + 1L - n2) / n2)) : (n2 * (m / n2));
				for (i = 0L; (i < SignalLength); i++) {
					j = -KernelLength;
					k = m;
					q = Kernel;
					Sum = 0.0;
					while (j < 0L) {
						p = InputData + (ptrdiff_t)k;
						kp = ((k - SignalLength) < j) ? (j) : (k - SignalLength);
						if (kp < 0L) {
							for (n = kp; (n < 0L); n++) {
								Sum += *p++ * *q--;
							}
							k -= kp;
							j -= kp;
						}
						p = InputData + (ptrdiff_t)(n2 - k - 1L);
						km = ((k - n2) < j) ? (j) : (k - n2);
						if (km < 0L) {
							for (n = km; (n < 0L); n++) {
								Sum += *p-- * *q--;
							}
							j -= km;
						}
						k = 0L;
					}
					if (++m == n2) {
						m = 0L;
					}
					*OutputData++ = Sum;
				}
			}
			break;
		case MirrorOnBounds:
			Kernel += (ptrdiff_t)(KernelLength - 1L);
			if (SignalLength == 1L) {
				Sum = 0.0;
				for (j = -KernelLength; (j < 0L); j++) {
					Sum += *Kernel--;
				}
				*OutputData = *InputData * Sum;
			}
			else {
				n2 = 2L * (SignalLength - 1L);
				m = 1L + KernelOrigin - KernelLength;
				m -= (m < 0L) ? (n2 * ((m + 1L - n2) / n2)) : (n2 * (m / n2));
				for (i = 0L; (i < SignalLength); i++) {
					j = -KernelLength;
					k = m;
					q = Kernel;
					Sum = 0.0;
					while (j < 0L) {
						p = InputData + (ptrdiff_t)k;
						kp = ((k - SignalLength) < j) ? (j) : (k - SignalLength);
						if (kp < 0L) {
							for (n = kp; (n < 0L); n++) {
								Sum += *p++ * *q--;
							}
							k -= kp;
							j -= kp;
						}
						p = InputData + (ptrdiff_t)(n2 - k);
						km = ((k - n2) < j) ? (j) : (k - n2);
						if (km < 0L) {
							for (n = km; (n < 0L); n++) {
								Sum += *p-- * *q--;
							}
							j -= km;
						}
						k = 0L;
					}
					if (++m == n2) {
						m = 0L;
					}
					*OutputData++ = Sum;
				}
			}
			break;
		case Periodic:
			Kernel += (ptrdiff_t)(KernelLength - 1L);
			k = KernelOrigin - (KernelLength - 1L);
			k = (k < 0L) ? (k + SignalLength * ((SignalLength - 1L - k) / SignalLength))
				: ((SignalLength <= k) ? (k - SignalLength * (k / SignalLength)) : (k));
			for (i = -SignalLength; (i < 0L); i++) {
				p = InputData + (ptrdiff_t)k;
				q = Kernel;
				Sum = 0.0;
				j = k - SignalLength;
				if ((j + KernelLength) <= 0L) {
					for (j = -KernelLength; (j < 0L); j++) {
						Sum += *p++ * *q--;
					}
				}
				else {
					while (j++ < 0L) {
						Sum += *p++ * *q--;
					}
					p = InputData;
					m = 1L - (KernelLength + k) / SignalLength;
					for (n = m; (n < 0L); n++) {
						for (j = -SignalLength; (j < 0L); j++) {
							Sum += *p++ * *q--;
						}
						p = InputData;
					}
					for (j = (1L - m) * SignalLength - KernelLength - k; (j < 0L); j++) {
						Sum += *p++ * *q--;
					}
				}
				if (++k == SignalLength) {
					k = 0L;
				}
				*OutputData++ = Sum;
			}
			break;
		default:
			Status = ERROR;
			WRITE_ERROR(FirConvolve, "Invalid boundary convention")
			break;
	}
/**/DEBUG_WRITE_LEAVING(FirConvolve, "Done")
	return(Status);
} /* end FirConvolve */

/*--------------------------------------------------------------------------*/
extern int		FirConvolveAntiSymmetric
				(
					double	InputData[],		/* data to process */
					double	OutputData[],		/* result */
					long	SignalLength,		/* length of the 1D data array */
					double	HalfKernel[],		/* causal part of the kernel */
					long	KernelHalfLength,	/* length of the causal part of the kernel */
					enum TBoundaryConvention
							Convention			/* boundary convention */
				)

/* the specified boundary convention applies to the input data only, not to the kernel */
/* the boundary convention applied to the kernel is FiniteDataSupport */
/* the input and the output have the same length */
/* the origin for the kernel is its leftmost sample; it corresponds to its anti-symmetry axis */
/* the value of the kernel at the origin is expected to be HalfKernel[0] = 0.0 */
/* the full length of the symmetric kernel is (2L * KernelHalfLength - 1L) */
/* success: return(!ERROR); failure: return(ERROR) */
/* general structure is as follows:
	for (i = 0L; (i < SignalLength); i++) {
		Sum = 0.0;
		for (j = 1L; (j < KernelHalfLength); j++) {
			Sum += (InputData[i + j] - InputData[i - j]) * HalfKernel[j];
		}
		OutputData[i] = Sum;
	}
*/

{ /* begin FirConvolveAntiSymmetric */

	double	*q, *u, *v;
	double	*u1, *u2, *v1, *v2;
	double	Sum;
	double	c1, c2;
	long	i, j;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(FirConvolveAntiSymmetric, InputData, Status,
/**/	"No input data")
/**/DEBUG_CHECK_NULL_POINTER(FirConvolveAntiSymmetric, OutputData, Status,
/**/	"No output data")
/**/DEBUG_CHECK_RANGE_LONG(FirConvolveAntiSymmetric, SignalLength, 1L, LONG_MAX, Status,
/**/	"Invalid signal length (should be strictly positive)")
/**/DEBUG_CHECK_NULL_POINTER(FirConvolveAntiSymmetric, HalfKernel, Status,
/**/	"No kernel")
/**/DEBUG_CHECK_RANGE_LONG(FirConvolveAntiSymmetric, KernelHalfLength, 1L, LONG_MAX, Status,
/**/	"Invalid kernel length (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(FirConvolveAntiSymmetric, Status)
/**/DEBUG_WRITE_ENTERING(FirConvolveAntiSymmetric,
/**/	"About to perform FIR symmetric convolution")

	if (*HalfKernel != 0.0) {
		Status = ERROR;
		WRITE_ERROR(FirConvolveAntiSymmetric,
			"Invalid anti-symmetric kernel (should be 0.0 at the origin)")
/**/	DEBUG_WRITE_LEAVING(FirConvolveAntiSymmetric, "Done")
		return(Status);
	}
	switch (Convention) {
		case AntiMirrorOnBounds:
			Status = ERROR;
			WRITE_ERROR(FirConvolveAntiSymmetric, "Not yet implemented")
			break;
		case FiniteDataSupport:
			switch (KernelHalfLength) {
				case 1L:
					for (i = -SignalLength; (i < 0L); i++) {
						*OutputData++ = 0.0;
					}
					break;
				case 2L:
					if (3L <= SignalLength) {
						c1 = *++HalfKernel;
						u = InputData;
						*OutputData++ = -*++u * c1;
						for (i = 2L - SignalLength; (i < 0L); i++) {
							v = InputData;
							InputData = u++;
							*OutputData++ = (*v - *u) * c1;
						}
						*OutputData = *InputData * c1;
					}
					else {
						switch (SignalLength) {
							case 2L:
								c1 = *++HalfKernel;
								v = InputData++;
								*OutputData++ = -*InputData * c1;
								*OutputData = *v * c1;
								break;
							case 1L:
								*OutputData = 0.0;
								break;
							default:
								Status = ERROR;
								WRITE_ERROR(FirConvolveAntiSymmetric,
									"Unexpected signal length")
								break;
						}
					}
					break;
				case 3L:
					if (5L <= SignalLength) {
						c1 = *++HalfKernel;
						c2 = *++HalfKernel;
						u2 = InputData;
						u2++;
						u1 = u2++;
						*OutputData++ = -*u1 * c1 - *u2 * c2;
						v1 = InputData;
						InputData = u1;
						u1 = u2++;
						*OutputData++ = (*v1 - *u1) * c1 - *u2 * c2;
						for (i = 4L - SignalLength; (i < 0L); i++) {
							v2 = v1;
							v1 = InputData;
							InputData = u1;
							u1 = u2++;
							*OutputData++ = (*v1 - *u1) * c1 + (*v2 - *u2) * c2;
						}
						*OutputData++ = (*InputData - *u2) * c1 + *v1 * c2;
						*OutputData = *u1 * c1 + *InputData * c2;
					}
					else {
						switch (SignalLength) {
							case 4L:
								c1 = *++HalfKernel;
								c2 = *++HalfKernel;
								u2 = InputData;
								u2++;
								u1 = u2++;
								*OutputData++ = -*u1 * c1 - *u2 * c2;
								v = u2++;
								*OutputData++ = (*InputData - *v) * c1 - *u2 * c2;
								*OutputData++ = (*u1 - *u2) * c1 + *InputData * c2;
								*OutputData = *v * c1 + *u1 * c2;
								break;
							case 3L:
								c1 = *++HalfKernel;
								c2 = *++HalfKernel;
								u2 = InputData;
								u2++;
								u1 = u2++;
								*OutputData++ = -*u1 * c1 - *u2 * c2;
								*OutputData++ = (*InputData - *u2) * c1;
								*OutputData = *u1 * c1 + *InputData * c2;
								break;
							case 2L:
								c1 = *++HalfKernel;
								v = InputData++;
								*OutputData++ = -*InputData * c1;
								*OutputData = *v * c1;
								break;
							case 1L:
								*OutputData = 0.0;
								break;
							default:
								Status = ERROR;
								WRITE_ERROR(FirConvolveAntiSymmetric,
									"Unexpected signal length")
								break;
						}
					}
					break;
				default:
					if ((2L * KernelHalfLength - 1L) <= SignalLength) {
						for (i = 1L - KernelHalfLength; (i < 0L); i++) {
							q = HalfKernel;
							u = v = InputData;
							InputData++;
							Sum = 0.0;
							for (j = 1L - KernelHalfLength - i; (j < 0L); j++) {
								Sum -= (*++u - *--v) * *++q;
							}
							for (j = i; (j < 0L); j++) {
								Sum -= *++u * *++q;
							}
							*OutputData++ = Sum;
						}
						for (i = 2L * KernelHalfLength - SignalLength - 2L; (i < 0L); i++) {
							q = HalfKernel;
							u = v = InputData;
							InputData++;
							Sum = 0.0;
							for (j = 1L - KernelHalfLength; (j < 0L); j++) {
								Sum -= (*++u - *--v) * *++q;
							}
							*OutputData++ = Sum;
						}
						for (i = 1L - KernelHalfLength; (i < 0L); i++) {
							q = HalfKernel;
							u = v = InputData;
							InputData++;
							Sum = 0.0;
							for (j = i + 1L; (j < 0L); j++) {
								Sum -= (*++u - *--v) * *++q;
							}
							for (j = -KernelHalfLength - i; (j < 0L); j++) {
								Sum += *--v * *++q;
							}
							*OutputData++ = Sum;
						}
					}
					else if (KernelHalfLength <= SignalLength) {
						for (i = KernelHalfLength - SignalLength - 1L; (i < 0L); i++) {
							q = HalfKernel;
							u = v = InputData;
							InputData++;
							Sum = 0.0;
							for (j = KernelHalfLength - SignalLength - i - 1L; (j < 0L); j++) {
								Sum -= (*++u - *--v) * *++q;
							}
							for (j = SignalLength - 2L * KernelHalfLength + i + 2L; (j < 0L); j++) {
								Sum -= *++u * *++q;
							}
							*OutputData++ = Sum;
						}
						for (i = SignalLength - (SignalLength / 2L) - KernelHalfLength + 1L; (i < 0L); i++) {
							q = HalfKernel;
							u = v = InputData;
							InputData++;
							Sum = 0.0;
							for (j = -(SignalLength / 2L) - i; (j < 0L); j++) {
								Sum -= (*++u - *--v) * *++q;
							}
							for (j = 2L * (SignalLength / 2L) - SignalLength + 2L * i + 1L; (j < 0L); j++) {
								Sum -= *++u * *++q;
							}
							*OutputData++ = Sum;
						}
						if ((2L * (SignalLength / 2L) + 1L) == SignalLength) {
							q = HalfKernel;
							u = v = InputData;
							InputData++;
							Sum = 0.0;
							for (j = -(SignalLength / 2L); (j < 0L); j++) {
								Sum -= (*++u - *--v) * *++q;
							}
							*OutputData++ = Sum;
							i = 1L;
						}
						for (i += SignalLength / 2L - KernelHalfLength + 1L; (i < 0L); i++) {
							q = HalfKernel;
							u = v = InputData;
							InputData++;
							Sum = 0.0;
							for (j = KernelHalfLength - SignalLength + i; (j < 0L); j++) {
								Sum -= (*++u - *--v) * *++q;
							}
							for (j = SignalLength - 2L * KernelHalfLength - 2L * i + 1L; (j < 0L); j++) {
								Sum += *--v * *++q;
							}
							*OutputData++ = Sum;
						}
						for (i = KernelHalfLength - SignalLength - 1L; (i < 0L); i++) {
							q = HalfKernel;
							u = v = InputData;
							InputData++;
							Sum = 0.0;
							for (j = i + 1L; (j < 0L); j++) {
								Sum -= (*++u - *--v) * *++q;
							}
							for (j = -KernelHalfLength - i; (j < 0L); j++) {
								Sum += *--v * *++q;
							}
							*OutputData++ = Sum;
						}
					}
					else {
						for (i = -(SignalLength / 2L); (i < 0L); i++) {
							q = HalfKernel;
							u = v = InputData;
							InputData++;
							Sum = 0.0;
							for (j = -(SignalLength / 2L) - i; (j < 0L); j++) {
								Sum -= (*++u - *--v) * *++q;
							}
							for (j = 2L * (SignalLength / 2L) - SignalLength + 2L * i + 1L; (j < 0L); j++) {
								Sum -= *++u * *++q;
							}
							*OutputData++ = Sum;
						}
						if ((2L * (SignalLength / 2L) + 1L) == SignalLength) {
							q = HalfKernel;
							u = v = InputData;
							InputData++;
							Sum = 0.0;
							for (j = -(SignalLength / 2L) - i; (j < 0L); j++) {
								Sum -= (*++u - *--v) * *++q;
							}
							*OutputData++ = Sum;
							i = 1L;
						}
						for (i += SignalLength / 2L - SignalLength; (i < 0L); i++) {
							q = HalfKernel;
							u = v = InputData;
							InputData++;
							Sum = 0.0;
							for (j = i + 1L; (j < 0L); j++) {
								Sum -= (*++u - *--v) * *++q;
							}
							for (j = -SignalLength - 2L * i - 1L; (j < 0L); j++) {
								Sum += *--v * *++q;
							}
							*OutputData++ = Sum;
						}
					}
					break;
			}
			break;
		case MirrorOffBounds:
			Status = ERROR;
			WRITE_ERROR(FirConvolveAntiSymmetric, "Not yet implemented")
			break;
		case MirrorOnBounds:
			switch (KernelHalfLength) {
				case 1L:
					for (i = -SignalLength; (i < 0L); i++) {
						*OutputData++ = 0.0;
					}
					break;
				case 2L:
					if (3L <= SignalLength) {
						c1 = *++HalfKernel;
						u = InputData;
						u++;
						*OutputData++ = 0.0;
						for (i = 2L - SignalLength; (i < 0L); i++) {
							v = InputData;
							InputData = u++;
							*OutputData++ = (*v - *u) * c1;
						}
						*OutputData = 0.0;
					}
					else {
						switch (SignalLength) {
							case 2L:
								*OutputData++ = 0.0;
								*OutputData = 0.0;
								break;
							case 1L:
								*OutputData = 0.0;
								break;
							default:
								Status = ERROR;
								WRITE_ERROR(FirConvolveAntiSymmetric,
									"Unexpected signal length")
								break;
						}
					}
					break;
				case 3L:
					if (5L <= SignalLength) {
						c1 = *++HalfKernel;
						c2 = *++HalfKernel;
						u2 = InputData;
						u2++;
						u1 = u2++;
						*OutputData++ = 0.0;
						v1 = InputData;
						InputData = u1;
						u1 = u2++;
						*OutputData++ = (*v1 - *u1) * c1 + (*InputData - *u2) * c2;
						for (i = 4L - SignalLength; (i < 0L); i++) {
							v2 = v1;
							v1 = InputData;
							InputData = u1;
							u1 = u2++;
							*OutputData++ = (*v1 - *u1) * c1 + (*v2 - *u2) * c2;
						}
						*OutputData++ = (*InputData - *u2) * c1 + (*v1 - *u1) * c2;
						*OutputData = 0.0;
					}
					else {
						switch (SignalLength) {
							case 4L:
								c1 = *++HalfKernel;
								c2 = *++HalfKernel;
								u2 = InputData;
								u2++;
								u1 = u2++;
								v = u2++;
								*OutputData++ = 0.0;
								*OutputData++ = (*InputData - *v) * c1 + (*u1 - *u2) * c2;
								*OutputData++ = (*u1 - *u2) * c1 + (*InputData - *v) * c2;
								*OutputData = 0.0;
								break;
							case 3L:
								c1 = *++HalfKernel;
								u2 = InputData + (ptrdiff_t)2;
								*OutputData++ = 0.0;
								*OutputData++ = (*InputData - *u2) * c1;
								*OutputData = 0.0;
								break;
							case 2L:
								*OutputData++ = 0.0;
								*OutputData = 0.0;
								break;
							case 1L:
								*OutputData = 0.0;
								break;
							default:
								Status = ERROR;
								WRITE_ERROR(FirConvolveAntiSymmetric,
									"Unexpected signal length")
								break;
						}
					}
					break;
				default:
					if ((2L * KernelHalfLength - 1L) <= SignalLength) {
						for (i = 1L - KernelHalfLength; (i < 0L); i++) {
							q = HalfKernel;
							u = v = InputData;
							InputData++;
							Sum = 0.0;
							for (j = 1L - KernelHalfLength - i; (j < 0L); j++) {
								Sum -= (*++u - *--v) * *++q;
							}
							for (j = i; (j < 0L); j++) {
								Sum -= (*++u - *++v) * *++q;
							}
							*OutputData++ = Sum;
						}
						for (i = 2L * KernelHalfLength - SignalLength - 2L; (i < 0L); i++) {
							q = HalfKernel;
							u = v = InputData;
							InputData++;
							Sum = 0.0;
							for (j = 1L - KernelHalfLength; (j < 0L); j++) {
								Sum -= (*++u - *--v) * *++q;
							}
							*OutputData++ = Sum;
						}
						for (i = 1L - KernelHalfLength; (i < 0L); i++) {
							q = HalfKernel;
							u = v = InputData;
							InputData++;
							Sum = 0.0;
							for (j = i + 1L; (j < 0L); j++) {
								Sum -= (*++u - *--v) * *++q;
							}
							for (j = -KernelHalfLength - i; (j < 0L); j++) {
								Sum -= (*--u - *--v) * *++q;
							}
							*OutputData++ = Sum;
						}
					}
					else if (KernelHalfLength <= SignalLength) {
						for (i = KernelHalfLength - SignalLength - 1L; (i < 0L); i++) {
							q = HalfKernel;
							u = v = InputData;
							InputData++;
							Sum = 0.0;
							for (j = KernelHalfLength - SignalLength - i - 1L; (j < 0L); j++) {
								Sum -= (*++u - *--v) * *++q;
							}
							for (j = SignalLength - 2L * KernelHalfLength + i + 2L; (j < 0L); j++) {
								Sum -= (*++u - *++v) * *++q;
							}
							*OutputData++ = Sum;
						}
						for (i = SignalLength - (SignalLength / 2L) - KernelHalfLength + 1L; (i < 0L); i++) {
							q = HalfKernel;
							u = v = InputData;
							InputData++;
							Sum = 0.0;
							for (j = -(SignalLength / 2L) - i; (j < 0L); j++) {
								Sum -= (*++u - *--v) * *++q;
							}
							for (j = 2L * (SignalLength / 2L) - SignalLength + 2L * i + 1L; (j < 0L); j++) {
								Sum -= (*++u - *++v) * *++q;
							}
							for (j = SignalLength - (SignalLength / 2L) - KernelHalfLength - i; (j < 0L); j++) {
								Sum -= (*--u - *++v) * *++q;
							}
							*OutputData++ = Sum;
						}
						if ((2L * (SignalLength / 2L) + 1L) == SignalLength) {
							q = HalfKernel;
							u = v = InputData;
							InputData++;
							Sum = 0.0;
							for (j = -(SignalLength / 2L); (j < 0L); j++) {
								Sum -= (*++u - *--v) * *++q;
							}
							for (j = SignalLength / 2L - KernelHalfLength + 1L; (j < 0L); j++) {
								Sum -= (*--u - *++v) * *++q;
							}
							*OutputData++ = Sum;
							i = 1L;
						}
						for (i += SignalLength / 2L - KernelHalfLength + 1L; (i < 0L); i++) {
							q = HalfKernel;
							u = v = InputData;
							InputData++;
							Sum = 0.0;
							for (j = KernelHalfLength - SignalLength + i; (j < 0L); j++) {
								Sum -= (*++u - *--v) * *++q;
							}
							for (j = SignalLength - 2L * KernelHalfLength - 2L * i + 1L; (j < 0L); j++) {
								Sum -= (*--u - *--v) * *++q;
							}
							for (j = i; (j < 0L); j++) {
								Sum -= (*--u - *++v) * *++q;
							}
							*OutputData++ = Sum;
						}
						for (i = KernelHalfLength - SignalLength - 1L; (i < 0L); i++) {
							q = HalfKernel;
							u = v = InputData;
							InputData++;
							Sum = 0.0;
							for (j = i + 1L; (j < 0L); j++) {
								Sum -= (*++u - *--v) * *++q;
							}
							for (j = -KernelHalfLength - i; (j < 0L); j++) {
								Sum -= (*--u - *--v) * *++q;
							}
							*OutputData++ = Sum;
						}
					}
					else {
						if (SignalLength == 1L) {
							*OutputData = 0.0;
						}
						else {
							v1 = InputData;
							u1 = InputData + (ptrdiff_t)(SignalLength - 1L);
							for (i = -SignalLength; (i < 0L); i++) {
								q = HalfKernel;
								u = v = InputData;
								InputData++;
								Sum = 0.0;
								j = 1L - KernelHalfLength;
								while ((j < 0L) && (u < u1) && (v1 < v)) {
									Sum -= (*++u - *--v) * *++q;
									j++;
								}
								if ((u == u1) && (v == v1)) {
									while (j < 0L) {
										while ((j < 0L) && (v < u1)) {
											Sum -= (*--u - *++v) * *++q;
											j++;
										}
										while ((j < 0L) && (u < u1)) {
											Sum -= (*++u - *--v) * *++q;
											j++;
										}
									}
								}
								else if (u == u1) {
									while (j < 0L) {
										while ((j < 0L) && (v1 < v)) {
											Sum -= (*--u - *--v) * *++q;
											j++;
										}
										while ((j < 0L) && (v1 < u)) {
											Sum -= (*--u - *++v) * *++q;
											j++;
										}
										while ((j < 0L) && (v < u1)) {
											Sum -= (*++u - *++v) * *++q;
											j++;
										}
										while ((j < 0L) && (u < u1)) {
											Sum -= (*++u - *--v) * *++q;
											j++;
										}
									}
								}
								else {
									while (j < 0L) {
										while ((j < 0L) && (u < u1)) {
											Sum -= (*++u - *++v) * *++q;
											j++;
										}
										while ((j < 0L) && (v < u1)) {
											Sum -= (*--u - *++v) * *++q;
											j++;
										}
										while ((j < 0L) && (v1 < u)) {
											Sum -= (*--u - *--v) * *++q;
											j++;
										}
										while ((j < 0L) && (v1 < v)) {
											Sum -= (*++u - *--v) * *++q;
											j++;
										}
									}
								}
								*OutputData++ = Sum;
							}
						}
					}
					break;
			}
			break;
		case Periodic:
			switch (KernelHalfLength) {
				case 1L:
					for (i = -SignalLength; (i < 0L); i++) {
						*OutputData++ = 0.0;
					}
					break;
				case 2L:
					if (3L <= SignalLength) {
						c1 = *++HalfKernel;
						u = u1 = InputData;
						v = InputData + (ptrdiff_t)SignalLength;
						*OutputData++ = (*--v - *++u) * c1;
						for (i = 2L - SignalLength; (i < 0L); i++) {
							v = InputData;
							InputData = u++;
							*OutputData++ = (*v - *u) * c1;
						}
						*OutputData = (*InputData - *u1) * c1;
					}
					else {
						switch (SignalLength) {
							case 2L:
								*OutputData++ = 0.0;
								*OutputData = 0.0;
								break;
							case 1L:
								*OutputData = 0.0;
								break;
							default:
								Status = ERROR;
								WRITE_ERROR(FirConvolveAntiSymmetric,
									"Unexpected signal length")
								break;
						}
					}
					break;
				case 3L:
					if (5L <= SignalLength) {
						c1 = *++HalfKernel;
						c2 = *++HalfKernel;
						u = u2 = InputData;
						u2++;
						u1 = u2++;
						v = u++;
						v1 = InputData + (ptrdiff_t)SignalLength;
						v2 = v1--;
						*OutputData++ = (*--v2 - *u1) * c1 + (*--v1 - *u2) * c2;
						v1 = InputData;
						InputData = u1;
						u1 = u2++;
						*OutputData++ = (*v1 - *u1) * c1 + (*v2 - *u2) * c2;
						for (i = 4L - SignalLength; (i < 0L); i++) {
							v2 = v1;
							v1 = InputData;
							InputData = u1;
							u1 = u2++;
							*OutputData++ = (*v1 - *u1) * c1 + (*v2 - *u2) * c2;
						}
						*OutputData++ = (*InputData - *u2) * c1 + (*v1 - *v) * c2;
						*OutputData = (*u1 - *v) * c1 + (*InputData - *u) * c2;
					}
					else {
						switch (SignalLength) {
							case 4L:
								c1 = *++HalfKernel;
								c2 = *++HalfKernel;
								u2 = InputData;
								u2++;
								u1 = u2++;
								v = u2++;
								c2 = (*u1 - *u2) * c1;
								c1 *= *v - *InputData;
								*OutputData++ = -c2;
								*OutputData++ = -c1;
								*OutputData++ = c2;
								*OutputData = c1;
								break;
							case 3L:
								c1 = *++HalfKernel;
								c2 = *++HalfKernel;
								u2 = InputData;
								u2++;
								u1 = u2++;
								c1 -= c2;
								*OutputData++ = (*u2 - *u1) * c1;
								*OutputData++ = (*InputData - *u2) * c1;
								*OutputData = (*u1 - *InputData) * c1;
								break;
							case 2L:
								*OutputData++ = 0.0;
								*OutputData = 0.0;
								break;
							case 1L:
								*OutputData = 0.0;
								break;
							default:
								Status = ERROR;
								WRITE_ERROR(FirConvolveAntiSymmetric,
									"Unexpected signal length")
								break;
						}
					}
					break;
				default:
					v2 = InputData;
					v1 = v2--;
					u1 = InputData + (ptrdiff_t)SignalLength;
					u2 = u1--;
					if ((2L * KernelHalfLength - 1L) <= SignalLength) {
						for (i = 1L - KernelHalfLength; (i < 0L); i++) {
							q = HalfKernel;
							u = v = InputData;
							InputData++;
							Sum = 0.0;
							for (j = 1L - KernelHalfLength - i; (j < 0L); j++) {
								Sum -= (*++u - *--v) * *++q;
							}
							v = u2;
							for (j = i; (j < 0L); j++) {
								Sum -= (*++u - *--v) * *++q;
							}
							*OutputData++ = Sum;
						}
						for (i = 2L * KernelHalfLength - SignalLength - 2L; (i < 0L); i++) {
							q = HalfKernel;
							u = v = InputData;
							InputData++;
							Sum = 0.0;
							for (j = 1L - KernelHalfLength; (j < 0L); j++) {
								Sum -= (*++u - *--v) * *++q;
							}
							*OutputData++ = Sum;
						}
						for (i = 1L - KernelHalfLength; (i < 0L); i++) {
							q = HalfKernel;
							u = v = InputData;
							InputData++;
							Sum = 0.0;
							for (j = i + 1L; (j < 0L); j++) {
								Sum -= (*++u - *--v) * *++q;
							}
							u = v2;
							for (j = -KernelHalfLength - i; (j < 0L); j++) {
								Sum -= (*++u - *--v) * *++q;
							}
							*OutputData++ = Sum;
						}
					}
					else if (KernelHalfLength <= SignalLength) {
						for (i = KernelHalfLength - SignalLength - 1L; (i < 0L); i++) {
							q = HalfKernel;
							u = v = InputData;
							InputData++;
							Sum = 0.0;
							for (j = KernelHalfLength - SignalLength - i - 1L; (j < 0L); j++) {
								Sum -= (*++u - *--v) * *++q;
							}
							v = u2;
							for (j = SignalLength - 2L * KernelHalfLength + i + 2L; (j < 0L); j++) {
								Sum -= (*++u - *--v) * *++q;
							}
							*OutputData++ = Sum;
						}
						for (i = SignalLength - (SignalLength / 2L) - KernelHalfLength + 1L; (i < 0L); i++) {
							q = HalfKernel;
							u = v = InputData;
							InputData++;
							Sum = 0.0;
							for (j = -(SignalLength / 2L) - i; (j < 0L); j++) {
								Sum -= (*++u - *--v) * *++q;
							}
							v = u2;
							for (j = 2L * (SignalLength / 2L) - SignalLength + 2L * i + 1L; (j < 0L); j++) {
								Sum -= (*++u - *--v) * *++q;
							}
							u = v2;
							for (j = SignalLength - (SignalLength / 2L) - KernelHalfLength - i; (j < 0L); j++) {
								Sum -= (*++u - *--v) * *++q;
							}
							*OutputData++ = Sum;
						}
						if ((2L * (SignalLength / 2L) + 1L) == SignalLength) {
							q = HalfKernel;
							u = v = InputData;
							InputData++;
							Sum = 0.0;
							for (j = -(SignalLength / 2L); (j < 0L); j++) {
								Sum -= (*++u - *--v) * *++q;
							}
							v = u2;
							u = v2;
							for (j = SignalLength / 2L - KernelHalfLength + 1L; (j < 0L); j++) {
								Sum -= (*++u - *--v) * *++q;
							}
							*OutputData++ = Sum;
							i = 1L;
						}
						for (i += SignalLength / 2L - KernelHalfLength + 1L; (i < 0L); i++) {
							q = HalfKernel;
							u = v = InputData;
							InputData++;
							Sum = 0.0;
							for (j = KernelHalfLength - SignalLength + i; (j < 0L); j++) {
								Sum -= (*++u - *--v) * *++q;
							}
							u = v2;
							for (j = SignalLength - 2L * KernelHalfLength - 2L * i + 1L; (j < 0L); j++) {
								Sum -= (*++u - *--v) * *++q;
							}
							v = u2;
							for (j = i; (j < 0L); j++) {
								Sum -= (*++u - *--v) * *++q;
							}
							*OutputData++ = Sum;
						}
						for (i = KernelHalfLength - SignalLength - 1L; (i < 0L); i++) {
							q = HalfKernel;
							u = v = InputData;
							InputData++;
							Sum = 0.0;
							for (j = i + 1L; (j < 0L); j++) {
								Sum -= (*++u - *--v) * *++q;
							}
							u = v2;
							for (j = -KernelHalfLength - i; (j < 0L); j++) {
								Sum -= (*++u - *--v) * *++q;
							}
							*OutputData++ = Sum;
						}
					}
					else {
						for (i = -SignalLength; (i < 0L); i++) {
							q = HalfKernel;
							u = v = InputData;
							InputData++;
							Sum = 0.0;
							for (j = 1L - KernelHalfLength; (j < 0L); j++) {
								if (u == u1) {
									u = v2;
								}
								if (v == v1) {
									v = u2;
								}
								Sum -= (*++u - *--v) * *++q;
							}
							*OutputData++ = Sum;
						}
					}
					break;
			}
			break;
		default:
			Status = ERROR;
			WRITE_ERROR(FirConvolveAntiSymmetric, "Invalid boundary convention")
			break;
	}
/**/DEBUG_WRITE_LEAVING(FirConvolveAntiSymmetric, "Done")
	return(Status);
} /* end FirConvolveAntiSymmetric */

/*--------------------------------------------------------------------------*/
extern int		FirConvolveAntiSymmetricVolume
				(
					float	*VolumeSource,		/* data to process */
					float	*VolumeDestination,	/* result */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					double	HalfKernel[],		/* causal part of the kernel */
					long	KernelHalfLength,	/* length of the causal part of the kernel */
					enum TBoundaryConvention
							Convention			/* boundary convention */
				)

/* the specified boundary convention applies to the input data only, not to the kernel */
/* the boundary convention applied to the kernel is FiniteDataSupport */
/* VolumeSource is a (float)volume of size (Nx x Ny x Nz) */
/* OutputData is a (float)volume of size (Nx x Ny x Nz) */
/* the origin for the kernel is its leftmost sample; it corresponds to its anti-symmetry axis */
/* the full length of the amti-symmetric kernel is (2L * KernelHalfLength - 1L) */
/* the 1D kernel is applied successively to each principal direction in a separable fashion */
/* success: return(!ERROR); failure: return(ERROR) */
/* general structure is as follows:
	for (i = 0L; (i < SignalLength); i++) {
		Sum = InputData[i] * HalfKernel[0];
		for (j = 1L; (j < KernelHalfLength); j++) {
			Sum += (InputData[i + j] - InputData[i - j]) * HalfKernel[j];
		}
		OutputData[i] = Sum;
	}
*/

{ /* begin FirConvolveAntiSymmetricVolume */

	double	*InBuffer = (double *)NULL, *OutBuffer = (double *)NULL;
	long	i, j, k;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(FirConvolveAntiSymmetricVolume, VolumeSource, Status,
/**/	"No VolumeSource")
/**/DEBUG_CHECK_NULL_POINTER(FirConvolveAntiSymmetricVolume, VolumeDestination, Status,
/**/	"No OutputData")
/**/DEBUG_CHECK_RANGE_LONG(FirConvolveAntiSymmetricVolume, Nx, 1L, LONG_MAX, Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(FirConvolveAntiSymmetricVolume, Ny, 1L, LONG_MAX, Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(FirConvolveAntiSymmetricVolume, Nz, 1L, LONG_MAX, Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_CHECK_NULL_POINTER(FirConvolveAntiSymmetricVolume, HalfKernel, Status,
/**/	"No kernel")
/**/DEBUG_CHECK_RANGE_LONG(FirConvolveAntiSymmetricVolume, KernelHalfLength, 1L, LONG_MAX, Status,
/**/	"Invalid kernel length (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(FirConvolveAntiSymmetricVolume, Status)
/**/DEBUG_WRITE_ENTERING(FirConvolveAntiSymmetricVolume,
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
			WRITE_ERROR(FirConvolveAntiSymmetricVolume, "Invalid boundary convention")
/**/		DEBUG_WRITE_LEAVING(FirConvolveAntiSymmetricVolume, "Done")
			return(Status);
	}
	if (1L < Nx) {
		AllocateLineDouble(&InBuffer, Nx, &Status);
		if (Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(FirConvolveAntiSymmetricVolume, "Done")
			return(Status);
		}
		AllocateLineDouble(&OutBuffer, Nx, &Status);
		if (Status == ERROR) {
			FreeLineDouble(&InBuffer);
/**/		DEBUG_WRITE_LEAVING(FirConvolveAntiSymmetricVolume, "Done")
			return(Status);
		}
		for (k = 0L; (k < Nz); k++) {
			for (j = 0L; (j < Ny); j++) {
				Status = GetxFloatToDouble(VolumeSource, Nx, Ny, Nz, 0L, j, k, InBuffer, Nx);
				if (Status == ERROR) {
					FreeLineDouble(&OutBuffer);
					FreeLineDouble(&InBuffer);
/**/				DEBUG_WRITE_LEAVING(FirConvolveAntiSymmetricVolume, "Done")
					return(Status);
				}
				Status = FirConvolveAntiSymmetric(InBuffer, OutBuffer, Nx,
					HalfKernel, KernelHalfLength, Convention);
				if (Status == ERROR) {
					FreeLineDouble(&OutBuffer);
					FreeLineDouble(&InBuffer);
/**/				DEBUG_WRITE_LEAVING(FirConvolveAntiSymmetricVolume, "Done")
					return(Status);
				}
				Status = PutxDoubleToFloat(VolumeDestination, Nx, Ny, Nz, 0L, j, k, OutBuffer, Nx);
				if (Status == ERROR) {
					FreeLineDouble(&OutBuffer);
					FreeLineDouble(&InBuffer);
/**/				DEBUG_WRITE_LEAVING(FirConvolveAntiSymmetricVolume, "Done")
					return(Status);
				}
			}
		}
		Status = FreeLineDouble(&OutBuffer);
		if (Status == ERROR) {
			FreeLineDouble(&InBuffer);
/**/		DEBUG_WRITE_LEAVING(FirConvolveAntiSymmetricVolume, "Done")
			return(Status);
		}
		Status = FreeLineDouble(&InBuffer);
		if (Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(FirConvolveAntiSymmetricVolume, "Done")
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
/**/		DEBUG_WRITE_LEAVING(FirConvolveAntiSymmetricVolume, "Done")
			return(Status);
		}
		AllocateLineDouble(&OutBuffer, Ny, &Status);
		if (Status == ERROR) {
			FreeLineDouble(&InBuffer);
/**/		DEBUG_WRITE_LEAVING(FirConvolveAntiSymmetricVolume, "Done")
			return(Status);
		}
		for (k = 0L; (k < Nz); k++) {
			for (i = 0L; (i < Nx); i++) {
				Status = GetyFloatToDouble(VolumeDestination, Nx, Ny, Nz, i, 0L, k, InBuffer, Ny);
				if (Status == ERROR) {
					FreeLineDouble(&OutBuffer);
					FreeLineDouble(&InBuffer);
/**/				DEBUG_WRITE_LEAVING(FirConvolveAntiSymmetricVolume, "Done")
					return(Status);
				}
				Status = FirConvolveAntiSymmetric(InBuffer, OutBuffer, Ny,
					HalfKernel, KernelHalfLength, Convention);
				if (Status == ERROR) {
					FreeLineDouble(&OutBuffer);
					FreeLineDouble(&InBuffer);
/**/				DEBUG_WRITE_LEAVING(FirConvolveAntiSymmetricVolume, "Done")
					return(Status);
				}
				Status = PutyDoubleToFloat(VolumeDestination, Nx, Ny, Nz, i, 0L, k, OutBuffer, Ny);
				if (Status == ERROR) {
					FreeLineDouble(&OutBuffer);
					FreeLineDouble(&InBuffer);
/**/				DEBUG_WRITE_LEAVING(FirConvolveAntiSymmetricVolume, "Done")
					return(Status);
				}
			}
		}
		Status = FreeLineDouble(&OutBuffer);
		if (Status == ERROR) {
			FreeLineDouble(&InBuffer);
/**/		DEBUG_WRITE_LEAVING(FirConvolveAntiSymmetricVolume, "Done")
			return(Status);
		}
		Status = FreeLineDouble(&InBuffer);
		if (Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(FirConvolveAntiSymmetricVolume, "Done")
			return(Status);
		}
	}
	if (1L < Nz) {
		AllocateLineDouble(&InBuffer, Nz, &Status);
		if (Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(FirConvolveAntiSymmetricVolume, "Done")
			return(Status);
		}
		AllocateLineDouble(&OutBuffer, Nz, &Status);
		if (Status == ERROR) {
			FreeLineDouble(&InBuffer);
/**/		DEBUG_WRITE_LEAVING(FirConvolveAntiSymmetricVolume, "Done")
			return(Status);
		}
		for (j = 0L; (j < Ny); j++) {
			for (i = 0L; (i < Nx); i++) {
				Status = GetzFloatToDouble(VolumeDestination, Nx, Ny, Nz, i, j, 0L, InBuffer, Nz);
				if (Status == ERROR) {
					FreeLineDouble(&OutBuffer);
					FreeLineDouble(&InBuffer);
/**/				DEBUG_WRITE_LEAVING(FirConvolveAntiSymmetricVolume, "Done")
					return(Status);
				}
				Status = FirConvolveAntiSymmetric(InBuffer, OutBuffer, Nz,
					HalfKernel, KernelHalfLength, Convention);
				if (Status == ERROR) {
					FreeLineDouble(&OutBuffer);
					FreeLineDouble(&InBuffer);
/**/				DEBUG_WRITE_LEAVING(FirConvolveAntiSymmetricVolume, "Done")
					return(Status);
				}
				Status = PutzDoubleToFloat(VolumeDestination, Nx, Ny, Nz, i, j, 0L, OutBuffer, Nz);
				if (Status == ERROR) {
					FreeLineDouble(&OutBuffer);
					FreeLineDouble(&InBuffer);
/**/				DEBUG_WRITE_LEAVING(FirConvolveAntiSymmetricVolume, "Done")
					return(Status);
				}
			}
		}
		Status = FreeLineDouble(&OutBuffer);
		if (Status == ERROR) {
			FreeLineDouble(&InBuffer);
/**/		DEBUG_WRITE_LEAVING(FirConvolveAntiSymmetricVolume, "Done")
			return(Status);
		}
		Status = FreeLineDouble(&InBuffer);
		if (Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(FirConvolveAntiSymmetricVolume, "Done")
			return(Status);
		}
	}
/**/DEBUG_WRITE_LEAVING(FirConvolveAntiSymmetricVolume, "Done")
	return(Status);
} /* end FirConvolveAntiSymmetricVolume */

/*--------------------------------------------------------------------------*/
extern int		FirConvolveSymmetric
				(
					double	InputData[],		/* data to process */
					double	OutputData[],		/* result */
					long	SignalLength,		/* length of the 1D data array */
					double	HalfKernel[],		/* causal part of the kernel */
					long	KernelHalfLength,	/* length of the causal part of the kernel */
					enum TBoundaryConvention
							Convention			/* boundary convention */
				)

/* the specified boundary convention applies to the input data only, not to the kernel */
/* the boundary convention applied to the kernel is FiniteDataSupport */
/* the input and the output have the same length */
/* the origin for the kernel is its leftmost sample; it corresponds to its symmetry axis */
/* the full length of the symmetric kernel is (2L * KernelHalfLength - 1L) */
/* success: return(!ERROR); failure: return(ERROR) */
/* general structure is as follows:
	for (i = 0L; (i < SignalLength); i++) {
		Sum = InputData[i] * HalfKernel[0];
		for (j = 1L; (j < KernelHalfLength); j++) {
			Sum += (InputData[i + j] + InputData[i - j]) * HalfKernel[j];
		}
		OutputData[i] = Sum;
	}
*/

{ /* begin FirConvolveSymmetric */

	double	*q, *u, *v;
	double	*u1, *u2, *v1, *v2;
	double	Sum;
	double	c0, c1, c2;
	long	i, j;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(FirConvolveSymmetric, InputData, Status,
/**/	"No input data")
/**/DEBUG_CHECK_NULL_POINTER(FirConvolveSymmetric, OutputData, Status,
/**/	"No output data")
/**/DEBUG_CHECK_RANGE_LONG(FirConvolveSymmetric, SignalLength, 1L, LONG_MAX, Status,
/**/	"Invalid signal length (should be strictly positive)")
/**/DEBUG_CHECK_NULL_POINTER(FirConvolveSymmetric, HalfKernel, Status,
/**/	"No kernel")
/**/DEBUG_CHECK_RANGE_LONG(FirConvolveSymmetric, KernelHalfLength, 1L, LONG_MAX, Status,
/**/	"Invalid kernel length (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(FirConvolveSymmetric, Status)
/**/DEBUG_WRITE_ENTERING(FirConvolveSymmetric,
/**/	"About to perform FIR symmetric convolution")

	switch (Convention) {
		case AntiMirrorOnBounds:
			Status = ERROR;
			WRITE_ERROR(FirConvolveSymmetric, "Not yet implemented")
			break;
		case FiniteDataSupport:
			switch (KernelHalfLength) {
				case 1L:
					c0 = *HalfKernel;
					for (i = -SignalLength; (i < 0L); i++) {
						*OutputData++ = *InputData++ * c0;
					}
					break;
				case 2L:
					c0 = *HalfKernel++;
					c1 = *HalfKernel;
					if (3L <= SignalLength) {
						u = InputData;
						*OutputData++ = *InputData * c0 + *++u * c1;
						for (i = 2L - SignalLength; (i < 0L); i++) {
							v = InputData;
							InputData = u++;
							*OutputData++ = *InputData * c0 + (*u + *v) * c1;
						}
						*OutputData = *u * c0 + *InputData * c1;
					}
					else {
						switch (SignalLength) {
							case 2L:
								v = InputData++;
								*OutputData++ = *v * c0 + *InputData * c1;
								*OutputData = *InputData * c0 + *v * c1;
								break;
							case 1L:
								*OutputData = *InputData * c0;
								break;
							default:
								Status = ERROR;
								WRITE_ERROR(FirConvolveSymmetric,
									"Unexpected signal length")
								break;
						}
					}
					break;
				case 3L:
					c0 = *HalfKernel++;
					c1 = *HalfKernel++;
					c2 = *HalfKernel;
					if (5L <= SignalLength) {
						u2 = InputData;
						u2++;
						u1 = u2++;
						*OutputData++ = *InputData * c0 + *u1 * c1 + *u2 * c2;
						v1 = InputData;
						InputData = u1;
						u1 = u2++;
						*OutputData++ = *InputData * c0 + (*u1 + *v1) * c1 + *u2 * c2;
						for (i = 4L - SignalLength; (i < 0L); i++) {
							v2 = v1;
							v1 = InputData;
							InputData = u1;
							u1 = u2++;
							*OutputData++ = *InputData * c0 + (*u1 + *v1) * c1 + (*u2 + *v2) * c2;
						}
						*OutputData++ = *u1 * c0 + (*u2 + *InputData) * c1 + *v1 * c2;
						*OutputData = *u2 * c0 + *u1 * c1 + *InputData * c2;
					}
					else {
						switch (SignalLength) {
							case 4L:
								u2 = InputData;
								u2++;
								u1 = u2++;
								*OutputData++ = *InputData * c0 + *u1 * c1 + *u2 * c2;
								v = u2++;
								*OutputData++ = *u1 * c0 + (*v + *InputData) * c1 + *u2 * c2;
								*OutputData++ = *v * c0 + (*u2 + *u1) * c1 + *InputData * c2;
								*OutputData = *u2 * c0 + *v * c1 + *u1 * c2;
								break;
							case 3L:
								u2 = InputData;
								u2++;
								u1 = u2++;
								*OutputData++ = *InputData * c0 + *u1 * c1 + *u2 * c2;
								*OutputData++ = *u1 * c0 + (*u2 + *InputData) * c1;
								*OutputData = *u2 * c0 + *u1 * c1 + *InputData * c2;
								break;
							case 2L:
								v = InputData++;
								*OutputData++ = *v * c0 + *InputData * c1;
								*OutputData = *InputData * c0 + *v * c1;
								break;
							case 1L:
								*OutputData = *InputData * c0;
								break;
							default:
								Status = ERROR;
								WRITE_ERROR(FirConvolveSymmetric,
									"Unexpected signal length")
								break;
						}
					}
					break;
				default:
					if ((2L * KernelHalfLength - 1L) <= SignalLength) {
						for (i = 1L - KernelHalfLength; (i < 0L); i++) {
							q = HalfKernel;
							u = v = InputData;
							Sum = *InputData++ * *q;
							for (j = 1L - KernelHalfLength - i; (j < 0L); j++) {
								Sum += (*++u + *--v) * *++q;
							}
							for (j = i; (j < 0L); j++) {
								Sum += *++u * *++q;
							}
							*OutputData++ = Sum;
						}
						for (i = 2L * KernelHalfLength - SignalLength - 2L; (i < 0L); i++) {
							q = HalfKernel;
							u = v = InputData;
							Sum = *InputData++ * *q;
							for (j = 1L - KernelHalfLength; (j < 0L); j++) {
								Sum += (*++u + *--v) * *++q;
							}
							*OutputData++ = Sum;
						}
						for (i = 1L - KernelHalfLength; (i < 0L); i++) {
							q = HalfKernel;
							u = v = InputData;
							Sum = *InputData++ * *q;
							for (j = i + 1L; (j < 0L); j++) {
								Sum += (*++u + *--v) * *++q;
							}
							for (j = -KernelHalfLength - i; (j < 0L); j++) {
								Sum += *--v * *++q;
							}
							*OutputData++ = Sum;
						}
					}
					else if (KernelHalfLength <= SignalLength) {
						for (i = KernelHalfLength - SignalLength - 1L; (i < 0L); i++) {
							q = HalfKernel;
							u = v = InputData;
							Sum = *InputData++ * *q;
							for (j = KernelHalfLength - SignalLength - i - 1L; (j < 0L); j++) {
								Sum += (*++u + *--v) * *++q;
							}
							for (j = SignalLength - 2L * KernelHalfLength + i + 2L; (j < 0L); j++) {
								Sum += *++u * *++q;
							}
							*OutputData++ = Sum;
						}
						for (i = SignalLength - (SignalLength / 2L) - KernelHalfLength + 1L; (i < 0L); i++) {
							q = HalfKernel;
							u = v = InputData;
							Sum = *InputData++ * *q;
							for (j = -(SignalLength / 2L) - i; (j < 0L); j++) {
								Sum += (*++u + *--v) * *++q;
							}
							for (j = 2L * (SignalLength / 2L) - SignalLength + 2L * i + 1L; (j < 0L); j++) {
								Sum += *++u * *++q;
							}
							*OutputData++ = Sum;
						}
						if ((2L * (SignalLength / 2L) + 1L) == SignalLength) {
							q = HalfKernel;
							u = v = InputData;
							Sum = *InputData++ * *q;
							for (j = -(SignalLength / 2L); (j < 0L); j++) {
								Sum += (*++u + *--v) * *++q;
							}
							*OutputData++ = Sum;
							i = 1L;
						}
						for (i += SignalLength / 2L - KernelHalfLength + 1L; (i < 0L); i++) {
							q = HalfKernel;
							u = v = InputData;
							Sum = *InputData++ * *q;
							for (j = KernelHalfLength - SignalLength + i; (j < 0L); j++) {
								Sum += (*++u + *--v) * *++q;
							}
							for (j = SignalLength - 2L * KernelHalfLength - 2L * i + 1L; (j < 0L); j++) {
								Sum += *--v * *++q;
							}
							*OutputData++ = Sum;
						}
						for (i = KernelHalfLength - SignalLength - 1L; (i < 0L); i++) {
							q = HalfKernel;
							u = v = InputData;
							Sum = *InputData++ * *q;
							for (j = i + 1L; (j < 0L); j++) {
								Sum += (*++u + *--v) * *++q;
							}
							for (j = -KernelHalfLength - i; (j < 0L); j++) {
								Sum += *--v * *++q;
							}
							*OutputData++ = Sum;
						}
					}
					else {
						for (i = -(SignalLength / 2L); (i < 0L); i++) {
							q = HalfKernel;
							u = v = InputData;
							Sum = *InputData++ * *q;
							for (j = -(SignalLength / 2L) - i; (j < 0L); j++) {
								Sum += (*++u + *--v) * *++q;
							}
							for (j = 2L * (SignalLength / 2L) - SignalLength + 2L * i + 1L; (j < 0L); j++) {
								Sum += *++u * *++q;
							}
							*OutputData++ = Sum;
						}
						if ((2L * (SignalLength / 2L) + 1L) == SignalLength) {
							q = HalfKernel;
							u = v = InputData;
							Sum = *InputData++ * *q;
							for (j = -(SignalLength / 2L) - i; (j < 0L); j++) {
								Sum += (*++u + *--v) * *++q;
							}
							*OutputData++ = Sum;
							i = 1L;
						}
						for (i += SignalLength / 2L - SignalLength; (i < 0L); i++) {
							q = HalfKernel;
							u = v = InputData;
							Sum = *InputData++ * *q;
							for (j = i + 1L; (j < 0L); j++) {
								Sum += (*++u + *--v) * *++q;
							}
							for (j = -SignalLength - 2L * i - 1L; (j < 0L); j++) {
								Sum += *--v * *++q;
							}
							*OutputData++ = Sum;
						}
					}
					break;
			}
			break;
		case MirrorOffBounds:
			Status = ERROR;
			WRITE_ERROR(FirConvolveSymmetric, "Not yet implemented")
			break;
		case MirrorOnBounds:
			switch (KernelHalfLength) {
				case 1L:
					c0 = *HalfKernel;
					for (i = -SignalLength; (i < 0L); i++) {
						*OutputData++ = *InputData++ * c0;
					}
					break;
				case 2L:
					c0 = *HalfKernel++;
					c1 = *HalfKernel;
					u = InputData;
					if (3L <= SignalLength) {
						*OutputData++ = *InputData * c0 + *++u * 2.0 * c1;
						for (i = 2L - SignalLength; (i < 0L); i++) {
							v = InputData;
							InputData = u++;
							*OutputData++ = *InputData * c0 + (*u + *v) * c1;
						}
						*OutputData = *u * c0 + *InputData * 2.0 * c1;
					}
					else {
						switch (SignalLength) {
							case 2L:
								c1 += c1;
								*OutputData++ = *InputData * c0 + *++u * c1;
								*OutputData = *u * c0 + *InputData * c1;
								break;
							case 1L:
								*OutputData = *InputData * (c0 + 2.0 * c1);
								break;
							default:
								Status = ERROR;
								WRITE_ERROR(FirConvolveSymmetric,
									"Unexpected signal length")
								break;
						}
					}
					break;
				case 3L:
					c0 = *HalfKernel++;
					c1 = *HalfKernel++;
					c2 = *HalfKernel;
					u2 = InputData;
					u2++;
					u1 = u2++;
					if (5L <= SignalLength) {
						*OutputData++ = *InputData * c0 + *u1 * 2.0 * c1 + *u2 * 2.0 * c2;
						v1 = InputData;
						InputData = u1;
						u1 = u2++;
						*OutputData++ = *InputData * c0 + (*u1 + *v1) * c1 + (*u2 + *InputData) * c2;
						for (i = 4L - SignalLength; (i < 0L); i++) {
							v2 = v1;
							v1 = InputData;
							InputData = u1;
							u1 = u2++;
							*OutputData++ = *InputData * c0 + (*u1 + *v1) * c1 + (*u2 + *v2) * c2;
						}
						*OutputData++ = *u1 * c0 + (*u2 + *InputData) * c1 + (*u1 + *v1) * c2;
						*OutputData = *u2 * c0 + *u1 * 2.0 * c1 + *InputData * 2.0 * c2;
					}
					else {
						switch (SignalLength) {
							case 4L:
								v = u2++;
								*OutputData++ = *InputData * c0 + *u1 * 2.0 * c1 + *v * 2.0 * c2;
								*OutputData++ = *u1 * c0 + (*v + *InputData) * c1 + (*u2 + *u1) * c2;
								*OutputData++ = *v * c0 + (*u2 + *u1) * c1 + (*v + *InputData) * c2;
								*OutputData = *u2 * c0 + *v * 2.0 * c1 + *u1 * 2.0 * c2;
								break;
							case 3L:
								c2 *= 2.0;
								*OutputData++ = *InputData * c0 + *u1 * 2.0 * c1 + *u2 * c2;
								*OutputData++ = *u1 * c0 + (*u2 + *InputData) * c1 + *u1 * c2;
								*OutputData = *u2 * c0 + *u1 * 2.0 * c1 + *InputData * c2;
								break;
							case 2L:
								c0 += 2.0 * c2;
								c1 *= 2.0;
								*OutputData++ = *InputData * c0 + *u1 * c1;
								*OutputData = *u1 * c0 + *InputData * c1;
								break;
							case 1L:
								*OutputData = *InputData * (c0 + 2.0 * c1 + 2.0 * c2);
								break;
							default:
								Status = ERROR;
								WRITE_ERROR(FirConvolveSymmetric,
									"Unexpected signal length")
								break;
						}
					}
					break;
				default:
					if ((2L * KernelHalfLength - 1L) <= SignalLength) {
						for (i = 1L - KernelHalfLength; (i < 0L); i++) {
							q = HalfKernel;
							u = v = InputData;
							Sum = *InputData++ * *q;
							for (j = 1L - KernelHalfLength - i; (j < 0L); j++) {
								Sum += (*++u + *--v) * *++q;
							}
							for (j = i; (j < 0L); j++) {
								Sum += (*++u + *++v) * *++q;
							}
							*OutputData++ = Sum;
						}
						for (i = 2L * KernelHalfLength - SignalLength - 2L; (i < 0L); i++) {
							q = HalfKernel;
							u = v = InputData;
							Sum = *InputData++ * *q;
							for (j = 1L - KernelHalfLength; (j < 0L); j++) {
								Sum += (*++u + *--v) * *++q;
							}
							*OutputData++ = Sum;
						}
						for (i = 1L - KernelHalfLength; (i < 0L); i++) {
							q = HalfKernel;
							u = v = InputData;
							Sum = *InputData++ * *q;
							for (j = i + 1L; (j < 0L); j++) {
								Sum += (*++u + *--v) * *++q;
							}
							for (j = -KernelHalfLength - i; (j < 0L); j++) {
								Sum += (*--u + *--v) * *++q;
							}
							*OutputData++ = Sum;
						}
					}
					else if (KernelHalfLength <= SignalLength) {
						for (i = KernelHalfLength - SignalLength - 1L; (i < 0L); i++) {
							q = HalfKernel;
							u = v = InputData;
							Sum = *InputData++ * *q;
							for (j = KernelHalfLength - SignalLength - i - 1L; (j < 0L); j++) {
								Sum += (*++u + *--v) * *++q;
							}
							for (j = SignalLength - 2L * KernelHalfLength + i + 2L; (j < 0L); j++) {
								Sum += (*++u + *++v) * *++q;
							}
							*OutputData++ = Sum;
						}
						for (i = SignalLength - (SignalLength / 2L) - KernelHalfLength + 1L; (i < 0L); i++) {
							q = HalfKernel;
							u = v = InputData;
							Sum = *InputData++ * *q;
							for (j = -(SignalLength / 2L) - i; (j < 0L); j++) {
								Sum += (*++u + *--v) * *++q;
							}
							for (j = 2L * (SignalLength / 2L) - SignalLength + 2L * i + 1L; (j < 0L); j++) {
								Sum += (*++u + *++v) * *++q;
							}
							for (j = SignalLength - (SignalLength / 2L) - KernelHalfLength - i; (j < 0L); j++) {
								Sum += (*--u + *++v) * *++q;
							}
							*OutputData++ = Sum;
						}
						if ((2L * (SignalLength / 2L) + 1L) == SignalLength) {
							q = HalfKernel;
							u = v = InputData;
							Sum = *InputData++ * *q;
							for (j = -(SignalLength / 2L); (j < 0L); j++) {
								Sum += (*++u + *--v) * *++q;
							}
							for (j = SignalLength / 2L - KernelHalfLength + 1L; (j < 0L); j++) {
								Sum += (*--u + *++v) * *++q;
							}
							*OutputData++ = Sum;
							i = 1L;
						}
						for (i += SignalLength / 2L - KernelHalfLength + 1L; (i < 0L); i++) {
							q = HalfKernel;
							u = v = InputData;
							Sum = *InputData++ * *q;
							for (j = KernelHalfLength - SignalLength + i; (j < 0L); j++) {
								Sum += (*++u + *--v) * *++q;
							}
							for (j = SignalLength - 2L * KernelHalfLength - 2L * i + 1L; (j < 0L); j++) {
								Sum += (*--u + *--v) * *++q;
							}
							for (j = i; (j < 0L); j++) {
								Sum += (*--u + *++v) * *++q;
							}
							*OutputData++ = Sum;
						}
						for (i = KernelHalfLength - SignalLength - 1L; (i < 0L); i++) {
							q = HalfKernel;
							u = v = InputData;
							Sum = *InputData++ * *q;
							for (j = i + 1L; (j < 0L); j++) {
								Sum += (*++u + *--v) * *++q;
							}
							for (j = -KernelHalfLength - i; (j < 0L); j++) {
								Sum += (*--u + *--v) * *++q;
							}
							*OutputData++ = Sum;
						}
					}
					else {
						if (SignalLength == 1L) {
							Sum = *HalfKernel++;
							for (j = 1L - KernelHalfLength; (j < 0L); j++) {
								Sum += 2.0 * *HalfKernel++;
							}
							*OutputData = Sum * *InputData;
						}
						else {
							v1 = InputData;
							u1 = InputData + (ptrdiff_t)(SignalLength - 1L);
							for (i = -SignalLength; (i < 0L); i++) {
								q = HalfKernel;
								u = v = InputData;
								Sum = *InputData++ * *q;
								j = 1L - KernelHalfLength;
								while ((j < 0L) && (u < u1) && (v1 < v)) {
									Sum += (*++u + *--v) * *++q;
									j++;
								}
								if ((u == u1) && (v == v1)) {
									while (j < 0L) {
										while ((j < 0L) && (v < u1)) {
											Sum += (*--u + *++v) * *++q;
											j++;
										}
										while ((j < 0L) && (u < u1)) {
											Sum += (*++u + *--v) * *++q;
											j++;
										}
									}
								}
								else if (u == u1) {
									while (j < 0L) {
										while ((j < 0L) && (v1 < v)) {
											Sum += (*--u + *--v) * *++q;
											j++;
										}
										while ((j < 0L) && (v1 < u)) {
											Sum += (*--u + *++v) * *++q;
											j++;
										}
										while ((j < 0L) && (v < u1)) {
											Sum += (*++u + *++v) * *++q;
											j++;
										}
										while ((j < 0L) && (u < u1)) {
											Sum += (*++u + *--v) * *++q;
											j++;
										}
									}
								}
								else {
									while (j < 0L) {
										while ((j < 0L) && (u < u1)) {
											Sum += (*++u + *++v) * *++q;
											j++;
										}
										while ((j < 0L) && (v < u1)) {
											Sum += (*--u + *++v) * *++q;
											j++;
										}
										while ((j < 0L) && (v1 < u)) {
											Sum += (*--u + *--v) * *++q;
											j++;
										}
										while ((j < 0L) && (v1 < v)) {
											Sum += (*++u + *--v) * *++q;
											j++;
										}
									}
								}
								*OutputData++ = Sum;
							}
						}
					}
					break;
			}
			break;
		case Periodic:
			switch (KernelHalfLength) {
				case 1L:
					c0 = *HalfKernel;
					for (i = -SignalLength; (i < 0L); i++) {
						*OutputData++ = *InputData++ * c0;
					}
					break;
				case 2L:
					c0 = *HalfKernel++;
					c1 = *HalfKernel;
					u = InputData;
					if (3L <= SignalLength) {
						u1 = InputData;
						v = InputData + (ptrdiff_t)SignalLength;
						*OutputData++ = *InputData * c0 + (*++u + *--v) * c1;
						for (i = 2L - SignalLength; (i < 0L); i++) {
							v = InputData;
							InputData = u++;
							*OutputData++ = *InputData * c0 + (*u + *v) * c1;
						}
						*OutputData = *u * c0 + (*u1 + *InputData) * c1;
					}
					else {
						switch (SignalLength) {
							case 2L:
								c1 += c1;
								*OutputData++ = *InputData * c0 + *++u * c1;
								*OutputData = *u * c0 + *InputData * c1;
								break;
							case 1L:
								*OutputData = *InputData * (c0 + 2.0 * c1);
								break;
							default:
								Status = ERROR;
								WRITE_ERROR(FirConvolveSymmetric,
									"Unexpected signal length")
								break;
						}
					}
					break;
				case 3L:
					c0 = *HalfKernel++;
					c1 = *HalfKernel++;
					c2 = *HalfKernel;
					u2 = InputData;
					u2++;
					u1 = u2++;
					if (5L <= SignalLength) {
						u = InputData;
						v = u++;
						v1 = InputData + (ptrdiff_t)SignalLength;
						v2 = v1--;
						*OutputData++ = *InputData * c0 + (*u1 + *--v2) * c1 + (*u2 + *--v1) * c2;
						v1 = InputData;
						InputData = u1;
						u1 = u2++;
						*OutputData++ = *InputData * c0 + (*u1 + *v1) * c1 + (*u2 + *v2) * c2;
						for (i = 4L - SignalLength; (i < 0L); i++) {
							v2 = v1;
							v1 = InputData;
							InputData = u1;
							u1 = u2++;
							*OutputData++ = *InputData * c0 + (*u1 + *v1) * c1 + (*u2 + *v2) * c2;
						}
						*OutputData++ = *u1 * c0 + (*u2 + *InputData) * c1 + (*v + *v1) * c2;
						*OutputData = *u2 * c0 + (*v + *u1) * c1 + (*u + *InputData) * c2;
					}
					else {
						switch (SignalLength) {
							case 4L:
								c2 *= 2.0;
								v = u2++;
								*OutputData++ = *InputData * c0 + (*u1 + *u2) * c1 + *v * c2;
								*OutputData++ = *u1 * c0 + (*v + *InputData) * c1 + *u2 * c2;
								*OutputData++ = *v * c0 + (*u2 + *u1) * c1 + *InputData * c2;
								*OutputData = *u2 * c0 + (*InputData + *v) * c1 + *u1 * c2;
								break;
							case 3L:
								c1 += c2;
								*OutputData++ = *InputData * c0 + (*u1 + *u2) * c1;
								*OutputData++ = *u1 * c0 + (*u2 + *InputData) * c1;
								*OutputData = *u2 * c0 + (*InputData + *u1) * c1;
								break;
							case 2L:
								c0 += 2.0 * c2;
								c1 += c1;
								*OutputData++ = *InputData * c0 + *u1 * c1;
								*OutputData = *u1 * c0 + *InputData * c1;
								break;
							case 1L:
								*OutputData = *InputData * (c0 + 2.0 * c1 + 2.0 * c2);
								break;
							default:
								Status = ERROR;
								WRITE_ERROR(FirConvolveSymmetric,
									"Unexpected signal length")
								break;
						}
					}
					break;
				default:
					v2 = InputData;
					v1 = v2--;
					u1 = InputData + (ptrdiff_t)SignalLength;
					u2 = u1--;
					if ((2L * KernelHalfLength - 1L) <= SignalLength) {
						for (i = 1L - KernelHalfLength; (i < 0L); i++) {
							q = HalfKernel;
							u = v = InputData;
							Sum = *InputData++ * *q;
							for (j = 1L - KernelHalfLength - i; (j < 0L); j++) {
								Sum += (*++u + *--v) * *++q;
							}
							v = u2;
							for (j = i; (j < 0L); j++) {
								Sum += (*++u + *--v) * *++q;
							}
							*OutputData++ = Sum;
						}
						for (i = 2L * KernelHalfLength - SignalLength - 2L; (i < 0L); i++) {
							q = HalfKernel;
							u = v = InputData;
							Sum = *InputData++ * *q;
							for (j = 1L - KernelHalfLength; (j < 0L); j++) {
								Sum += (*++u + *--v) * *++q;
							}
							*OutputData++ = Sum;
						}
						for (i = 1L - KernelHalfLength; (i < 0L); i++) {
							q = HalfKernel;
							u = v = InputData;
							Sum = *InputData++ * *q;
							for (j = i + 1L; (j < 0L); j++) {
								Sum += (*++u + *--v) * *++q;
							}
							u = v2;
							for (j = -KernelHalfLength - i; (j < 0L); j++) {
								Sum += (*++u + *--v) * *++q;
							}
							*OutputData++ = Sum;
						}
					}
					else if (KernelHalfLength <= SignalLength) {
						for (i = KernelHalfLength - SignalLength - 1L; (i < 0L); i++) {
							q = HalfKernel;
							u = v = InputData;
							Sum = *InputData++ * *q;
							for (j = KernelHalfLength - SignalLength - i - 1L; (j < 0L); j++) {
								Sum += (*++u + *--v) * *++q;
							}
							v = u2;
							for (j = SignalLength - 2L * KernelHalfLength + i + 2L; (j < 0L); j++) {
								Sum += (*++u + *--v) * *++q;
							}
							*OutputData++ = Sum;
						}
						for (i = SignalLength - (SignalLength / 2L) - KernelHalfLength + 1L; (i < 0L); i++) {
							q = HalfKernel;
							u = v = InputData;
							Sum = *InputData++ * *q;
							for (j = -(SignalLength / 2L) - i; (j < 0L); j++) {
								Sum += (*++u + *--v) * *++q;
							}
							v = u2;
							for (j = 2L * (SignalLength / 2L) - SignalLength + 2L * i + 1L; (j < 0L); j++) {
								Sum += (*++u + *--v) * *++q;
							}
							u = v2;
							for (j = SignalLength - (SignalLength / 2L) - KernelHalfLength - i; (j < 0L); j++) {
								Sum += (*++u + *--v) * *++q;
							}
							*OutputData++ = Sum;
						}
						if ((2L * (SignalLength / 2L) + 1L) == SignalLength) {
							q = HalfKernel;
							u = v = InputData;
							Sum = *InputData++ * *q;
							for (j = -(SignalLength / 2L); (j < 0L); j++) {
								Sum += (*++u + *--v) * *++q;
							}
							v = u2;
							u = v2;
							for (j = SignalLength / 2L - KernelHalfLength + 1L; (j < 0L); j++) {
								Sum += (*++u + *--v) * *++q;
							}
							*OutputData++ = Sum;
							i = 1L;
						}
						for (i += SignalLength / 2L - KernelHalfLength + 1L; (i < 0L); i++) {
							q = HalfKernel;
							u = v = InputData;
							Sum = *InputData++ * *q;
							for (j = KernelHalfLength - SignalLength + i; (j < 0L); j++) {
								Sum += (*++u + *--v) * *++q;
							}
							u = v2;
							for (j = SignalLength - 2L * KernelHalfLength - 2L * i + 1L; (j < 0L); j++) {
								Sum += (*++u + *--v) * *++q;
							}
							v = u2;
							for (j = i; (j < 0L); j++) {
								Sum += (*++u + *--v) * *++q;
							}
							*OutputData++ = Sum;
						}
						for (i = KernelHalfLength - SignalLength - 1L; (i < 0L); i++) {
							q = HalfKernel;
							u = v = InputData;
							Sum = *InputData++ * *q;
							for (j = i + 1L; (j < 0L); j++) {
								Sum += (*++u + *--v) * *++q;
							}
							u = v2;
							for (j = -KernelHalfLength - i; (j < 0L); j++) {
								Sum += (*++u + *--v) * *++q;
							}
							*OutputData++ = Sum;
						}
					}
					else {
						for (i = -SignalLength; (i < 0L); i++) {
							q = HalfKernel;
							u = v = InputData;
							Sum = *InputData++ * *q;
							for (j = 1L - KernelHalfLength; (j < 0L); j++) {
								if (u == u1) {
									u = v2;
								}
								if (v == v1) {
									v = u2;
								}
								Sum += (*++u + *--v) * *++q;
							}
							*OutputData++ = Sum;
						}
					}
					break;
			}
			break;
		default:
			Status = ERROR;
			WRITE_ERROR(FirConvolveSymmetric, "Invalid boundary convention")
			break;
	}
/**/DEBUG_WRITE_LEAVING(FirConvolveSymmetric, "Done")
	return(Status);
} /* end FirConvolveSymmetric */

/*--------------------------------------------------------------------------*/
extern int		FirConvolveSymmetricVolume
				(
					float	*VolumeSource,		/* data to process */
					float	*VolumeDestination,	/* result */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					double	HalfKernel[],		/* causal part of the kernel */
					long	KernelHalfLength,	/* length of the causal part of the kernel */
					enum TBoundaryConvention
							Convention			/* boundary convention */
				)

/* the specified boundary convention applies to the input data only, not to the kernel */
/* the boundary convention applied to the kernel is FiniteDataSupport */
/* VolumeSource is a (float)volume of size (Nx x Ny x Nz) */
/* OutputData is a (float)volume of size (Nx x Ny x Nz) */
/* the origin for the kernel is its leftmost sample; it corresponds to its symmetry axis */
/* the full length of the symmetric kernel is (2L * KernelHalfLength - 1L) */
/* the 1D kernel is applied successively to each principal direction in a separable fashion */
/* success: return(!ERROR); failure: return(ERROR) */
/* general structure is as follows:
	for (i = 0L; (i < SignalLength); i++) {
		Sum = InputData[i] * HalfKernel[0];
		for (j = 1L; (j < KernelHalfLength); j++) {
			Sum += (InputData[i + j] + InputData[i - j]) * HalfKernel[j];
		}
		OutputData[i] = Sum;
	}
*/

{ /* begin FirConvolveSymmetricVolume */

	double	*InBuffer = (double *)NULL, *OutBuffer = (double *)NULL;
	long	i, j, k;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(FirConvolveSymmetricVolume, VolumeSource, Status,
/**/	"No VolumeSource")
/**/DEBUG_CHECK_NULL_POINTER(FirConvolveSymmetricVolume, VolumeDestination, Status,
/**/	"No OutputData")
/**/DEBUG_CHECK_RANGE_LONG(FirConvolveSymmetricVolume, Nx, 1L, LONG_MAX, Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(FirConvolveSymmetricVolume, Ny, 1L, LONG_MAX, Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(FirConvolveSymmetricVolume, Nz, 1L, LONG_MAX, Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_CHECK_NULL_POINTER(FirConvolveSymmetricVolume, HalfKernel, Status,
/**/	"No kernel")
/**/DEBUG_CHECK_RANGE_LONG(FirConvolveSymmetricVolume, KernelHalfLength, 1L, LONG_MAX, Status,
/**/	"Invalid kernel length (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(FirConvolveSymmetricVolume, Status)
/**/DEBUG_WRITE_ENTERING(FirConvolveSymmetricVolume,
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
			WRITE_ERROR(FirConvolveSymmetricVolume, "Invalid boundary convention")
/**/		DEBUG_WRITE_LEAVING(FirConvolveSymmetricVolume, "Done")
			return(Status);
	}
	if (1L < Nx) {
		AllocateLineDouble(&InBuffer, Nx, &Status);
		if (Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(FirConvolveSymmetricVolume, "Done")
			return(Status);
		}
		AllocateLineDouble(&OutBuffer, Nx, &Status);
		if (Status == ERROR) {
			FreeLineDouble(&InBuffer);
/**/		DEBUG_WRITE_LEAVING(FirConvolveSymmetricVolume, "Done")
			return(Status);
		}
		for (k = 0L; (k < Nz); k++) {
			for (j = 0L; (j < Ny); j++) {
				Status = GetxFloatToDouble(VolumeSource, Nx, Ny, Nz, 0L, j, k, InBuffer, Nx);
				if (Status == ERROR) {
					FreeLineDouble(&OutBuffer);
					FreeLineDouble(&InBuffer);
/**/				DEBUG_WRITE_LEAVING(FirConvolveSymmetricVolume, "Done")
					return(Status);
				}
				Status = FirConvolveSymmetric(InBuffer, OutBuffer, Nx,
					HalfKernel, KernelHalfLength, Convention);
				if (Status == ERROR) {
					FreeLineDouble(&OutBuffer);
					FreeLineDouble(&InBuffer);
/**/				DEBUG_WRITE_LEAVING(FirConvolveSymmetricVolume, "Done")
					return(Status);
				}
				Status = PutxDoubleToFloat(VolumeDestination, Nx, Ny, Nz, 0L, j, k, OutBuffer, Nx);
				if (Status == ERROR) {
					FreeLineDouble(&OutBuffer);
					FreeLineDouble(&InBuffer);
/**/				DEBUG_WRITE_LEAVING(FirConvolveSymmetricVolume, "Done")
					return(Status);
				}
			}
		}
		Status = FreeLineDouble(&OutBuffer);
		if (Status == ERROR) {
			FreeLineDouble(&InBuffer);
/**/		DEBUG_WRITE_LEAVING(FirConvolveSymmetricVolume, "Done")
			return(Status);
		}
		Status = FreeLineDouble(&InBuffer);
		if (Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(FirConvolveSymmetricVolume, "Done")
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
/**/		DEBUG_WRITE_LEAVING(FirConvolveSymmetricVolume, "Done")
			return(Status);
		}
		AllocateLineDouble(&OutBuffer, Ny, &Status);
		if (Status == ERROR) {
			FreeLineDouble(&InBuffer);
/**/		DEBUG_WRITE_LEAVING(FirConvolveSymmetricVolume, "Done")
			return(Status);
		}
		for (k = 0L; (k < Nz); k++) {
			for (i = 0L; (i < Nx); i++) {
				Status = GetyFloatToDouble(VolumeDestination, Nx, Ny, Nz, i, 0L, k, InBuffer, Ny);
				if (Status == ERROR) {
					FreeLineDouble(&OutBuffer);
					FreeLineDouble(&InBuffer);
/**/				DEBUG_WRITE_LEAVING(FirConvolveSymmetricVolume, "Done")
					return(Status);
				}
				Status = FirConvolveSymmetric(InBuffer, OutBuffer, Ny,
					HalfKernel, KernelHalfLength, Convention);
				if (Status == ERROR) {
					FreeLineDouble(&OutBuffer);
					FreeLineDouble(&InBuffer);
/**/				DEBUG_WRITE_LEAVING(FirConvolveSymmetricVolume, "Done")
					return(Status);
				}
				Status = PutyDoubleToFloat(VolumeDestination, Nx, Ny, Nz, i, 0L, k, OutBuffer, Ny);
				if (Status == ERROR) {
					FreeLineDouble(&OutBuffer);
					FreeLineDouble(&InBuffer);
/**/				DEBUG_WRITE_LEAVING(FirConvolveSymmetricVolume, "Done")
					return(Status);
				}
			}
		}
		Status = FreeLineDouble(&OutBuffer);
		if (Status == ERROR) {
			FreeLineDouble(&InBuffer);
/**/		DEBUG_WRITE_LEAVING(FirConvolveSymmetricVolume, "Done")
			return(Status);
		}
		Status = FreeLineDouble(&InBuffer);
		if (Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(FirConvolveSymmetricVolume, "Done")
			return(Status);
		}
	}
	if (1L < Nz) {
		AllocateLineDouble(&InBuffer, Nz, &Status);
		if (Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(FirConvolveSymmetricVolume, "Done")
			return(Status);
		}
		AllocateLineDouble(&OutBuffer, Nz, &Status);
		if (Status == ERROR) {
			FreeLineDouble(&InBuffer);
/**/		DEBUG_WRITE_LEAVING(FirConvolveSymmetricVolume, "Done")
			return(Status);
		}
		for (j = 0L; (j < Ny); j++) {
			for (i = 0L; (i < Nx); i++) {
				Status = GetzFloatToDouble(VolumeDestination, Nx, Ny, Nz, i, j, 0L, InBuffer, Nz);
				if (Status == ERROR) {
					FreeLineDouble(&OutBuffer);
					FreeLineDouble(&InBuffer);
/**/				DEBUG_WRITE_LEAVING(FirConvolveSymmetricVolume, "Done")
					return(Status);
				}
				Status = FirConvolveSymmetric(InBuffer, OutBuffer, Nz,
					HalfKernel, KernelHalfLength, Convention);
				if (Status == ERROR) {
					FreeLineDouble(&OutBuffer);
					FreeLineDouble(&InBuffer);
/**/				DEBUG_WRITE_LEAVING(FirConvolveSymmetricVolume, "Done")
					return(Status);
				}
				Status = PutzDoubleToFloat(VolumeDestination, Nx, Ny, Nz, i, j, 0L, OutBuffer, Nz);
				if (Status == ERROR) {
					FreeLineDouble(&OutBuffer);
					FreeLineDouble(&InBuffer);
/**/				DEBUG_WRITE_LEAVING(FirConvolveSymmetricVolume, "Done")
					return(Status);
				}
			}
		}
		Status = FreeLineDouble(&OutBuffer);
		if (Status == ERROR) {
			FreeLineDouble(&InBuffer);
/**/		DEBUG_WRITE_LEAVING(FirConvolveSymmetricVolume, "Done")
			return(Status);
		}
		Status = FreeLineDouble(&InBuffer);
		if (Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(FirConvolveSymmetricVolume, "Done")
			return(Status);
		}
	}
/**/DEBUG_WRITE_LEAVING(FirConvolveSymmetricVolume, "Done")
	return(Status);
} /* end FirConvolveSymmetricVolume */

