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

/*****************************************************************************
 *	Toolbox includes
 ****************************************************************************/
#include	"dft.h"
#include	"dht.h"
#include	"fourierconvolve.h"
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
extern int		ManyConvolveFourier
				(
					double	Data[],				/* data for in-place processing */
					double	KernelDht[],		/* discrete Hartley transform of the kernel */
					double	CaS[],				/* Hartley transform coefficients */
					double	*ScratchBuffer,		/* scratch buffer */
					long	SignalLength,		/* signal length */
					double	Shift				/* additional translation */
				)

/* Fourier convolution with any kernel */
/* same conventions as for OneConvolveFourier */
/* preprocessing steps have to be carrried out elsewhere */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin ManyConvolveFourier */

	double	*p, *q, *u, *v, *t;
	double	c, s, ct, st;
	double	butterfly;
	long	i;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(ManyConvolveFourier, Data, Status,
/**/	"No in-place data")
/**/DEBUG_CHECK_NULL_POINTER(ManyConvolveFourier, KernelDht, Status,
/**/	"No kernel")
/**/DEBUG_CHECK_NULL_POINTER(ManyConvolveFourier, CaS, Status,
/**/	"No Hartley-coefficients")
/**/DEBUG_CHECK_NULL_POINTER(ManyConvolveFourier, ScratchBuffer, Status,
/**/	"No buffer")
/**/DEBUG_CHECK_RANGE_LONG(ManyConvolveFourier, SignalLength, 1L, LONG_MAX, Status,
/**/	"Invalid length (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(ManyConvolveFourier, Status)
/**/DEBUG_WRITE_ENTERING(ManyConvolveFourier,
/**/	"About to perform a preprocessed Fourier convolution")

	Status = DiscreteHartleyTransform(Data, ScratchBuffer, CaS, SignalLength, TRUE);
	if (Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(ManyConvolveFourier, "Done")
		return(Status);
	}
	p = Data;
	q = Data + (ptrdiff_t)SignalLength;
	u = KernelDht;
	v = KernelDht + (ptrdiff_t)SignalLength;
	t = ScratchBuffer;
	*t++ = *p++ * *u++;
	for (q--, v--, i = 1L; (i < SignalLength); u++, v--, i++) {
		*t++ = (*p++ * (*u + *v) + *q-- * (*u - *v)) / 2.0;
	}
	if (Shift == 0.0) {
		Data = (double *)memcpy(Data, ScratchBuffer,
			(size_t)(SignalLength * (long)sizeof(double)));
	}
	else {
		p = Data;
		q = Data + (ptrdiff_t)SignalLength;
		u = ScratchBuffer;
		v = ScratchBuffer + (ptrdiff_t)SignalLength;
		c = ct = cos(2.0 * PI * Shift / (double)SignalLength);
		s = st = -sin(2.0 * PI * Shift / (double)SignalLength);
		*p++ = *u++;
		q--;
		v--;
		while (p < q) {
			*p++ = *u * c - *v * s;
			*q-- = *v-- * c + *u++ * s;
			butterfly = c * ct - s * st;
			s = s * ct + c * st;
			c = butterfly;
		}
		if (p == q) {
			*p = *u * c;
		}
	}
	Status = DiscreteHartleyTransform(Data, ScratchBuffer, CaS, SignalLength, FALSE);
/**/DEBUG_WRITE_LEAVING(ManyConvolveFourier, "Done")
	return(Status);
} /* end ManyConvolveFourier */

/*--------------------------------------------------------------------------*/
extern int		ManyConvolveFourierSymmetricKernel
				(
					double	Data[],				/* data for in-place processing */
					double	SymmetricKernelDht[],
												/* discrete Hartley transform of the kernel */
					double	CaS[],				/* Hartley transform coefficients */
					double	*ScratchBuffer,		/* scratch buffer */
					long	SignalLength,		/* signal length */
					double	Shift				/* additional translation */
				)

/* Fourier convolution with a symmetric kernel */
/* same conventions as for OneConvolveFourierSymmetricKernel */
/* preprocessing steps have to be carrried out elsewhere */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin ManyConvolveFourierSymmetricKernel */

	double	*p, *q, *u, *v, *t;
	double	c, s, ct, st;
	double	butterfly;
	long	i;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(ManyConvolveFourierSymmetricKernel, Data, Status,
/**/	"No in-place data")
/**/DEBUG_CHECK_NULL_POINTER(ManyConvolveFourierSymmetricKernel, SymmetricKernelDht, Status,
/**/	"No kernel")
/**/DEBUG_CHECK_NULL_POINTER(ManyConvolveFourierSymmetricKernel, CaS, Status,
/**/	"No Hartley-coefficients")
/**/DEBUG_CHECK_NULL_POINTER(ManyConvolveFourierSymmetricKernel, ScratchBuffer, Status,
/**/	"No buffer")
/**/DEBUG_CHECK_RANGE_LONG(ManyConvolveFourierSymmetricKernel, SignalLength, 1L, LONG_MAX, Status,
/**/	"Invalid length (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(ManyConvolveFourierSymmetricKernel, Status)
/**/DEBUG_WRITE_ENTERING(ManyConvolveFourierSymmetricKernel,
/**/	"About to perform a preprocessed Fourier convolution with a symmetric kernel")

	Status = DiscreteHartleyTransform(Data, ScratchBuffer, CaS, SignalLength, TRUE);
	if (Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(ManyConvolveFourierSymmetricKernel, "Done")
		return(Status);
	}
	p = SymmetricKernelDht;
	q = Data;
	t = ScratchBuffer;
	for (i = 0L; (i < SignalLength); i++) {
		*t++ = *q++ * *p++;
	}
	if (Shift == 0.0) {
		Data = (double *)memcpy(Data, ScratchBuffer,
			(size_t)(SignalLength * (long)sizeof(double)));
	}
	else {
		p = Data;
		q = Data + (ptrdiff_t)SignalLength;
		u = ScratchBuffer;
		v = ScratchBuffer + (ptrdiff_t)SignalLength;
		c = ct = cos(2.0 * PI * Shift / (double)SignalLength);
		s = st = -sin(2.0 * PI * Shift / (double)SignalLength);
		*p++ = *u++;
		q--;
		v--;
		while (p < q) {
			*p++ = *u * c - *v * s;
			*q-- = *v-- * c + *u++ * s;
			butterfly = c * ct - s * st;
			s = s * ct + c * st;
			c = butterfly;
		}
		if (p == q) {
			*p = *u * c;
		}
	}
	Status = DiscreteHartleyTransform(Data, ScratchBuffer, CaS, SignalLength, FALSE);
/**/DEBUG_WRITE_LEAVING(ManyConvolveFourierSymmetricKernel, "Done")
	return(Status);
} /* end ManyConvolveFourierSymmetricKernel */

/*--------------------------------------------------------------------------*/
extern int		OneConvolveFourier
				(
					double	Data[],				/* data for in-place processing */
					double	Kernel[],			/* kernel for in-place processing */
					long	SignalLength,		/* signal length */
					double	Shift,				/* additional translation */
					int		*Status				/* error management */
				)

/* Fourier convolution with any kernel */
/* the kernel has an infinite, periodic (SignalLength) support */
/* the kernel origin (hot spot) is at index [0] */
/* the highest coordinate (SignalLength-2)/2 is at index [(SignalLength-2)/2] for SignalLength even */
/* the highest coordinate (SignalLength-1)/2 is at index [(SignalLength-1)/2] for SignalLength odd */
/* the lowest coordinate -SignalLength/2 is at index [SignalLength/2] for SignalLength even */
/* the lowest coordinate -(SignalLength-1)/2 is at index [(SignalLength+1)/2] for SignalLength odd */
/* the coordinate -1 is at index [SignalLength-1] for SignalLength even or odd */
/* the signal has an infinite, periodic (SignalLength) support */
/* the output has same length SignalLength than the input */
/* the output has already been allocated */
/* the result of the convolution overwrites the input */
/* the inverse DHT of the kernel overwrites the kernel */
/* if the kernel is finite support, don't forget to pad! */
/* success: return(!ERROR); failure: return(ERROR) */
/* the returned value is duplicated in Status */

{ /* begin OneConvolveFourier */

	double	*CaS = (double *)NULL, *ScratchBuffer = (double *)NULL;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(OneConvolveFourier, Data, *Status,
/**/	"No in-place data")
/**/DEBUG_CHECK_NULL_POINTER(OneConvolveFourier, Kernel, *Status,
/**/	"No kernel")
/**/DEBUG_CHECK_RANGE_LONG(OneConvolveFourier, SignalLength, 1L, LONG_MAX, *Status,
/**/	"Invalid length (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(OneConvolveFourier, *Status)
/**/DEBUG_WRITE_ENTERING(OneConvolveFourier,
/**/	"About to perform Fourier convolution")

	AllocateLineDouble(&ScratchBuffer, SignalLength, Status);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(OneConvolveFourier, "Done")
		return(*Status);
	}
	AllocateLineDouble(&CaS, SignalLength, Status);
	if (*Status == ERROR) {
		FreeLineDouble(&ScratchBuffer);
/**/	DEBUG_WRITE_LEAVING(OneConvolveFourier, "Done")
		return(*Status);
	}
	*Status = GetCaS(CaS, SignalLength);
	if (*Status == ERROR) {
		FreeLineDouble(&CaS);
		FreeLineDouble(&ScratchBuffer);
/**/	DEBUG_WRITE_LEAVING(OneConvolveFourier, "Done")
		return(*Status);
	}
	*Status = DiscreteHartleyTransform(Kernel, ScratchBuffer, CaS, SignalLength, FALSE);
	if (*Status == ERROR) {
		FreeLineDouble(&CaS);
		FreeLineDouble(&ScratchBuffer);
/**/	DEBUG_WRITE_LEAVING(OneConvolveFourier, "Done")
		return(*Status);
	}
	*Status = ManyConvolveFourier(Data, Kernel, CaS, ScratchBuffer, SignalLength, Shift);
	if (*Status == ERROR) {
		FreeLineDouble(&CaS);
		FreeLineDouble(&ScratchBuffer);
/**/	DEBUG_WRITE_LEAVING(OneConvolveFourier, "Done")
		return(*Status);
	}
	*Status = FreeLineDouble(&CaS);
	if (*Status == ERROR) {
		FreeLineDouble(&ScratchBuffer);
/**/	DEBUG_WRITE_LEAVING(OneConvolveFourier, "Done")
		return(*Status);
	}
	*Status = FreeLineDouble(&ScratchBuffer);
/**/DEBUG_WRITE_LEAVING(OneConvolveFourier, "Done")
	return(*Status);
} /* end oneConvolveFourier */

/*--------------------------------------------------------------------------*/
extern int		OneConvolveFourierSymmetricKernel
				(
					double	Data[],				/* data for in-place processing */
					double	SymmetricKernel[],	/* symmetric kernel for in-place processing */
					long	SignalLength,		/* signal length */
					double	Shift,				/* additional translation */
					int		*Status				/* error management */
				)

/* Fourier convolution with a symmetric kernel */
/* the kernel has an infinite, periodic (SignalLength) support */
/* the kernel origin (hot spot) is at index [0] */
/* the highest coordinate (SignalLength-2)/2 is at index [(SignalLength-2)/2] for SignalLength even */
/* for symmetry, kernel[(SignalLength-2)/2] = 0 for SignalLength even */
/* the highest coordinate (SignalLength-1)/2 is at index [(SignalLength-1)/2] for SignalLength odd */
/* no special symmetry requirement for SignalLength odd */
/* the lowest coordinate -SignalLength/2 is at index [SignalLength/2] for SignalLength even */
/* the lowest coordinate -(SignalLength-1)/2 is at index [(SignalLength+1)/2] for SignalLength odd */
/* the coordinate -1 is at index [SignalLength-1] for SignalLength even or odd */
/* the signal has an infinite, periodic (SignalLength) support */
/* the output has the same length SignalLength than the input */
/* the output has already been allocated */
/* the kernel is even-symmetric around 0 (k[x] = k[-x]) */
/* only the coord. >= 0 of the kernel need be given */
/* the result of the convolution overwrites the input */
/* the inverse DHT of the kernel overwrites the kernel */
/* if the kernel is finite support, don't forget to pad! */
/* observe the separate symmetries for SignalLength odd/even! */
/* for a symmetric kernel, DHT <=> DFT */
/* success: return(!ERROR); failure: return(ERROR) */
/* the returned value is duplicated in Status */

{ /* begin OneConvolveFourierSymmetricKernel */

	double	*CaS = (double *)NULL, *ScratchBuffer = (double *)NULL;
	double	*p, *q;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(OneConvolveFourierSymmetricKernel, Data, *Status,
/**/	"No in-place data")
/**/DEBUG_CHECK_NULL_POINTER(OneConvolveFourierSymmetricKernel, SymmetricKernel, *Status,
/**/	"No kernel")
/**/DEBUG_CHECK_RANGE_LONG(OneConvolveFourierSymmetricKernel, SignalLength, 1L, LONG_MAX, *Status,
/**/	"Invalid length (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(OneConvolveFourierSymmetricKernel, *Status)
/**/DEBUG_WRITE_ENTERING(OneConvolveFourierSymmetricKernel,
/**/	"About to perform Fourier symmetric convolution")

	AllocateLineDouble(&ScratchBuffer, SignalLength, Status);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(OneConvolveFourierSymmetricKernel, "Done")
		return(*Status);
	}
	AllocateLineDouble(&CaS, SignalLength, Status);
	if (*Status == ERROR) {
		FreeLineDouble(&ScratchBuffer);
/**/	DEBUG_WRITE_LEAVING(OneConvolveFourierSymmetricKernel, "Done")
		return(*Status);
	}
	*Status = GetCaS(CaS, SignalLength);
	if (*Status == ERROR) {
		FreeLineDouble(&CaS);
		FreeLineDouble(&ScratchBuffer);
/**/	DEBUG_WRITE_LEAVING(OneConvolveFourierSymmetricKernel, "Done")
		return(*Status);
	}
	p = SymmetricKernel;
	q = SymmetricKernel + (ptrdiff_t)SignalLength;
	while (++p < --q) {
		*q = *p;
	}
	if (p == q) {
		*q = 0.0;
	}
	*Status = DiscreteHartleyTransform(SymmetricKernel, ScratchBuffer, CaS, SignalLength, FALSE);
	if (*Status == ERROR) {
		FreeLineDouble(&CaS);
		FreeLineDouble(&ScratchBuffer);
/**/	DEBUG_WRITE_LEAVING(OneConvolveFourierSymmetricKernel, "Done")
		return(*Status);
	}
	*Status = ManyConvolveFourierSymmetricKernel(Data, SymmetricKernel, CaS,
		ScratchBuffer, SignalLength, Shift);
	if (*Status == ERROR) {
		FreeLineDouble(&CaS);
		FreeLineDouble(&ScratchBuffer);
/**/	DEBUG_WRITE_LEAVING(OneConvolveFourierSymmetricKernel, "Done")
		return(*Status);
	}
	*Status = FreeLineDouble(&CaS);
	if (*Status == ERROR) {
		FreeLineDouble(&ScratchBuffer);
/**/	DEBUG_WRITE_LEAVING(OneConvolveFourierSymmetricKernel, "Done")
		return(*Status);
	}
	*Status = FreeLineDouble(&ScratchBuffer);
/**/DEBUG_WRITE_LEAVING(OneConvolveFourierSymmetricKernel, "Done")
	return(*Status);
} /* end OneConvolveFourierSymmetricKernel */

