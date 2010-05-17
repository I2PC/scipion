/*****************************************************************************
 *	System includes
 ****************************************************************************/
#include	<stdio.h>

/*****************************************************************************
 *	Toolbox defines
 ****************************************************************************/
#include	"configs.h"
#include	"debug.h"
#include	"error.h"
#include	"tboundaryconvention.h"
#include	"tsplinebasis.h"

/*****************************************************************************
 *	Toolbox includes
 ****************************************************************************/
#include	"changebasis.h"
#include	"convert.h"
#include	"firconvolve.h"
#include	"getput.h"
#include	"gradient.h"
#include	"kernel.h"
#include	"messagedisplay.h"

/*****************************************************************************
 *	Conditional includes
 ****************************************************************************/
#ifdef		DEBUG
#include	<float.h>
#include	<limits.h>
#include	<stddef.h>
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
extern int		LinearGradient
				(
					double	InputData[],		/* input 1D data */
					double	OutputGradient[],	/* output 1D gradient */
					long	SignalLength,		/* signal length */
					long	Degree,				/* degree of the spline space */
					enum TBoundaryConvention
							Convention,			/* boundary convention */
					double	Tolerance,			/* admissible relative error */
					int		*Status				/* error management */
				)

/* compute the gradient of 1D data */
/* InputData and OutputGradient have size SignalLength */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin LinearGradient */

	double	*Bcoeff = (double *)NULL, *HalfKernel = (double *)NULL;
	long	HalfDegree = (1L < Degree) ? (Degree / 2L) : (1L);
	long	i;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(LinearGradient, InputData, *Status,
/**/	"Missing input")
/**/DEBUG_CHECK_NULL_POINTER(LinearGradient, OutputGradient, *Status,
/**/	"Missing output")
/**/DEBUG_CHECK_RANGE_LONG(LinearGradient, SignalLength, 1L, LONG_MAX, *Status,
/**/	"Invalid length (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(LinearGradient, Degree, 1L, LONG_MAX, *Status,
/**/	"Invalid degree (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_DOUBLE(LinearGradient, Tolerance, 0.0, DBL_MAX, *Status,
/**/	"Invalid Tolerance (should be positive)")
/**/DEBUG_RETURN_ON_ERROR(LinearGradient, *Status)
/**/DEBUG_WRITE_ENTERING(LinearGradient,
/**/	"About to compute a 1D gradient")

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
			WRITE_ERROR(LinearGradient, "Invalid boundary convention")
/**/		DEBUG_WRITE_LEAVING(LinearGradient, "Done")
			return(*Status);
	}
	AllocateLineDouble(&Bcoeff, SignalLength, Status);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(LinearGradient, "Done")
		return(*Status);
	}
	AllocateLineDouble(&HalfKernel, HalfDegree + 1L, Status);
	if (*Status == ERROR) {
		FreeLineDouble(&Bcoeff);
/**/	DEBUG_WRITE_LEAVING(LinearGradient, "Done")
		return(*Status);
	}
	ChangeBasis(InputData, Bcoeff, SignalLength, CardinalSpline, BasicSpline,
		Degree, Convention, Tolerance, Status);
	if (*Status == ERROR) {
		FreeLineDouble(&HalfKernel);
		FreeLineDouble(&Bcoeff);
/**/	DEBUG_WRITE_LEAVING(LinearGradient, "Done")
		return(*Status);
	}
	for (i = 0L; (i < HalfDegree); i++) {
		*Status = Bspline(Degree - 1L, 0.5 + (double)i, &(HalfKernel[i]));
		if (*Status == ERROR) {
			FreeLineDouble(&HalfKernel);
			FreeLineDouble(&Bcoeff);
/**/		DEBUG_WRITE_LEAVING(LinearGradient, "Done")
			return(*Status);
		}
	}
	HalfKernel[i] = 0.0;
	do {
		HalfKernel[i] -= HalfKernel[i - 1L];
	} while (0L < --i);
	HalfKernel[i] = 0.0;
	*Status = FirConvolveAntiSymmetric(Bcoeff, OutputGradient, SignalLength,
		HalfKernel, HalfDegree + 1L, Convention);
	if (*Status == ERROR) {
		FreeLineDouble(&HalfKernel);
		FreeLineDouble(&Bcoeff);
/**/	DEBUG_WRITE_LEAVING(LinearGradient, "Done")
		return(*Status);
	}
	*Status = FreeLineDouble(&HalfKernel);
	if (*Status == ERROR) {
		FreeLineDouble(&Bcoeff);
/**/	DEBUG_WRITE_LEAVING(LinearGradient, "Done")
		return(*Status);
	}
	*Status = FreeLineDouble(&Bcoeff);
/**/DEBUG_WRITE_LEAVING(LinearGradient, "Done")
	return(*Status);
} /* end LinearGradient */

/*--------------------------------------------------------------------------*/
extern int		PlanarGradient
				(
					float	*InputImage,		/* input 2D data */
					float	*OutputGradient[],	/* output 2D gradient array [x, y] */
					long	Nx,					/* width of the image */
					long	Ny,					/* height of the image */
					long	Degree,				/* degree of the spline space */
					enum TBoundaryConvention
							Convention,			/* boundary convention */
					double	Tolerance,			/* admissible relative error */
					int		*Status				/* error management */
				)

/* compute the gradient of 2D data */
/* InputImage has size (Nx x Ny) */
/* OutputGradient is an array of two elements */
/* each element of OutputGradient is a (float *) pointer to an image */
/* each element has size (Nx x Ny) */
/* the first element is the gradient along x */
/* the second element is the gradient along y */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin PlanarGradient */

	float	*p;
	double	*Line1 = (double *)NULL, *Line2 = (double *)NULL, *HalfKernel = (double *)NULL;
	long	HalfDegree = (1L < Degree) ? (Degree / 2L) : (1L);
	long	x, y;
	long	i;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(PlanarGradient, InputImage, *Status,
/**/	"Missing input")
/**/DEBUG_CHECK_NULL_POINTER(PlanarGradient, OutputGradient, *Status,
/**/	"Missing output")
/**/DEBUG_CHECK_NULL_POINTER(PlanarGradient, OutputGradient[0], *Status,
/**/	"Missing output x gradient")
/**/DEBUG_CHECK_NULL_POINTER(PlanarGradient, OutputGradient[1], *Status,
/**/	"Missing output y gradient")
/**/DEBUG_CHECK_RANGE_LONG(PlanarGradient, Nx, 1L, LONG_MAX, *Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(PlanarGradient, Ny, 1L, LONG_MAX, *Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(PlanarGradient, Degree, 1L, LONG_MAX, *Status,
/**/	"Invalid degree (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_DOUBLE(PlanarGradient, Tolerance, 0.0, DBL_MAX, *Status,
/**/	"Invalid Tolerance (should be positive)")
/**/DEBUG_RETURN_ON_ERROR(PlanarGradient, *Status)
/**/DEBUG_WRITE_ENTERING(PlanarGradient,
/**/	"About to compute a 2D gradient")

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
			WRITE_ERROR(PlanarGradient, "Invalid boundary convention")
/**/		DEBUG_WRITE_LEAVING(PlanarGradient, "Done")
			return(*Status);
	}
	AllocateLineDouble(&HalfKernel, HalfDegree + 1L, Status);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(PlanarGradient, "Done")
		return(*Status);
	}
	for (i = 0L; (i < HalfDegree); i++) {
		*Status = Bspline(Degree - 1L, 0.5 + (double)i, &(HalfKernel[i]));
		if (*Status == ERROR) {
			FreeLineDouble(&HalfKernel);
/**/		DEBUG_WRITE_LEAVING(PlanarGradient, "Done")
			return(*Status);
		}
	}
	HalfKernel[i] = 0.0;
	do {
		HalfKernel[i] -= HalfKernel[i - 1L];
	} while (0L < --i);
	HalfKernel[i] = 0.0;
	if (1L < Nx) {
		AllocateLineDouble(&Line1, Nx, Status);
		if (*Status == ERROR) {
			FreeLineDouble(&HalfKernel);
/**/		DEBUG_WRITE_LEAVING(PlanarGradient, "Done")
			return(*Status);
		}
		AllocateLineDouble(&Line2, Nx, Status);
		if (*Status == ERROR) {
			FreeLineDouble(&Line1);
			FreeLineDouble(&HalfKernel);
/**/		DEBUG_WRITE_LEAVING(PlanarGradient, "Done")
			return(*Status);
		}
		for (y = 0L; (y < Ny); y++) {
			*Status = GetxFloatToDouble(InputImage, Nx, Ny, 1L, 0L, y, 0L, Line1, Nx);
			if (*Status == ERROR) {
				FreeLineDouble(&Line2);
				FreeLineDouble(&Line1);
				FreeLineDouble(&HalfKernel);
/**/			DEBUG_WRITE_LEAVING(PlanarGradient, "Done")
				return(*Status);
			}
			ChangeBasis(Line1, Line2, Nx, CardinalSpline, BasicSpline,
				Degree, Convention, Tolerance, Status);
			if (*Status == ERROR) {
				FreeLineDouble(&Line2);
				FreeLineDouble(&Line1);
				FreeLineDouble(&HalfKernel);
/**/			DEBUG_WRITE_LEAVING(PlanarGradient, "Done")
				return(*Status);
			}
			*Status = FirConvolveAntiSymmetric(Line2, Line1, Nx,
				HalfKernel, HalfDegree + 1L, Convention);
			if (*Status == ERROR) {
				FreeLineDouble(&Line2);
				FreeLineDouble(&Line1);
				FreeLineDouble(&HalfKernel);
/**/			DEBUG_WRITE_LEAVING(PlanarGradient, "Done")
				return(*Status);
			}
			*Status = PutxDoubleToFloat(OutputGradient[0], Nx, Ny, 1L, 0L, y, 0L, Line1, Nx);
			if (*Status == ERROR) {
				FreeLineDouble(&Line2);
				FreeLineDouble(&Line1);
				FreeLineDouble(&HalfKernel);
/**/			DEBUG_WRITE_LEAVING(PlanarGradient, "Done")
				return(*Status);
			}
		}
		*Status = FreeLineDouble(&Line2);
		if (*Status == ERROR) {
			FreeLineDouble(&Line1);
			FreeLineDouble(&HalfKernel);
/**/		DEBUG_WRITE_LEAVING(PlanarGradient, "Done")
			return(*Status);
		}
		*Status = FreeLineDouble(&Line1);
		if (*Status == ERROR) {
			FreeLineDouble(&HalfKernel);
/**/		DEBUG_WRITE_LEAVING(PlanarGradient, "Done")
			return(*Status);
		}
	}
	else {
		p = OutputGradient[0];
		for (y = Ny; (0L < y); y--) {
			*p++ = 0.0F;
		}
	}
	if (1L < Ny) {
		AllocateLineDouble(&Line1, Ny, Status);
		if (*Status == ERROR) {
			FreeLineDouble(&HalfKernel);
/**/		DEBUG_WRITE_LEAVING(PlanarGradient, "Done")
			return(*Status);
		}
		AllocateLineDouble(&Line2, Ny, Status);
		if (*Status == ERROR) {
			FreeLineDouble(&Line1);
			FreeLineDouble(&HalfKernel);
/**/		DEBUG_WRITE_LEAVING(PlanarGradient, "Done")
			return(*Status);
		}
		for (x = 0L; (x < Nx); x++) {
			*Status = GetyFloatToDouble(InputImage, Nx, Ny, 1L, x, 0L, 0L, Line1, Ny);
			if (*Status == ERROR) {
				FreeLineDouble(&Line2);
				FreeLineDouble(&Line1);
				FreeLineDouble(&HalfKernel);
/**/			DEBUG_WRITE_LEAVING(PlanarGradient, "Done")
				return(*Status);
			}
			ChangeBasis(Line1, Line2, Ny, CardinalSpline, BasicSpline,
				Degree, Convention, Tolerance, Status);
			if (*Status == ERROR) {
				FreeLineDouble(&Line2);
				FreeLineDouble(&Line1);
				FreeLineDouble(&HalfKernel);
/**/			DEBUG_WRITE_LEAVING(PlanarGradient, "Done")
				return(*Status);
			}
			*Status = FirConvolveAntiSymmetric(Line2, Line1, Ny,
				HalfKernel, HalfDegree + 1L, Convention);
			if (*Status == ERROR) {
				FreeLineDouble(&Line2);
				FreeLineDouble(&Line1);
				FreeLineDouble(&HalfKernel);
/**/			DEBUG_WRITE_LEAVING(PlanarGradient, "Done")
				return(*Status);
			}
			*Status = PutyDoubleToFloat(OutputGradient[1], Nx, Ny, 1L, x, 0L, 0L, Line1, Ny);
			if (*Status == ERROR) {
				FreeLineDouble(&Line2);
				FreeLineDouble(&Line1);
				FreeLineDouble(&HalfKernel);
/**/			DEBUG_WRITE_LEAVING(PlanarGradient, "Done")
				return(*Status);
			}
		}
		*Status = FreeLineDouble(&Line2);
		if (*Status == ERROR) {
			FreeLineDouble(&Line1);
			FreeLineDouble(&HalfKernel);
/**/		DEBUG_WRITE_LEAVING(PlanarGradient, "Done")
			return(*Status);
		}
		*Status = FreeLineDouble(&Line1);
		if (*Status == ERROR) {
			FreeLineDouble(&HalfKernel);
/**/		DEBUG_WRITE_LEAVING(PlanarGradient, "Done")
			return(*Status);
		}
	}
	else {
		p = OutputGradient[1];
		for (x = Nx; (0L < x); x--) {
			*p++ = 0.0F;
		}
	}
	*Status = FreeLineDouble(&HalfKernel);
/**/DEBUG_WRITE_LEAVING(PlanarGradient, "Done")
	return(*Status);
} /* end PlanarGradient */

/*--------------------------------------------------------------------------*/
extern int		VolumetricGradient
				(
					float	*InputVolume,		/* input 3D data */
					float	*OutputGradient[],	/* output 3D gradient array [x, y, z] */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					long	Degree,				/* degree of the spline space */
					enum TBoundaryConvention
							Convention,			/* boundary convention */
					double	Tolerance,			/* admissible relative error */
					int		*Status				/* error management */
				)

/* compute the gradient of 3D data */
/* InputVolume has size (Nx x Ny x Nz) */
/* OutputGradient is an array of three elements */
/* each element of OutputGradient is a (float *) pointer to a volume */
/* each element of OutputGradient has size (Nx x Ny x Nz) */
/* the first element is the gradient along x */
/* the second element is the gradient along y */
/* the third element is the gradient along z */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin VolumetricGradient */

	float	*p;
	double	*Line1 = (double *)NULL, *Line2 = (double *)NULL, *HalfKernel = (double *)NULL;
	long	HalfDegree = (1L < Degree) ? (Degree / 2L) : (1L);
	long	x, y, z;
	long	i;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(VolumetricGradient, InputVolume, *Status,
/**/	"Missing input")
/**/DEBUG_CHECK_NULL_POINTER(VolumetricGradient, OutputGradient, *Status,
/**/	"Missing output")
/**/DEBUG_CHECK_NULL_POINTER(VolumetricGradient, OutputGradient[0], *Status,
/**/	"Missing output x gradient")
/**/DEBUG_CHECK_NULL_POINTER(VolumetricGradient, OutputGradient[1], *Status,
/**/	"Missing output y gradient")
/**/DEBUG_CHECK_NULL_POINTER(VolumetricGradient, OutputGradient[2], *Status,
/**/	"Missing output z gradient")
/**/DEBUG_CHECK_RANGE_LONG(VolumetricGradient, Nx, 1L, LONG_MAX, *Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(VolumetricGradient, Ny, 1L, LONG_MAX, *Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(VolumetricGradient, Nz, 1L, LONG_MAX, *Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(VolumetricGradient, Degree, 1L, LONG_MAX, *Status,
/**/	"Invalid degree (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_DOUBLE(VolumetricGradient, Tolerance, 0.0, DBL_MAX, *Status,
/**/	"Invalid Tolerance (should be positive)")
/**/DEBUG_RETURN_ON_ERROR(VolumetricGradient, *Status)
/**/DEBUG_WRITE_ENTERING(VolumetricGradient,
/**/	"About to compute a 3D gradient")

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
			WRITE_ERROR(VolumetricGradient, "Invalid boundary convention")
/**/		DEBUG_WRITE_LEAVING(VolumetricGradient, "Done")
			return(*Status);
	}
	AllocateLineDouble(&HalfKernel, HalfDegree + 1L, Status);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(VolumetricGradient, "Done")
		return(*Status);
	}
	for (i = 0L; (i < HalfDegree); i++) {
		*Status = Bspline(Degree - 1L, 0.5 + (double)i, &(HalfKernel[i]));
		if (*Status == ERROR) {
			FreeLineDouble(&HalfKernel);
/**/		DEBUG_WRITE_LEAVING(VolumetricGradient, "Done")
			return(*Status);
		}
	}
	HalfKernel[i] = 0.0;
	do {
		HalfKernel[i] -= HalfKernel[i - 1L];
	} while (0L < --i);
	HalfKernel[i] = 0.0;
	if (1L < Nx) {
		AllocateLineDouble(&Line1, Nx, Status);
		if (*Status == ERROR) {
			FreeLineDouble(&HalfKernel);
/**/		DEBUG_WRITE_LEAVING(VolumetricGradient, "Done")
			return(*Status);
		}
		AllocateLineDouble(&Line2, Nx, Status);
		if (*Status == ERROR) {
			FreeLineDouble(&Line1);
			FreeLineDouble(&HalfKernel);
/**/		DEBUG_WRITE_LEAVING(VolumetricGradient, "Done")
			return(*Status);
		}
		for (z = 0L; (z < Nz); z++) {
			for (y = 0L; (y < Ny); y++) {
				*Status = GetxFloatToDouble(InputVolume, Nx, Ny, Nz, 0L, y, z, Line1, Nx);
				if (*Status == ERROR) {
					FreeLineDouble(&Line2);
					FreeLineDouble(&Line1);
					FreeLineDouble(&HalfKernel);
/**/				DEBUG_WRITE_LEAVING(VolumetricGradient, "Done")
					return(*Status);
				}
				ChangeBasis(Line1, Line2, Nx, CardinalSpline, BasicSpline,
					Degree, Convention, Tolerance, Status);
				if (*Status == ERROR) {
					FreeLineDouble(&Line2);
					FreeLineDouble(&Line1);
					FreeLineDouble(&HalfKernel);
/**/				DEBUG_WRITE_LEAVING(VolumetricGradient, "Done")
					return(*Status);
				}
				*Status = FirConvolveAntiSymmetric(Line2, Line1, Nx,
					HalfKernel, HalfDegree + 1L, Convention);
				if (*Status == ERROR) {
					FreeLineDouble(&Line2);
					FreeLineDouble(&Line1);
					FreeLineDouble(&HalfKernel);
/**/				DEBUG_WRITE_LEAVING(VolumetricGradient, "Done")
					return(*Status);
				}
				*Status = PutxDoubleToFloat(OutputGradient[0], Nx, Ny, Nz, 0L, y, z, Line1, Nx);
				if (*Status == ERROR) {
					FreeLineDouble(&Line2);
					FreeLineDouble(&Line1);
					FreeLineDouble(&HalfKernel);
/**/				DEBUG_WRITE_LEAVING(VolumetricGradient, "Done")
					return(*Status);
				}
			}
		}
		*Status = FreeLineDouble(&Line2);
		if (*Status == ERROR) {
			FreeLineDouble(&Line1);
			FreeLineDouble(&HalfKernel);
/**/		DEBUG_WRITE_LEAVING(VolumetricGradient, "Done")
			return(*Status);
		}
		*Status = FreeLineDouble(&Line1);
		if (*Status == ERROR) {
			FreeLineDouble(&HalfKernel);
/**/		DEBUG_WRITE_LEAVING(VolumetricGradient, "Done")
			return(*Status);
		}
	}
	else {
		p = OutputGradient[0];
		for (z = Nz; (0L < z); z--) {
			for (y = Ny; (0L < y); y--) {
				*p++ = 0.0F;
			}
		}
	}
	if (1L < Ny) {
		AllocateLineDouble(&Line1, Ny, Status);
		if (*Status == ERROR) {
			FreeLineDouble(&HalfKernel);
/**/		DEBUG_WRITE_LEAVING(VolumetricGradient, "Done")
			return(*Status);
		}
		AllocateLineDouble(&Line2, Ny, Status);
		if (*Status == ERROR) {
			FreeLineDouble(&Line1);
			FreeLineDouble(&HalfKernel);
/**/		DEBUG_WRITE_LEAVING(VolumetricGradient, "Done")
			return(*Status);
		}
		for (z = 0L; (z < Nz); z++) {
			for (x = 0L; (x < Nx); x++) {
				*Status = GetyFloatToDouble(InputVolume, Nx, Ny, Nz, x, 0L, z, Line1, Ny);
				if (*Status == ERROR) {
					FreeLineDouble(&Line2);
					FreeLineDouble(&Line1);
					FreeLineDouble(&HalfKernel);
/**/				DEBUG_WRITE_LEAVING(VolumetricGradient, "Done")
					return(*Status);
				}
				ChangeBasis(Line1, Line2, Ny, CardinalSpline, BasicSpline,
					Degree, Convention, Tolerance, Status);
				if (*Status == ERROR) {
					FreeLineDouble(&Line2);
					FreeLineDouble(&Line1);
					FreeLineDouble(&HalfKernel);
/**/				DEBUG_WRITE_LEAVING(VolumetricGradient, "Done")
					return(*Status);
				}
				*Status = FirConvolveAntiSymmetric(Line2, Line1, Ny,
					HalfKernel, HalfDegree + 1L, Convention);
				if (*Status == ERROR) {
					FreeLineDouble(&Line2);
					FreeLineDouble(&Line1);
					FreeLineDouble(&HalfKernel);
/**/				DEBUG_WRITE_LEAVING(VolumetricGradient, "Done")
					return(*Status);
				}
				*Status = PutyDoubleToFloat(OutputGradient[1], Nx, Ny, Nz, x, 0L, z, Line1, Ny);
				if (*Status == ERROR) {
					FreeLineDouble(&Line2);
					FreeLineDouble(&Line1);
					FreeLineDouble(&HalfKernel);
/**/				DEBUG_WRITE_LEAVING(VolumetricGradient, "Done")
					return(*Status);
				}
			}
		}
		*Status = FreeLineDouble(&Line2);
		if (*Status == ERROR) {
			FreeLineDouble(&Line1);
			FreeLineDouble(&HalfKernel);
/**/		DEBUG_WRITE_LEAVING(VolumetricGradient, "Done")
			return(*Status);
		}
		*Status = FreeLineDouble(&Line1);
		if (*Status == ERROR) {
			FreeLineDouble(&HalfKernel);
/**/		DEBUG_WRITE_LEAVING(VolumetricGradient, "Done")
			return(*Status);
		}
	}
	else {
		p = OutputGradient[1];
		for (z = Nz; (0L < z); z--) {
			for (x = Nx; (0L < x); x--) {
				*p++ = 0.0F;
			}
		}
	}
	if (1L < Nz) {
		AllocateLineDouble(&Line1, Nz, Status);
		if (*Status == ERROR) {
			FreeLineDouble(&HalfKernel);
/**/		DEBUG_WRITE_LEAVING(VolumetricGradient, "Done")
			return(*Status);
		}
		AllocateLineDouble(&Line2, Nz, Status);
		if (*Status == ERROR) {
			FreeLineDouble(&Line1);
			FreeLineDouble(&HalfKernel);
/**/		DEBUG_WRITE_LEAVING(VolumetricGradient, "Done")
			return(*Status);
		}
		for (y = 0L; (y < Ny); y++) {
			for (x = 0L; (x < Nx); x++) {
				*Status = GetzFloatToDouble(InputVolume, Nx, Ny, Nz, x, y, 0L, Line1, Nz);
				if (*Status == ERROR) {
					FreeLineDouble(&Line2);
					FreeLineDouble(&Line1);
					FreeLineDouble(&HalfKernel);
/**/				DEBUG_WRITE_LEAVING(VolumetricGradient, "Done")
					return(*Status);
				}
				ChangeBasis(Line1, Line2, Nz, CardinalSpline, BasicSpline,
					Degree, Convention, Tolerance, Status);
				if (*Status == ERROR) {
					FreeLineDouble(&Line2);
					FreeLineDouble(&Line1);
					FreeLineDouble(&HalfKernel);
/**/				DEBUG_WRITE_LEAVING(VolumetricGradient, "Done")
					return(*Status);
				}
				*Status = FirConvolveAntiSymmetric(Line2, Line1, Nz,
					HalfKernel, HalfDegree + 1L, Convention);
				if (*Status == ERROR) {
					FreeLineDouble(&Line2);
					FreeLineDouble(&Line1);
					FreeLineDouble(&HalfKernel);
/**/				DEBUG_WRITE_LEAVING(VolumetricGradient, "Done")
					return(*Status);
				}
				*Status = PutzDoubleToFloat(OutputGradient[2], Nx, Ny, Nz, x, y, 0L, Line1, Nz);
				if (*Status == ERROR) {
					FreeLineDouble(&Line2);
					FreeLineDouble(&Line1);
					FreeLineDouble(&HalfKernel);
/**/				DEBUG_WRITE_LEAVING(VolumetricGradient, "Done")
					return(*Status);
				}
			}
		}
		*Status = FreeLineDouble(&Line2);
		if (*Status == ERROR) {
			FreeLineDouble(&Line1);
			FreeLineDouble(&HalfKernel);
/**/		DEBUG_WRITE_LEAVING(VolumetricGradient, "Done")
			return(*Status);
		}
		*Status = FreeLineDouble(&Line1);
		if (*Status == ERROR) {
			FreeLineDouble(&HalfKernel);
/**/		DEBUG_WRITE_LEAVING(VolumetricGradient, "Done")
			return(*Status);
		}
	}
	else {
		p = OutputGradient[2];
		for (y = Ny; (0L < y); y--) {
			for (x = Nx; (0L < x); x--) {
				*p++ = 0.0F;
			}
		}
	}
	*Status = FreeLineDouble(&HalfKernel);
/**/DEBUG_WRITE_LEAVING(VolumetricGradient, "Done")
	return(*Status);
} /* end VolumetricGradient */

