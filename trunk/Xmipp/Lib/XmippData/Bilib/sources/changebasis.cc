/*****************************************************************************
 *	System includes
 ****************************************************************************/
#include	<stdio.h>
#include	<string.h>

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
#include	"firconvolve.h"
#include	"getpoles.h"
#include	"getput.h"
#include	"iirconvolve.h"
#include	"kernel.h"
#include	"messagedisplay.h"

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
extern int		ChangeBasis
				(
					double	InputData[],		/* data to process */
					double	OutputData[],		/* result */
					long	SignalLength,		/* signal length */
					enum TSplineBasis
							FromBasis,			/* input basis */
					enum TSplineBasis
							ToBasis,			/* output basis */
					long	Degree,				/* degree of the representation space */
					enum TBoundaryConvention
							Convention,			/* boundary convention */
					double	Tolerance,			/* admissible relative error */
					int		*Status				/* error management */
				)

/* change spline coefficients from a source basis into a destination basis */
/* InputData is a (double)vector array of size SignalLength */
/* OutputData is a (double)vector array of size SignalLength */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin ChangeBasis */

	double	*p;
	double	*HalfKernel = (double *)NULL, *RealPoles = (double *)NULL, *Buffer = (double *)NULL;
	long	KernelHalfLength, PoleNumber;
	long	i;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(ChangeBasis, InputData, *Status,
/**/	"No InputData")
/**/DEBUG_CHECK_NULL_POINTER(ChangeBasis, OutputData, *Status,
/**/	"No OutputData")
/**/DEBUG_CHECK_RANGE_LONG(ChangeBasis, SignalLength, 1L, LONG_MAX, *Status,
/**/	"Invalid length (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(ChangeBasis, Degree, 0L, LONG_MAX, *Status,
/**/	"Invalid degree (should be positive)")
/**/DEBUG_CHECK_RANGE_DOUBLE(ChangeBasis, Tolerance, 0.0, DBL_MAX, *Status,
/**/	"Invalid Tolerance (should be positive)")
/**/DEBUG_RETURN_ON_ERROR(ChangeBasis, *Status)
/**/DEBUG_WRITE_ENTERING(ChangeBasis,
/**/	"About to perform a change of basis for spline coefficients")

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
			WRITE_ERROR(ChangeBasis, "Invalid boundary convention")
/**/		DEBUG_WRITE_LEAVING(ChangeBasis, "Done")
			return(*Status);
	}
	switch (FromBasis) {
		case BasicSpline:
			switch (ToBasis) {
				case BasicSpline:
					OutputData = (double *)memcpy(OutputData, InputData,
						(size_t)(SignalLength * (long)sizeof(double)));
					break;
				case CardinalSpline:
					KernelHalfLength = Degree / 2L + 1L;
					AllocateLineDouble(&HalfKernel, KernelHalfLength, Status);
					if (*Status == ERROR) {
/**/					DEBUG_WRITE_LEAVING(ChangeBasis, "Done")
						return(*Status);
					}
					p = HalfKernel;
					for (i = 0L; (i < KernelHalfLength); i++) {
						*Status = Bspline(Degree, (double)i, p++);
						if (*Status == ERROR) {
							FreeLineDouble(&HalfKernel);
/**/						DEBUG_WRITE_LEAVING(ChangeBasis, "Done")
							return(*Status);
						}
					}
					*Status = FirConvolveSymmetric(InputData, OutputData, SignalLength,
						HalfKernel, KernelHalfLength, Convention);
					if (*Status == ERROR) {
						FreeLineDouble(&HalfKernel);
/**/					DEBUG_WRITE_LEAVING(ChangeBasis, "Done")
						return(*Status);
					}
					*Status = FreeLineDouble(&HalfKernel);
					break;
				case DualSpline:
					ChangeBasis(InputData, OutputData, SignalLength,
						BasicSpline, CardinalSpline, 2L * Degree + 1L,
						Convention, Tolerance, Status);
					break;
				case OrthogonalSpline:
					*Status = ERROR;
					WRITE_ERROR(ChangeBasis, "Not yet implemented")
					break;
				default:
					*Status = ERROR;
					WRITE_ERROR(ChangeBasis, "Invalid destination basis")
					break;
			}
			break;
		case CardinalSpline:
			switch (ToBasis) {
				case BasicSpline:
					switch (Degree) {
						case 0L:
						case 1L:
							OutputData = (double *)memcpy(OutputData, InputData,
								(size_t)(SignalLength * (long)sizeof(double)));
							break;
						default:
							PoleNumber = Degree / 2L;
							AllocateLineDouble(&RealPoles, PoleNumber, Status);
							if (*Status == ERROR) {
/**/							DEBUG_WRITE_LEAVING(ChangeBasis, "Done")
								return(*Status);
							}
							GetBsplinePoles(RealPoles, Degree, Tolerance, Status);
							if (*Status == ERROR) {
								FreeLineDouble(&RealPoles);
/**/							DEBUG_WRITE_LEAVING(ChangeBasis, "Done")
								return(*Status);
							}
							*Status = IirConvolvePoles(InputData, OutputData, SignalLength,
								RealPoles, PoleNumber, Convention, Tolerance);
							if (*Status == ERROR) {
								FreeLineDouble(&RealPoles);
/**/							DEBUG_WRITE_LEAVING(ChangeBasis, "Done")
								return(*Status);
							}
							*Status = FreeLineDouble(&RealPoles);
							break;
					}
					break;
				case CardinalSpline:
					OutputData = (double *)memcpy(OutputData, InputData,
						(size_t)(SignalLength * (long)sizeof(double)));
					break;
				case DualSpline:
					AllocateLineDouble(&Buffer, SignalLength, Status);
					if (*Status == ERROR) {
/**/					DEBUG_WRITE_LEAVING(ChangeBasis, "Done")
						return(*Status);
					}
					ChangeBasis(InputData, Buffer, SignalLength,
						CardinalSpline, BasicSpline, Degree,
						Convention, Tolerance, Status);
					if (*Status == ERROR) {
						FreeLineDouble(&Buffer);
/**/					DEBUG_WRITE_LEAVING(ChangeBasis, "Done")
						return(*Status);
					}
					ChangeBasis(Buffer, OutputData, SignalLength,
						BasicSpline, CardinalSpline, 2L * Degree + 1L,
						Convention, Tolerance, Status);
					if (*Status == ERROR) {
						FreeLineDouble(&Buffer);
/**/					DEBUG_WRITE_LEAVING(ChangeBasis, "Done")
						return(*Status);
					}
					*Status = FreeLineDouble(&Buffer);
					break;
				case OrthogonalSpline:
					AllocateLineDouble(&Buffer, SignalLength, Status);
					if (*Status == ERROR) {
/**/					DEBUG_WRITE_LEAVING(ChangeBasis, "Done")
						return(*Status);
					}
					ChangeBasis(InputData, Buffer, SignalLength,
						CardinalSpline, BasicSpline, Degree,
						Convention, Tolerance, Status);
					if (*Status == ERROR) {
						FreeLineDouble(&Buffer);
/**/					DEBUG_WRITE_LEAVING(ChangeBasis, "Done")
						return(*Status);
					}
					ChangeBasis(Buffer, OutputData, SignalLength,
						BasicSpline, OrthogonalSpline, Degree,
						Convention, Tolerance, Status);
					if (*Status == ERROR) {
						FreeLineDouble(&Buffer);
/**/					DEBUG_WRITE_LEAVING(ChangeBasis, "Done")
						return(*Status);
					}
					*Status = FreeLineDouble(&Buffer);
					break;
				default:
					*Status = ERROR;
					WRITE_ERROR(ChangeBasis, "Invalid destination basis")
					break;
			}
			break;
		case DualSpline:
			switch (ToBasis) {
				case BasicSpline:
					ChangeBasis(InputData, OutputData, SignalLength,
						CardinalSpline, BasicSpline, 2L * Degree + 1L,
						Convention, Tolerance, Status);
					break;
				case CardinalSpline:
					AllocateLineDouble(&Buffer, SignalLength, Status);
					if (*Status == ERROR) {
/**/					DEBUG_WRITE_LEAVING(ChangeBasis, "Done")
						return(*Status);
					}
					ChangeBasis(InputData, Buffer, SignalLength,
						CardinalSpline, BasicSpline, 2L * Degree + 1L,
						Convention, Tolerance, Status);
					if (*Status == ERROR) {
						FreeLineDouble(&Buffer);
/**/					DEBUG_WRITE_LEAVING(ChangeBasis, "Done")
						return(*Status);
					}
					ChangeBasis(Buffer, OutputData, SignalLength,
						BasicSpline, CardinalSpline, Degree,
						Convention, Tolerance, Status);
					if (*Status == ERROR) {
						FreeLineDouble(&Buffer);
/**/					DEBUG_WRITE_LEAVING(ChangeBasis, "Done")
						return(*Status);
					}
					*Status = FreeLineDouble(&Buffer);
					break;
				case DualSpline:
					OutputData = (double *)memcpy(OutputData, InputData,
						(size_t)(SignalLength * (long)sizeof(double)));
					break;
				case OrthogonalSpline:
					ChangeBasis(InputData, OutputData, SignalLength,
						OrthogonalSpline, BasicSpline, Degree,
						Convention, Tolerance, Status);
					break;
				default:
					*Status = ERROR;
					WRITE_ERROR(ChangeBasis, "Invalid destination basis")
					break;
			}
			break;
		case OrthogonalSpline:
			switch (ToBasis) {
				case BasicSpline:
					*Status = ERROR;
					WRITE_ERROR(ChangeBasis, "Not yet implemented")
					break;
				case CardinalSpline:
					AllocateLineDouble(&Buffer, SignalLength, Status);
					if (*Status == ERROR) {
/**/					DEBUG_WRITE_LEAVING(ChangeBasis, "Done")
						return(*Status);
					}
					ChangeBasis(InputData, Buffer, SignalLength,
						OrthogonalSpline, BasicSpline, Degree,
						Convention, Tolerance, Status);
					if (*Status == ERROR) {
						FreeLineDouble(&Buffer);
/**/					DEBUG_WRITE_LEAVING(ChangeBasis, "Done")
						return(*Status);
					}
					ChangeBasis(Buffer, OutputData, SignalLength,
						BasicSpline, CardinalSpline, Degree,
						Convention, Tolerance, Status);
					break;
				case DualSpline:
					ChangeBasis(InputData, OutputData, SignalLength,
						BasicSpline, OrthogonalSpline, Degree,
						Convention, Tolerance, Status);
					break;
				case OrthogonalSpline:
					OutputData = (double *)memcpy(OutputData, InputData,
						(size_t)(SignalLength * (long)sizeof(double)));
					break;
				default:
					*Status = ERROR;
					WRITE_ERROR(ChangeBasis, "Invalid destination basis")
					break;
			}
			break;
		default:
			*Status = ERROR;
			WRITE_ERROR(ChangeBasis, "Invalid source basis")
			break;
	}
/**/DEBUG_WRITE_LEAVING(ChangeBasis, "Done")
	return(*Status);
} /* end ChangeBasis */

/*--------------------------------------------------------------------------*/
extern int		ChangeBasisVolume
				(
					float	*VolumeSource,		/* data to process */
					float	*VolumeDestination,	/* result */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					enum TSplineBasis
							FromBasis,			/* input basis */
					enum TSplineBasis
							ToBasis,			/* output basis */
					long	Degree,				/* degree of the representation space */
					enum TBoundaryConvention
							Convention,			/* boundary convention */
					double	Tolerance,			/* admissible relative error */
					int		*Status				/* error management */
				)

/* change a volume of spline coefficients from a source basis into a destination basis */
/* VolumeSource is a (float)volume of size (Nx x Ny x Nz) */
/* OutputData is a (float)volume of size (Nx x Ny x Nz) */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin ChangeBasisVolume */

	double	*p;
	double	*HalfKernel = (double *)NULL, *RealPoles = (double *)NULL;
	float	*Buffer;
	long	KernelHalfLength, PoleNumber;
	long	i;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(ChangeBasisVolume, VolumeSource, *Status,
/**/	"No VolumeSource")
/**/DEBUG_CHECK_NULL_POINTER(ChangeBasisVolume, VolumeDestination, *Status,
/**/	"No OutputData")
/**/DEBUG_CHECK_RANGE_LONG(ChangeBasisVolume, Nx, 1L, LONG_MAX, *Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(ChangeBasisVolume, Ny, 1L, LONG_MAX, *Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(ChangeBasisVolume, Nz, 1L, LONG_MAX, *Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(ChangeBasisVolume, Degree, 0L, LONG_MAX, *Status,
/**/	"Invalid degree (should be positive)")
/**/DEBUG_CHECK_RANGE_DOUBLE(ChangeBasisVolume, Tolerance, 0.0, DBL_MAX, *Status,
/**/	"Invalid Tolerance (should be positive)")
/**/DEBUG_RETURN_ON_ERROR(ChangeBasisVolume, *Status)
/**/DEBUG_WRITE_ENTERING(ChangeBasisVolume,
/**/	"About to perform a change of basis for a volume")

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
			WRITE_ERROR(ChangeBasisVolume, "Invalid boundary convention")
/**/		DEBUG_WRITE_LEAVING(ChangeBasisVolume, "Done")
			return(*Status);
	}
	switch (FromBasis) {
		case BasicSpline:
			switch (ToBasis) {
				case BasicSpline:
					VolumeDestination = (float *)memcpy(VolumeDestination, VolumeSource,
						(size_t)(Nx * Ny * Nz * (long)sizeof(float)));
					break;
				case CardinalSpline:
					KernelHalfLength = Degree / 2L + 1L;
					AllocateLineDouble(&HalfKernel, KernelHalfLength, Status);
					if (*Status == ERROR) {
/**/					DEBUG_WRITE_LEAVING(ChangeBasisVolume, "Done")
						return(*Status);
					}
					p = HalfKernel;
					for (i = 0L; (i < KernelHalfLength); i++) {
						*Status = Bspline(Degree, (double)i, p++);
						if (*Status == ERROR) {
							FreeLineDouble(&HalfKernel);
/**/						DEBUG_WRITE_LEAVING(ChangeBasisVolume, "Done")
							return(*Status);
						}
					}
					*Status = FirConvolveSymmetricVolume(VolumeSource, VolumeDestination,
						Nx, Ny, Nz, HalfKernel, KernelHalfLength, Convention);
					if (*Status == ERROR) {
						FreeLineDouble(&HalfKernel);
/**/					DEBUG_WRITE_LEAVING(ChangeBasisVolume, "Done")
						return(*Status);
					}
					*Status = FreeLineDouble(&HalfKernel);
					break;
				case DualSpline:
					ChangeBasisVolume(VolumeSource, VolumeDestination, Nx, Ny, Nz,
						BasicSpline, CardinalSpline, 2L * Degree + 1L,
						Convention, Tolerance, Status);
					break;
				case OrthogonalSpline:
					*Status = ERROR;
					WRITE_ERROR(ChangeBasisVolume, "Not yet implemented")
					break;
				default:
					*Status = ERROR;
					WRITE_ERROR(ChangeBasisVolume, "Invalid destination basis")
					break;
			}
			break;
		case CardinalSpline:
			switch (ToBasis) {
				case BasicSpline:
					switch (Degree) {
						case 0L:
						case 1L:
							VolumeDestination = (float *)memcpy(VolumeDestination, VolumeSource,
								(size_t)(Nx * Ny * Nz * (long)sizeof(float)));
							break;
						default:
							PoleNumber = Degree / 2L;
							AllocateLineDouble(&RealPoles, PoleNumber, Status);
							if (*Status == ERROR) {
/**/							DEBUG_WRITE_LEAVING(ChangeBasisVolume, "Done")
								return(*Status);
							}
							GetBsplinePoles(RealPoles, Degree, Tolerance, Status);
							if (*Status == ERROR) {
								FreeLineDouble(&RealPoles);
/**/							DEBUG_WRITE_LEAVING(ChangeBasisVolume, "Done")
								return(*Status);
							}
							IirConvolvePolesVolume(VolumeSource, VolumeDestination,
								Nx, Ny, Nz, RealPoles, PoleNumber, Convention, Tolerance, Status);
							if (*Status == ERROR) {
								FreeLineDouble(&RealPoles);
/**/							DEBUG_WRITE_LEAVING(ChangeBasisVolume, "Done")
								return(*Status);
							}
							*Status = FreeLineDouble(&RealPoles);
							break;
					}
					break;
				case CardinalSpline:
					VolumeDestination = (float *)memcpy(VolumeDestination, VolumeSource,
						(size_t)(Nx * Ny * Nz * (long)sizeof(float)));
					break;
				case DualSpline:
					AllocateVolumeFloat(&Buffer, Nx, Ny, Nz, Status);
					if (*Status == ERROR) {
/**/					DEBUG_WRITE_LEAVING(ChangeBasisVolume, "Done")
						return(*Status);
					}
					ChangeBasisVolume(VolumeSource, Buffer, Nx, Ny, Nz,
						CardinalSpline, BasicSpline, Degree,
						Convention, Tolerance, Status);
					if (*Status == ERROR) {
						FreeVolumeFloat(&Buffer);
/**/					DEBUG_WRITE_LEAVING(ChangeBasisVolume, "Done")
						return(*Status);
					}
					ChangeBasisVolume(Buffer, VolumeDestination, Nx, Ny, Nz,
						BasicSpline, CardinalSpline, 2L * Degree + 1L,
						Convention, Tolerance, Status);
					if (*Status == ERROR) {
						FreeVolumeFloat(&Buffer);
/**/					DEBUG_WRITE_LEAVING(ChangeBasisVolume, "Done")
						return(*Status);
					}
					*Status = FreeVolumeFloat(&Buffer);
					break;
				case OrthogonalSpline:
					AllocateVolumeFloat(&Buffer, Nx, Ny, Nz, Status);
					if (*Status == ERROR) {
/**/					DEBUG_WRITE_LEAVING(ChangeBasisVolume, "Done")
						return(*Status);
					}
					ChangeBasisVolume(VolumeSource, Buffer, Nx, Ny, Nz,
						CardinalSpline, BasicSpline, Degree,
						Convention, Tolerance, Status);
					if (*Status == ERROR) {
						FreeVolumeFloat(&Buffer);
/**/					DEBUG_WRITE_LEAVING(ChangeBasisVolume, "Done")
						return(*Status);
					}
					ChangeBasisVolume(Buffer, VolumeDestination, Nx, Ny, Nz,
						BasicSpline, OrthogonalSpline, Degree,
						Convention, Tolerance, Status);
					if (*Status == ERROR) {
						FreeVolumeFloat(&Buffer);
/**/					DEBUG_WRITE_LEAVING(ChangeBasisVolume, "Done")
						return(*Status);
					}
					*Status = FreeVolumeFloat(&Buffer);
					break;
				default:
					*Status = ERROR;
					WRITE_ERROR(ChangeBasisVolume, "Invalid destination basis")
					break;
			}
			break;
		case DualSpline:
			switch (ToBasis) {
				case BasicSpline:
					ChangeBasisVolume(VolumeSource, VolumeDestination, Nx, Ny, Nz,
						CardinalSpline, BasicSpline, 2L * Degree + 1L,
						Convention, Tolerance, Status);
					break;
				case CardinalSpline:
					AllocateVolumeFloat(&Buffer, Nx, Ny, Nz, Status);
					if (*Status == ERROR) {
/**/					DEBUG_WRITE_LEAVING(ChangeBasisVolume, "Done")
						return(*Status);
					}
					ChangeBasisVolume(VolumeSource, Buffer, Nx, Ny, Nz,
						CardinalSpline, BasicSpline, 2L * Degree + 1L,
						Convention, Tolerance, Status);
					if (*Status == ERROR) {
						FreeVolumeFloat(&Buffer);
/**/					DEBUG_WRITE_LEAVING(ChangeBasisVolume, "Done")
						return(*Status);
					}
					ChangeBasisVolume(Buffer, VolumeDestination, Nx, Ny, Nz,
						BasicSpline, CardinalSpline, Degree,
						Convention, Tolerance, Status);
					if (*Status == ERROR) {
						FreeVolumeFloat(&Buffer);
/**/					DEBUG_WRITE_LEAVING(ChangeBasisVolume, "Done")
						return(*Status);
					}
					*Status = FreeVolumeFloat(&Buffer);
					break;
				case DualSpline:
					VolumeDestination = (float *)memcpy(VolumeDestination, VolumeSource,
						(size_t)(Nx * Ny * Nz * (long)sizeof(float)));
					break;
				case OrthogonalSpline:
					ChangeBasisVolume(VolumeSource, VolumeDestination, Nx, Ny, Nz,
						OrthogonalSpline, BasicSpline, Degree,
						Convention, Tolerance, Status);
					break;
				default:
					*Status = ERROR;
					WRITE_ERROR(ChangeBasisVolume, "Invalid destination basis")
					break;
			}
			break;
		case OrthogonalSpline:
			switch (ToBasis) {
				case BasicSpline:
					*Status = ERROR;
					WRITE_ERROR(ChangeBasisVolume, "Not yet implemented")
					break;
				case CardinalSpline:
					AllocateVolumeFloat(&Buffer, Nx, Ny, Nz, Status);
					if (*Status == ERROR) {
/**/					DEBUG_WRITE_LEAVING(ChangeBasisVolume, "Done")
						return(*Status);
					}
					ChangeBasisVolume(VolumeSource, Buffer, Nx, Ny, Nz,
						OrthogonalSpline, BasicSpline, Degree,
						Convention, Tolerance, Status);
					if (*Status == ERROR) {
						FreeVolumeFloat(&Buffer);
/**/					DEBUG_WRITE_LEAVING(ChangeBasisVolume, "Done")
						return(*Status);
					}
					ChangeBasisVolume(Buffer, VolumeDestination, Nx, Ny, Nz,
						BasicSpline, CardinalSpline, Degree,
						Convention, Tolerance, Status);
					break;
				case DualSpline:
					ChangeBasisVolume(VolumeSource, VolumeDestination, Nx, Ny, Nz,
						BasicSpline, OrthogonalSpline, Degree,
						Convention, Tolerance, Status);
					break;
				case OrthogonalSpline:
					VolumeDestination = (float *)memcpy(VolumeDestination, VolumeSource,
						(size_t)(Nx * Ny * Nz * (long)sizeof(float)));
					break;
				default:
					*Status = ERROR;
					WRITE_ERROR(ChangeBasisVolume, "Invalid destination basis")
					break;
			}
			break;
		default:
			*Status = ERROR;
			WRITE_ERROR(ChangeBasisVolume, "Invalid source basis")
			break;
	}
/**/DEBUG_WRITE_LEAVING(ChangeBasisVolume, "Done")
	return(*Status);
} /* end ChangeBasisVolume */

