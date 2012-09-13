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

/*****************************************************************************
 *	Toolbox includes
 ****************************************************************************/
#include	"../headers/convert.h"
#include	"fold.h"
#include	"getput.h"
#include	"messagedisplay.h"
#include	"morphology.h"

/*****************************************************************************
 *	Conditional includes
 ****************************************************************************/
#ifdef		DEBUG
#include	<limits.h>
#endif

/*****************************************************************************
 *	Declaration of static procedures
 ****************************************************************************/
/*--------------------------------------------------------------------------*/
static int		MaxFilterLine
				(
					double	*BufferIn,
					double	*BufferOut,
					long	Nx,
					long	Kx,
					enum TBoundaryConvention
							Convention
				);

/*--------------------------------------------------------------------------*/
static int		MinFilterLine
				(
					double	*BufferIn,
					double	*BufferOut,
					long	Nx,
					long	Kx,
					enum TBoundaryConvention
							Convention
				);

/*****************************************************************************
 *	Definition of static procedures
 ****************************************************************************/
/*--------------------------------------------------------------------------*/
static int		MaxFilterLine
				(
					double	*BufferIn,
					double	*BufferOut,
					long	N,
					long	K,
					enum TBoundaryConvention
							Convention
				)

{ /* begin MaxFilterLine */

	double	f, fmax;
	long	Ok = (K - 1L) / 2L;
	long	i, n, k;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(MaxFilterLine, BufferIn, Status,
/**/	"Missing BufferIn")
/**/DEBUG_CHECK_NULL_POINTER(MaxFilterLine, BufferOut, Status,
/**/	"Missing BufferOut")
/**/DEBUG_CHECK_RANGE_LONG(MaxFilterLine, N, 1L, LONG_MAX, Status,
/**/	"Invalid data width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MaxFilterLine, K, 1L, LONG_MAX, Status,
/**/	"Invalid kernel width (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(MaxFilterLine, Status)
/**/DEBUG_WRITE_ENTERING(MaxFilterLine,
/**/	"About to execute MaxFilterLine")

	switch (Convention) {
		case FiniteDataSupport:
		case MirrorOffBounds:
		case MirrorOnBounds:
		case Periodic:
			break;
		case AntiMirrorOnBounds:
		case FiniteCoefficientSupport:
		default:
			Status = ERROR;
			WRITE_ERROR(MaxFilterLine, "Invalid boundary convention")
/**/		DEBUG_WRITE_LEAVING(MaxFilterLine, "Done")
			return(Status);
	}
	for (n = Ok; (n < (N + Ok)); n++) {
		Status = GetFoldedIndex(n, &i, N, Convention);
		if (Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(MaxFilterLine, "Done")
			return(Status);
		}
		fmax = BufferIn[i];
		for (k = n + 1L - K; (k < n); k++) {
			Status = GetFoldedIndex(k, &i, N, Convention);
			if (Status == ERROR) {
/**/			DEBUG_WRITE_LEAVING(MaxFilterLine, "Done")
				return(Status);
			}
			f = BufferIn[i];
			fmax = (f < fmax) ? (fmax) : (f);
		}
		*BufferOut++ = fmax;
	}
/**/DEBUG_WRITE_LEAVING(MaxFilterLine, "Done")
	return(Status);
} /* end MaxFilterLine */

/*--------------------------------------------------------------------------*/
static int		MinFilterLine
				(
					double	*BufferIn,
					double	*BufferOut,
					long	N,
					long	K,
					enum TBoundaryConvention
							Convention
				)

{ /* begin MinFilterLine */

	double	f, fmin;
	long	Ok = (K - 1L) / 2L;
	long	i, n, k;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(MinFilterLine, BufferIn, Status,
/**/	"Missing BufferIn")
/**/DEBUG_CHECK_NULL_POINTER(MinFilterLine, BufferOut, Status,
/**/	"Missing BufferOut")
/**/DEBUG_CHECK_RANGE_LONG(MinFilterLine, N, 1L, LONG_MAX, Status,
/**/	"Invalid data width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MinFilterLine, K, 1L, LONG_MAX, Status,
/**/	"Invalid kernel width (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(MinFilterLine, Status)
/**/DEBUG_WRITE_ENTERING(MinFilterLine,
/**/	"About to execute MinFilterLine")

	switch (Convention) {
		case FiniteDataSupport:
		case MirrorOffBounds:
		case MirrorOnBounds:
		case Periodic:
			break;
		case AntiMirrorOnBounds:
		case FiniteCoefficientSupport:
		default:
			Status = ERROR;
			WRITE_ERROR(MinFilterLine, "Invalid boundary convention")
/**/		DEBUG_WRITE_LEAVING(MinFilterLine, "Done")
			return(Status);
	}
	for (n = -Ok; (n < (N - Ok)); n++) {
		Status = GetFoldedIndex(n, &i, N, Convention);
		if (Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(MinFilterLine, "Done")
			return(Status);
		}
		fmin = BufferIn[i];
		for (k = n + 1L; (k < (n + K)); k++) {
			Status = GetFoldedIndex(k, &i, N, Convention);
			if (Status == ERROR) {
/**/			DEBUG_WRITE_LEAVING(MinFilterLine, "Done")
				return(Status);
			}
			f = BufferIn[i];
			fmin = (f < fmin) ? (f) : (fmin);
		}
		*BufferOut++ = fmin;
	}
/**/DEBUG_WRITE_LEAVING(MinFilterLine, "Done")
	return(Status);
} /* end MinFilterLine */

/*****************************************************************************
 *	Definition of extern procedures
 ****************************************************************************/
/*--------------------------------------------------------------------------*/
extern int		BrightTopHatFloat
				(
					float	*VolumeSource,		/* data to process */
					float	*VolumeDestination,	/* result */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					float	*Kernel,			/* structural element */
					long	Kx,					/* width of the kernel */
					long	Ky,					/* height of the kernel */
					long	Kz,					/* depth of the kernel */
					long	Ox,					/* kernel X origin */
					long	Oy,					/* kernel Y origin */
					long	Oz,					/* kernel Z origin */
					enum TBoundaryConvention
							Convention,			/* boundary convention */
					int		*Status				/* error management */
				)

/* perform the grey-scale morphological operation known as bright top hat */
/* VolumeSource is the grey-scale volume to process */
/* VolumeDestination is the resulting grey-scale volume */
/* both VolumeSource and VolumeDestination have size (Nx, Ny, Nz) */
/* Kernel is the grey-scale structural element by which VolumeSource is processed */
/* (Kx, Ky, Kz) is the size (bounding box) of the structural element */
/* (Ox, Oy, Oz) is the origin of the structural element */
/* Convention is the boundary convention applied to VolumeSource before processing */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin BrightTopHatFloat */

	float	*p;
	float	*Opened = (float *)NULL;
	long	x, y, z;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(BrightTopHatFloat, VolumeSource, *Status,
/**/	"Missing VolumeSource")
/**/DEBUG_CHECK_NULL_POINTER(BrightTopHatFloat, VolumeDestination, *Status,
/**/	"Missing VolumeDestination")
/**/DEBUG_CHECK_NULL_POINTER(BrightTopHatFloat, Kernel, *Status,
/**/	"Missing Kernel")
/**/DEBUG_CHECK_RANGE_LONG(BrightTopHatFloat, Nx, 1L, LONG_MAX, *Status,
/**/	"Invalid data width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(BrightTopHatFloat, Ny, 1L, LONG_MAX, *Status,
/**/	"Invalid data height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(BrightTopHatFloat, Nz, 1L, LONG_MAX, *Status,
/**/	"Invalid data depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(BrightTopHatFloat, Kx, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(BrightTopHatFloat, Ky, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(BrightTopHatFloat, Kz, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(BrightTopHatFloat, *Status)
/**/DEBUG_WRITE_ENTERING(BrightTopHatFloat,
/**/	"About to execute BrightTopHatFloat")

	switch (Convention) {
		case FiniteDataSupport:
		case MirrorOffBounds:
		case MirrorOnBounds:
		case Periodic:
			break;
		case AntiMirrorOnBounds:
		case FiniteCoefficientSupport:
		default:
			*Status = ERROR;
			WRITE_ERROR(BrightTopHatFloat, "Invalid boundary convention")
/**/		DEBUG_WRITE_LEAVING(BrightTopHatFloat, "Done")
			return(*Status);
	}
	AllocateVolumeFloat(&Opened, Nx, Ny, Nz, Status);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(BrightTopHatFloat, "Done")
		return(*Status);
	}
	OpeningFloat(VolumeSource, Opened, Nx, Ny, Nz, Kernel, Kx, Ky, Kz,
		Ox, Oy, Oz, Convention, Status);
	if (*Status == ERROR) {
		FreeVolumeFloat(&Opened);
/**/	DEBUG_WRITE_LEAVING(BrightTopHatFloat, "Done")
		return(*Status);
	}
	p = Opened;
	for (z = -Nz; (z < 0L); z++) {
		for (y = -Ny; (y < 0L); y++) {
			for (x = -Nx; (x < 0L); x++) {
				*VolumeDestination++ = *VolumeSource++ - *p++;
			}
		}
	}
	*Status = FreeVolumeFloat(&Opened);
/**/DEBUG_WRITE_LEAVING(BrightTopHatFloat, "Done")
	return(*Status);
} /* end BrightTopHatFloat */

/*--------------------------------------------------------------------------*/
extern int		BrightTopHatFloatCuboid
				(
					float	*Volume,			/* data to process in-place */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					long	Kx,					/* width of the kernel */
					long	Ky,					/* height of the kernel */
					long	Kz,					/* depth of the kernel */
					enum TBoundaryConvention
							Convention,			/* boundary convention */
					int		*Status				/* error management */
				)

/* perform a morphological grey-scale bright top hat filter using a cuboid as structural element */
/* Volume is the grey-scale volume to process in-place */
/* Volume has size (Nx, Ny, Nz) */
/* (Kx, Ky, Kz) is the size of the structural element that is filled with the value 1.0F */
/* the origin of the structural element is ((Kx - 1L) / 2L, (Ky - 1L) / 2L, (Kz - 1L) / 2L) */
/* Convention is the boundary convention applied to Volume before processing */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin BrightTopHatFloatCuboid */

	float	*p;
	float	*Opened = (float *)NULL;
	long	x, y, z;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(BrightTopHatFloatCuboid, Volume, *Status,
/**/	"Missing Volume")
/**/DEBUG_CHECK_RANGE_LONG(BrightTopHatFloatCuboid, Nx, 1L, LONG_MAX, *Status,
/**/	"Invalid data width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(BrightTopHatFloatCuboid, Ny, 1L, LONG_MAX, *Status,
/**/	"Invalid data height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(BrightTopHatFloatCuboid, Nz, 1L, LONG_MAX, *Status,
/**/	"Invalid data depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(BrightTopHatFloatCuboid, Kx, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(BrightTopHatFloatCuboid, Ky, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(BrightTopHatFloatCuboid, Kz, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(BrightTopHatFloatCuboid, *Status)
/**/DEBUG_WRITE_ENTERING(BrightTopHatFloatCuboid,
/**/	"About to execute BrightTopHatFloatCuboid")

	switch (Convention) {
		case FiniteDataSupport:
		case MirrorOffBounds:
		case MirrorOnBounds:
		case Periodic:
			break;
		case AntiMirrorOnBounds:
		case FiniteCoefficientSupport:
		default:
			*Status = ERROR;
			WRITE_ERROR(BrightTopHatFloatCuboid, "Invalid boundary convention")
/**/		DEBUG_WRITE_LEAVING(BrightTopHatFloatCuboid, "Done")
			return(*Status);
	}
	AllocateVolumeFloat(&Opened, Nx, Ny, Nz, Status);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(BrightTopHatFloatCuboid, "Done")
		return(*Status);
	}
	memcpy(Opened, Volume, (size_t)(Nx * Ny * Nz * (long)sizeof(float)));
	OpeningFloatCuboid(Opened, Nx, Ny, Nz, Kx, Ky, Kz, Convention, Status);
	if (*Status == ERROR) {
		FreeVolumeFloat(&Opened);
/**/	DEBUG_WRITE_LEAVING(BrightTopHatFloatCuboid, "Done")
		return(*Status);
	}
	p = Opened;
	for (z = -Nz; (z < 0L); z++) {
		for (y = -Ny; (y < 0L); y++) {
			for (x = -Nx; (x < 0L); x++) {
				*Volume++ -= *p++;
			}
		}
	}
	*Status = FreeVolumeFloat(&Opened);
/**/DEBUG_WRITE_LEAVING(BrightTopHatFloatCuboid, "Done")
	return(*Status);
} /* end BrightTopHatFloatCuboid */

/*--------------------------------------------------------------------------*/
extern int		BrightTopHatShort
				(
					short	*VolumeSource,		/* data to process */
					short	*VolumeDestination,	/* result */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					short	*Kernel,			/* structural element */
					long	Kx,					/* width of the kernel */
					long	Ky,					/* height of the kernel */
					long	Kz,					/* depth of the kernel */
					long	Ox,					/* kernel X origin */
					long	Oy,					/* kernel Y origin */
					long	Oz,					/* kernel Z origin */
					enum TBoundaryConvention
							Convention,			/* boundary convention */
					int		*Status				/* error management */
				)

/* perform the grey-scale morphological operation known as bright top hat */
/* VolumeSource is the grey-scale volume to process */
/* VolumeDestination is the resulting grey-scale volume */
/* both VolumeSource and VolumeDestination have size (Nx, Ny, Nz) */
/* Kernel is the grey-scale structural element by which VolumeSource is processed */
/* (Kx, Ky, Kz) is the size (bounding box) of the structural element */
/* (Ox, Oy, Oz) is the origin of the structural element */
/* Convention is the boundary convention applied to VolumeSource before processing */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin BrightTopHatShort */

	short	*p;
	short	*Opened = (short *)NULL;
	long	x, y, z;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(BrightTopHatShort, VolumeSource, *Status,
/**/	"Missing VolumeSource")
/**/DEBUG_CHECK_NULL_POINTER(BrightTopHatShort, VolumeDestination, *Status,
/**/	"Missing VolumeDestination")
/**/DEBUG_CHECK_NULL_POINTER(BrightTopHatShort, Kernel, *Status,
/**/	"Missing Kernel")
/**/DEBUG_CHECK_RANGE_LONG(BrightTopHatShort, Nx, 1L, LONG_MAX, *Status,
/**/	"Invalid data width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(BrightTopHatShort, Ny, 1L, LONG_MAX, *Status,
/**/	"Invalid data height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(BrightTopHatShort, Nz, 1L, LONG_MAX, *Status,
/**/	"Invalid data depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(BrightTopHatShort, Kx, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(BrightTopHatShort, Ky, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(BrightTopHatShort, Kz, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(BrightTopHatShort, *Status)
/**/DEBUG_WRITE_ENTERING(BrightTopHatShort,
/**/	"About to execute BrightTopHatShort")

	switch (Convention) {
		case FiniteDataSupport:
		case MirrorOffBounds:
		case MirrorOnBounds:
		case Periodic:
			break;
		case AntiMirrorOnBounds:
		case FiniteCoefficientSupport:
		default:
			*Status = ERROR;
			WRITE_ERROR(BrightTopHatShort, "Invalid boundary convention")
/**/		DEBUG_WRITE_LEAVING(BrightTopHatShort, "Done")
			return(*Status);
	}
	AllocateVolumeShort(&Opened, Nx, Ny, Nz, Status);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(BrightTopHatShort, "Done")
		return(*Status);
	}
	OpeningShort(VolumeSource, Opened, Nx, Ny, Nz, Kernel, Kx, Ky, Kz,
		Ox, Oy, Oz, Convention, Status);
	if (*Status == ERROR) {
		FreeVolumeShort(&Opened);
/**/	DEBUG_WRITE_LEAVING(BrightTopHatShort, "Done")
		return(*Status);
	}
	p = Opened;
	for (z = -Nz; (z < 0L); z++) {
		for (y = -Ny; (y < 0L); y++) {
			for (x = -Nx; (x < 0L); x++) {
				*VolumeDestination++ = ConvertIntToShort((int)*VolumeSource++ - (int)*p++);
			}
		}
	}
	*Status = FreeVolumeShort(&Opened);
/**/DEBUG_WRITE_LEAVING(BrightTopHatShort, "Done")
	return(*Status);
} /* end BrightTopHatShort */

/*--------------------------------------------------------------------------*/
extern int		BrightTopHatShortCuboid
				(
					short	*Volume,			/* data to process in-place */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					long	Kx,					/* width of the kernel */
					long	Ky,					/* height of the kernel */
					long	Kz,					/* depth of the kernel */
					enum TBoundaryConvention
							Convention,			/* boundary convention */
					int		*Status				/* error management */
				)

/* perform a morphological grey-scale bright top hat filter using a cuboid as structural element */
/* Volume is the grey-scale volume to process in-place */
/* Volume has size (Nx, Ny, Nz) */
/* (Kx, Ky, Kz) is the size of the structural element that is filled with the value 1.0F */
/* the origin of the structural element is ((Kx - 1L) / 2L, (Ky - 1L) / 2L, (Kz - 1L) / 2L) */
/* Convention is the boundary convention applied to Volume before processing */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin BrightTopHatShortCuboid */

	short	*p;
	short	*Opened = (short *)NULL;
	long	x, y, z;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(BrightTopHatShortCuboid, Volume, *Status,
/**/	"Missing Volume")
/**/DEBUG_CHECK_RANGE_LONG(BrightTopHatShortCuboid, Nx, 1L, LONG_MAX, *Status,
/**/	"Invalid data width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(BrightTopHatShortCuboid, Ny, 1L, LONG_MAX, *Status,
/**/	"Invalid data height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(BrightTopHatShortCuboid, Nz, 1L, LONG_MAX, *Status,
/**/	"Invalid data depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(BrightTopHatShortCuboid, Kx, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(BrightTopHatShortCuboid, Ky, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(BrightTopHatShortCuboid, Kz, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(BrightTopHatShortCuboid, *Status)
/**/DEBUG_WRITE_ENTERING(BrightTopHatShortCuboid,
/**/	"About to execute BrightTopHatShortCuboid")

	switch (Convention) {
		case FiniteDataSupport:
		case MirrorOffBounds:
		case MirrorOnBounds:
		case Periodic:
			break;
		case AntiMirrorOnBounds:
		case FiniteCoefficientSupport:
		default:
			*Status = ERROR;
			WRITE_ERROR(BrightTopHatShortCuboid, "Invalid boundary convention")
/**/		DEBUG_WRITE_LEAVING(BrightTopHatShortCuboid, "Done")
			return(*Status);
	}
	AllocateVolumeShort(&Opened, Nx, Ny, Nz, Status);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(BrightTopHatShortCuboid, "Done")
		return(*Status);
	}
	memcpy(Opened, Volume, (size_t)(Nx * Ny * Nz * (long)sizeof(short)));
	OpeningShortCuboid(Opened, Nx, Ny, Nz, Kx, Ky, Kz, Convention, Status);
	if (*Status == ERROR) {
		FreeVolumeShort(&Opened);
/**/	DEBUG_WRITE_LEAVING(BrightTopHatShortCuboid, "Done")
		return(*Status);
	}
	p = Opened;
	for (z = -Nz; (z < 0L); z++) {
		for (y = -Ny; (y < 0L); y++) {
			for (x = -Nx; (x < 0L); x++) {
				*Volume++ -= *p++;
			}
		}
	}
	*Status = FreeVolumeShort(&Opened);
/**/DEBUG_WRITE_LEAVING(BrightTopHatShortCuboid, "Done")
	return(*Status);
} /* end BrightTopHatShortCuboid */

/*--------------------------------------------------------------------------*/
extern int		ClosingFloat
				(
					float	*VolumeSource,		/* data to process */
					float	*VolumeDestination,	/* result */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					float	*Kernel,			/* structural element */
					long	Kx,					/* width of the kernel */
					long	Ky,					/* height of the kernel */
					long	Kz,					/* depth of the kernel */
					long	Ox,					/* kernel X origin */
					long	Oy,					/* kernel Y origin */
					long	Oz,					/* kernel Z origin */
					enum TBoundaryConvention
							Convention,			/* boundary convention */
					int		*Status				/* error management */
				)

/* perform the grey-scale morphological operation known as closing */
/* VolumeSource is the grey-scale volume to close */
/* VolumeDestination is the resulting grey-scale volume after closing */
/* both VolumeSource and VolumeDestination have size (Nx, Ny, Nz) */
/* Kernel is the grey-scale structural element by which VolumeSource is closed */
/* (Kx, Ky, Kz) is the size (bounding box) of the structural element */
/* (Ox, Oy, Oz) is the origin of the structural element */
/* Convention is the boundary convention applied to VolumeSource before closing */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin ClosingFloat */

	float	*Buffer = (float *)NULL;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(ClosingFloat, VolumeSource, *Status,
/**/	"Missing VolumeSource")
/**/DEBUG_CHECK_NULL_POINTER(ClosingFloat, VolumeDestination, *Status,
/**/	"Missing VolumeDestination")
/**/DEBUG_CHECK_NULL_POINTER(ClosingFloat, Kernel, *Status,
/**/	"Missing Kernel")
/**/DEBUG_CHECK_RANGE_LONG(ClosingFloat, Nx, 1L, LONG_MAX, *Status,
/**/	"Invalid data width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(ClosingFloat, Ny, 1L, LONG_MAX, *Status,
/**/	"Invalid data height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(ClosingFloat, Nz, 1L, LONG_MAX, *Status,
/**/	"Invalid data depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(ClosingFloat, Kx, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(ClosingFloat, Ky, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(ClosingFloat, Kz, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(ClosingFloat, *Status)
/**/DEBUG_WRITE_ENTERING(ClosingFloat,
/**/	"About to execute ClosingFloat")

	switch (Convention) {
		case FiniteDataSupport:
		case MirrorOffBounds:
		case MirrorOnBounds:
		case Periodic:
			break;
		case AntiMirrorOnBounds:
		case FiniteCoefficientSupport:
		default:
			*Status = ERROR;
			WRITE_ERROR(ClosingFloat, "Invalid boundary convention")
/**/		DEBUG_WRITE_LEAVING(ClosingFloat, "Done")
			return(*Status);
	}
	AllocateVolumeFloat(&Buffer, Nx, Ny, Nz, Status);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(ClosingFloat, "Done")
		return(*Status);
	}
	*Status = DilateFloat(VolumeSource, Buffer, Nx, Ny, Nz, Kernel, Kx, Ky, Kz,
		Ox, Oy, Oz, Convention);
	if (*Status == ERROR) {
		FreeVolumeFloat(&Buffer);
/**/	DEBUG_WRITE_LEAVING(ClosingFloat, "Done")
		return(*Status);
	}
	*Status = ErodeFloat(Buffer, VolumeDestination, Nx, Ny, Nz, Kernel, Kx, Ky, Kz,
		Ox, Oy, Oz, Convention);
	if (*Status == ERROR) {
		FreeVolumeFloat(&Buffer);
/**/	DEBUG_WRITE_LEAVING(ClosingFloat, "Done")
		return(*Status);
	}
	*Status = FreeVolumeFloat(&Buffer);
/**/DEBUG_WRITE_LEAVING(ClosingFloat, "Done")
	return(*Status);
} /* end ClosingFloat */

/*--------------------------------------------------------------------------*/
extern int		ClosingFloatCuboid
				(
					float	*Volume,			/* data to process in-place */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					long	Kx,					/* width of the kernel */
					long	Ky,					/* height of the kernel */
					long	Kz,					/* depth of the kernel */
					enum TBoundaryConvention
							Convention,			/* boundary convention */
					int		*Status				/* error management */
				)

/* perform a morphological grey-scale closing filter using a cuboid as structural element */
/* Volume is the grey-scale volume to close in-place */
/* Volume has size (Nx, Ny, Nz) */
/* (Kx, Ky, Kz) is the size of the structural element that is filled with the value 1.0F */
/* the origin of the structural element is ((Kx - 1L) / 2L, (Ky - 1L) / 2L, (Kz - 1L) / 2L) */
/* Convention is the boundary convention applied to Volume before closing */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin ClosingFloatCuboid */

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(ClosingFloatCuboid, Volume, *Status,
/**/	"Missing Volume")
/**/DEBUG_CHECK_RANGE_LONG(ClosingFloatCuboid, Nx, 1L, LONG_MAX, *Status,
/**/	"Invalid data width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(ClosingFloatCuboid, Ny, 1L, LONG_MAX, *Status,
/**/	"Invalid data height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(ClosingFloatCuboid, Nz, 1L, LONG_MAX, *Status,
/**/	"Invalid data depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(ClosingFloatCuboid, Kx, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(ClosingFloatCuboid, Ky, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(ClosingFloatCuboid, Kz, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(ClosingFloatCuboid, *Status)
/**/DEBUG_WRITE_ENTERING(ClosingFloatCuboid,
/**/	"About to execute ClosingFloatCuboid")

	switch (Convention) {
		case FiniteDataSupport:
		case MirrorOffBounds:
		case MirrorOnBounds:
		case Periodic:
			break;
		case AntiMirrorOnBounds:
		case FiniteCoefficientSupport:
		default:
			*Status = ERROR;
			WRITE_ERROR(ClosingFloatCuboid, "Invalid boundary convention")
/**/		DEBUG_WRITE_LEAVING(ClosingFloatCuboid, "Done")
			return(*Status);
	}
	MaxFilterFloatCuboid(Volume, Nx, Ny, Nz, Kx, Ky, Kz, Convention, Status);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(ClosingFloatCuboid, "Done")
		return(*Status);
	}
	MinFilterFloatCuboid(Volume, Nx, Ny, Nz, Kx, Ky, Kz, Convention, Status);
/**/DEBUG_WRITE_LEAVING(ClosingFloatCuboid, "Done")
	return(*Status);
} /* end ClosingFloatCuboid */

/*--------------------------------------------------------------------------*/
extern int		ClosingShort
				(
					short	*VolumeSource,		/* data to process */
					short	*VolumeDestination,	/* result */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					short	*Kernel,			/* structural element */
					long	Kx,					/* width of the kernel */
					long	Ky,					/* height of the kernel */
					long	Kz,					/* depth of the kernel */
					long	Ox,					/* kernel X origin */
					long	Oy,					/* kernel Y origin */
					long	Oz,					/* kernel Z origin */
					enum TBoundaryConvention
							Convention,			/* boundary convention */
					int		*Status				/* error management */
				)

/* perform the grey-scale morphological operation known as closing */
/* VolumeSource is the grey-scale volume to close */
/* VolumeDestination is the resulting grey-scale volume after closing */
/* both VolumeSource and VolumeDestination have size (Nx, Ny, Nz) */
/* Kernel is the grey-scale structural element by which VolumeSource is closed */
/* (Kx, Ky, Kz) is the size (bounding box) of the structural element */
/* (Ox, Oy, Oz) is the origin of the structural element */
/* Convention is the boundary convention applied to VolumeSource before closing */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin ClosingShort */

	short	*Buffer = (short *)NULL;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(ClosingShort, VolumeSource, *Status,
/**/	"Missing VolumeSource")
/**/DEBUG_CHECK_NULL_POINTER(ClosingShort, VolumeDestination, *Status,
/**/	"Missing VolumeDestination")
/**/DEBUG_CHECK_NULL_POINTER(ClosingShort, Kernel, *Status,
/**/	"Missing Kernel")
/**/DEBUG_CHECK_RANGE_LONG(ClosingShort, Nx, 1L, LONG_MAX, *Status,
/**/	"Invalid data width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(ClosingShort, Ny, 1L, LONG_MAX, *Status,
/**/	"Invalid data height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(ClosingShort, Nz, 1L, LONG_MAX, *Status,
/**/	"Invalid data depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(ClosingShort, Kx, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(ClosingShort, Ky, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(ClosingShort, Kz, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(ClosingShort, *Status)
/**/DEBUG_WRITE_ENTERING(ClosingShort,
/**/	"About to execute ClosingShort")

	switch (Convention) {
		case FiniteDataSupport:
		case MirrorOffBounds:
		case MirrorOnBounds:
		case Periodic:
			break;
		case AntiMirrorOnBounds:
		case FiniteCoefficientSupport:
		default:
			*Status = ERROR;
			WRITE_ERROR(ClosingShort, "Invalid boundary convention")
/**/		DEBUG_WRITE_LEAVING(ClosingShort, "Done")
			return(*Status);
	}
	AllocateVolumeShort(&Buffer, Nx, Ny, Nz, Status);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(ClosingShort, "Done")
		return(*Status);
	}
	*Status = DilateShort(VolumeSource, Buffer, Nx, Ny, Nz, Kernel, Kx, Ky, Kz,
		Ox, Oy, Oz, Convention);
	if (*Status == ERROR) {
		FreeVolumeShort(&Buffer);
/**/	DEBUG_WRITE_LEAVING(ClosingShort, "Done")
		return(*Status);
	}
	*Status = ErodeShort(Buffer, VolumeDestination, Nx, Ny, Nz, Kernel, Kx, Ky, Kz,
		Ox, Oy, Oz, Convention);
	if (*Status == ERROR) {
		FreeVolumeShort(&Buffer);
/**/	DEBUG_WRITE_LEAVING(ClosingShort, "Done")
		return(*Status);
	}
	*Status = FreeVolumeShort(&Buffer);
/**/DEBUG_WRITE_LEAVING(ClosingShort, "Done")
	return(*Status);
} /* end ClosingShort */

/*--------------------------------------------------------------------------*/
extern int		ClosingShortCuboid
				(
					short	*Volume,			/* data to process in-place */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					long	Kx,					/* width of the kernel */
					long	Ky,					/* height of the kernel */
					long	Kz,					/* depth of the kernel */
					enum TBoundaryConvention
							Convention,			/* boundary convention */
					int		*Status				/* error management */
				)

/* perform a morphological grey-scale closing filter using a cuboid as structural element */
/* Volume is the grey-scale volume to close in-place */
/* Volume has size (Nx, Ny, Nz) */
/* (Kx, Ky, Kz) is the size of the structural element that is filled with the value 1.0F */
/* the origin of the structural element is ((Kx - 1L) / 2L, (Ky - 1L) / 2L, (Kz - 1L) / 2L) */
/* Convention is the boundary convention applied to Volume before closing */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin ClosingShortCuboid */

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(ClosingShortCuboid, Volume, *Status,
/**/	"Missing Volume")
/**/DEBUG_CHECK_RANGE_LONG(ClosingShortCuboid, Nx, 1L, LONG_MAX, *Status,
/**/	"Invalid data width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(ClosingShortCuboid, Ny, 1L, LONG_MAX, *Status,
/**/	"Invalid data height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(ClosingShortCuboid, Nz, 1L, LONG_MAX, *Status,
/**/	"Invalid data depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(ClosingShortCuboid, Kx, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(ClosingShortCuboid, Ky, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(ClosingShortCuboid, Kz, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(ClosingShortCuboid, *Status)
/**/DEBUG_WRITE_ENTERING(ClosingShortCuboid,
/**/	"About to execute ClosingShortCuboid")

	switch (Convention) {
		case FiniteDataSupport:
		case MirrorOffBounds:
		case MirrorOnBounds:
		case Periodic:
			break;
		case AntiMirrorOnBounds:
		case FiniteCoefficientSupport:
		default:
			*Status = ERROR;
			WRITE_ERROR(ClosingShortCuboid, "Invalid boundary convention")
/**/		DEBUG_WRITE_LEAVING(ClosingShortCuboid, "Done")
			return(*Status);
	}
	MaxFilterShortCuboid(Volume, Nx, Ny, Nz, Kx, Ky, Kz, Convention, Status);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(ClosingShortCuboid, "Done")
		return(*Status);
	}
	MinFilterShortCuboid(Volume, Nx, Ny, Nz, Kx, Ky, Kz, Convention, Status);
/**/DEBUG_WRITE_LEAVING(ClosingShortCuboid, "Done")
	return(*Status);
} /* end ClosingShortCuboid */

/*--------------------------------------------------------------------------*/
extern int		DarkTopHatFloat
				(
					float	*VolumeSource,		/* data to process */
					float	*VolumeDestination,	/* result */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					float	*Kernel,			/* structural element */
					long	Kx,					/* width of the kernel */
					long	Ky,					/* height of the kernel */
					long	Kz,					/* depth of the kernel */
					long	Ox,					/* kernel X origin */
					long	Oy,					/* kernel Y origin */
					long	Oz,					/* kernel Z origin */
					enum TBoundaryConvention
							Convention,			/* boundary convention */
					int		*Status				/* error management */
				)

/* perform the grey-scale morphological operation known as dark top hat */
/* VolumeSource is the grey-scale volume to process */
/* VolumeDestination is the resulting grey-scale volume */
/* both VolumeSource and VolumeDestination have size (Nx, Ny, Nz) */
/* Kernel is the grey-scale structural element by which VolumeSource is processed */
/* (Kx, Ky, Kz) is the size (bounding box) of the structural element */
/* (Ox, Oy, Oz) is the origin of the structural element */
/* Convention is the boundary convention applied to VolumeSource before processing */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin DarkTopHatFloat */

	float	*p;
	float	*Closed = (float *)NULL;
	long	x, y, z;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(DarkTopHatFloat, VolumeSource, *Status,
/**/	"Missing VolumeSource")
/**/DEBUG_CHECK_NULL_POINTER(DarkTopHatFloat, VolumeDestination, *Status,
/**/	"Missing VolumeDestination")
/**/DEBUG_CHECK_NULL_POINTER(DarkTopHatFloat, Kernel, *Status,
/**/	"Missing Kernel")
/**/DEBUG_CHECK_RANGE_LONG(DarkTopHatFloat, Nx, 1L, LONG_MAX, *Status,
/**/	"Invalid data width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(DarkTopHatFloat, Ny, 1L, LONG_MAX, *Status,
/**/	"Invalid data height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(DarkTopHatFloat, Nz, 1L, LONG_MAX, *Status,
/**/	"Invalid data depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(DarkTopHatFloat, Kx, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(DarkTopHatFloat, Ky, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(DarkTopHatFloat, Kz, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(DarkTopHatFloat, *Status)
/**/DEBUG_WRITE_ENTERING(DarkTopHatFloat,
/**/	"About to execute DarkTopHatFloat")

	switch (Convention) {
		case FiniteDataSupport:
		case MirrorOffBounds:
		case MirrorOnBounds:
		case Periodic:
			break;
		case AntiMirrorOnBounds:
		case FiniteCoefficientSupport:
		default:
			*Status = ERROR;
			WRITE_ERROR(DarkTopHatFloat, "Invalid boundary convention")
/**/		DEBUG_WRITE_LEAVING(DarkTopHatFloat, "Done")
			return(*Status);
	}
	AllocateVolumeFloat(&Closed, Nx, Ny, Nz, Status);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(DarkTopHatFloat, "Done")
		return(*Status);
	}
	ClosingFloat(VolumeSource, Closed, Nx, Ny, Nz, Kernel, Kx, Ky, Kz,
		Ox, Oy, Oz, Convention, Status);
	if (*Status == ERROR) {
		FreeVolumeFloat(&Closed);
/**/	DEBUG_WRITE_LEAVING(DarkTopHatFloat, "Done")
		return(*Status);
	}
	p = Closed;
	for (z = -Nz; (z < 0L); z++) {
		for (y = -Ny; (y < 0L); y++) {
			for (x = -Nx; (x < 0L); x++) {
				*VolumeDestination++ = *VolumeSource++ - *p++;
			}
		}
	}
	*Status = FreeVolumeFloat(&Closed);
/**/DEBUG_WRITE_LEAVING(DarkTopHatFloat, "Done")
	return(*Status);
} /* end DarkTopHatFloat */

/*--------------------------------------------------------------------------*/
extern int		DarkTopHatFloatCuboid
				(
					float	*Volume,			/* data to process in-place */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					long	Kx,					/* width of the kernel */
					long	Ky,					/* height of the kernel */
					long	Kz,					/* depth of the kernel */
					enum TBoundaryConvention
							Convention,			/* boundary convention */
					int		*Status				/* error management */
				)

/* perform a morphological grey-scale dark top hat filter using a cuboid as structural element */
/* Volume is the grey-scale volume to process in-place */
/* Volume has size (Nx, Ny, Nz) */
/* (Kx, Ky, Kz) is the size of the structural element that is filled with the value 1.0F */
/* the origin of the structural element is ((Kx - 1L) / 2L, (Ky - 1L) / 2L, (Kz - 1L) / 2L) */
/* Convention is the boundary convention applied to Volume before processing */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin DarkTopHatFloatCuboid */

	float	*p;
	float	*Closed = (float *)NULL;
	long	x, y, z;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(DarkTopHatFloatCuboid, Volume, *Status,
/**/	"Missing Volume")
/**/DEBUG_CHECK_RANGE_LONG(DarkTopHatFloatCuboid, Nx, 1L, LONG_MAX, *Status,
/**/	"Invalid data width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(DarkTopHatFloatCuboid, Ny, 1L, LONG_MAX, *Status,
/**/	"Invalid data height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(DarkTopHatFloatCuboid, Nz, 1L, LONG_MAX, *Status,
/**/	"Invalid data depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(DarkTopHatFloatCuboid, Kx, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(DarkTopHatFloatCuboid, Ky, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(DarkTopHatFloatCuboid, Kz, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(DarkTopHatFloatCuboid, *Status)
/**/DEBUG_WRITE_ENTERING(DarkTopHatFloatCuboid,
/**/	"About to execute DarkTopHatFloatCuboid")

	switch (Convention) {
		case FiniteDataSupport:
		case MirrorOffBounds:
		case MirrorOnBounds:
		case Periodic:
			break;
		case AntiMirrorOnBounds:
		case FiniteCoefficientSupport:
		default:
			*Status = ERROR;
			WRITE_ERROR(DarkTopHatFloatCuboid, "Invalid boundary convention")
/**/		DEBUG_WRITE_LEAVING(DarkTopHatFloatCuboid, "Done")
			return(*Status);
	}
	AllocateVolumeFloat(&Closed, Nx, Ny, Nz, Status);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(DarkTopHatFloatCuboid, "Done")
		return(*Status);
	}
	memcpy(Closed, Volume, (size_t)(Nx * Ny * Nz * (long)sizeof(float)));
	ClosingFloatCuboid(Closed, Nx, Ny, Nz, Kx, Ky, Kz, Convention, Status);
	if (*Status == ERROR) {
		FreeVolumeFloat(&Closed);
/**/	DEBUG_WRITE_LEAVING(DarkTopHatFloatCuboid, "Done")
		return(*Status);
	}
	p = Closed;
	for (z = -Nz; (z < 0L); z++) {
		for (y = -Ny; (y < 0L); y++) {
			for (x = -Nx; (x < 0L); x++) {
				*Volume++ -= *p++;
			}
		}
	}
	*Status = FreeVolumeFloat(&Closed);
/**/DEBUG_WRITE_LEAVING(DarkTopHatFloatCuboid, "Done")
	return(*Status);
} /* end DarkTopHatFloatCuboid */

/*--------------------------------------------------------------------------*/
extern int		DarkTopHatShort
				(
					short	*VolumeSource,		/* data to process */
					short	*VolumeDestination,	/* result */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					short	*Kernel,			/* structural element */
					long	Kx,					/* width of the kernel */
					long	Ky,					/* height of the kernel */
					long	Kz,					/* depth of the kernel */
					long	Ox,					/* kernel X origin */
					long	Oy,					/* kernel Y origin */
					long	Oz,					/* kernel Z origin */
					enum TBoundaryConvention
							Convention,			/* boundary convention */
					int		*Status				/* error management */
				)

/* perform the grey-scale morphological operation known as dark top hat */
/* VolumeSource is the grey-scale volume to process */
/* VolumeDestination is the resulting grey-scale volume */
/* both VolumeSource and VolumeDestination have size (Nx, Ny, Nz) */
/* Kernel is the grey-scale structural element by which VolumeSource is processed */
/* (Kx, Ky, Kz) is the size (bounding box) of the structural element */
/* (Ox, Oy, Oz) is the origin of the structural element */
/* Convention is the boundary convention applied to VolumeSource before processing */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin DarkTopHatShort */

	short	*p;
	short	*Closed = (short *)NULL;
	long	x, y, z;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(DarkTopHatShort, VolumeSource, *Status,
/**/	"Missing VolumeSource")
/**/DEBUG_CHECK_NULL_POINTER(DarkTopHatShort, VolumeDestination, *Status,
/**/	"Missing VolumeDestination")
/**/DEBUG_CHECK_NULL_POINTER(DarkTopHatShort, Kernel, *Status,
/**/	"Missing Kernel")
/**/DEBUG_CHECK_RANGE_LONG(DarkTopHatShort, Nx, 1L, LONG_MAX, *Status,
/**/	"Invalid data width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(DarkTopHatShort, Ny, 1L, LONG_MAX, *Status,
/**/	"Invalid data height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(DarkTopHatShort, Nz, 1L, LONG_MAX, *Status,
/**/	"Invalid data depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(DarkTopHatShort, Kx, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(DarkTopHatShort, Ky, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(DarkTopHatShort, Kz, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(DarkTopHatShort, *Status)
/**/DEBUG_WRITE_ENTERING(DarkTopHatShort,
/**/	"About to execute DarkTopHatShort")

	switch (Convention) {
		case FiniteDataSupport:
		case MirrorOffBounds:
		case MirrorOnBounds:
		case Periodic:
			break;
		case AntiMirrorOnBounds:
		case FiniteCoefficientSupport:
		default:
			*Status = ERROR;
			WRITE_ERROR(DarkTopHatShort, "Invalid boundary convention")
/**/		DEBUG_WRITE_LEAVING(DarkTopHatShort, "Done")
			return(*Status);
	}
	AllocateVolumeShort(&Closed, Nx, Ny, Nz, Status);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(DarkTopHatShort, "Done")
		return(*Status);
	}
	ClosingShort(VolumeSource, Closed, Nx, Ny, Nz, Kernel, Kx, Ky, Kz,
		Ox, Oy, Oz, Convention, Status);
	if (*Status == ERROR) {
		FreeVolumeShort(&Closed);
/**/	DEBUG_WRITE_LEAVING(DarkTopHatShort, "Done")
		return(*Status);
	}
	p = Closed;
	for (z = -Nz; (z < 0L); z++) {
		for (y = -Ny; (y < 0L); y++) {
			for (x = -Nx; (x < 0L); x++) {
				*VolumeDestination++ = ConvertIntToShort((int)*VolumeSource++ - (int)*p++);
			}
		}
	}
	*Status = FreeVolumeShort(&Closed);
/**/DEBUG_WRITE_LEAVING(DarkTopHatShort, "Done")
	return(*Status);
} /* end DarkTopHatShort */

/*--------------------------------------------------------------------------*/
extern int		DarkTopHatShortCuboid
				(
					short	*Volume,			/* data to process in-place */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					long	Kx,					/* width of the kernel */
					long	Ky,					/* height of the kernel */
					long	Kz,					/* depth of the kernel */
					enum TBoundaryConvention
							Convention,			/* boundary convention */
					int		*Status				/* error management */
				)

/* perform a morphological grey-scale dark top hat filter using a cuboid as structural element */
/* Volume is the grey-scale volume to process in-place */
/* Volume has size (Nx, Ny, Nz) */
/* (Kx, Ky, Kz) is the size of the structural element that is filled with the value 1.0F */
/* the origin of the structural element is ((Kx - 1L) / 2L, (Ky - 1L) / 2L, (Kz - 1L) / 2L) */
/* Convention is the boundary convention applied to Volume before processing */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin DarkTopHatShortCuboid */

	short	*p;
	short	*Closed = (short *)NULL;
	long	x, y, z;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(DarkTopHatShortCuboid, Volume, *Status,
/**/	"Missing Volume")
/**/DEBUG_CHECK_RANGE_LONG(DarkTopHatShortCuboid, Nx, 1L, LONG_MAX, *Status,
/**/	"Invalid data width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(DarkTopHatShortCuboid, Ny, 1L, LONG_MAX, *Status,
/**/	"Invalid data height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(DarkTopHatShortCuboid, Nz, 1L, LONG_MAX, *Status,
/**/	"Invalid data depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(DarkTopHatShortCuboid, Kx, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(DarkTopHatShortCuboid, Ky, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(DarkTopHatShortCuboid, Kz, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(DarkTopHatShortCuboid, *Status)
/**/DEBUG_WRITE_ENTERING(DarkTopHatShortCuboid,
/**/	"About to execute DarkTopHatShortCuboid")

	switch (Convention) {
		case FiniteDataSupport:
		case MirrorOffBounds:
		case MirrorOnBounds:
		case Periodic:
			break;
		case AntiMirrorOnBounds:
		case FiniteCoefficientSupport:
		default:
			*Status = ERROR;
			WRITE_ERROR(DarkTopHatShortCuboid, "Invalid boundary convention")
/**/		DEBUG_WRITE_LEAVING(DarkTopHatShortCuboid, "Done")
			return(*Status);
	}
	AllocateVolumeShort(&Closed, Nx, Ny, Nz, Status);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(DarkTopHatShortCuboid, "Done")
		return(*Status);
	}
	memcpy(Closed, Volume, (size_t)(Nx * Ny * Nz * (long)sizeof(short)));
	ClosingShortCuboid(Closed, Nx, Ny, Nz, Kx, Ky, Kz, Convention, Status);
	if (*Status == ERROR) {
		FreeVolumeShort(&Closed);
/**/	DEBUG_WRITE_LEAVING(DarkTopHatShortCuboid, "Done")
		return(*Status);
	}
	p = Closed;
	for (z = -Nz; (z < 0L); z++) {
		for (y = -Ny; (y < 0L); y++) {
			for (x = -Nx; (x < 0L); x++) {
				*Volume++ -= *p++;
			}
		}
	}
	*Status = FreeVolumeShort(&Closed);
/**/DEBUG_WRITE_LEAVING(DarkTopHatShortCuboid, "Done")
	return(*Status);
} /* end DarkTopHatShortCuboid */

/*--------------------------------------------------------------------------*/
extern int		DilateFloat
				(
					float	*VolumeSource,		/* data to process */
					float	*VolumeDestination,	/* result */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					float	*Kernel,			/* structural element */
					long	Kx,					/* width of the kernel */
					long	Ky,					/* height of the kernel */
					long	Kz,					/* depth of the kernel */
					long	Ox,					/* kernel X origin */
					long	Oy,					/* kernel Y origin */
					long	Oz,					/* kernel Z origin */
					enum TBoundaryConvention
							Convention			/* boundary convention */
				)

/* perform the grey-scale morphological operation known as dilatation */
/* VolumeSource is the grey-scale volume to dilate */
/* VolumeDestination is the resulting grey-scale volume after dilatation */
/* both VolumeSource and VolumeDestination have size (Nx, Ny, Nz) */
/* Kernel is the grey-scale structural element by which VolumeSource is dilated */
/* (Kx, Ky, Kz) is the size (bounding box) of the structural element */
/* (Ox, Oy, Oz) is the origin of the structural element */
/* Convention is the boundary convention applied to VolumeSource before dilatation */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin DilateFloat */

	float	*p, *q;
	float	f, fmax;
	long	x, y, z;
	long	x0, y0, z0;
	long	i, j, k;
	const long
			Nxy = Nx * Ny;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(DilateFloat, VolumeSource, Status,
/**/	"Missing VolumeSource")
/**/DEBUG_CHECK_NULL_POINTER(DilateFloat, VolumeDestination, Status,
/**/	"Missing VolumeDestination")
/**/DEBUG_CHECK_NULL_POINTER(DilateFloat, Kernel, Status,
/**/	"Missing Kernel")
/**/DEBUG_CHECK_RANGE_LONG(DilateFloat, Nx, 1L, LONG_MAX, Status,
/**/	"Invalid data width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(DilateFloat, Ny, 1L, LONG_MAX, Status,
/**/	"Invalid data height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(DilateFloat, Nz, 1L, LONG_MAX, Status,
/**/	"Invalid data depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(DilateFloat, Kx, 1L, LONG_MAX, Status,
/**/	"Invalid kernel width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(DilateFloat, Ky, 1L, LONG_MAX, Status,
/**/	"Invalid kernel height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(DilateFloat, Kz, 1L, LONG_MAX, Status,
/**/	"Invalid kernel depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(DilateFloat, Status)
/**/DEBUG_WRITE_ENTERING(DilateFloat,
/**/	"About to execute DilateFloat")

	switch (Convention) {
		case FiniteDataSupport:
		case MirrorOffBounds:
		case MirrorOnBounds:
		case Periodic:
			break;
		case AntiMirrorOnBounds:
		case FiniteCoefficientSupport:
		default:
			Status = ERROR;
			WRITE_ERROR(DilateFloat, "Invalid boundary convention")
/**/		DEBUG_WRITE_LEAVING(DilateFloat, "Done")
			return(Status);
	}
	p = VolumeDestination;
	for (z = Oz; (z < (Nz + Oz)); z++) {
		for (y = Oy; (y < (Ny + Oy)); y++) {
			for (x = Ox; (x < (Nx + Ox)); x++) {
				q = Kernel;
				Status = GetFoldedIndex(z, &z0, Nz, Convention);
				if (Status == ERROR) {
/**/				DEBUG_WRITE_LEAVING(DilateFloat, "Done")
					return(Status);
				}
				Status = GetFoldedIndex(y, &y0, Ny, Convention);
				if (Status == ERROR) {
/**/				DEBUG_WRITE_LEAVING(DilateFloat, "Done")
					return(Status);
				}
				Status = GetFoldedIndex(x, &x0, Nx, Convention);
				if (Status == ERROR) {
/**/				DEBUG_WRITE_LEAVING(DilateFloat, "Done")
					return(Status);
				}
				fmax = VolumeSource[Nxy * z0 + Nx * y0 + x0] + *q;
				for (k = z + 1L - Kz; (k <= 0L); k++) {
					Status = GetFoldedIndex(k, &z0, Nz, Convention);
					if (Status == ERROR) {
/**/					DEBUG_WRITE_LEAVING(DilateFloat, "Done")
						return(Status);
					}
					z0 *= Nxy;
					for (j = y + 1L - Ky; (j <= 0L); j++) {
						Status = GetFoldedIndex(j, &y0, Ny, Convention);
						if (Status == ERROR) {
/**/						DEBUG_WRITE_LEAVING(DilateFloat, "Done")
							return(Status);
						}
						y0 *= Nx;
						y0 += z0;
						for (i = x + 1L - Kx; (i <= 0L); i++) {
							Status = GetFoldedIndex(i, &x0, Nx, Convention);
							if (Status == ERROR) {
/**/							DEBUG_WRITE_LEAVING(DilateFloat, "Done")
								return(Status);
							}
							f = VolumeSource[y0 + x0] + *q++;
							fmax = (f < fmax) ? (fmax) : (f);
						}
					}
				}
				*p++ = fmax;
			}
		}
	}
/**/DEBUG_WRITE_LEAVING(DilateFloat, "Done")
	return(Status);
} /* end DilateFloat */

/*--------------------------------------------------------------------------*/
extern int		DilateShort
				(
					short	*VolumeSource,		/* data to process */
					short	*VolumeDestination,	/* result */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					short	*Kernel,			/* structural element */
					long	Kx,					/* width of the kernel */
					long	Ky,					/* height of the kernel */
					long	Kz,					/* depth of the kernel */
					long	Ox,					/* kernel X origin */
					long	Oy,					/* kernel Y origin */
					long	Oz,					/* kernel Z origin */
					enum TBoundaryConvention
							Convention			/* boundary convention */
				)

/* perform the grey-scale morphological operation known as dilatation */
/* VolumeSource is the grey-scale volume to dilate */
/* VolumeDestination is the resulting grey-scale volume after dilatation */
/* both VolumeSource and VolumeDestination have size (Nx, Ny, Nz) */
/* Kernel is the grey-scale structural element by which VolumeSource is dilated */
/* (Kx, Ky, Kz) is the size (bounding box) of the structural element */
/* (Ox, Oy, Oz) is the origin of the structural element */
/* Convention is the boundary convention applied to VolumeSource before dilatation */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin DilateShort */

	short	*p, *q;
	short	f, fmax;
	long	x, y, z;
	long	x0, y0, z0;
	long	i, j, k;
	const long
			Nxy = Nx * Ny;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(DilateShort, VolumeSource, Status,
/**/	"Missing VolumeSource")
/**/DEBUG_CHECK_NULL_POINTER(DilateShort, VolumeDestination, Status,
/**/	"Missing VolumeDestination")
/**/DEBUG_CHECK_NULL_POINTER(DilateShort, Kernel, Status,
/**/	"Missing Kernel")
/**/DEBUG_CHECK_RANGE_LONG(DilateShort, Nx, 1L, LONG_MAX, Status,
/**/	"Invalid data width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(DilateShort, Ny, 1L, LONG_MAX, Status,
/**/	"Invalid data height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(DilateShort, Nz, 1L, LONG_MAX, Status,
/**/	"Invalid data depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(DilateShort, Kx, 1L, LONG_MAX, Status,
/**/	"Invalid kernel width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(DilateShort, Ky, 1L, LONG_MAX, Status,
/**/	"Invalid kernel height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(DilateShort, Kz, 1L, LONG_MAX, Status,
/**/	"Invalid kernel depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(DilateShort, Status)
/**/DEBUG_WRITE_ENTERING(DilateShort,
/**/	"About to execute DilateShort")

	switch (Convention) {
		case FiniteDataSupport:
		case MirrorOffBounds:
		case MirrorOnBounds:
		case Periodic:
			break;
		case AntiMirrorOnBounds:
		case FiniteCoefficientSupport:
		default:
			Status = ERROR;
			WRITE_ERROR(DilateShort, "Invalid boundary convention")
/**/		DEBUG_WRITE_LEAVING(DilateShort, "Done")
			return(Status);
	}
	p = VolumeDestination;
	for (z = Oz; (z < (Nz + Oz)); z++) {
		for (y = Oy; (y < (Ny + Oy)); y++) {
			for (x = Ox; (x < (Nx + Ox)); x++) {
				q = Kernel;
				Status = GetFoldedIndex(z, &z0, Nz, Convention);
				if (Status == ERROR) {
/**/				DEBUG_WRITE_LEAVING(DilateShort, "Done")
					return(Status);
				}
				Status = GetFoldedIndex(y, &y0, Ny, Convention);
				if (Status == ERROR) {
/**/				DEBUG_WRITE_LEAVING(DilateShort, "Done")
					return(Status);
				}
				Status = GetFoldedIndex(x, &x0, Nx, Convention);
				if (Status == ERROR) {
/**/				DEBUG_WRITE_LEAVING(DilateShort, "Done")
					return(Status);
				}
				fmax = ConvertIntToShort((int)VolumeSource[Nxy * z0 + Nx * y0 + x0] + (int)*q);
				for (k = z + 1L - Kz; (k <= 0L); k++) {
					Status = GetFoldedIndex(k, &z0, Nz, Convention);
					if (Status == ERROR) {
/**/					DEBUG_WRITE_LEAVING(DilateShort, "Done")
						return(Status);
					}
					z0 *= Nxy;
					for (j = y + 1L - Ky; (j <= 0L); j++) {
						Status = GetFoldedIndex(j, &y0, Ny, Convention);
						if (Status == ERROR) {
/**/						DEBUG_WRITE_LEAVING(DilateShort, "Done")
							return(Status);
						}
						y0 *= Nx;
						y0 += z0;
						for (i = x + 1L - Kx; (i <= 0L); i++) {
							Status = GetFoldedIndex(i, &x0, Nx, Convention);
							if (Status == ERROR) {
/**/							DEBUG_WRITE_LEAVING(DilateShort, "Done")
								return(Status);
							}
							f = ConvertIntToShort((int)VolumeSource[y0 + x0] + (int)*q++);
							fmax = (f < fmax) ? (fmax) : (f);
						}
					}
				}
				*p++ = fmax;
			}
		}
	}
/**/DEBUG_WRITE_LEAVING(DilateShort, "Done")
	return(Status);
} /* end DilateShort */

/*--------------------------------------------------------------------------*/
extern int		ErodeFloat
				(
					float	*VolumeSource,		/* data to process */
					float	*VolumeDestination,	/* result */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					float	*Kernel,			/* structural element */
					long	Kx,					/* width of the kernel */
					long	Ky,					/* height of the kernel */
					long	Kz,					/* depth of the kernel */
					long	Ox,					/* kernel X origin */
					long	Oy,					/* kernel Y origin */
					long	Oz,					/* kernel Z origin */
					enum TBoundaryConvention
							Convention			/* boundary convention */
				)

/* perform the grey-scale morphological operation known as erosion */
/* VolumeSource is the grey-scale volume to erode */
/* VolumeDestination is the resulting grey-scale volume after erosion */
/* both VolumeSource and VolumeDestination have size (Nx, Ny, Nz) */
/* Kernel is the grey-scale structural element by which VolumeSource is eroded */
/* (Kx, Ky, Kz) is the size (bounding box) of the structural element */
/* (Ox, Oy, Oz) is the origin of the structural element */
/* Convention is the boundary convention applied to VolumeSource before erosion */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin ErodeFloat */

	float	*p, *q;
	float	f, fmin;
	long	x, y, z;
	long	x0, y0, z0;
	long	i, j, k;
	const long
			Nxy = Nx * Ny;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(ErodeFloat, VolumeSource, Status,
/**/	"Missing VolumeSource")
/**/DEBUG_CHECK_NULL_POINTER(ErodeFloat, VolumeDestination, Status,
/**/	"Missing VolumeDestination")
/**/DEBUG_CHECK_NULL_POINTER(ErodeFloat, Kernel, Status,
/**/	"Missing Kernel")
/**/DEBUG_CHECK_RANGE_LONG(ErodeFloat, Nx, 1L, LONG_MAX, Status,
/**/	"Invalid data width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(ErodeFloat, Ny, 1L, LONG_MAX, Status,
/**/	"Invalid data height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(ErodeFloat, Nz, 1L, LONG_MAX, Status,
/**/	"Invalid data depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(ErodeFloat, Kx, 1L, LONG_MAX, Status,
/**/	"Invalid kernel width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(ErodeFloat, Ky, 1L, LONG_MAX, Status,
/**/	"Invalid kernel height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(ErodeFloat, Kz, 1L, LONG_MAX, Status,
/**/	"Invalid kernel depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(ErodeFloat, Status)
/**/DEBUG_WRITE_ENTERING(ErodeFloat,
/**/	"About to execute ErodeFloat")

	switch (Convention) {
		case FiniteDataSupport:
		case MirrorOffBounds:
		case MirrorOnBounds:
		case Periodic:
			break;
		case AntiMirrorOnBounds:
		case FiniteCoefficientSupport:
		default:
			Status = ERROR;
			WRITE_ERROR(ErodeFloat, "Invalid boundary convention")
/**/		DEBUG_WRITE_LEAVING(ErodeFloat, "Done")
			return(Status);
	}
	p = VolumeDestination;
	for (z = -Oz; (z < (Nz - Oz)); z++) {
		for (y = -Oy; (y < (Ny - Oy)); y++) {
			for (x = -Ox; (x < (Nx - Ox)); x++) {
				q = Kernel;
				Status = GetFoldedIndex(z, &z0, Nz, Convention);
				if (Status == ERROR) {
/**/				DEBUG_WRITE_LEAVING(ErodeFloat, "Done")
					return(Status);
				}
				Status = GetFoldedIndex(y, &y0, Ny, Convention);
				if (Status == ERROR) {
/**/				DEBUG_WRITE_LEAVING(ErodeFloat, "Done")
					return(Status);
				}
				Status = GetFoldedIndex(x, &x0, Nx, Convention);
				if (Status == ERROR) {
/**/				DEBUG_WRITE_LEAVING(ErodeFloat, "Done")
					return(Status);
				}
				fmin = VolumeSource[Nxy * z0 + Nx * y0 + x0] - *q;
				for (k = z; (k < (z + Kz)); k++) {
					Status = GetFoldedIndex(k, &z0, Nz, Convention);
					if (Status == ERROR) {
/**/					DEBUG_WRITE_LEAVING(ErodeFloat, "Done")
						return(Status);
					}
					z0 *= Nxy;
					for (j = y; (j < (y + Ky)); j++) {
						Status = GetFoldedIndex(j, &y0, Ny, Convention);
						if (Status == ERROR) {
/**/						DEBUG_WRITE_LEAVING(ErodeFloat, "Done")
							return(Status);
						}
						y0 *= Nx;
						y0 += z0;
						for (i = x; (i < (x + Kx)); i++) {
							Status = GetFoldedIndex(i, &x0, Nx, Convention);
							if (Status == ERROR) {
/**/							DEBUG_WRITE_LEAVING(ErodeFloat, "Done")
								return(Status);
							}
							f = VolumeSource[y0 + x0] - *q++;
							fmin = (f < fmin) ? (f) : (fmin);
						}
					}
				}
				*p++ = fmin;
			}
		}
	}
/**/DEBUG_WRITE_LEAVING(ErodeFloat, "Done")
	return(Status);
} /* end ErodeFloat */

/*--------------------------------------------------------------------------*/
extern int		ErodeShort
				(
					short	*VolumeSource,		/* data to process */
					short	*VolumeDestination,	/* result */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					short	*Kernel,			/* structural element */
					long	Kx,					/* width of the kernel */
					long	Ky,					/* height of the kernel */
					long	Kz,					/* depth of the kernel */
					long	Ox,					/* kernel X origin */
					long	Oy,					/* kernel Y origin */
					long	Oz,					/* kernel Z origin */
					enum TBoundaryConvention
							Convention			/* boundary convention */
				)

/* perform the grey-scale morphological operation known as erosion */
/* VolumeSource is the grey-scale volume to erode */
/* VolumeDestination is the resulting grey-scale volume after erosion */
/* both VolumeSource and VolumeDestination have size (Nx, Ny, Nz) */
/* Kernel is the grey-scale structural element by which VolumeSource is eroded */
/* (Kx, Ky, Kz) is the size (bounding box) of the structural element */
/* (Ox, Oy, Oz) is the origin of the structural element */
/* Convention is the boundary convention applied to VolumeSource before erosion */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin ErodeShort */

	short	*p, *q;
	short	f, fmin;
	long	x, y, z;
	long	x0, y0, z0;
	long	i, j, k;
	const long
			Nxy = Nx * Ny;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(ErodeShort, VolumeSource, Status,
/**/	"Missing VolumeSource")
/**/DEBUG_CHECK_NULL_POINTER(ErodeShort, VolumeDestination, Status,
/**/	"Missing VolumeDestination")
/**/DEBUG_CHECK_NULL_POINTER(ErodeShort, Kernel, Status,
/**/	"Missing Kernel")
/**/DEBUG_CHECK_RANGE_LONG(ErodeShort, Nx, 1L, LONG_MAX, Status,
/**/	"Invalid data width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(ErodeShort, Ny, 1L, LONG_MAX, Status,
/**/	"Invalid data height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(ErodeShort, Nz, 1L, LONG_MAX, Status,
/**/	"Invalid data depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(ErodeShort, Kx, 1L, LONG_MAX, Status,
/**/	"Invalid kernel width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(ErodeShort, Ky, 1L, LONG_MAX, Status,
/**/	"Invalid kernel height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(ErodeShort, Kz, 1L, LONG_MAX, Status,
/**/	"Invalid kernel depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(ErodeShort, Status)
/**/DEBUG_WRITE_ENTERING(ErodeShort,
/**/	"About to execute ErodeShort")

	switch (Convention) {
		case FiniteDataSupport:
		case MirrorOffBounds:
		case MirrorOnBounds:
		case Periodic:
			break;
		case AntiMirrorOnBounds:
		case FiniteCoefficientSupport:
		default:
			Status = ERROR;
			WRITE_ERROR(ErodeShort, "Invalid boundary convention")
/**/		DEBUG_WRITE_LEAVING(ErodeShort, "Done")
			return(Status);
	}
	p = VolumeDestination;
	for (z = -Oz; (z < (Nz - Oz)); z++) {
		for (y = -Oy; (y < (Ny - Oy)); y++) {
			for (x = -Ox; (x < (Nx - Ox)); x++) {
				q = Kernel;
				Status = GetFoldedIndex(z, &z0, Nz, Convention);
				if (Status == ERROR) {
/**/				DEBUG_WRITE_LEAVING(ErodeShort, "Done")
					return(Status);
				}
				Status = GetFoldedIndex(y, &y0, Ny, Convention);
				if (Status == ERROR) {
/**/				DEBUG_WRITE_LEAVING(ErodeShort, "Done")
					return(Status);
				}
				Status = GetFoldedIndex(x, &x0, Nx, Convention);
				if (Status == ERROR) {
/**/				DEBUG_WRITE_LEAVING(ErodeShort, "Done")
					return(Status);
				}
				fmin = ConvertIntToShort((int)VolumeSource[Nxy * z0 + Nx * y0 + x0] - (int)*q);
				for (k = z; (k < (z + Kz)); k++) {
					Status = GetFoldedIndex(k, &z0, Nz, Convention);
					if (Status == ERROR) {
/**/					DEBUG_WRITE_LEAVING(ErodeShort, "Done")
						return(Status);
					}
					z0 *= Nxy;
					for (j = y; (j < (y + Ky)); j++) {
						Status = GetFoldedIndex(j, &y0, Ny, Convention);
						if (Status == ERROR) {
/**/						DEBUG_WRITE_LEAVING(ErodeShort, "Done")
							return(Status);
						}
						y0 *= Nx;
						y0 += z0;
						for (i = x; (i < (x + Kx)); i++) {
							Status = GetFoldedIndex(i, &x0, Nx, Convention);
							if (Status == ERROR) {
/**/							DEBUG_WRITE_LEAVING(ErodeShort, "Done")
								return(Status);
							}
							f = ConvertIntToShort((int)VolumeSource[y0 + x0] - (int)*q++);
							fmin = (f < fmin) ? (f) : (fmin);
						}
					}
				}
				*p++ = fmin;
			}
		}
	}
/**/DEBUG_WRITE_LEAVING(ErodeShort, "Done")
	return(Status);
} /* end ErodeShort */

/*--------------------------------------------------------------------------*/
extern int		MaxFilterFloatCuboid
				(
					float	*Volume,			/* data to process in-place */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					long	Kx,					/* width of the kernel */
					long	Ky,					/* height of the kernel */
					long	Kz,					/* depth of the kernel */
					enum TBoundaryConvention
							Convention,			/* boundary convention */
					int		*Status				/* error management */
				)

/* perform a morphological grey-scale max filter using a cuboid as structural element */
/* Volume is the grey-scale volume to dilate in-place */
/* Volume has size (Nx, Ny, Nz) */
/* (Kx, Ky, Kz) is the size of the structural element that is filled with the value 1.0F */
/* the origin of the structural element is ((Kx - 1L) / 2L, (Ky - 1L) / 2L, (Kz - 1L) / 2L) */
/* Convention is the boundary convention applied to Volume before dilatation */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin MaxFilterFloatCuboid */

	double	*BufferIn = (double *)NULL, *BufferOut = (double *)NULL;
	long	i, j, k;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(MaxFilterFloatCuboid, Volume, *Status,
/**/	"Missing Volume")
/**/DEBUG_CHECK_RANGE_LONG(MaxFilterFloatCuboid, Nx, 1L, LONG_MAX, *Status,
/**/	"Invalid data width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MaxFilterFloatCuboid, Ny, 1L, LONG_MAX, *Status,
/**/	"Invalid data height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MaxFilterFloatCuboid, Nz, 1L, LONG_MAX, *Status,
/**/	"Invalid data depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MaxFilterFloatCuboid, Kx, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MaxFilterFloatCuboid, Ky, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MaxFilterFloatCuboid, Kz, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(MaxFilterFloatCuboid, *Status)
/**/DEBUG_WRITE_ENTERING(MaxFilterFloatCuboid,
/**/	"About to execute MaxFilterFloatCuboid")

	switch (Convention) {
		case FiniteDataSupport:
		case MirrorOffBounds:
		case MirrorOnBounds:
		case Periodic:
			break;
		case AntiMirrorOnBounds:
		case FiniteCoefficientSupport:
		default:
			*Status = ERROR;
			WRITE_ERROR(MaxFilterFloatCuboid, "Invalid boundary convention")
/**/		DEBUG_WRITE_LEAVING(MaxFilterFloatCuboid, "Done")
			return(*Status);
	}
	if (1L < Nx) {
		AllocateLineDouble(&BufferIn, Nx, Status);
		if (*Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(MaxFilterFloatCuboid, "Done")
			return(*Status);
		}
		AllocateLineDouble(&BufferOut, Nx, Status);
		if (*Status == ERROR) {
			FreeLineDouble(&BufferIn);
/**/		DEBUG_WRITE_LEAVING(MaxFilterFloatCuboid, "Done")
			return(*Status);
		}
		for (k = 0L; (k < Nz); k++) {
			for (j = 0L; (j < Ny); j++) {
				*Status = GetxFloatToDouble(Volume, Nx, Ny, Nz, 0L, j, k, BufferIn, Nx);
				if (*Status == ERROR) {
					FreeLineDouble(&BufferOut);
					FreeLineDouble(&BufferIn);
/**/				DEBUG_WRITE_LEAVING(MaxFilterFloatCuboid, "Done")
					return(*Status);
				}
				*Status = MaxFilterLine(BufferIn, BufferOut, Nx, Kx, Convention);
				if (*Status == ERROR) {
					FreeLineDouble(&BufferOut);
					FreeLineDouble(&BufferIn);
/**/				DEBUG_WRITE_LEAVING(MaxFilterFloatCuboid, "Done")
					return(*Status);
				}
				*Status = PutxDoubleToFloat(Volume, Nx, Ny, Nz, 0L, j, k, BufferOut, Nx);
				if (*Status == ERROR) {
					FreeLineDouble(&BufferOut);
					FreeLineDouble(&BufferIn);
/**/				DEBUG_WRITE_LEAVING(MaxFilterFloatCuboid, "Done")
					return(*Status);
				}
			}
		}
		*Status = FreeLineDouble(&BufferOut);
		if (*Status == ERROR) {
			FreeLineDouble(&BufferIn);
/**/		DEBUG_WRITE_LEAVING(MaxFilterFloatCuboid, "Done")
			return(*Status);
		}
		*Status = FreeLineDouble(&BufferIn);
		if (*Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(MaxFilterFloatCuboid, "Done")
			return(*Status);
		}
	}
	if (1L < Ny) {
		AllocateLineDouble(&BufferIn, Ny, Status);
		if (*Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(MaxFilterFloatCuboid, "Done")
			return(*Status);
		}
		AllocateLineDouble(&BufferOut, Ny, Status);
		if (*Status == ERROR) {
			FreeLineDouble(&BufferIn);
/**/		DEBUG_WRITE_LEAVING(MaxFilterFloatCuboid, "Done")
			return(*Status);
		}
		for (k = 0L; (k < Nz); k++) {
			for (i = 0L; (i < Nx); i++) {
				*Status = GetyFloatToDouble(Volume, Nx, Ny, Nz, i, 0L, k, BufferIn, Ny);
				if (*Status == ERROR) {
					FreeLineDouble(&BufferOut);
					FreeLineDouble(&BufferIn);
/**/				DEBUG_WRITE_LEAVING(MaxFilterFloatCuboid, "Done")
					return(*Status);
				}
				*Status = MaxFilterLine(BufferIn, BufferOut, Ny, Ky, Convention);
				if (*Status == ERROR) {
					FreeLineDouble(&BufferOut);
					FreeLineDouble(&BufferIn);
/**/				DEBUG_WRITE_LEAVING(MaxFilterFloatCuboid, "Done")
					return(*Status);
				}
				*Status = PutyDoubleToFloat(Volume, Nx, Ny, Nz, i, 0L, k, BufferOut, Ny);
				if (*Status == ERROR) {
					FreeLineDouble(&BufferOut);
					FreeLineDouble(&BufferIn);
/**/				DEBUG_WRITE_LEAVING(MaxFilterFloatCuboid, "Done")
					return(*Status);
				}
			}
		}
		*Status = FreeLineDouble(&BufferOut);
		if (*Status == ERROR) {
			FreeLineDouble(&BufferIn);
/**/		DEBUG_WRITE_LEAVING(MaxFilterFloatCuboid, "Done")
			return(*Status);
		}
		*Status = FreeLineDouble(&BufferIn);
		if (*Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(MaxFilterFloatCuboid, "Done")
			return(*Status);
		}
	}
	if (1L < Nz) {
		AllocateLineDouble(&BufferIn, Nz, Status);
		if (*Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(MaxFilterFloatCuboid, "Done")
			return(*Status);
		}
		AllocateLineDouble(&BufferOut, Nz, Status);
		if (*Status == ERROR) {
			FreeLineDouble(&BufferIn);
/**/		DEBUG_WRITE_LEAVING(MaxFilterFloatCuboid, "Done")
			return(*Status);
		}
		for (j = 0L; (j < Ny); j++) {
			for (i = 0L; (i < Nx); i++) {
				*Status = GetzFloatToDouble(Volume, Nx, Ny, Nz, i, j, 0L, BufferIn, Nz);
				if (*Status == ERROR) {
					FreeLineDouble(&BufferOut);
					FreeLineDouble(&BufferIn);
/**/				DEBUG_WRITE_LEAVING(MaxFilterFloatCuboid, "Done")
					return(*Status);
				}
				*Status = MaxFilterLine(BufferIn, BufferOut, Nz, Kz, Convention);
				if (*Status == ERROR) {
					FreeLineDouble(&BufferOut);
					FreeLineDouble(&BufferIn);
/**/				DEBUG_WRITE_LEAVING(MaxFilterFloatCuboid, "Done")
					return(*Status);
				}
				*Status = PutzDoubleToFloat(Volume, Nx, Ny, Nz, i, j, 0L, BufferOut, Nz);
				if (*Status == ERROR) {
					FreeLineDouble(&BufferOut);
					FreeLineDouble(&BufferIn);
/**/				DEBUG_WRITE_LEAVING(MaxFilterFloatCuboid, "Done")
					return(*Status);
				}
			}
		}
		*Status = FreeLineDouble(&BufferOut);
		if (*Status == ERROR) {
			FreeLineDouble(&BufferIn);
/**/		DEBUG_WRITE_LEAVING(MaxFilterFloatCuboid, "Done")
			return(*Status);
		}
		*Status = FreeLineDouble(&BufferIn);
		if (*Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(MaxFilterFloatCuboid, "Done")
			return(*Status);
		}
	}
/**/DEBUG_WRITE_LEAVING(MaxFilterFloatCuboid, "Done")
	return(*Status);
} /* end MaxFilterFloatCuboid */

/*--------------------------------------------------------------------------*/
extern int		MaxFilterShortCuboid
				(
					short	*Volume,			/* data to process in-place */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					long	Kx,					/* width of the kernel */
					long	Ky,					/* height of the kernel */
					long	Kz,					/* depth of the kernel */
					enum TBoundaryConvention
							Convention,			/* boundary convention */
					int		*Status				/* error management */
				)

/* perform a morphological grey-scale max filter using a cuboid as structural element */
/* Volume is the grey-scale volume to dilate in-place */
/* Volume has size (Nx, Ny, Nz) */
/* (Kx, Ky, Kz) is the size of the structural element that is filled with the value 1.0F */
/* the origin of the structural element is ((Kx - 1L) / 2L, (Ky - 1L) / 2L, (Kz - 1L) / 2L) */
/* Convention is the boundary convention applied to Volume before dilatation */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin MaxFilterShortCuboid */

	double	*BufferIn = (double *)NULL, *BufferOut = (double *)NULL;
	long	i, j, k;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(MaxFilterShortCuboid, Volume, *Status,
/**/	"Missing Volume")
/**/DEBUG_CHECK_RANGE_LONG(MaxFilterShortCuboid, Nx, 1L, LONG_MAX, *Status,
/**/	"Invalid data width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MaxFilterShortCuboid, Ny, 1L, LONG_MAX, *Status,
/**/	"Invalid data height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MaxFilterShortCuboid, Nz, 1L, LONG_MAX, *Status,
/**/	"Invalid data depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MaxFilterShortCuboid, Kx, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MaxFilterShortCuboid, Ky, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MaxFilterShortCuboid, Kz, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(MaxFilterShortCuboid, *Status)
/**/DEBUG_WRITE_ENTERING(MaxFilterShortCuboid,
/**/	"About to execute MaxFilterShortCuboid")

	switch (Convention) {
		case FiniteDataSupport:
		case MirrorOffBounds:
		case MirrorOnBounds:
		case Periodic:
			break;
		case AntiMirrorOnBounds:
		case FiniteCoefficientSupport:
		default:
			*Status = ERROR;
			WRITE_ERROR(MaxFilterShortCuboid, "Invalid boundary convention")
/**/		DEBUG_WRITE_LEAVING(MaxFilterShortCuboid, "Done")
			return(*Status);
	}
	if (1L < Nx) {
		AllocateLineDouble(&BufferIn, Nx, Status);
		if (*Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(MaxFilterShortCuboid, "Done")
			return(*Status);
		}
		AllocateLineDouble(&BufferOut, Nx, Status);
		if (*Status == ERROR) {
			FreeLineDouble(&BufferIn);
/**/		DEBUG_WRITE_LEAVING(MaxFilterShortCuboid, "Done")
			return(*Status);
		}
		for (k = 0L; (k < Nz); k++) {
			for (j = 0L; (j < Ny); j++) {
				*Status = GetxShortToDouble(Volume, Nx, Ny, Nz, 0L, j, k, BufferIn, Nx);
				if (*Status == ERROR) {
					FreeLineDouble(&BufferOut);
					FreeLineDouble(&BufferIn);
/**/				DEBUG_WRITE_LEAVING(MaxFilterShortCuboid, "Done")
					return(*Status);
				}
				*Status = MaxFilterLine(BufferIn, BufferOut, Nx, Kx, Convention);
				if (*Status == ERROR) {
					FreeLineDouble(&BufferOut);
					FreeLineDouble(&BufferIn);
/**/				DEBUG_WRITE_LEAVING(MaxFilterShortCuboid, "Done")
					return(*Status);
				}
				*Status = PutxDoubleToShort(Volume, Nx, Ny, Nz, 0L, j, k, BufferOut, Nx);
				if (*Status == ERROR) {
					FreeLineDouble(&BufferOut);
					FreeLineDouble(&BufferIn);
/**/				DEBUG_WRITE_LEAVING(MaxFilterShortCuboid, "Done")
					return(*Status);
				}
			}
		}
		*Status = FreeLineDouble(&BufferOut);
		if (*Status == ERROR) {
			FreeLineDouble(&BufferIn);
/**/		DEBUG_WRITE_LEAVING(MaxFilterShortCuboid, "Done")
			return(*Status);
		}
		*Status = FreeLineDouble(&BufferIn);
		if (*Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(MaxFilterShortCuboid, "Done")
			return(*Status);
		}
	}
	if (1L < Ny) {
		AllocateLineDouble(&BufferIn, Ny, Status);
		if (*Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(MaxFilterShortCuboid, "Done")
			return(*Status);
		}
		AllocateLineDouble(&BufferOut, Ny, Status);
		if (*Status == ERROR) {
			FreeLineDouble(&BufferIn);
/**/		DEBUG_WRITE_LEAVING(MaxFilterShortCuboid, "Done")
			return(*Status);
		}
		for (k = 0L; (k < Nz); k++) {
			for (i = 0L; (i < Nx); i++) {
				*Status = GetyShortToDouble(Volume, Nx, Ny, Nz, i, 0L, k, BufferIn, Ny);
				if (*Status == ERROR) {
					FreeLineDouble(&BufferOut);
					FreeLineDouble(&BufferIn);
/**/				DEBUG_WRITE_LEAVING(MaxFilterShortCuboid, "Done")
					return(*Status);
				}
				*Status = MaxFilterLine(BufferIn, BufferOut, Ny, Ky, Convention);
				if (*Status == ERROR) {
					FreeLineDouble(&BufferOut);
					FreeLineDouble(&BufferIn);
/**/				DEBUG_WRITE_LEAVING(MaxFilterShortCuboid, "Done")
					return(*Status);
				}
				*Status = PutyDoubleToShort(Volume, Nx, Ny, Nz, i, 0L, k, BufferOut, Ny);
				if (*Status == ERROR) {
					FreeLineDouble(&BufferOut);
					FreeLineDouble(&BufferIn);
/**/				DEBUG_WRITE_LEAVING(MaxFilterShortCuboid, "Done")
					return(*Status);
				}
			}
		}
		*Status = FreeLineDouble(&BufferOut);
		if (*Status == ERROR) {
			FreeLineDouble(&BufferIn);
/**/		DEBUG_WRITE_LEAVING(MaxFilterShortCuboid, "Done")
			return(*Status);
		}
		*Status = FreeLineDouble(&BufferIn);
		if (*Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(MaxFilterShortCuboid, "Done")
			return(*Status);
		}
	}
	if (1L < Nz) {
		AllocateLineDouble(&BufferIn, Nz, Status);
		if (*Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(MaxFilterShortCuboid, "Done")
			return(*Status);
		}
		AllocateLineDouble(&BufferOut, Nz, Status);
		if (*Status == ERROR) {
			FreeLineDouble(&BufferIn);
/**/		DEBUG_WRITE_LEAVING(MaxFilterShortCuboid, "Done")
			return(*Status);
		}
		for (j = 0L; (j < Ny); j++) {
			for (i = 0L; (i < Nx); i++) {
				*Status = GetzShortToDouble(Volume, Nx, Ny, Nz, i, j, 0L, BufferIn, Nz);
				if (*Status == ERROR) {
					FreeLineDouble(&BufferOut);
					FreeLineDouble(&BufferIn);
/**/				DEBUG_WRITE_LEAVING(MaxFilterShortCuboid, "Done")
					return(*Status);
				}
				*Status = MaxFilterLine(BufferIn, BufferOut, Nz, Kz, Convention);
				if (*Status == ERROR) {
					FreeLineDouble(&BufferOut);
					FreeLineDouble(&BufferIn);
/**/				DEBUG_WRITE_LEAVING(MaxFilterShortCuboid, "Done")
					return(*Status);
				}
				*Status = PutzDoubleToShort(Volume, Nx, Ny, Nz, i, j, 0L, BufferOut, Nz);
				if (*Status == ERROR) {
					FreeLineDouble(&BufferOut);
					FreeLineDouble(&BufferIn);
/**/				DEBUG_WRITE_LEAVING(MaxFilterShortCuboid, "Done")
					return(*Status);
				}
			}
		}
		*Status = FreeLineDouble(&BufferOut);
		if (*Status == ERROR) {
			FreeLineDouble(&BufferIn);
/**/		DEBUG_WRITE_LEAVING(MaxFilterShortCuboid, "Done")
			return(*Status);
		}
		*Status = FreeLineDouble(&BufferIn);
		if (*Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(MaxFilterShortCuboid, "Done")
			return(*Status);
		}
	}
/**/DEBUG_WRITE_LEAVING(MaxFilterShortCuboid, "Done")
	return(*Status);
} /* end MaxFilterShortCuboid */

/*--------------------------------------------------------------------------*/
extern int		MinFilterFloatCuboid
				(
					float	*Volume,			/* data to process in-place */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					long	Kx,					/* width of the kernel */
					long	Ky,					/* height of the kernel */
					long	Kz,					/* depth of the kernel */
					enum TBoundaryConvention
							Convention,			/* boundary convention */
					int		*Status				/* error management */
				)

/* perform a morphological grey-scale min filter using a cuboid as structural element */
/* Volume is the grey-scale volume to erode in-place */
/* Volume has size (Nx, Ny, Nz) */
/* (Kx, Ky, Kz) is the size of the structural element that is filled with the value 1.0F */
/* the origin of the structural element is ((Kx - 1L) / 2L, (Ky - 1L) / 2L, (Kz - 1L) / 2L) */
/* Convention is the boundary convention applied to Volume before erosion */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin MinFilterFloatCuboid */

	double	*BufferIn = (double *)NULL, *BufferOut = (double *)NULL;
	long	i, j, k;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(MinFilterFloatCuboid, Volume, *Status,
/**/	"Missing Volume")
/**/DEBUG_CHECK_RANGE_LONG(MinFilterFloatCuboid, Nx, 1L, LONG_MAX, *Status,
/**/	"Invalid data width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MinFilterFloatCuboid, Ny, 1L, LONG_MAX, *Status,
/**/	"Invalid data height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MinFilterFloatCuboid, Nz, 1L, LONG_MAX, *Status,
/**/	"Invalid data depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MinFilterFloatCuboid, Kx, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MinFilterFloatCuboid, Ky, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MinFilterFloatCuboid, Kz, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(MinFilterFloatCuboid, *Status)
/**/DEBUG_WRITE_ENTERING(MinFilterFloatCuboid,
/**/	"About to execute MinFilterFloatCuboid")

	switch (Convention) {
		case FiniteDataSupport:
		case MirrorOffBounds:
		case MirrorOnBounds:
		case Periodic:
			break;
		case AntiMirrorOnBounds:
		case FiniteCoefficientSupport:
		default:
			*Status = ERROR;
			WRITE_ERROR(MinFilterFloatCuboid, "Invalid boundary convention")
/**/		DEBUG_WRITE_LEAVING(MinFilterFloatCuboid, "Done")
			return(*Status);
	}
	if (1L < Nx) {
		AllocateLineDouble(&BufferIn, Nx, Status);
		if (*Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(MinFilterFloatCuboid, "Done")
			return(*Status);
		}
		AllocateLineDouble(&BufferOut, Nx, Status);
		if (*Status == ERROR) {
			FreeLineDouble(&BufferIn);
/**/		DEBUG_WRITE_LEAVING(MinFilterFloatCuboid, "Done")
			return(*Status);
		}
		for (k = 0L; (k < Nz); k++) {
			for (j = 0L; (j < Ny); j++) {
				*Status = GetxFloatToDouble(Volume, Nx, Ny, Nz, 0L, j, k, BufferIn, Nx);
				if (*Status == ERROR) {
					FreeLineDouble(&BufferOut);
					FreeLineDouble(&BufferIn);
/**/				DEBUG_WRITE_LEAVING(MinFilterFloatCuboid, "Done")
					return(*Status);
				}
				*Status = MinFilterLine(BufferIn, BufferOut, Nx, Kx, Convention);
				if (*Status == ERROR) {
					FreeLineDouble(&BufferOut);
					FreeLineDouble(&BufferIn);
/**/				DEBUG_WRITE_LEAVING(MinFilterFloatCuboid, "Done")
					return(*Status);
				}
				*Status = PutxDoubleToFloat(Volume, Nx, Ny, Nz, 0L, j, k, BufferOut, Nx);
				if (*Status == ERROR) {
					FreeLineDouble(&BufferOut);
					FreeLineDouble(&BufferIn);
/**/				DEBUG_WRITE_LEAVING(MinFilterFloatCuboid, "Done")
					return(*Status);
				}
			}
		}
		*Status = FreeLineDouble(&BufferOut);
		if (*Status == ERROR) {
			FreeLineDouble(&BufferIn);
/**/		DEBUG_WRITE_LEAVING(MinFilterFloatCuboid, "Done")
			return(*Status);
		}
		*Status = FreeLineDouble(&BufferIn);
		if (*Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(MinFilterFloatCuboid, "Done")
			return(*Status);
		}
	}
	if (1L < Ny) {
		AllocateLineDouble(&BufferIn, Ny, Status);
		if (*Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(MinFilterFloatCuboid, "Done")
			return(*Status);
		}
		AllocateLineDouble(&BufferOut, Ny, Status);
		if (*Status == ERROR) {
			FreeLineDouble(&BufferIn);
/**/		DEBUG_WRITE_LEAVING(MinFilterFloatCuboid, "Done")
			return(*Status);
		}
		for (k = 0L; (k < Nz); k++) {
			for (i = 0L; (i < Nx); i++) {
				*Status = GetyFloatToDouble(Volume, Nx, Ny, Nz, i, 0L, k, BufferIn, Ny);
				if (*Status == ERROR) {
					FreeLineDouble(&BufferOut);
					FreeLineDouble(&BufferIn);
/**/				DEBUG_WRITE_LEAVING(MinFilterFloatCuboid, "Done")
					return(*Status);
				}
				*Status = MinFilterLine(BufferIn, BufferOut, Ny, Ky, Convention);
				if (*Status == ERROR) {
					FreeLineDouble(&BufferOut);
					FreeLineDouble(&BufferIn);
/**/				DEBUG_WRITE_LEAVING(MinFilterFloatCuboid, "Done")
					return(*Status);
				}
				*Status = PutyDoubleToFloat(Volume, Nx, Ny, Nz, i, 0L, k, BufferOut, Ny);
				if (*Status == ERROR) {
					FreeLineDouble(&BufferOut);
					FreeLineDouble(&BufferIn);
/**/				DEBUG_WRITE_LEAVING(MinFilterFloatCuboid, "Done")
					return(*Status);
				}
			}
		}
		*Status = FreeLineDouble(&BufferOut);
		if (*Status == ERROR) {
			FreeLineDouble(&BufferIn);
/**/		DEBUG_WRITE_LEAVING(MinFilterFloatCuboid, "Done")
			return(*Status);
		}
		*Status = FreeLineDouble(&BufferIn);
		if (*Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(MinFilterFloatCuboid, "Done")
			return(*Status);
		}
	}
	if (1L < Nz) {
		AllocateLineDouble(&BufferIn, Nz, Status);
		if (*Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(MinFilterFloatCuboid, "Done")
			return(*Status);
		}
		AllocateLineDouble(&BufferOut, Nz, Status);
		if (*Status == ERROR) {
			FreeLineDouble(&BufferIn);
/**/		DEBUG_WRITE_LEAVING(MinFilterFloatCuboid, "Done")
			return(*Status);
		}
		for (j = 0L; (j < Ny); j++) {
			for (i = 0L; (i < Nx); i++) {
				*Status = GetzFloatToDouble(Volume, Nx, Ny, Nz, i, j, 0L, BufferIn, Nz);
				if (*Status == ERROR) {
					FreeLineDouble(&BufferOut);
					FreeLineDouble(&BufferIn);
/**/				DEBUG_WRITE_LEAVING(MinFilterFloatCuboid, "Done")
					return(*Status);
				}
				*Status = MinFilterLine(BufferIn, BufferOut, Nz, Kz, Convention);
				if (*Status == ERROR) {
					FreeLineDouble(&BufferOut);
					FreeLineDouble(&BufferIn);
/**/				DEBUG_WRITE_LEAVING(MinFilterFloatCuboid, "Done")
					return(*Status);
				}
				*Status = PutzDoubleToFloat(Volume, Nx, Ny, Nz, i, j, 0L, BufferOut, Nz);
				if (*Status == ERROR) {
					FreeLineDouble(&BufferOut);
					FreeLineDouble(&BufferIn);
/**/				DEBUG_WRITE_LEAVING(MinFilterFloatCuboid, "Done")
					return(*Status);
				}
			}
		}
		*Status = FreeLineDouble(&BufferOut);
		if (*Status == ERROR) {
			FreeLineDouble(&BufferIn);
/**/		DEBUG_WRITE_LEAVING(MinFilterFloatCuboid, "Done")
			return(*Status);
		}
		*Status = FreeLineDouble(&BufferIn);
		if (*Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(MinFilterFloatCuboid, "Done")
			return(*Status);
		}
	}
/**/DEBUG_WRITE_LEAVING(MinFilterFloatCuboid, "Done")
	return(*Status);
} /* end MinFilterFloatCuboid */

/*--------------------------------------------------------------------------*/
extern int		MinFilterShortCuboid
				(
					short	*Volume,			/* data to process in-place */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					long	Kx,					/* width of the kernel */
					long	Ky,					/* height of the kernel */
					long	Kz,					/* depth of the kernel */
					enum TBoundaryConvention
							Convention,			/* boundary convention */
					int		*Status				/* error management */
				)

/* perform a morphological grey-scale min filter using a cuboid as structural element */
/* Volume is the grey-scale volume to erode in-place */
/* Volume has size (Nx, Ny, Nz) */
/* (Kx, Ky, Kz) is the size of the structural element that is filled with the value 1.0F */
/* the origin of the structural element is ((Kx - 1L) / 2L, (Ky - 1L) / 2L, (Kz - 1L) / 2L) */
/* Convention is the boundary convention applied to Volume before erosion */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin MinFilterShortCuboid */

	double	*BufferIn = (double *)NULL, *BufferOut = (double *)NULL;
	long	i, j, k;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(MinFilterShortCuboid, Volume, *Status,
/**/	"Missing Volume")
/**/DEBUG_CHECK_RANGE_LONG(MinFilterShortCuboid, Nx, 1L, LONG_MAX, *Status,
/**/	"Invalid data width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MinFilterShortCuboid, Ny, 1L, LONG_MAX, *Status,
/**/	"Invalid data height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MinFilterShortCuboid, Nz, 1L, LONG_MAX, *Status,
/**/	"Invalid data depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MinFilterShortCuboid, Kx, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MinFilterShortCuboid, Ky, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MinFilterShortCuboid, Kz, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(MinFilterShortCuboid, *Status)
/**/DEBUG_WRITE_ENTERING(MinFilterShortCuboid,
/**/	"About to execute MinFilterShortCuboid")

	switch (Convention) {
		case FiniteDataSupport:
		case MirrorOffBounds:
		case MirrorOnBounds:
		case Periodic:
			break;
		case AntiMirrorOnBounds:
		case FiniteCoefficientSupport:
		default:
			*Status = ERROR;
			WRITE_ERROR(MinFilterShortCuboid, "Invalid boundary convention")
/**/		DEBUG_WRITE_LEAVING(MinFilterShortCuboid, "Done")
			return(*Status);
	}
	if (1L < Nx) {
		AllocateLineDouble(&BufferIn, Nx, Status);
		if (*Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(MinFilterShortCuboid, "Done")
			return(*Status);
		}
		AllocateLineDouble(&BufferOut, Nx, Status);
		if (*Status == ERROR) {
			FreeLineDouble(&BufferIn);
/**/		DEBUG_WRITE_LEAVING(MinFilterShortCuboid, "Done")
			return(*Status);
		}
		for (k = 0L; (k < Nz); k++) {
			for (j = 0L; (j < Ny); j++) {
				*Status = GetxShortToDouble(Volume, Nx, Ny, Nz, 0L, j, k, BufferIn, Nx);
				if (*Status == ERROR) {
					FreeLineDouble(&BufferOut);
					FreeLineDouble(&BufferIn);
/**/				DEBUG_WRITE_LEAVING(MinFilterShortCuboid, "Done")
					return(*Status);
				}
				*Status = MinFilterLine(BufferIn, BufferOut, Nx, Kx, Convention);
				if (*Status == ERROR) {
					FreeLineDouble(&BufferOut);
					FreeLineDouble(&BufferIn);
/**/				DEBUG_WRITE_LEAVING(MinFilterShortCuboid, "Done")
					return(*Status);
				}
				*Status = PutxDoubleToShort(Volume, Nx, Ny, Nz, 0L, j, k, BufferOut, Nx);
				if (*Status == ERROR) {
					FreeLineDouble(&BufferOut);
					FreeLineDouble(&BufferIn);
/**/				DEBUG_WRITE_LEAVING(MinFilterShortCuboid, "Done")
					return(*Status);
				}
			}
		}
		*Status = FreeLineDouble(&BufferOut);
		if (*Status == ERROR) {
			FreeLineDouble(&BufferIn);
/**/		DEBUG_WRITE_LEAVING(MinFilterShortCuboid, "Done")
			return(*Status);
		}
		*Status = FreeLineDouble(&BufferIn);
		if (*Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(MinFilterShortCuboid, "Done")
			return(*Status);
		}
	}
	if (1L < Ny) {
		AllocateLineDouble(&BufferIn, Ny, Status);
		if (*Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(MinFilterShortCuboid, "Done")
			return(*Status);
		}
		AllocateLineDouble(&BufferOut, Ny, Status);
		if (*Status == ERROR) {
			FreeLineDouble(&BufferIn);
/**/		DEBUG_WRITE_LEAVING(MinFilterShortCuboid, "Done")
			return(*Status);
		}
		for (k = 0L; (k < Nz); k++) {
			for (i = 0L; (i < Nx); i++) {
				*Status = GetyShortToDouble(Volume, Nx, Ny, Nz, i, 0L, k, BufferIn, Ny);
				if (*Status == ERROR) {
					FreeLineDouble(&BufferOut);
					FreeLineDouble(&BufferIn);
/**/				DEBUG_WRITE_LEAVING(MinFilterShortCuboid, "Done")
					return(*Status);
				}
				*Status = MinFilterLine(BufferIn, BufferOut, Ny, Ky, Convention);
				if (*Status == ERROR) {
					FreeLineDouble(&BufferOut);
					FreeLineDouble(&BufferIn);
/**/				DEBUG_WRITE_LEAVING(MinFilterShortCuboid, "Done")
					return(*Status);
				}
				*Status = PutyDoubleToShort(Volume, Nx, Ny, Nz, i, 0L, k, BufferOut, Ny);
				if (*Status == ERROR) {
					FreeLineDouble(&BufferOut);
					FreeLineDouble(&BufferIn);
/**/				DEBUG_WRITE_LEAVING(MinFilterShortCuboid, "Done")
					return(*Status);
				}
			}
		}
		*Status = FreeLineDouble(&BufferOut);
		if (*Status == ERROR) {
			FreeLineDouble(&BufferIn);
/**/		DEBUG_WRITE_LEAVING(MinFilterShortCuboid, "Done")
			return(*Status);
		}
		*Status = FreeLineDouble(&BufferIn);
		if (*Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(MinFilterShortCuboid, "Done")
			return(*Status);
		}
	}
	if (1L < Nz) {
		AllocateLineDouble(&BufferIn, Nz, Status);
		if (*Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(MinFilterShortCuboid, "Done")
			return(*Status);
		}
		AllocateLineDouble(&BufferOut, Nz, Status);
		if (*Status == ERROR) {
			FreeLineDouble(&BufferIn);
/**/		DEBUG_WRITE_LEAVING(MinFilterShortCuboid, "Done")
			return(*Status);
		}
		for (j = 0L; (j < Ny); j++) {
			for (i = 0L; (i < Nx); i++) {
				*Status = GetzShortToDouble(Volume, Nx, Ny, Nz, i, j, 0L, BufferIn, Nz);
				if (*Status == ERROR) {
					FreeLineDouble(&BufferOut);
					FreeLineDouble(&BufferIn);
/**/				DEBUG_WRITE_LEAVING(MinFilterShortCuboid, "Done")
					return(*Status);
				}
				*Status = MinFilterLine(BufferIn, BufferOut, Nz, Kz, Convention);
				if (*Status == ERROR) {
					FreeLineDouble(&BufferOut);
					FreeLineDouble(&BufferIn);
/**/				DEBUG_WRITE_LEAVING(MinFilterShortCuboid, "Done")
					return(*Status);
				}
				*Status = PutzDoubleToShort(Volume, Nx, Ny, Nz, i, j, 0L, BufferOut, Nz);
				if (*Status == ERROR) {
					FreeLineDouble(&BufferOut);
					FreeLineDouble(&BufferIn);
/**/				DEBUG_WRITE_LEAVING(MinFilterShortCuboid, "Done")
					return(*Status);
				}
			}
		}
		*Status = FreeLineDouble(&BufferOut);
		if (*Status == ERROR) {
			FreeLineDouble(&BufferIn);
/**/		DEBUG_WRITE_LEAVING(MinFilterShortCuboid, "Done")
			return(*Status);
		}
		*Status = FreeLineDouble(&BufferIn);
		if (*Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(MinFilterShortCuboid, "Done")
			return(*Status);
		}
	}
/**/DEBUG_WRITE_LEAVING(MinFilterShortCuboid, "Done")
	return(*Status);
} /* end MinFilterShortCuboid */

/*--------------------------------------------------------------------------*/
extern int		MorphologicalGradientFloat
				(
					float	*VolumeSource,		/* data to process */
					float	*VolumeDestination,	/* result */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					float	*Kernel,			/* structural element */
					long	Kx,					/* width of the kernel */
					long	Ky,					/* height of the kernel */
					long	Kz,					/* depth of the kernel */
					long	Ox,					/* kernel X origin */
					long	Oy,					/* kernel Y origin */
					long	Oz,					/* kernel Z origin */
					enum TBoundaryConvention
							Convention,			/* boundary convention */
					int		*Status				/* error management */
				)

/* perform the grey-scale morphological operation known as morphological gradient */
/* VolumeSource is the grey-scale volume to process */
/* VolumeDestination is the resulting grey-scale volume after the gradient operation */
/* both VolumeSource and VolumeDestination have size (Nx, Ny, Nz) */
/* Kernel is the grey-scale structural element by which VolumeSource is processed */
/* (Kx, Ky, Kz) is the size (bounding box) of the structural element */
/* (Ox, Oy, Oz) is the origin of the structural element */
/* Convention is the boundary convention applied to VolumeSource before the gradient operation */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin MorphologicalGradientFloat */

	float	*p, *q;
	float	*Dilated = (float *)NULL, *Eroded = (float *)NULL;
	long	x, y, z;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(MorphologicalGradientFloat, VolumeSource, *Status,
/**/	"Missing VolumeSource")
/**/DEBUG_CHECK_NULL_POINTER(MorphologicalGradientFloat, VolumeDestination, *Status,
/**/	"Missing VolumeDestination")
/**/DEBUG_CHECK_NULL_POINTER(MorphologicalGradientFloat, Kernel, *Status,
/**/	"Missing Kernel")
/**/DEBUG_CHECK_RANGE_LONG(MorphologicalGradientFloat, Nx, 1L, LONG_MAX, *Status,
/**/	"Invalid data width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MorphologicalGradientFloat, Ny, 1L, LONG_MAX, *Status,
/**/	"Invalid data height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MorphologicalGradientFloat, Nz, 1L, LONG_MAX, *Status,
/**/	"Invalid data depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MorphologicalGradientFloat, Kx, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MorphologicalGradientFloat, Ky, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MorphologicalGradientFloat, Kz, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(MorphologicalGradientFloat, *Status)
/**/DEBUG_WRITE_ENTERING(MorphologicalGradientFloat,
/**/	"About to execute MorphologicalGradientFloat")

	switch (Convention) {
		case FiniteDataSupport:
		case MirrorOffBounds:
		case MirrorOnBounds:
		case Periodic:
			break;
		case AntiMirrorOnBounds:
		case FiniteCoefficientSupport:
		default:
			*Status = ERROR;
			WRITE_ERROR(MorphologicalGradientFloat, "Invalid boundary convention")
/**/		DEBUG_WRITE_LEAVING(MorphologicalGradientFloat, "Done")
			return(*Status);
	}
	AllocateVolumeFloat(&Dilated, Nx, Ny, Nz, Status);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(MorphologicalGradientFloat, "Done")
		return(*Status);
	}
	*Status = DilateFloat(VolumeSource, Dilated, Nx, Ny, Nz, Kernel, Kx, Ky, Kz,
		Ox, Oy, Oz, Convention);
	if (*Status == ERROR) {
		FreeVolumeFloat(&Dilated);
/**/	DEBUG_WRITE_LEAVING(MorphologicalGradientFloat, "Done")
		return(*Status);
	}
	AllocateVolumeFloat(&Eroded, Nx, Ny, Nz, Status);
	if (*Status == ERROR) {
		FreeVolumeFloat(&Dilated);
/**/	DEBUG_WRITE_LEAVING(MorphologicalGradientFloat, "Done")
		return(*Status);
	}
	*Status = ErodeFloat(VolumeSource, Eroded, Nx, Ny, Nz, Kernel, Kx, Ky, Kz,
		Ox, Oy, Oz, Convention);
	if (*Status == ERROR) {
		FreeVolumeFloat(&Eroded);
		FreeVolumeFloat(&Dilated);
/**/	DEBUG_WRITE_LEAVING(MorphologicalGradientFloat, "Done")
		return(*Status);
	}
	p = Dilated;
	q = Eroded;
	for (z = -Nz; (z < 0L); z++) {
		for (y = -Ny; (y < 0L); y++) {
			for (x = -Nx; (x < 0L); x++) {
				*VolumeDestination++ = *p++ - *q++;
			}
		}
	}
	*Status = FreeVolumeFloat(&Eroded);
	if (*Status == ERROR) {
		FreeVolumeFloat(&Dilated);
/**/	DEBUG_WRITE_LEAVING(MorphologicalGradientFloat, "Done")
		return(*Status);
	}
	*Status = FreeVolumeFloat(&Dilated);
/**/DEBUG_WRITE_LEAVING(MorphologicalGradientFloat, "Done")
	return(*Status);
} /* end MorphologicalGradientFloat */

/*--------------------------------------------------------------------------*/
extern int		MorphologicalGradientFloatCuboid
				(
					float	*Volume,			/* data to process in-place */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					long	Kx,					/* width of the kernel */
					long	Ky,					/* height of the kernel */
					long	Kz,					/* depth of the kernel */
					enum TBoundaryConvention
							Convention,			/* boundary convention */
					int		*Status				/* error management */
				)

/* perform a morphological grey-scale gradient filter using a cuboid as structural element */
/* Volume is the grey-scale volume to process in-place */
/* Volume has size (Nx, Ny, Nz) */
/* (Kx, Ky, Kz) is the size of the structural element that is filled with the value 1.0F */
/* the origin of the structural element is ((Kx - 1L) / 2L, (Ky - 1L) / 2L, (Kz - 1L) / 2L) */
/* Convention is the boundary convention applied to Volume before processing */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin MorphologicalGradientFloatCuboid */

	float	*p, *q;
	float	*Dilated = (float *)NULL, *Eroded = (float *)NULL;
	long	x, y, z;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(MorphologicalGradientFloatCuboid, Volume, *Status,
/**/	"Missing Volume")
/**/DEBUG_CHECK_RANGE_LONG(MorphologicalGradientFloatCuboid, Nx, 1L, LONG_MAX, *Status,
/**/	"Invalid data width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MorphologicalGradientFloatCuboid, Ny, 1L, LONG_MAX, *Status,
/**/	"Invalid data height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MorphologicalGradientFloatCuboid, Nz, 1L, LONG_MAX, *Status,
/**/	"Invalid data depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MorphologicalGradientFloatCuboid, Kx, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MorphologicalGradientFloatCuboid, Ky, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MorphologicalGradientFloatCuboid, Kz, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(MorphologicalGradientFloatCuboid, *Status)
/**/DEBUG_WRITE_ENTERING(MorphologicalGradientFloatCuboid,
/**/	"About to execute MorphologicalGradientFloatCuboid")

	switch (Convention) {
		case FiniteDataSupport:
		case MirrorOffBounds:
		case MirrorOnBounds:
		case Periodic:
			break;
		case AntiMirrorOnBounds:
		case FiniteCoefficientSupport:
		default:
			*Status = ERROR;
			WRITE_ERROR(MorphologicalGradientFloatCuboid, "Invalid boundary convention")
/**/		DEBUG_WRITE_LEAVING(MorphologicalGradientFloatCuboid, "Done")
			return(*Status);
	}
	AllocateVolumeFloat(&Dilated, Nx, Ny, Nz, Status);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(MorphologicalGradientFloatCuboid, "Done")
		return(*Status);
	}
	memcpy(Dilated, Volume, (size_t)(Nx * Ny * Nz * (long)sizeof(float)));
	MaxFilterFloatCuboid(Dilated, Nx, Ny, Nz, Kx, Ky, Kz, Convention, Status);
	if (*Status == ERROR) {
		FreeVolumeFloat(&Dilated);
/**/	DEBUG_WRITE_LEAVING(MorphologicalGradientFloatCuboid, "Done")
		return(*Status);
	}
	AllocateVolumeFloat(&Eroded, Nx, Ny, Nz, Status);
	if (*Status == ERROR) {
		FreeVolumeFloat(&Dilated);
/**/	DEBUG_WRITE_LEAVING(MorphologicalGradientFloatCuboid, "Done")
		return(*Status);
	}
	memcpy(Eroded, Volume, (size_t)(Nx * Ny * Nz * (long)sizeof(float)));
	MinFilterFloatCuboid(Eroded, Nx, Ny, Nz, Kx, Ky, Kz, Convention, Status);
	if (*Status == ERROR) {
		FreeVolumeFloat(&Eroded);
		FreeVolumeFloat(&Dilated);
/**/	DEBUG_WRITE_LEAVING(MorphologicalGradientFloatCuboid, "Done")
		return(*Status);
	}
	p = Dilated;
	q = Eroded;
	for (z = -Nz; (z < 0L); z++) {
		for (y = -Ny; (y < 0L); y++) {
			for (x = -Nx; (x < 0L); x++) {
				*Volume++ = *p++ - *q++;
			}
		}
	}
	*Status = FreeVolumeFloat(&Eroded);
	if (*Status == ERROR) {
		FreeVolumeFloat(&Dilated);
/**/	DEBUG_WRITE_LEAVING(MorphologicalGradientFloatCuboid, "Done")
		return(*Status);
	}
	*Status = FreeVolumeFloat(&Dilated);
/**/DEBUG_WRITE_LEAVING(MorphologicalGradientFloatCuboid, "Done")
	return(*Status);
} /* end MorphologicalGradientFloatCuboid */

/*--------------------------------------------------------------------------*/
extern int		MorphologicalGradientShort
				(
					short	*VolumeSource,		/* data to process */
					short	*VolumeDestination,	/* result */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					short	*Kernel,			/* structural element */
					long	Kx,					/* width of the kernel */
					long	Ky,					/* height of the kernel */
					long	Kz,					/* depth of the kernel */
					long	Ox,					/* kernel X origin */
					long	Oy,					/* kernel Y origin */
					long	Oz,					/* kernel Z origin */
					enum TBoundaryConvention
							Convention,			/* boundary convention */
					int		*Status				/* error management */
				)

/* perform the grey-scale morphological operation known as morphological gradient */
/* VolumeSource is the grey-scale volume to process */
/* VolumeDestination is the resulting grey-scale volume after the gradient operation */
/* both VolumeSource and VolumeDestination have size (Nx, Ny, Nz) */
/* Kernel is the grey-scale structural element by which VolumeSource is processed */
/* (Kx, Ky, Kz) is the size (bounding box) of the structural element */
/* (Ox, Oy, Oz) is the origin of the structural element */
/* Convention is the boundary convention applied to VolumeSource before the gradient operation */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin MorphologicalGradientShort */

	short	*p, *q;
	short	*Dilated = (short *)NULL, *Eroded = (short *)NULL;
	long	x, y, z;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(MorphologicalGradientShort, VolumeSource, *Status,
/**/	"Missing VolumeSource")
/**/DEBUG_CHECK_NULL_POINTER(MorphologicalGradientShort, VolumeDestination, *Status,
/**/	"Missing VolumeDestination")
/**/DEBUG_CHECK_NULL_POINTER(MorphologicalGradientShort, Kernel, *Status,
/**/	"Missing Kernel")
/**/DEBUG_CHECK_RANGE_LONG(MorphologicalGradientShort, Nx, 1L, LONG_MAX, *Status,
/**/	"Invalid data width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MorphologicalGradientShort, Ny, 1L, LONG_MAX, *Status,
/**/	"Invalid data height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MorphologicalGradientShort, Nz, 1L, LONG_MAX, *Status,
/**/	"Invalid data depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MorphologicalGradientShort, Kx, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MorphologicalGradientShort, Ky, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MorphologicalGradientShort, Kz, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(MorphologicalGradientShort, *Status)
/**/DEBUG_WRITE_ENTERING(MorphologicalGradientShort,
/**/	"About to execute MorphologicalGradientShort")

	switch (Convention) {
		case FiniteDataSupport:
		case MirrorOffBounds:
		case MirrorOnBounds:
		case Periodic:
			break;
		case AntiMirrorOnBounds:
		case FiniteCoefficientSupport:
		default:
			*Status = ERROR;
			WRITE_ERROR(MorphologicalGradientShort, "Invalid boundary convention")
/**/		DEBUG_WRITE_LEAVING(MorphologicalGradientShort, "Done")
			return(*Status);
	}
	AllocateVolumeShort(&Dilated, Nx, Ny, Nz, Status);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(MorphologicalGradientShort, "Done")
		return(*Status);
	}
	*Status = DilateShort(VolumeSource, Dilated, Nx, Ny, Nz, Kernel, Kx, Ky, Kz,
		Ox, Oy, Oz, Convention);
	if (*Status == ERROR) {
		FreeVolumeShort(&Dilated);
/**/	DEBUG_WRITE_LEAVING(MorphologicalGradientShort, "Done")
		return(*Status);
	}
	AllocateVolumeShort(&Eroded, Nx, Ny, Nz, Status);
	if (*Status == ERROR) {
		FreeVolumeShort(&Dilated);
/**/	DEBUG_WRITE_LEAVING(MorphologicalGradientShort, "Done")
		return(*Status);
	}
	*Status = ErodeShort(VolumeSource, Eroded, Nx, Ny, Nz, Kernel, Kx, Ky, Kz,
		Ox, Oy, Oz, Convention);
	if (*Status == ERROR) {
		FreeVolumeShort(&Eroded);
		FreeVolumeShort(&Dilated);
/**/	DEBUG_WRITE_LEAVING(MorphologicalGradientShort, "Done")
		return(*Status);
	}
	p = Dilated;
	q = Eroded;
	for (z = -Nz; (z < 0L); z++) {
		for (y = -Ny; (y < 0L); y++) {
			for (x = -Nx; (x < 0L); x++) {
				*VolumeDestination++ = ConvertIntToShort((int)*p++ - (int)*q++);
			}
		}
	}
	*Status = FreeVolumeShort(&Eroded);
	if (*Status == ERROR) {
		FreeVolumeShort(&Dilated);
/**/	DEBUG_WRITE_LEAVING(MorphologicalGradientShort, "Done")
		return(*Status);
	}
	*Status = FreeVolumeShort(&Dilated);
/**/DEBUG_WRITE_LEAVING(MorphologicalGradientShort, "Done")
	return(*Status);
} /* end MorphologicalGradientShort */

/*--------------------------------------------------------------------------*/
extern int		MorphologicalGradientShortCuboid
				(
					short	*Volume,			/* data to process in-place */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					long	Kx,					/* width of the kernel */
					long	Ky,					/* height of the kernel */
					long	Kz,					/* depth of the kernel */
					enum TBoundaryConvention
							Convention,			/* boundary convention */
					int		*Status				/* error management */
				)

/* perform a morphological grey-scale gradient filter using a cuboid as structural element */
/* Volume is the grey-scale volume to process in-place */
/* Volume has size (Nx, Ny, Nz) */
/* (Kx, Ky, Kz) is the size of the structural element that is filled with the value 1.0F */
/* the origin of the structural element is ((Kx - 1L) / 2L, (Ky - 1L) / 2L, (Kz - 1L) / 2L) */
/* Convention is the boundary convention applied to Volume before processing */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin MorphologicalGradientShortCuboid */

	short	*p, *q;
	short	*Dilated = (short *)NULL, *Eroded = (short *)NULL;
	long	x, y, z;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(MorphologicalGradientShortCuboid, Volume, *Status,
/**/	"Missing Volume")
/**/DEBUG_CHECK_RANGE_LONG(MorphologicalGradientShortCuboid, Nx, 1L, LONG_MAX, *Status,
/**/	"Invalid data width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MorphologicalGradientShortCuboid, Ny, 1L, LONG_MAX, *Status,
/**/	"Invalid data height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MorphologicalGradientShortCuboid, Nz, 1L, LONG_MAX, *Status,
/**/	"Invalid data depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MorphologicalGradientShortCuboid, Kx, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MorphologicalGradientShortCuboid, Ky, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MorphologicalGradientShortCuboid, Kz, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(MorphologicalGradientShortCuboid, *Status)
/**/DEBUG_WRITE_ENTERING(MorphologicalGradientShortCuboid,
/**/	"About to execute MorphologicalGradientShortCuboid")

	switch (Convention) {
		case FiniteDataSupport:
		case MirrorOffBounds:
		case MirrorOnBounds:
		case Periodic:
			break;
		case AntiMirrorOnBounds:
		case FiniteCoefficientSupport:
		default:
			*Status = ERROR;
			WRITE_ERROR(MorphologicalGradientShortCuboid, "Invalid boundary convention")
/**/		DEBUG_WRITE_LEAVING(MorphologicalGradientShortCuboid, "Done")
			return(*Status);
	}
	AllocateVolumeShort(&Dilated, Nx, Ny, Nz, Status);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(MorphologicalGradientShortCuboid, "Done")
		return(*Status);
	}
	memcpy(Dilated, Volume, (size_t)(Nx * Ny * Nz * (long)sizeof(short)));
	MaxFilterShortCuboid(Dilated, Nx, Ny, Nz, Kx, Ky, Kz, Convention, Status);
	if (*Status == ERROR) {
		FreeVolumeShort(&Dilated);
/**/	DEBUG_WRITE_LEAVING(MorphologicalGradientShortCuboid, "Done")
		return(*Status);
	}
	AllocateVolumeShort(&Eroded, Nx, Ny, Nz, Status);
	if (*Status == ERROR) {
		FreeVolumeShort(&Dilated);
/**/	DEBUG_WRITE_LEAVING(MorphologicalGradientShortCuboid, "Done")
		return(*Status);
	}
	memcpy(Eroded, Volume, (size_t)(Nx * Ny * Nz * (long)sizeof(short)));
	MinFilterShortCuboid(Eroded, Nx, Ny, Nz, Kx, Ky, Kz, Convention, Status);
	if (*Status == ERROR) {
		FreeVolumeShort(&Eroded);
		FreeVolumeShort(&Dilated);
/**/	DEBUG_WRITE_LEAVING(MorphologicalGradientShortCuboid, "Done")
		return(*Status);
	}
	p = Dilated;
	q = Eroded;
	for (z = -Nz; (z < 0L); z++) {
		for (y = -Ny; (y < 0L); y++) {
			for (x = -Nx; (x < 0L); x++) {
				*Volume++ = ConvertIntToShort((int)*p++ - (int)*q++);
			}
		}
	}
	*Status = FreeVolumeShort(&Eroded);
	if (*Status == ERROR) {
		FreeVolumeShort(&Dilated);
/**/	DEBUG_WRITE_LEAVING(MorphologicalGradientShortCuboid, "Done")
		return(*Status);
	}
	*Status = FreeVolumeShort(&Dilated);
/**/DEBUG_WRITE_LEAVING(MorphologicalGradientShortCuboid, "Done")
	return(*Status);
} /* end MorphologicalGradientShortCuboid */

/*--------------------------------------------------------------------------*/
extern int		OpeningFloat
				(
					float	*VolumeSource,		/* data to process */
					float	*VolumeDestination,	/* result */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					float	*Kernel,			/* structural element */
					long	Kx,					/* width of the kernel */
					long	Ky,					/* height of the kernel */
					long	Kz,					/* depth of the kernel */
					long	Ox,					/* kernel X origin */
					long	Oy,					/* kernel Y origin */
					long	Oz,					/* kernel Z origin */
					enum TBoundaryConvention
							Convention,			/* boundary convention */
					int		*Status				/* error management */
				)

/* perform the grey-scale morphological operation known as opening */
/* VolumeSource is the grey-scale volume to open */
/* VolumeDestination is the resulting grey-scale volume after opening */
/* both VolumeSource and VolumeDestination have size (Nx, Ny, Nz) */
/* Kernel is the grey-scale structural element by which VolumeSource is opened */
/* (Kx, Ky, Kz) is the size (bounding box) of the structural element */
/* (Ox, Oy, Oz) is the origin of the structural element */
/* Convention is the boundary convention applied to VolumeSource before opening */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin OpeningFloat */

	float	*Buffer = (float *)NULL;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(OpeningFloat, VolumeSource, *Status,
/**/	"Missing VolumeSource")
/**/DEBUG_CHECK_NULL_POINTER(OpeningFloat, VolumeDestination, *Status,
/**/	"Missing VolumeDestination")
/**/DEBUG_CHECK_NULL_POINTER(OpeningFloat, Kernel, *Status,
/**/	"Missing Kernel")
/**/DEBUG_CHECK_RANGE_LONG(OpeningFloat, Nx, 1L, LONG_MAX, *Status,
/**/	"Invalid data width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(OpeningFloat, Ny, 1L, LONG_MAX, *Status,
/**/	"Invalid data height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(OpeningFloat, Nz, 1L, LONG_MAX, *Status,
/**/	"Invalid data depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(OpeningFloat, Kx, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(OpeningFloat, Ky, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(OpeningFloat, Kz, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(OpeningFloat, *Status)
/**/DEBUG_WRITE_ENTERING(OpeningFloat,
/**/	"About to execute OpeningFloat")

	switch (Convention) {
		case FiniteDataSupport:
		case MirrorOffBounds:
		case MirrorOnBounds:
		case Periodic:
			break;
		case AntiMirrorOnBounds:
		case FiniteCoefficientSupport:
		default:
			*Status = ERROR;
			WRITE_ERROR(OpeningFloat, "Invalid boundary convention")
/**/		DEBUG_WRITE_LEAVING(OpeningFloat, "Done")
			return(*Status);
	}
	AllocateVolumeFloat(&Buffer, Nx, Ny, Nz, Status);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(OpeningFloat, "Done")
		return(*Status);
	}
	*Status = ErodeFloat(VolumeSource, Buffer, Nx, Ny, Nz, Kernel, Kx, Ky, Kz,
		Ox, Oy, Oz, Convention);
	if (*Status == ERROR) {
		FreeVolumeFloat(&Buffer);
/**/	DEBUG_WRITE_LEAVING(OpeningFloat, "Done")
		return(*Status);
	}
	*Status = DilateFloat(Buffer, VolumeDestination, Nx, Ny, Nz, Kernel, Kx, Ky, Kz,
		Ox, Oy, Oz, Convention);
	if (*Status == ERROR) {
		FreeVolumeFloat(&Buffer);
/**/	DEBUG_WRITE_LEAVING(OpeningFloat, "Done")
		return(*Status);
	}
	*Status = FreeVolumeFloat(&Buffer);
/**/DEBUG_WRITE_LEAVING(OpeningFloat, "Done")
	return(*Status);
} /* end OpeningFloat */

/*--------------------------------------------------------------------------*/
extern int		OpeningFloatCuboid
				(
					float	*Volume,			/* data to process in-place */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					long	Kx,					/* width of the kernel */
					long	Ky,					/* height of the kernel */
					long	Kz,					/* depth of the kernel */
					enum TBoundaryConvention
							Convention,			/* boundary convention */
					int		*Status				/* error management */
				)

/* perform a morphological grey-scale opening filter using a cuboid as structural element */
/* Volume is the grey-scale volume to open in-place */
/* Volume has size (Nx, Ny, Nz) */
/* (Kx, Ky, Kz) is the size of the structural element that is filled with the value 1.0F */
/* the origin of the structural element is ((Kx - 1L) / 2L, (Ky - 1L) / 2L, (Kz - 1L) / 2L) */
/* Convention is the boundary convention applied to Volume before opening */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin OpeningFloatCuboid */

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(OpeningFloatCuboid, Volume, *Status,
/**/	"Missing Volume")
/**/DEBUG_CHECK_RANGE_LONG(OpeningFloatCuboid, Nx, 1L, LONG_MAX, *Status,
/**/	"Invalid data width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(OpeningFloatCuboid, Ny, 1L, LONG_MAX, *Status,
/**/	"Invalid data height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(OpeningFloatCuboid, Nz, 1L, LONG_MAX, *Status,
/**/	"Invalid data depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(OpeningFloatCuboid, Kx, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(OpeningFloatCuboid, Ky, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(OpeningFloatCuboid, Kz, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(OpeningFloatCuboid, *Status)
/**/DEBUG_WRITE_ENTERING(OpeningFloatCuboid,
/**/	"About to execute OpeningFloatCuboid")

	switch (Convention) {
		case FiniteDataSupport:
		case MirrorOffBounds:
		case MirrorOnBounds:
		case Periodic:
			break;
		case AntiMirrorOnBounds:
		case FiniteCoefficientSupport:
		default:
			*Status = ERROR;
			WRITE_ERROR(OpeningFloatCuboid, "Invalid boundary convention")
/**/		DEBUG_WRITE_LEAVING(OpeningFloatCuboid, "Done")
			return(*Status);
	}
	MinFilterFloatCuboid(Volume, Nx, Ny, Nz, Kx, Ky, Kz, Convention, Status);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(OpeningFloatCuboid, "Done")
		return(*Status);
	}
	MaxFilterFloatCuboid(Volume, Nx, Ny, Nz, Kx, Ky, Kz, Convention, Status);
/**/DEBUG_WRITE_LEAVING(OpeningFloatCuboid, "Done")
	return(*Status);
} /* end OpeningFloatCuboid */

/*--------------------------------------------------------------------------*/
extern int		OpeningShort
				(
					short	*VolumeSource,		/* data to process */
					short	*VolumeDestination,	/* result */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					short	*Kernel,			/* structural element */
					long	Kx,					/* width of the kernel */
					long	Ky,					/* height of the kernel */
					long	Kz,					/* depth of the kernel */
					long	Ox,					/* kernel X origin */
					long	Oy,					/* kernel Y origin */
					long	Oz,					/* kernel Z origin */
					enum TBoundaryConvention
							Convention,			/* boundary convention */
					int		*Status				/* error management */
				)

/* perform the grey-scale morphological operation known as opening */
/* VolumeSource is the grey-scale volume to open */
/* VolumeDestination is the resulting grey-scale volume after opening */
/* both VolumeSource and VolumeDestination have size (Nx, Ny, Nz) */
/* Kernel is the grey-scale structural element by which VolumeSource is opened */
/* (Kx, Ky, Kz) is the size (bounding box) of the structural element */
/* (Ox, Oy, Oz) is the origin of the structural element */
/* Convention is the boundary convention applied to VolumeSource before opening */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin OpeningShort */

	short	*Buffer = (short *)NULL;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(OpeningShort, VolumeSource, *Status,
/**/	"Missing VolumeSource")
/**/DEBUG_CHECK_NULL_POINTER(OpeningShort, VolumeDestination, *Status,
/**/	"Missing VolumeDestination")
/**/DEBUG_CHECK_NULL_POINTER(OpeningShort, Kernel, *Status,
/**/	"Missing Kernel")
/**/DEBUG_CHECK_RANGE_LONG(OpeningShort, Nx, 1L, LONG_MAX, *Status,
/**/	"Invalid data width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(OpeningShort, Ny, 1L, LONG_MAX, *Status,
/**/	"Invalid data height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(OpeningShort, Nz, 1L, LONG_MAX, *Status,
/**/	"Invalid data depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(OpeningShort, Kx, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(OpeningShort, Ky, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(OpeningShort, Kz, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(OpeningShort, *Status)
/**/DEBUG_WRITE_ENTERING(OpeningShort,
/**/	"About to execute OpeningShort")

	switch (Convention) {
		case FiniteDataSupport:
		case MirrorOffBounds:
		case MirrorOnBounds:
		case Periodic:
			break;
		case AntiMirrorOnBounds:
		case FiniteCoefficientSupport:
		default:
			*Status = ERROR;
			WRITE_ERROR(OpeningShort, "Invalid boundary convention")
/**/		DEBUG_WRITE_LEAVING(OpeningShort, "Done")
			return(*Status);
	}
	AllocateVolumeShort(&Buffer, Nx, Ny, Nz, Status);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(OpeningShort, "Done")
		return(*Status);
	}
	*Status = ErodeShort(VolumeSource, Buffer, Nx, Ny, Nz, Kernel, Kx, Ky, Kz,
		Ox, Oy, Oz, Convention);
	if (*Status == ERROR) {
		FreeVolumeShort(&Buffer);
/**/	DEBUG_WRITE_LEAVING(OpeningShort, "Done")
		return(*Status);
	}
	*Status = DilateShort(Buffer, VolumeDestination, Nx, Ny, Nz, Kernel, Kx, Ky, Kz,
		Ox, Oy, Oz, Convention);
	if (*Status == ERROR) {
		FreeVolumeShort(&Buffer);
/**/	DEBUG_WRITE_LEAVING(OpeningShort, "Done")
		return(*Status);
	}
	*Status = FreeVolumeShort(&Buffer);
/**/DEBUG_WRITE_LEAVING(OpeningShort, "Done")
	return(*Status);
} /* end OpeningShort */

/*--------------------------------------------------------------------------*/
extern int		OpeningShortCuboid
				(
					short	*Volume,			/* data to process in-place */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					long	Kx,					/* width of the kernel */
					long	Ky,					/* height of the kernel */
					long	Kz,					/* depth of the kernel */
					enum TBoundaryConvention
							Convention,			/* boundary convention */
					int		*Status				/* error management */
				)

/* perform a morphological grey-scale opening filter using a cuboid as structural element */
/* Volume is the grey-scale volume to open in-place */
/* Volume has size (Nx, Ny, Nz) */
/* (Kx, Ky, Kz) is the size of the structural element that is filled with the value 1.0F */
/* the origin of the structural element is ((Kx - 1L) / 2L, (Ky - 1L) / 2L, (Kz - 1L) / 2L) */
/* Convention is the boundary convention applied to Volume before opening */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin OpeningShortCuboid */

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(OpeningShortCuboid, Volume, *Status,
/**/	"Missing Volume")
/**/DEBUG_CHECK_RANGE_LONG(OpeningShortCuboid, Nx, 1L, LONG_MAX, *Status,
/**/	"Invalid data width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(OpeningShortCuboid, Ny, 1L, LONG_MAX, *Status,
/**/	"Invalid data height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(OpeningShortCuboid, Nz, 1L, LONG_MAX, *Status,
/**/	"Invalid data depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(OpeningShortCuboid, Kx, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(OpeningShortCuboid, Ky, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(OpeningShortCuboid, Kz, 1L, LONG_MAX, *Status,
/**/	"Invalid kernel depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(OpeningShortCuboid, *Status)
/**/DEBUG_WRITE_ENTERING(OpeningShortCuboid,
/**/	"About to execute OpeningShortCuboid")

	switch (Convention) {
		case FiniteDataSupport:
		case MirrorOffBounds:
		case MirrorOnBounds:
		case Periodic:
			break;
		case AntiMirrorOnBounds:
		case FiniteCoefficientSupport:
		default:
			*Status = ERROR;
			WRITE_ERROR(OpeningShortCuboid, "Invalid boundary convention")
/**/		DEBUG_WRITE_LEAVING(OpeningShortCuboid, "Done")
			return(*Status);
	}
	MinFilterShortCuboid(Volume, Nx, Ny, Nz, Kx, Ky, Kz, Convention, Status);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(OpeningShortCuboid, "Done")
		return(*Status);
	}
	MaxFilterShortCuboid(Volume, Nx, Ny, Nz, Kx, Ky, Kz, Convention, Status);
/**/DEBUG_WRITE_LEAVING(OpeningShortCuboid, "Done")
	return(*Status);
} /* end OpeningShortCuboid */

