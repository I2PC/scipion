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

/*****************************************************************************
 *	Toolbox includes
 ****************************************************************************/
#include	"messagedisplay.h"
#include	"flip.h"

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
extern int		FlipXvolumeFloat
				(
					float	*VolumeSource,		/* float input data */
					float	*VolumeDestination,	/* float output data */
					long	Nx,					/* width */
					long	Ny,					/* height */
					long	Nz					/* depth */
				)

/* flip the x-axis of a volume */
/* input VolumeSource is a (float)volume of size (Nx x Ny x Nz) */
/* output VolumeDestination is a (float)volume of size (Nx x Ny x Nz) */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin FlipXvolumeFloat */

	float	*p, *q;
	long	i, j, k;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(FlipXvolumeFloat, VolumeSource, Status,
/**/	"No VolumeSource")
/**/DEBUG_CHECK_NULL_POINTER(FlipXvolumeFloat, VolumeDestination, Status,
/**/	"No VolumeDestination")
/**/DEBUG_CHECK_RANGE_LONG(FlipXvolumeFloat, Nx, 1L, LONG_MAX, Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(FlipXvolumeFloat, Ny, 1L, LONG_MAX, Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(FlipXvolumeFloat, Nz, 1L, LONG_MAX, Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(FlipXvolumeFloat, Status)
/**/DEBUG_WRITE_ENTERING(FlipXvolumeFloat,
/**/	"About to flip the x-axis of a float volume")

	p = VolumeSource;
	q = VolumeDestination + (ptrdiff_t)Nx;
	for (k = 0L; (k < Nz); k++) {
		for (j = 0L; (j < Ny); j++) {
			for (i = 0L; (i < Nx); i++) {
				*--q = *p++;
			}
			q += (ptrdiff_t)(2L * Nx);
		}
	}
/**/DEBUG_WRITE_LEAVING(FlipXvolumeFloat, "Done")
	return(Status);
} /* end FlipXvolumeFloat */

/*--------------------------------------------------------------------------*/
extern int		FlipYvolumeFloat
				(
					float	*VolumeSource,		/* float input data */
					float	*VolumeDestination,	/* float output data */
					long	Nx,					/* width */
					long	Ny,					/* height */
					long	Nz					/* depth */
				)

/* flip the y-axis of a volume */
/* input VolumeSource is a (float)volume of size (Nx x Ny x Nz) */
/* output VolumeDestination is a (float)volume of size (Nx x Ny x Nz) */

{ /* begin FlipYvolumeFloat */

	float	*p, *q;
	long	i, j, k;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(FlipYvolumeFloat, VolumeSource, Status,
/**/	"No VolumeSource")
/**/DEBUG_CHECK_NULL_POINTER(FlipYvolumeFloat, VolumeDestination, Status,
/**/	"No VolumeDestination")
/**/DEBUG_CHECK_RANGE_LONG(FlipYvolumeFloat, Nx, 1L, LONG_MAX, Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(FlipYvolumeFloat, Ny, 1L, LONG_MAX, Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(FlipYvolumeFloat, Nz, 1L, LONG_MAX, Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(FlipYvolumeFloat, Status)
/**/DEBUG_WRITE_ENTERING(FlipYvolumeFloat,
/**/	"About to flip the y-axis of a float volume")

	VolumeDestination += (ptrdiff_t)(Nx * (Ny - 1L));
	for (k = 0L; (k < Nz); k++) {
		for (i = 0L; (i < Nx); i++) {
			p = VolumeSource++;
			q = VolumeDestination++;
			for (j = 0L; (j < Ny); i++) {
				*q = *p;
				p += (ptrdiff_t)Nx;
				q -= (ptrdiff_t)Nx;
			}
		}
		VolumeSource += (ptrdiff_t)(Nx * (Ny - 1L));
		VolumeDestination += (ptrdiff_t)(Nx * (Ny - 1L));
	}
/**/DEBUG_WRITE_LEAVING(FlipYvolumeFloat, "Done")
	return(Status);
} /* end FlipYvolumeFloat */

/*--------------------------------------------------------------------------*/
extern int		FlipZvolumeFloat
				(
					float	*VolumeSource,		/* float input data */
					float	*VolumeDestination,	/* float output data */
					long	Nx,					/* width */
					long	Ny,					/* height */
					long	Nz					/* depth */
				)

/* flip the z-axis of a volume */
/* input VolumeSource is a (float)volume of size (Nx x Ny x Nz) */
/* output VolumeDestination is a (float)volume of size (Nx x Ny x Nz) */

{ /* begin FlipZvolumeFloat */

	float	*p, *q;
	long	i, j, k;
	long	Nxy = Nx * Ny;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(FlipZvolumeFloat, VolumeSource, Status,
/**/	"No VolumeSource")
/**/DEBUG_CHECK_NULL_POINTER(FlipZvolumeFloat, VolumeDestination, Status,
/**/	"No VolumeDestination")
/**/DEBUG_CHECK_RANGE_LONG(FlipZvolumeFloat, Nx, 1L, LONG_MAX, Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(FlipZvolumeFloat, Ny, 1L, LONG_MAX, Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(FlipZvolumeFloat, Nz, 1L, LONG_MAX, Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(FlipZvolumeFloat, Status)
/**/DEBUG_WRITE_ENTERING(FlipZvolumeFloat,
/**/	"About to flip the z-axis of a float volume")

	VolumeDestination += (ptrdiff_t)(Nxy * (Nz - 1L));
	for (j = 0L; (j < Ny); i++) {
		for (i = 0L; (i < Nx); i++) {
			p = VolumeSource++;
			q = VolumeDestination++;
			for (k = 0L; (k < Nz); k++) {
				*q = *p;
				p += (ptrdiff_t)Nxy;
				q -= (ptrdiff_t)Nxy;
			}
		}
	}
/**/DEBUG_WRITE_LEAVING(FlipZvolumeFloat, "Done")
	return(Status);
} /* end FlipZvolumeFloat */

