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
#include	"swap.h"

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
extern int		SwapXyVolumeFloat
				(
					float	*VolumeSource,		/* float input data */
					float	*VolumeDestination,	/* float output data */
					long	Nxy,				/* input width, output height */
					long	Nyx,				/* input height, output width */
					long	Nz					/* depth */
				)

/* swap the x-axis and the y-axis of a volume */
/* input VolumeSource is a (float)volume of size (Nxy x Nyx x Nz) */
/* output VolumeDestination is a (float)volume of size (Nyx x Nxy x Nz) */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin SwapXyVolumeFloat */

	const long
			N = Nxy * Nyx;
	long	k, i, j;
	long	n;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(SwapXyVolumeFloat, VolumeSource, Status,
/**/	"No VolumeSource")
/**/DEBUG_CHECK_NULL_POINTER(SwapXyVolumeFloat, VolumeDestination, Status,
/**/	"No VolumeDestination")
/**/DEBUG_CHECK_RANGE_LONG(SwapXyVolumeFloat, Nxy, 1L, LONG_MAX, Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(SwapXyVolumeFloat, Nyx, 1L, LONG_MAX, Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(SwapXyVolumeFloat, Nz, 1L, LONG_MAX, Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(SwapXyVolumeFloat, Status)
/**/DEBUG_WRITE_ENTERING(SwapXyVolumeFloat,
/**/	"About to swap x<->y axis of a float volume")

	for (k = 0L; (k < Nz); k++) {
		for (j = 0L; (j < Nyx); j++) {
			n = 0L;
			for (i = 0L; (i < Nxy); i++) {
				VolumeDestination[n] = *VolumeSource++;
				n += Nyx;
			}
			VolumeDestination++;
		}
		VolumeDestination += (ptrdiff_t)N;
	}
/**/DEBUG_WRITE_LEAVING(SwapXyVolumeFloat, "Done")
	return(Status);
} /* end SwapXyVolumeFloat */

/*--------------------------------------------------------------------------*/
extern int		SwapYzVolumeFloat
				(
					float	*VolumeSource,		/* float input data */
					float	*VolumeDestination,	/* float output data */
					long	Nx,					/* width */
					long	Nyz,				/* input height, output depth */
					long	Nzy					/* input depth, output height */
				)

/* swap the y-axis and the z-axis of a volume */
/* input VolumeSource is a (float)volume of size (Nx x Nyz x Nzy) */
/* output VolumeDestination is a (float)volume of size (Nx x Nzy x Nyz) */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin SwapYzVolumeFloat */

	float	*p;
	const long
			N = Nx * Nzy;
	long	k, j;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(SwapYzVolumeFloat, VolumeSource, Status,
/**/	"No VolumeSource")
/**/DEBUG_CHECK_NULL_POINTER(SwapYzVolumeFloat, VolumeDestination, Status,
/**/	"No VolumeDestination")
/**/DEBUG_CHECK_RANGE_LONG(SwapYzVolumeFloat, Nx, 1L, LONG_MAX, Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(SwapYzVolumeFloat, Nyz, 1L, LONG_MAX, Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(SwapYzVolumeFloat, Nzy, 1L, LONG_MAX, Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(SwapYzVolumeFloat, Status)
/**/DEBUG_WRITE_ENTERING(SwapYzVolumeFloat,
/**/	"About to swap y<->z axis of a float volume")

	for (k = 0L; (k < Nzy); k++) {
		p = VolumeDestination;
		for (j = 0L; (j < Nyz); j++) {
			p = (float *)memcpy(p, VolumeSource, (size_t)(Nx * (long)sizeof(float)));
			VolumeSource += (ptrdiff_t)Nx;
			p += (ptrdiff_t)N;
		}
		VolumeDestination += (ptrdiff_t)Nx;
	}
/**/DEBUG_WRITE_LEAVING(SwapYzVolumeFloat, "Done")
	return(Status);
} /* end SwapYzVolumeFloat */

/*--------------------------------------------------------------------------*/
extern int		SwapZxVolumeFloat
				(
					float	*VolumeSource,		/* float input data */
					float	*VolumeDestination,	/* float output data */
					long	Nxz,				/* input width, output depth */
					long	Ny,					/* height */
					long	Nzx					/* input depth, output width */
				)

/* swap the z-axis and the x-axis of a volume */
/* input VolumeSource is a (float)volume of size (Nxz x Ny x Nzx) */
/* output VolumeDestination is a (float)volume of size (Nzx x Ny x Nxz) */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin SwapZxVolumeFloat */

	float	*p, *q;
	const long
			N = Nzx * Ny;
	long	k, j, i;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(SwapZxVolumeFloat, VolumeSource, Status,
/**/	"No VolumeSource")
/**/DEBUG_CHECK_NULL_POINTER(SwapZxVolumeFloat, VolumeDestination, Status,
/**/	"No VolumeDestination")
/**/DEBUG_CHECK_RANGE_LONG(SwapZxVolumeFloat, Nxz, 1L, LONG_MAX, Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(SwapZxVolumeFloat, Ny, 1L, LONG_MAX, Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(SwapZxVolumeFloat, Nzx, 1L, LONG_MAX, Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(SwapZxVolumeFloat, Status)
/**/DEBUG_WRITE_ENTERING(SwapZxVolumeFloat,
/**/	"About to swap z<->x axis of a float volume")

	for (k = 0L; (k < Nzx); k++) {
		q = VolumeDestination++;
		for (j = 0L; (j < Ny); j++) {
			p = q;
			for (i = 0L; (i < Nxz); i++) {
				*p = *VolumeSource++;
				p += (ptrdiff_t)N;
			}
			q += (ptrdiff_t)Nzx;
		}
	}
/**/DEBUG_WRITE_LEAVING(SwapZxVolumeFloat, "Done")
	return(Status);
} /* end SwapZxVolumeFloat */

