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
#include	"../headers/convert.h"
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
extern int		AllocateLineDouble
				(
					double	*(Line[]),			/* double output pointer */
					long	LineLength,			/* length of the line */
					int		*Status				/* error management */
				)

/* Allocates a 1D array of double */
/* success: return(!ERROR); failure: return(ERROR) and set *Line to NULL */
/* the returned value is duplicated in Status */

{ /* begin AllocateLineDouble */

	*Status = !ERROR;

/**/DEBUG_CHECK_RANGE_LONG(AllocateLineDouble, LineLength, 1L, LONG_MAX, *Status,
/**/	"Invalid length (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(AllocateLineDouble, *Status)
/**/DEBUG_WRITE_ENTERING(AllocateLineDouble,
/**/	"About to allocate a 1D array of double")
#ifdef DEBUG
/**/if (*Line != (double *)NULL) {
/**/	WRITE_WARNING(AllocateLineDouble, "Line may have been previously allocated")
/**/}
#endif

	*Line = (double *)malloc((size_t)(LineLength * (long)sizeof(double)));
	if (*Line == (double *)NULL) {
		*Status = ERROR;
		WRITE_ERROR(AllocateLineDouble, "Unable to perform allocation")
	}
/**/DEBUG_WRITE_LEAVING(AllocateLineDouble, "Done")
	return(*Status);
} /* end AllocateLineDouble */

/*--------------------------------------------------------------------------*/
extern int		AllocateVolumeFloat
				(
					float	**Volume,			/* float output pointer */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					int		*Status				/* error management */
				)

/* Allocates a 3D array of float */
/* success: return(!ERROR); failure: return(ERROR) and set *Volume to NULL */

{ /* begin AllocateVolumeFloat */

	*Status = !ERROR;

/**/DEBUG_CHECK_RANGE_LONG(AllocateVolumeFloat, Nx, 1L, LONG_MAX, *Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(AllocateVolumeFloat, Ny, 1L, LONG_MAX, *Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(AllocateVolumeFloat, Nz, 1L, LONG_MAX, *Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(AllocateVolumeFloat, *Status)
/**/DEBUG_WRITE_ENTERING(AllocateVolumeFloat,
/**/	"About to allocate a 3D array of float")
#ifdef DEBUG
/**/if (*Volume != (float *)NULL) {
/**/	WRITE_WARNING(AllocateVolumeFloat, "Volume may have been previously allocated")
/**/}
#endif

	*Volume = (float *)malloc((size_t)(Nx * Ny * Nz * (long)sizeof(float)));
	if (*Volume == (float *)NULL) {
		*Status = ERROR;
		WRITE_ERROR(AllocateVolumeFloat, "Unable to perform allocation")
	}
/**/DEBUG_WRITE_LEAVING(AllocateVolumeFloat, "Done")
	return(*Status);
} /* end AllocateVolumeFloat */

/*--------------------------------------------------------------------------*/
extern int		AllocateVolumeShort
				(
					short	**Volume,			/* short output pointer */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					int		*Status				/* error management */
				)

/* Allocates a 3D array of short */
/* success: return(!ERROR); failure: return(ERROR) and set Volume to NULL */

{ /* begin AllocateVolumeShort */

	*Status = !ERROR;

/**/DEBUG_CHECK_RANGE_LONG(AllocateVolumeShort, Nx, 1L, LONG_MAX, *Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(AllocateVolumeShort, Ny, 1L, LONG_MAX, *Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(AllocateVolumeShort, Nz, 1L, LONG_MAX, *Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(AllocateVolumeShort, *Status)
/**/DEBUG_WRITE_ENTERING(AllocateVolumeShort,
/**/	"About to allocate a 3D array of short")
#ifdef DEBUG
/**/if (*Volume != (short *)NULL) {
/**/	WRITE_WARNING(AllocateVolumeShort, "Volume may have been previously allocated")
/**/}
#endif

	*Volume = (short *)malloc((size_t)(Nx * Ny * Nz * (long)sizeof(short)));
	if (*Volume == (short *)NULL) {
		*Status = ERROR;
		WRITE_ERROR(AllocateVolumeShort, "Unable to perform allocation")
	}
/**/DEBUG_WRITE_LEAVING(AllocateVolumeShort, "Done")
	return(*Status);
} /* end AllocateVolumeShort */

/*--------------------------------------------------------------------------*/
extern int		CopyFloatToFloat
				(
					float	*VolumeSource,		/* float input data */
					long	NxSource,			/* width of the input */
					long	NySource,			/* height of the input */
					long	NzSource,			/* depth of the input */
					long	XSource,			/* x coordinate to get from */
					long	YSource,			/* y coordinate to get from */
					long	ZSource,			/* z coordinate to get from */
					float	*VolumeDestination,	/* float output data */
					long	NxDestination,		/* width of the output */
					long	NyDestination,		/* height of the output */
					long	NzDestination,		/* depth of the output */
					long	XDestination,		/* x coordinate to put into */
					long	YDestination,		/* y coordinate to put into */
					long	ZDestination,		/* z coordinate to put into */
					long	NxCopy,				/* width of the block to copy */
					long	NyCopy,				/* height of the block to copy */
					long	NzCopy				/* depth of the block to copy */
				)

/* copies a sub-volume from a float input volume into a portion of a float output volume */
/* input and output may not share their memory space */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin CopyFloatToFloat */

	float	*P, *Q;
	const ptrdiff_t
			NxNySource = (ptrdiff_t)(NxSource * NySource),
			NxNyDestination = (ptrdiff_t)(NxDestination * NyDestination);
	long	Y, Z;
	int		Status = !ERROR;

#if defined(CODEWARRIOR) && !defined(DEBUG)
#pragma unused(NzSource, NzDestination)
#endif
/**/DEBUG_CHECK_NULL_POINTER(CopyFloatToFloat, VolumeSource, Status,
/**/	"No input")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToFloat, NxSource, 1L, LONG_MAX, Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToFloat, NySource, 1L, LONG_MAX, Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToFloat, NzSource, 1L, LONG_MAX, Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToFloat, XSource, 0L, NxSource - 1L, Status,
/**/	"Invalid X coordinate")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToFloat, YSource, 0L, NySource - 1L, Status,
/**/	"Invalid Y coordinate")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToFloat, ZSource, 0L, NzSource - 1L, Status,
/**/	"Invalid Z coordinate")
/**/DEBUG_CHECK_NULL_POINTER(CopyFloatToFloat, VolumeDestination, Status,
/**/	"No output")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToFloat, NxDestination, 1L, LONG_MAX, Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToFloat, NyDestination, 1L, LONG_MAX, Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToFloat, NzDestination, 1L, LONG_MAX, Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToFloat, XDestination, 0L, NxDestination - 1L, Status,
/**/	"Invalid X coordinate")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToFloat, YDestination, 0L, NyDestination - 1L, Status,
/**/	"Invalid Y coordinate")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToFloat, ZDestination, 0L, NzDestination - 1L, Status,
/**/	"Invalid Z coordinate")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToFloat, NxCopy, 1L, NxSource - XSource, Status,
/**/	"Invalid width (negative or excessive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToFloat, NyCopy, 1L, NySource - YSource, Status,
/**/	"Invalid height (negative or excessive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToFloat, NzCopy, 1L, NzSource - ZSource, Status,
/**/	"Invalid depth (negative or excessive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToFloat, NxCopy, 1L, NxDestination - XDestination, Status,
/**/	"Invalid width (negative or excessive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToFloat, NyCopy, 1L, NyDestination - YDestination, Status,
/**/	"Invalid height (negative or excessive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToFloat, NzCopy, 1L, NzDestination - ZDestination, Status,
/**/	"Invalid depth (negative or excessive)")
/**/DEBUG_RETURN_ON_ERROR(CopyFloatToFloat, Status)
/**/DEBUG_WRITE_ENTERING(CopyFloatToFloat,
/**/	"About to copy a float subvolume into a float subvolume")
#ifdef DEBUG
/**/if (((VolumeSource + (ptrdiff_t)(NxSource * NySource * NzSource))
/**/	> (VolumeDestination + (ptrdiff_t)(XDestination + NxDestination
/**/	* (YDestination + NyDestination * ZDestination))))
/**/	&& ((VolumeDestination + (ptrdiff_t)(NxDestination * NyDestination
/**/	* NzDestination)) > (VolumeSource + (ptrdiff_t)(XSource + NxSource
/**/	* (YSource + NySource * ZSource))))) {
/**/	WRITE_WARNING(CopyFloatToFloat, "Data overlap: memcpy may fail")
/**/}
#endif

	VolumeSource += (ptrdiff_t)(XSource + NxSource * (YSource + NySource * ZSource));
	VolumeDestination += (ptrdiff_t)(XDestination + NxDestination * (YDestination
		+ NyDestination * ZDestination));
	for (Z = -NzCopy; (Z < 0L); Z++) {
		P = VolumeSource;
		Q = VolumeDestination;
		for (Y = -NyCopy; (Y < 0L); Y++) {
			Q = (float *)memcpy(Q, P, (size_t)(NxCopy * (long)sizeof(float)));
			P += (ptrdiff_t)NxSource;
			Q += (ptrdiff_t)NxDestination;
		}
		VolumeSource += NxNySource;
		VolumeDestination += NxNyDestination;
	}
/**/DEBUG_WRITE_LEAVING(CopyFloatToFloat, "Done")
	return(Status);
} /* end CopyFloatToFloat */

/*--------------------------------------------------------------------------*/
extern int		CopyFloatToShort
				(
					float	*VolumeSource,		/* float input data */
					long	NxSource,			/* width of the input */
					long	NySource,			/* height of the input */
					long	NzSource,			/* depth of the input */
					long	XSource,			/* x coordinate to get from */
					long	YSource,			/* y coordinate to get from */
					long	ZSource,			/* z coordinate to get from */
					short	*VolumeDestination,	/* short output data */
					long	NxDestination,		/* width of the output */
					long	NyDestination,		/* height of the output */
					long	NzDestination,		/* depth of the output */
					long	XDestination,		/* x coordinate to put into */
					long	YDestination,		/* y coordinate to put into */
					long	ZDestination,		/* z coordinate to put into */
					long	NxCopy,				/* width of the block to copy */
					long	NyCopy,				/* height of the block to copy */
					long	NzCopy				/* depth of the block to copy */
				)

/* copies a sub-volume from a float input volume into a portion of a short output volume */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin CopyFloatToShort */

	float	*P;
	short	*Q;
	const ptrdiff_t
			IgnoredDestination = (ptrdiff_t)(NxDestination - NxCopy),
			IgnoredSource = (ptrdiff_t)(NxSource - NxCopy),
			NxNyDestination = (ptrdiff_t)(NxDestination * NyDestination),
			NxNySource = (ptrdiff_t)(NxSource * NySource);
	long	X, Y, Z;
	int		Status = !ERROR;

#if defined(CODEWARRIOR) && !defined(DEBUG)
#pragma unused(NzSource, NzDestination)
#endif
/**/DEBUG_CHECK_NULL_POINTER(CopyFloatToShort, VolumeSource, Status,
/**/	"No input")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToShort, NxSource, 1L, LONG_MAX, Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToShort, NySource, 1L, LONG_MAX, Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToShort, NzSource, 1L, LONG_MAX, Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToShort, XSource, 0L, NxSource - 1L, Status,
/**/	"Invalid X coordinate")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToShort, YSource, 0L, NySource - 1L, Status,
/**/	"Invalid Y coordinate")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToShort, ZSource, 0L, NzSource - 1L, Status,
/**/	"Invalid Z coordinate")
/**/DEBUG_CHECK_NULL_POINTER(CopyFloatToShort, VolumeDestination, Status,
/**/	"No output")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToShort, NxDestination, 1L, LONG_MAX, Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToShort, NyDestination, 1L, LONG_MAX, Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToShort, NzDestination, 1L, LONG_MAX, Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToShort, XDestination, 0L, NxDestination - 1L, Status,
/**/	"Invalid X coordinate")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToShort, YDestination, 0L, NyDestination - 1L, Status,
/**/	"Invalid Y coordinate")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToShort, ZDestination, 0L, NzDestination - 1L, Status,
/**/	"Invalid Z coordinate")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToShort, NxCopy, 1L, NxSource - XSource, Status,
/**/	"Invalid width (negative or excessive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToShort, NyCopy, 1L, NySource - YSource, Status,
/**/	"Invalid height (negative or excessive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToShort, NzCopy, 1L, NzSource - ZSource, Status,
/**/	"Invalid depth (negative or excessive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToShort, NxCopy, 1L, NxDestination - XDestination, Status,
/**/	"Invalid width (negative or excessive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToShort, NyCopy, 1L, NyDestination - YDestination, Status,
/**/	"Invalid height (negative or excessive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToShort, NzCopy, 1L, NzDestination - ZDestination, Status,
/**/	"Invalid depth (negative or excessive)")
/**/DEBUG_RETURN_ON_ERROR(CopyFloatToShort, Status)
/**/DEBUG_WRITE_ENTERING(CopyFloatToShort,
/**/	"About to copy a float subvolume into a short subvolume")

	VolumeSource += (ptrdiff_t)(XSource + NxSource * (YSource + NySource * ZSource));
	VolumeDestination += (ptrdiff_t)(XDestination + NxDestination * (YDestination
		+ NyDestination * ZDestination));
	for (Z = -NzCopy; (Z < 0L); Z++) {
		P = VolumeSource;
		Q = VolumeDestination;
		for (Y = -NyCopy; (Y < 0L); Y++) {
			for (X = -NxCopy; (X < 0L); X++) {
				*Q++ = ConvertFloatToShort(*P++);
			}
			P += IgnoredSource;
			Q += IgnoredDestination;
		}
		VolumeSource += NxNySource;
		VolumeDestination += NxNyDestination;
	}
/**/DEBUG_WRITE_LEAVING(CopyFloatToShort, "Done")
	return(Status);
} /* end CopyFloatToShort */

/*--------------------------------------------------------------------------*/
extern int		CopyShortToFloat
				(
					short	*VolumeSource,		/* short input data */
					long	NxSource,			/* width of the input */
					long	NySource,			/* height of the input */
					long	NzSource,			/* depth of the input */
					long	XSource,			/* x coordinate to get from */
					long	YSource,			/* y coordinate to get from */
					long	ZSource,			/* z coordinate to get from */
					float	*VolumeDestination,	/* float output data */
					long	NxDestination,		/* width of the output */
					long	NyDestination,		/* height of the output */
					long	NzDestination,		/* depth of the output */
					long	XDestination,		/* x coordinate to put into */
					long	YDestination,		/* y coordinate to put into */
					long	ZDestination,		/* z coordinate to put into */
					long	NxCopy,				/* width of the block to copy */
					long	NyCopy,				/* height of the block to copy */
					long	NzCopy				/* depth of the block to copy */
				)

/* copies a sub-volume from a short input volume into a portion of a float output volume */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin CopyShortToFloat */

	float	*Q;
	short	*P;
	const ptrdiff_t
			IgnoredDestination = (ptrdiff_t)(NxDestination - NxCopy),
			IgnoredSource = (ptrdiff_t)(NxSource - NxCopy),
			NxNyDestination = (ptrdiff_t)(NxDestination * NyDestination),
			NxNySource = (ptrdiff_t)(NxSource * NySource);
	long	X, Y, Z;
	int		Status = !ERROR;

#if defined(CODEWARRIOR) && !defined(DEBUG)
#pragma unused(NzSource, NzDestination)
#endif
/**/DEBUG_CHECK_NULL_POINTER(CopyShortToFloat, VolumeSource, Status,
/**/	"No input")
/**/DEBUG_CHECK_RANGE_LONG(CopyShortToFloat, NxSource, 1L, LONG_MAX, Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyShortToFloat, NySource, 1L, LONG_MAX, Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyShortToFloat, NzSource, 1L, LONG_MAX, Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyShortToFloat, XSource, 0L, NxSource - 1L, Status,
/**/	"Invalid X coordinate")
/**/DEBUG_CHECK_RANGE_LONG(CopyShortToFloat, YSource, 0L, NySource - 1L, Status,
/**/	"Invalid Y coordinate")
/**/DEBUG_CHECK_RANGE_LONG(CopyShortToFloat, ZSource, 0L, NzSource - 1L, Status,
/**/	"Invalid Z coordinate")
/**/DEBUG_CHECK_NULL_POINTER(CopyShortToFloat, VolumeDestination, Status,
/**/	"No output")
/**/DEBUG_CHECK_RANGE_LONG(CopyShortToFloat, NxDestination, 1L, LONG_MAX, Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyShortToFloat, NyDestination, 1L, LONG_MAX, Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyShortToFloat, NzDestination, 1L, LONG_MAX, Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyShortToFloat, XDestination, 0L, NxDestination - 1L, Status,
/**/	"Invalid X coordinate")
/**/DEBUG_CHECK_RANGE_LONG(CopyShortToFloat, YDestination, 0L, NyDestination - 1L, Status,
/**/	"Invalid Y coordinate")
/**/DEBUG_CHECK_RANGE_LONG(CopyShortToFloat, ZDestination, 0L, NzDestination - 1L, Status,
/**/	"Invalid Z coordinate")
/**/DEBUG_CHECK_RANGE_LONG(CopyShortToFloat, NxCopy, 1L, NxSource - XSource, Status,
/**/	"Invalid width (negative or excessive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyShortToFloat, NyCopy, 1L, NySource - YSource, Status,
/**/	"Invalid height (negative or excessive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyShortToFloat, NzCopy, 1L, NzSource - ZSource, Status,
/**/	"Invalid depth (negative or excessive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyShortToFloat, NxCopy, 1L, NxDestination - XDestination, Status,
/**/	"Invalid width (negative or excessive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyShortToFloat, NyCopy, 1L, NyDestination - YDestination, Status,
/**/	"Invalid height (negative or excessive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyShortToFloat, NzCopy, 1L, NzDestination - ZDestination, Status,
/**/	"Invalid depth (negative or excessive)")
/**/DEBUG_RETURN_ON_ERROR(CopyShortToFloat, Status)
/**/DEBUG_WRITE_ENTERING(CopyShortToFloat,
/**/	"About to copy a short subvolume into a float subvolume")

	VolumeSource += (ptrdiff_t)(XSource + NxSource * (YSource + NySource * ZSource));
	VolumeDestination += (ptrdiff_t)(XDestination + NxDestination * (YDestination
		+ NyDestination * ZDestination));
	for (Z = -NzCopy; (Z < 0L); Z++) {
		P = VolumeSource;
		Q = VolumeDestination;
		for (Y = -NyCopy; (Y < 0L); Y++) {
			for (X = -NxCopy; (X < 0L); X++) {
				*Q++ = (float)*P++;
			}
			P += IgnoredSource;
			Q += IgnoredDestination;
		}
		VolumeSource += NxNySource;
		VolumeDestination += NxNyDestination;
	}
/**/DEBUG_WRITE_LEAVING(CopyShortToFloat, "Done")
	return(Status);
} /* end CopyShortToFloat */

/*--------------------------------------------------------------------------*/
extern int		CopyShortToShort
				(
					short	*VolumeSource,		/* short input data */
					long	NxSource,			/* width of the input */
					long	NySource,			/* height of the input */
					long	NzSource,			/* depth of the input */
					long	XSource,			/* x coordinate to get from */
					long	YSource,			/* y coordinate to get from */
					long	ZSource,			/* z coordinate to get from */
					short	*VolumeDestination,	/* short output data */
					long	NxDestination,		/* width of the output */
					long	NyDestination,		/* height of the output */
					long	NzDestination,		/* depth of the output */
					long	XDestination,		/* x coordinate to put into */
					long	YDestination,		/* y coordinate to put into */
					long	ZDestination,		/* z coordinate to put into */
					long	NxCopy,				/* width of the block to copy */
					long	NyCopy,				/* height of the block to copy */
					long	NzCopy				/* depth of the block to copy */
				)

/* copies a sub-volume from a short input volume into a portion of a short output volume */
/* input and output may not share their memory space */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin CopyShortToShort */

	short	*P, *Q;
	const ptrdiff_t
			NxNySource = (ptrdiff_t)(NxSource * NySource),
			NxNyDestination = (ptrdiff_t)(NxDestination * NyDestination);
	long	Y, Z;
	int		Status = !ERROR;

#if defined(CODEWARRIOR) && !defined(DEBUG)
#pragma unused(NzSource, NzDestination)
#endif
/**/DEBUG_CHECK_NULL_POINTER(CopyShortToShort, VolumeSource, Status,
/**/	"No input")
/**/DEBUG_CHECK_RANGE_LONG(CopyShortToShort, NxSource, 1L, LONG_MAX, Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyShortToShort, NySource, 1L, LONG_MAX, Status,
/**/	"Invalid heigth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyShortToShort, NzSource, 1L, LONG_MAX, Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyShortToShort, XSource, 0L, NxSource - 1L, Status,
/**/	"Invalid X coordinate")
/**/DEBUG_CHECK_RANGE_LONG(CopyShortToShort, YSource, 0L, NySource - 1L, Status,
/**/	"Invalid Y coordinate")
/**/DEBUG_CHECK_RANGE_LONG(CopyShortToShort, ZSource, 0L, NzSource - 1L, Status,
/**/	"Invalid Z coordinate")
/**/DEBUG_CHECK_NULL_POINTER(CopyShortToShort, VolumeDestination, Status,
/**/	"No output")
/**/DEBUG_CHECK_RANGE_LONG(CopyShortToShort, NxDestination, 1L, LONG_MAX, Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyShortToShort, NyDestination, 1L, LONG_MAX, Status,
/**/	"Invalid heigth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyShortToShort, NzDestination, 1L, LONG_MAX, Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyShortToShort, XDestination, 0L, NxDestination - 1L, Status,
/**/	"Invalid X coordinate")
/**/DEBUG_CHECK_RANGE_LONG(CopyShortToShort, YDestination, 0L, NyDestination - 1L, Status,
/**/	"Invalid Y coordinate")
/**/DEBUG_CHECK_RANGE_LONG(CopyShortToShort, ZDestination, 0L, NzDestination - 1L, Status,
/**/	"Invalid Z coordinate")
/**/DEBUG_CHECK_RANGE_LONG(CopyShortToShort, NxCopy, 1L, NxSource - XSource, Status,
/**/	"Invalid width (negative or excessive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyShortToShort, NyCopy, 1L, NySource - YSource, Status,
/**/	"Invalid height (negative or excessive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyShortToShort, NzCopy, 1L, NzSource - ZSource, Status,
/**/	"Invalid depth (negative or excessive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyShortToShort, NxCopy, 1L, NxDestination - XDestination, Status,
/**/	"Invalid width (negative or excessive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyShortToShort, NyCopy, 1L, NyDestination - YDestination, Status,
/**/	"Invalid height (negative or excessive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyShortToShort, NzCopy, 1L, NzDestination - ZDestination, Status,
/**/	"Invalid depth (negative or excessive)")
/**/DEBUG_RETURN_ON_ERROR(CopyShortToShort, Status)
/**/DEBUG_WRITE_ENTERING(CopyShortToShort,
/**/	"About to copy a short subvolume into a short subvolume")
#ifdef DEBUG
/**/if (((VolumeSource + (ptrdiff_t)(NxSource * NySource * NzSource))
/**/	> (VolumeDestination + (ptrdiff_t)(XDestination + NxDestination
/**/	* (YDestination + NyDestination * ZDestination))))
/**/	&& ((VolumeDestination + (ptrdiff_t)(NxDestination * NyDestination
/**/	* NzDestination)) > (VolumeSource + (ptrdiff_t)(XSource + NxSource
/**/	* (YSource + NySource * ZSource))))) {
/**/	WRITE_WARNING(CopyShortToShort, "Data overlap: memcpy may fail")
/**/}
#endif

	VolumeSource += (ptrdiff_t)(XSource + NxSource * (YSource + NySource * ZSource));
	VolumeDestination += (ptrdiff_t)(XDestination + NxDestination * (YDestination
		+ NyDestination * ZDestination));
	for (Z = -NzCopy; (Z < 0L); Z++) {
		P = VolumeSource;
		Q = VolumeDestination;
		for (Y = -NyCopy; (Y < 0L); Y++) {
			Q = (short *)memcpy(Q, P, (size_t)(NxCopy * (long)sizeof(short)));
			P += (ptrdiff_t)NxSource;
			Q += (ptrdiff_t)NxDestination;
		}
		VolumeSource += NxNySource;
		VolumeDestination += NxNyDestination;
	}
/**/DEBUG_WRITE_LEAVING(CopyShortToShort, "Done")
	return(Status);
} /* end CopyShortToShort */

/*--------------------------------------------------------------------------*/
extern int		FreeLineDouble
				(
					double	*(Line[])			/* 1D double array */
				)

/* Frees a 1D array of double */
/* success: return(!ERROR) and set *Line to NULL; failure: return(ERROR); */

{ /* begin FreeLineDouble */

	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(FreeLineDouble, *Line, Status,
/**/	"Nothing to free")
/**/DEBUG_RETURN_ON_ERROR(FreeLineDouble, Status)
/**/DEBUG_WRITE_ENTERING(FreeLineDouble,
/**/	"About to free a 1D array of double")

	free(*Line);
	*Line = (double *)NULL;
/**/DEBUG_WRITE_LEAVING(FreeLineDouble, "Done")
	return(Status);
} /* end FreeLineDouble */

/*--------------------------------------------------------------------------*/
extern int		FreeVolumeFloat
				(
					float	**Volume			/* 3D float array */
				)

/* Frees a 3D array of float */
/* success: return(!ERROR) and set *Volume to NULL; failure: return(ERROR); */

{ /* begin FreeVolumeFloat */

	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(FreeVolumeFloat, *Volume, Status,
/**/	"Nothing to free")
/**/DEBUG_RETURN_ON_ERROR(FreeVolumeFloat, Status)
/**/DEBUG_WRITE_ENTERING(FreeVolumeFloat,
/**/	"About to free a 3D array of float")

	free(*Volume);
	*Volume = (float *)NULL;
/**/DEBUG_WRITE_LEAVING(FreeVolumeFloat, "Done")
	return(Status);
} /* end FreeVolumeFloat */

/*--------------------------------------------------------------------------*/
extern int		FreeVolumeShort
				(
					short	**Volume			/* 3D short array */
				)

/* Frees a 3D array of short */
/* success: return(!ERROR) and set Volume to NULL; failure: return(ERROR); */

{ /* begin FreeVolumeShort */

	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(FreeVolumeShort, *Volume, Status,
/**/	"Nothing to free")
/**/DEBUG_RETURN_ON_ERROR(FreeVolumeShort, Status)
/**/DEBUG_WRITE_ENTERING(FreeVolumeShort,
/**/	"About to free a 3D array of short")

	free(*Volume);
	*Volume = (short *)NULL;
/**/DEBUG_WRITE_LEAVING(FreeVolumeShort, "Done")
	return(Status);
} /* end FreeVolumeShort */

/*--------------------------------------------------------------------------*/
extern int		GetxFloatToDouble
				(
					float	*VolumeSource,		/* float input data */
					long	NxSource,			/* width of the input */
					long	NySource,			/* height of the input */
					long	NzSource,			/* depth of the input */
					long	XSource,			/* x coordinate to get from */
					long	YSource,			/* y coordinate to get from */
					long	ZSource,			/* z coordinate to get from */
					double	RowDestination[],	/* double output data */
					long	NxCopy				/* length of the output */
				)

/* copies a horizontal line from a float input volume into a 1D double output array */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin GetxFloatToDouble */

	double	*End;
	int		Status = !ERROR;

#if defined(CODEWARRIOR) && !defined(DEBUG)
#pragma unused(NzSource)
#endif
/**/DEBUG_CHECK_NULL_POINTER(GetxFloatToDouble, VolumeSource, Status,
/**/	"No input volume")
/**/DEBUG_CHECK_RANGE_LONG(GetxFloatToDouble, NxSource, 1L, LONG_MAX, Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(GetxFloatToDouble, NySource, 1L, LONG_MAX, Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(GetxFloatToDouble, NzSource, 1L, LONG_MAX, Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(GetxFloatToDouble, XSource, 0L, NxSource - 1L, Status,
/**/	"Invalid X coordinate")
/**/DEBUG_CHECK_RANGE_LONG(GetxFloatToDouble, YSource, 0L, NySource - 1L, Status,
/**/	"Invalid Y coordinate")
/**/DEBUG_CHECK_RANGE_LONG(GetxFloatToDouble, ZSource, 0L, NzSource - 1L, Status,
/**/	"Invalid Z coordinate")
/**/DEBUG_CHECK_NULL_POINTER(GetxFloatToDouble, RowDestination, Status,
/**/	"No output line")
/**/DEBUG_CHECK_RANGE_LONG(GetxFloatToDouble, NxCopy, 1L, NxSource - XSource, Status,
/**/	"Invalid width (negative or excessive)")
/**/DEBUG_RETURN_ON_ERROR(GetxFloatToDouble, Status)
/**/DEBUG_WRITE_ENTERING(GetxFloatToDouble,
/**/	"About to get a (double)row from a (float)volume")

	VolumeSource += (ptrdiff_t)(XSource + NxSource * (YSource + NySource * ZSource));
	End = RowDestination + (ptrdiff_t)NxCopy;
	while (RowDestination < End) {
		*RowDestination++ = (double)*VolumeSource++;
	}
/**/DEBUG_WRITE_LEAVING(GetxFloatToDouble, "Done")
	return(Status);
} /* end GetxFloatToDouble */

/*--------------------------------------------------------------------------*/
extern int		GetxShortToDouble
				(
					short	*VolumeSource,		/* short input data */
					long	NxSource,			/* width of the input */
					long	NySource,			/* height of the input */
					long	NzSource,			/* depth of the input */
					long	XSource,			/* x coordinate to get from */
					long	YSource,			/* y coordinate to get from */
					long	ZSource,			/* z coordinate to get from */
					double	RowDestination[],	/* double output data */
					long	NxCopy				/* length of the output */
				)

/* copies a horizontal line from a short input volume into a 1D double output array */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin GetxShortToDouble */

	double	*End;
	int		Status = !ERROR;

#if defined(CODEWARRIOR) && !defined(DEBUG)
#pragma unused(NzSource)
#endif
/**/DEBUG_CHECK_NULL_POINTER(GetxShortToDouble, VolumeSource, Status,
/**/	"No input volume")
/**/DEBUG_CHECK_RANGE_LONG(GetxShortToDouble, NxSource, 1L, LONG_MAX, Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(GetxShortToDouble, NySource, 1L, LONG_MAX, Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(GetxShortToDouble, NzSource, 1L, LONG_MAX, Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(GetxShortToDouble, XSource, 0L, NxSource - 1L, Status,
/**/	"Invalid X coordinate")
/**/DEBUG_CHECK_RANGE_LONG(GetxShortToDouble, YSource, 0L, NySource - 1L, Status,
/**/	"Invalid Y coordinate")
/**/DEBUG_CHECK_RANGE_LONG(GetxShortToDouble, ZSource, 0L, NzSource - 1L, Status,
/**/	"Invalid Z coordinate")
/**/DEBUG_CHECK_NULL_POINTER(GetxShortToDouble, RowDestination, Status,
/**/	"No output line")
/**/DEBUG_CHECK_RANGE_LONG(GetxShortToDouble, NxCopy, 1L, NxSource - XSource, Status,
/**/	"Invalid width (negative or excessive)")
/**/DEBUG_RETURN_ON_ERROR(GetxShortToDouble, Status)
/**/DEBUG_WRITE_ENTERING(GetxShortToDouble,
/**/	"About to get a (double)row from a (short)volume")

	VolumeSource += (ptrdiff_t)(XSource + NxSource * (YSource + NySource * ZSource));
	End = RowDestination + (ptrdiff_t)NxCopy;
	while (RowDestination < End) {
		*RowDestination++ = (double)*VolumeSource++;
	}
/**/DEBUG_WRITE_LEAVING(GetxShortToDouble, "Done")
	return(Status);
} /* end GetxShortToDouble */

/*--------------------------------------------------------------------------*/
extern int		GetyFloatToDouble
				(
					float	*VolumeSource,		/* float input data */
					long	NxSource,			/* width of the input */
					long	NySource,			/* height of the input */
					long	NzSource,			/* depth of the input */
					long	XSource,			/* x coordinate to get from */
					long	YSource,			/* y coordinate to get from */
					long	ZSource,			/* z coordinate to get from */
					double	ColumnDestination[],/* double output data */
					long	NyCopy				/* length of the output */
				)

/* copies a vertical line from a float input volume into a 1D double output array */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin GetyFloatToDouble */

	double	*End;
	int		Status = !ERROR;

#if defined(CODEWARRIOR) && !defined(DEBUG)
#pragma unused(NzSource)
#endif
/**/DEBUG_CHECK_NULL_POINTER(GetyFloatToDouble, VolumeSource, Status,
/**/	"No input volume")
/**/DEBUG_CHECK_RANGE_LONG(GetyFloatToDouble, NxSource, 1L, LONG_MAX, Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(GetyFloatToDouble, NySource, 1L, LONG_MAX, Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(GetyFloatToDouble, NzSource, 1L, LONG_MAX, Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(GetyFloatToDouble, XSource, 0L, NxSource - 1L, Status,
/**/	"Invalid X coordinate")
/**/DEBUG_CHECK_RANGE_LONG(GetyFloatToDouble, YSource, 0L, NySource - 1L, Status,
/**/	"Invalid Y coordinate")
/**/DEBUG_CHECK_RANGE_LONG(GetyFloatToDouble, ZSource, 0L, NzSource - 1L, Status,
/**/	"Invalid Z coordinate")
/**/DEBUG_CHECK_NULL_POINTER(GetyFloatToDouble, ColumnDestination, Status,
/**/	"No output line")
/**/DEBUG_CHECK_RANGE_LONG(GetyFloatToDouble, NyCopy, 1L, NySource - YSource, Status,
/**/	"Invalid height (negative or excessive)")
/**/DEBUG_RETURN_ON_ERROR(GetyFloatToDouble, Status)
/**/DEBUG_WRITE_ENTERING(GetyFloatToDouble,
/**/	"About to get a (double)column from a (float)volume")

	VolumeSource += (ptrdiff_t)(XSource + NxSource * (YSource + NySource * ZSource));
	End = ColumnDestination + (ptrdiff_t)NyCopy;
	while (ColumnDestination < End) {
		*ColumnDestination++ = (double)*VolumeSource;
		VolumeSource += NxSource;
	}
/**/DEBUG_WRITE_LEAVING(GetyFloatToDouble, "Done")
	return(Status);
} /* end GetyFloatToDouble */

/*--------------------------------------------------------------------------*/
extern int		GetyShortToDouble
				(
					short	*VolumeSource,		/* short input data */
					long	NxSource,			/* width of the input */
					long	NySource,			/* height of the input */
					long	NzSource,			/* depth of the input */
					long	XSource,			/* x coordinate to get from */
					long	YSource,			/* y coordinate to get from */
					long	ZSource,			/* z coordinate to get from */
					double	ColumnDestination[],/* double output data */
					long	NyCopy				/* length of the output */
				)

/* copies a vertical line from a short input volume into a 1D double output array */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin GetyShortToDouble */

	double	*End;
	int		Status = !ERROR;

#if defined(CODEWARRIOR) && !defined(DEBUG)
#pragma unused(NzSource)
#endif
/**/DEBUG_CHECK_NULL_POINTER(GetyShortToDouble, VolumeSource, Status,
/**/	"No input volume")
/**/DEBUG_CHECK_RANGE_LONG(GetyShortToDouble, NxSource, 1L, LONG_MAX, Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(GetyShortToDouble, NySource, 1L, LONG_MAX, Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(GetyShortToDouble, NzSource, 1L, LONG_MAX, Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(GetyShortToDouble, XSource, 0L, NxSource - 1L, Status,
/**/	"Invalid X coordinate")
/**/DEBUG_CHECK_RANGE_LONG(GetyShortToDouble, YSource, 0L, NySource - 1L, Status,
/**/	"Invalid Y coordinate")
/**/DEBUG_CHECK_RANGE_LONG(GetyShortToDouble, ZSource, 0L, NzSource - 1L, Status,
/**/	"Invalid Z coordinate")
/**/DEBUG_CHECK_NULL_POINTER(GetyShortToDouble, ColumnDestination, Status,
/**/	"No output line")
/**/DEBUG_CHECK_RANGE_LONG(GetyShortToDouble, NyCopy, 1L, NySource - YSource, Status,
/**/	"Invalid height (negative or excessive)")
/**/DEBUG_RETURN_ON_ERROR(GetyShortToDouble, Status)
/**/DEBUG_WRITE_ENTERING(GetyShortToDouble,
/**/	"About to get a (double)column from a (short)volume")

	VolumeSource += (ptrdiff_t)(XSource + NxSource * (YSource + NySource * ZSource));
	End = ColumnDestination + (ptrdiff_t)NyCopy;
	while (ColumnDestination < End) {
		*ColumnDestination++ = (double)*VolumeSource;
		VolumeSource += NxSource;
	}
/**/DEBUG_WRITE_LEAVING(GetyShortToDouble, "Done")
	return(Status);
} /* end GetyShortToDouble */

/*--------------------------------------------------------------------------*/
extern int		GetzFloatToDouble
				(
					float	*VolumeSource,		/* float input data */
					long	NxSource,			/* width of the input */
					long	NySource,			/* height of the input */
					long	NzSource,			/* depth of the input */
					long	XSource,			/* x coordinate to get from */
					long	YSource,			/* y coordinate to get from */
					long	ZSource,			/* z coordinate to get from */
					double	PillarDestination[],/* double output data */
					long	NzCopy				/* length of the output */
				)

/* copies a line perpendicular to the screen from a float input volume into a 1D double output array */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin GetzFloatToDouble */

	double	*End;
	const ptrdiff_t
			NxNySource = (ptrdiff_t)(NxSource * NySource);
	int		Status = !ERROR;

#if defined(CODEWARRIOR) && !defined(DEBUG)
#pragma unused(NzSource)
#endif
/**/DEBUG_CHECK_NULL_POINTER(GetzFloatToDouble, VolumeSource, Status,
/**/	"No input volume")
/**/DEBUG_CHECK_RANGE_LONG(GetzFloatToDouble, NxSource, 1L, LONG_MAX, Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(GetzFloatToDouble, NySource, 1L, LONG_MAX, Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(GetzFloatToDouble, NzSource, 1L, LONG_MAX, Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(GetzFloatToDouble, XSource, 0L, NxSource - 1L, Status,
/**/	"Invalid X coordinate")
/**/DEBUG_CHECK_RANGE_LONG(GetzFloatToDouble, YSource, 0L, NySource - 1L, Status,
/**/	"Invalid Y coordinate")
/**/DEBUG_CHECK_RANGE_LONG(GetzFloatToDouble, ZSource, 0L, NzSource - 1L, Status,
/**/	"Invalid Z coordinate")
/**/DEBUG_CHECK_NULL_POINTER(GetzFloatToDouble, PillarDestination, Status,
/**/	"No output line")
/**/DEBUG_CHECK_RANGE_LONG(GetzFloatToDouble, NzCopy, 1L, NzSource - ZSource, Status,
/**/	"Invalid depth (negative or excessive)")
/**/DEBUG_RETURN_ON_ERROR(GetzFloatToDouble, Status)
/**/DEBUG_WRITE_ENTERING(GetzFloatToDouble,
/**/	"About to get a (double)pillar from a (float)volume")

	VolumeSource += (ptrdiff_t)(XSource + NxSource * (YSource + NySource * ZSource));
	End = PillarDestination + (ptrdiff_t)NzCopy;
	while (PillarDestination < End) {
		*PillarDestination++ = (double)*VolumeSource;
		VolumeSource += NxNySource;
	}
/**/DEBUG_WRITE_LEAVING(GetzFloatToDouble, "Done")
	return(Status);
} /* end GetzFloatToDouble */

/*--------------------------------------------------------------------------*/
extern int		GetzShortToDouble
				(
					short	*VolumeSource,		/* short input data */
					long	NxSource,			/* width of the input */
					long	NySource,			/* height of the input */
					long	NzSource,			/* depth of the input */
					long	XSource,			/* x coordinate to get from */
					long	YSource,			/* y coordinate to get from */
					long	ZSource,			/* z coordinate to get from */
					double	PillarDestination[],/* double output data */
					long	NzCopy				/* length of the output */
				)

/* copies a line perpendicular to the screen from a short input volume into a 1D double output array */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin GetzShortToDouble */

	double	*End;
	const ptrdiff_t
			NxNySource = (ptrdiff_t)(NxSource * NySource);
	int		Status = !ERROR;

#if defined(CODEWARRIOR) && !defined(DEBUG)
#pragma unused(NzSource)
#endif
/**/DEBUG_CHECK_NULL_POINTER(GetzShortToDouble, VolumeSource, Status,
/**/	"No input volume")
/**/DEBUG_CHECK_RANGE_LONG(GetzShortToDouble, NxSource, 1L, LONG_MAX, Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(GetzShortToDouble, NySource, 1L, LONG_MAX, Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(GetzShortToDouble, NzSource, 1L, LONG_MAX, Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(GetzShortToDouble, XSource, 0L, NxSource - 1L, Status,
/**/	"Invalid X coordinate")
/**/DEBUG_CHECK_RANGE_LONG(GetzShortToDouble, YSource, 0L, NySource - 1L, Status,
/**/	"Invalid Y coordinate")
/**/DEBUG_CHECK_RANGE_LONG(GetzShortToDouble, ZSource, 0L, NzSource - 1L, Status,
/**/	"Invalid Z coordinate")
/**/DEBUG_CHECK_NULL_POINTER(GetzShortToDouble, PillarDestination, Status,
/**/	"No output line")
/**/DEBUG_CHECK_RANGE_LONG(GetzShortToDouble, NzCopy, 1L, NzSource - ZSource, Status,
/**/	"Invalid depth (negative or excessive)")
/**/DEBUG_RETURN_ON_ERROR(GetzShortToDouble, Status)
/**/DEBUG_WRITE_ENTERING(GetzShortToDouble,
/**/	"About to get a (double)pillar from a (short)volume")

	VolumeSource += (ptrdiff_t)(XSource + NxSource * (YSource + NySource * ZSource));
	End = PillarDestination + (ptrdiff_t)NzCopy;
	while (PillarDestination < End) {
		*PillarDestination++ = (double)*VolumeSource;
		VolumeSource += NxNySource;
	}
/**/DEBUG_WRITE_LEAVING(GetzShortToDouble, "Done")
	return(Status);
} /* end GetzShortToDouble */

/*--------------------------------------------------------------------------*/
extern int		PutxDoubleToFloat
				(
					float	*VolumeDestination,	/* float output data */
					long	NxDestination,		/* width of the output */
					long	NyDestination,		/* height of the output */
					long	NzDestination,		/* depth of the output */
					long	XDestination,		/* x coordinate to put into */
					long	YDestination,		/* y coordinate to put into */
					long	ZDestination,		/* z coordinate to put into */
					double	RowSource[],		/* double input data */
					long	NxCopy				/* length of the input */
				)

/* copies a 1D double input array into a horizontal line of a float output volume */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin PutxDoubleToFloat */

	double	*End;
	int		Status = !ERROR;

#if defined(CODEWARRIOR) && !defined(DEBUG)
#pragma unused(NzDestination)
#endif
/**/DEBUG_CHECK_NULL_POINTER(PutxDoubleToFloat, VolumeDestination, Status,
/**/	"No output volume")
/**/DEBUG_CHECK_RANGE_LONG(PutxDoubleToFloat, NxDestination, 1L, LONG_MAX, Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(PutxDoubleToFloat, NyDestination, 1L, LONG_MAX, Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(PutxDoubleToFloat, NzDestination, 1L, LONG_MAX, Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(PutxDoubleToFloat, XDestination, 0L, NxDestination - 1L, Status,
/**/	"Invalid X coordinate")
/**/DEBUG_CHECK_RANGE_LONG(PutxDoubleToFloat, YDestination, 0L, NyDestination - 1L, Status,
/**/	"Invalid Y coordinate")
/**/DEBUG_CHECK_RANGE_LONG(PutxDoubleToFloat, ZDestination, 0L, NzDestination - 1L, Status,
/**/	"Invalid Z coordinate")
/**/DEBUG_CHECK_NULL_POINTER(PutxDoubleToFloat, RowSource, Status,
/**/	"No output line")
/**/DEBUG_CHECK_RANGE_LONG(PutxDoubleToFloat, NxCopy, 1L, NxDestination - XDestination, Status,
/**/	"Invalid width (negative or excessive)")
/**/DEBUG_RETURN_ON_ERROR(PutxDoubleToFloat, Status)
/**/DEBUG_WRITE_ENTERING(PutxDoubleToFloat,
/**/	"About to put a (double)row into a (float)volume")

	VolumeDestination += (ptrdiff_t)(XDestination + NxDestination * (YDestination
		+ NyDestination * ZDestination));
	End = RowSource + (ptrdiff_t)NxCopy;
	while (RowSource < End) {
		*VolumeDestination++ = ConvertDoubleToFloat(*RowSource++);
	}
/**/DEBUG_WRITE_LEAVING(PutxDoubleToFloat, "Done")
	return(Status);
} /* end PutxDoubleToFloat */

/*--------------------------------------------------------------------------*/
extern int		PutxDoubleToShort
				(
					short	*VolumeDestination,	/* short output data */
					long	NxDestination,		/* width of the output */
					long	NyDestination,		/* height of the output */
					long	NzDestination,		/* depth of the output */
					long	XDestination,		/* x coordinate to put into */
					long	YDestination,		/* y coordinate to put into */
					long	ZDestination,		/* z coordinate to put into */
					double	RowSource[],		/* double input data */
					long	NxCopy				/* length of the input */
				)

/* copies a 1D double input array into a horizontal line of a short output volume */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin PutxDoubleToShort */

	double	*End;
	int		Status = !ERROR;

#if defined(CODEWARRIOR) && !defined(DEBUG)
#pragma unused(NzDestination)
#endif
/**/DEBUG_CHECK_NULL_POINTER(PutxDoubleToShort, VolumeDestination, Status,
/**/	"No output volume")
/**/DEBUG_CHECK_RANGE_LONG(PutxDoubleToShort, NxDestination, 1L, LONG_MAX, Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(PutxDoubleToShort, NyDestination, 1L, LONG_MAX, Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(PutxDoubleToShort, NzDestination, 1L, LONG_MAX, Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(PutxDoubleToShort, XDestination, 0L, NxDestination - 1L, Status,
/**/	"Invalid X coordinate")
/**/DEBUG_CHECK_RANGE_LONG(PutxDoubleToShort, YDestination, 0L, NyDestination - 1L, Status,
/**/	"Invalid Y coordinate")
/**/DEBUG_CHECK_RANGE_LONG(PutxDoubleToShort, ZDestination, 0L, NzDestination - 1L, Status,
/**/	"Invalid Z coordinate")
/**/DEBUG_CHECK_NULL_POINTER(PutxDoubleToShort, RowSource, Status,
/**/	"No output line")
/**/DEBUG_CHECK_RANGE_LONG(PutxDoubleToShort, NxCopy, 1L, NxDestination - XDestination, Status,
/**/	"Invalid width (negative or excessive)")
/**/DEBUG_RETURN_ON_ERROR(PutxDoubleToShort, Status)
/**/DEBUG_WRITE_ENTERING(PutxDoubleToShort,
/**/	"About to put a (double)row into a (short)volume")

	VolumeDestination += (ptrdiff_t)(XDestination + NxDestination * (YDestination
		+ NyDestination * ZDestination));
	End = RowSource + (ptrdiff_t)NxCopy;
	while (RowSource < End) {
		*VolumeDestination++ = ConvertDoubleToShort(*RowSource++);
	}
/**/DEBUG_WRITE_LEAVING(PutxDoubleToShort, "Done")
	return(Status);
} /* end PutxDoubleToShort */

/*--------------------------------------------------------------------------*/
extern int		PutyDoubleToFloat
				(
					float	*VolumeDestination,	/* float output data */
					long	NxDestination,		/* width of the output */
					long	NyDestination,		/* height of the output */
					long	NzDestination,		/* depth of the output */
					long	XDestination,		/* x coordinate to put into */
					long	YDestination,		/* y coordinate to put into */
					long	ZDestination,		/* z coordinate to put into */
					double	ColumnSource[],		/* double input data */
					long	NyCopy				/* length of the input */
				)

/* copies a 1D double input array into a vertical line of a float output volume */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin PutyDoubleToFloat */

	double	*End;
	int		Status = !ERROR;

#if defined(CODEWARRIOR) && !defined(DEBUG)
#pragma unused(NzDestination)
#endif
/**/DEBUG_CHECK_NULL_POINTER(PutyDoubleToFloat, VolumeDestination, Status,
/**/	"No output volume")
/**/DEBUG_CHECK_RANGE_LONG(PutyDoubleToFloat, NxDestination, 1L, LONG_MAX, Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(PutyDoubleToFloat, NyDestination, 1L, LONG_MAX, Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(PutyDoubleToFloat, NzDestination, 1L, LONG_MAX, Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(PutyDoubleToFloat, XDestination, 0L, NxDestination - 1L, Status,
/**/	"Invalid X coordinate")
/**/DEBUG_CHECK_RANGE_LONG(PutyDoubleToFloat, YDestination, 0L, NyDestination - 1L, Status,
/**/	"Invalid Y coordinate")
/**/DEBUG_CHECK_RANGE_LONG(PutyDoubleToFloat, ZDestination, 0L, NzDestination - 1L, Status,
/**/	"Invalid Z coordinate")
/**/DEBUG_CHECK_NULL_POINTER(PutyDoubleToFloat, ColumnSource, Status,
/**/	"No output line")
/**/DEBUG_CHECK_RANGE_LONG(PutyDoubleToFloat, NyCopy, 1L, NyDestination - YDestination, Status,
/**/	"Invalid height (negative or excessive)")
/**/DEBUG_RETURN_ON_ERROR(PutyDoubleToFloat, Status)
/**/DEBUG_WRITE_ENTERING(PutyDoubleToFloat,
/**/	"About to put a (double)column into a (float)volume")

	VolumeDestination += (ptrdiff_t)(XDestination + NxDestination * (YDestination
		+ NyDestination * ZDestination));
	End = ColumnSource + (ptrdiff_t)NyCopy;
	while (ColumnSource < End) {
		*VolumeDestination = ConvertDoubleToFloat(*ColumnSource++);
		VolumeDestination += NxDestination;
	}
/**/DEBUG_WRITE_LEAVING(PutyDoubleToFloat, "Done")
	return(Status);
} /* end PutyDoubleToFloat */

/*--------------------------------------------------------------------------*/
extern int		PutyDoubleToShort
				(
					short	*VolumeDestination,	/* short output data */
					long	NxDestination,		/* width of the output */
					long	NyDestination,		/* height of the output */
					long	NzDestination,		/* depth of the output */
					long	XDestination,		/* x coordinate to put into */
					long	YDestination,		/* y coordinate to put into */
					long	ZDestination,		/* z coordinate to put into */
					double	ColumnSource[],		/* double input data */
					long	NyCopy				/* length of the input */
				)

/* copies a 1D double input array into a vertical line of a short output volume */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin PutyDoubleToShort */

	double	*End;
	int		Status = !ERROR;

#if defined(CODEWARRIOR) && !defined(DEBUG)
#pragma unused(NzDestination)
#endif
/**/DEBUG_CHECK_NULL_POINTER(PutyDoubleToShort, VolumeDestination, Status,
/**/	"No output volume")
/**/DEBUG_CHECK_RANGE_LONG(PutyDoubleToShort, NxDestination, 1L, LONG_MAX, Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(PutyDoubleToShort, NyDestination, 1L, LONG_MAX, Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(PutyDoubleToShort, NzDestination, 1L, LONG_MAX, Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(PutyDoubleToShort, XDestination, 0L, NxDestination - 1L, Status,
/**/	"Invalid X coordinate")
/**/DEBUG_CHECK_RANGE_LONG(PutyDoubleToShort, YDestination, 0L, NyDestination - 1L, Status,
/**/	"Invalid Y coordinate")
/**/DEBUG_CHECK_RANGE_LONG(PutyDoubleToShort, ZDestination, 0L, NzDestination - 1L, Status,
/**/	"Invalid Z coordinate")
/**/DEBUG_CHECK_NULL_POINTER(PutyDoubleToShort, ColumnSource, Status,
/**/	"No output line")
/**/DEBUG_CHECK_RANGE_LONG(PutyDoubleToShort, NyCopy, 1L, NyDestination - YDestination, Status,
/**/	"Invalid height (negative or excessive)")
/**/DEBUG_RETURN_ON_ERROR(PutyDoubleToShort, Status)
/**/DEBUG_WRITE_ENTERING(PutyDoubleToShort,
/**/	"About to put a (double)column into a (short)volume")

	VolumeDestination += (ptrdiff_t)(XDestination + NxDestination * (YDestination
		+ NyDestination * ZDestination));
	End = ColumnSource + (ptrdiff_t)NyCopy;
	while (ColumnSource < End) {
		*VolumeDestination = ConvertDoubleToShort(*ColumnSource++);
		VolumeDestination += NxDestination;
	}
/**/DEBUG_WRITE_LEAVING(PutyDoubleToShort, "Done")
	return(Status);
} /* end PutyDoubleToShort */

/*--------------------------------------------------------------------------*/
extern int		PutzDoubleToFloat
				(
					float	*VolumeDestination,	/* float output data */
					long	NxDestination,		/* width of the output */
					long	NyDestination,		/* height of the output */
					long	NzDestination,		/* depth of the output */
					long	XDestination,		/* x coordinate to put into */
					long	YDestination,		/* y coordinate to put into */
					long	ZDestination,		/* z coordinate to put into */
					double	PillarSource[],		/* double input data */
					long	NzCopy				/* length of the input */
				)

/* copies a 1D double input array into a line perpendicular to the screen of a float output volume */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin PutzDoubleToFloat */

	double	*End;
	const ptrdiff_t
			NxNyDestination = (ptrdiff_t)(NxDestination * NyDestination);
	int		Status = !ERROR;

#if defined(CODEWARRIOR) && !defined(DEBUG)
#pragma unused(NzDestination)
#endif
/**/DEBUG_CHECK_NULL_POINTER(PutzDoubleToFloat, VolumeDestination, Status,
/**/	"No output volume")
/**/DEBUG_CHECK_RANGE_LONG(PutzDoubleToFloat, NxDestination, 1L, LONG_MAX, Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(PutzDoubleToFloat, NyDestination, 1L, LONG_MAX, Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(PutzDoubleToFloat, NzDestination, 1L, LONG_MAX, Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(PutzDoubleToFloat, XDestination, 0L, NxDestination - 1L, Status,
/**/	"Invalid X coordinate")
/**/DEBUG_CHECK_RANGE_LONG(PutzDoubleToFloat, YDestination, 0L, NyDestination - 1L, Status,
/**/	"Invalid Y coordinate")
/**/DEBUG_CHECK_RANGE_LONG(PutzDoubleToFloat, ZDestination, 0L, NzDestination - 1L, Status,
/**/	"Invalid Z coordinate")
/**/DEBUG_CHECK_NULL_POINTER(PutzDoubleToFloat, PillarSource, Status,
/**/	"No output line")
/**/DEBUG_CHECK_RANGE_LONG(PutzDoubleToFloat, NzCopy, 1L, NzDestination - ZDestination, Status,
/**/	"Invalid depth (negative or excessive)")
/**/DEBUG_RETURN_ON_ERROR(PutzDoubleToFloat, Status)
/**/DEBUG_WRITE_ENTERING(PutzDoubleToFloat,
/**/	"About to put a (double)pillar into a (float)volume")

	VolumeDestination += (ptrdiff_t)(XDestination + NxDestination * (YDestination
		+ NyDestination * ZDestination));
	End = PillarSource + (ptrdiff_t)NzCopy;
	while (PillarSource < End) {
		*VolumeDestination = ConvertDoubleToFloat(*PillarSource++);
		VolumeDestination += NxNyDestination;
	}
/**/DEBUG_WRITE_LEAVING(PutzDoubleToFloat, "Done")
	return(Status);
} /* end PutzDoubleToFloat */

/*--------------------------------------------------------------------------*/
extern int		PutzDoubleToShort
				(
					short	*VolumeDestination,	/* short output data */
					long	NxDestination,		/* width of the output */
					long	NyDestination,		/* height of the output */
					long	NzDestination,		/* depth of the output */
					long	XDestination,		/* x coordinate to put into */
					long	YDestination,		/* y coordinate to put into */
					long	ZDestination,		/* z coordinate to put into */
					double	PillarSource[],		/* double input data */
					long	NzCopy				/* length of the input */
				)

/* copies a 1D double input array into a line perpendicular to the screen of a short output volume */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin PutzDoubleToShort */

	double	*End;
	const ptrdiff_t
			NxNyDestination = (ptrdiff_t)(NxDestination * NyDestination);
	int		Status = !ERROR;

#if defined(CODEWARRIOR) && !defined(DEBUG)
#pragma unused(NzDestination)
#endif
/**/DEBUG_CHECK_NULL_POINTER(PutzDoubleToShort, VolumeDestination, Status,
/**/	"No output volume")
/**/DEBUG_CHECK_RANGE_LONG(PutzDoubleToShort, NxDestination, 1L, LONG_MAX, Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(PutzDoubleToShort, NyDestination, 1L, LONG_MAX, Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(PutzDoubleToShort, NzDestination, 1L, LONG_MAX, Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(PutzDoubleToShort, XDestination, 0L, NxDestination - 1L, Status,
/**/	"Invalid X coordinate")
/**/DEBUG_CHECK_RANGE_LONG(PutzDoubleToShort, YDestination, 0L, NyDestination - 1L, Status,
/**/	"Invalid Y coordinate")
/**/DEBUG_CHECK_RANGE_LONG(PutzDoubleToShort, ZDestination, 0L, NzDestination - 1L, Status,
/**/	"Invalid Z coordinate")
/**/DEBUG_CHECK_NULL_POINTER(PutzDoubleToShort, PillarSource, Status,
/**/	"No output line")
/**/DEBUG_CHECK_RANGE_LONG(PutzDoubleToShort, NzCopy, 1L, NzDestination - ZDestination, Status,
/**/	"Invalid depth (negative or excessive)")
/**/DEBUG_RETURN_ON_ERROR(PutzDoubleToShort, Status)
/**/DEBUG_WRITE_ENTERING(PutzDoubleToShort,
/**/	"About to put a (double)pillar into a (short)volume")

	VolumeDestination += (ptrdiff_t)(XDestination + NxDestination * (YDestination
		+ NyDestination * ZDestination));
	End = PillarSource + (ptrdiff_t)NzCopy;
	while (PillarSource < End) {
		*VolumeDestination = ConvertDoubleToShort(*PillarSource++);
		VolumeDestination += NxNyDestination;
	}
/**/DEBUG_WRITE_LEAVING(PutzDoubleToShort, "Done")
	return(Status);
} /* end PutzDoubleToShort */

