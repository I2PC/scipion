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
#include	"getputd.h"
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
extern int		AllocateVolumeDouble
				(
					double	**Volume,			/* double output pointer */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					int		*Status				/* error management */
				)

/* Allocates a 3D array of double */
/* success: return(!ERROR); failure: return(ERROR) and set *Volume to NULL */

{ /* begin AllocateVolumeDouble */

	*Status = !ERROR;

/**/DEBUG_CHECK_RANGE_LONG(AllocateVolumeDouble, Nx, 1L, LONG_MAX, *Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(AllocateVolumeDouble, Ny, 1L, LONG_MAX, *Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(AllocateVolumeDouble, Nz, 1L, LONG_MAX, *Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(AllocateVolumeDouble, *Status)
/**/DEBUG_WRITE_ENTERING(AllocateVolumeDouble,
/**/	"About to allocate a 3D array of double")
#ifdef DEBUG
/**/if (*Volume != (double *)NULL) {
/**/	WRITE_WARNING(AllocateVolumeDouble, "Volume may have been previously allocated")
/**/}
#endif

	*Volume = (double *)malloc((size_t)(Nx * Ny * Nz * (long)sizeof(double)));
	if (*Volume == (double *)NULL) {
		*Status = ERROR;
		WRITE_ERROR(AllocateVolumeDouble, "Unable to perform allocation")
	}
/**/DEBUG_WRITE_LEAVING(AllocateVolumeDouble, "Done")
	return(*Status);
} /* end AllocateVolumeDouble */

/*--------------------------------------------------------------------------*/
extern int		CopyDoubleToFloat
				(
					double	*VolumeSource,		/* double input data */
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

/* copies a sub-volume from a double input volume into a portion of a float output volume */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin CopyDoubleToFloat */

	double	*P;
	float	*Q;
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
/**/DEBUG_CHECK_NULL_POINTER(CopyDoubleToFloat, VolumeSource, Status,
/**/	"No input")
/**/DEBUG_CHECK_RANGE_LONG(CopyDoubleToFloat, NxSource, 1L, LONG_MAX, Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyDoubleToFloat, NySource, 1L, LONG_MAX, Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyDoubleToFloat, NzSource, 1L, LONG_MAX, Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyDoubleToFloat, XSource, 0L, NxSource - 1L, Status,
/**/	"Invalid X coordinate")
/**/DEBUG_CHECK_RANGE_LONG(CopyDoubleToFloat, YSource, 0L, NySource - 1L, Status,
/**/	"Invalid Y coordinate")
/**/DEBUG_CHECK_RANGE_LONG(CopyDoubleToFloat, ZSource, 0L, NzSource - 1L, Status,
/**/	"Invalid Z coordinate")
/**/DEBUG_CHECK_NULL_POINTER(CopyDoubleToFloat, VolumeDestination, Status,
/**/	"No output")
/**/DEBUG_CHECK_RANGE_LONG(CopyDoubleToFloat, NxDestination, 1L, LONG_MAX, Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyDoubleToFloat, NyDestination, 1L, LONG_MAX, Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyDoubleToFloat, NzDestination, 1L, LONG_MAX, Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyDoubleToFloat, XDestination, 0L, NxDestination - 1L, Status,
/**/	"Invalid X coordinate")
/**/DEBUG_CHECK_RANGE_LONG(CopyDoubleToFloat, YDestination, 0L, NyDestination - 1L, Status,
/**/	"Invalid Y coordinate")
/**/DEBUG_CHECK_RANGE_LONG(CopyDoubleToFloat, ZDestination, 0L, NzDestination - 1L, Status,
/**/	"Invalid Z coordinate")
/**/DEBUG_CHECK_RANGE_LONG(CopyDoubleToFloat, NxCopy, 1L, NxSource - XSource, Status,
/**/	"Invalid width (negative or excessive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyDoubleToFloat, NyCopy, 1L, NySource - YSource, Status,
/**/	"Invalid height (negative or excessive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyDoubleToFloat, NzCopy, 1L, NzSource - ZSource, Status,
/**/	"Invalid depth (negative or excessive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyDoubleToFloat, NxCopy, 1L, NxDestination - XDestination, Status,
/**/	"Invalid width (negative or excessive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyDoubleToFloat, NyCopy, 1L, NyDestination - YDestination, Status,
/**/	"Invalid height (negative or excessive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyDoubleToFloat, NzCopy, 1L, NzDestination - ZDestination, Status,
/**/	"Invalid depth (negative or excessive)")
/**/DEBUG_RETURN_ON_ERROR(CopyDoubleToFloat, Status)
/**/DEBUG_WRITE_ENTERING(CopyDoubleToFloat,
/**/	"About to copy a double subvolume into a float subvolume")

	VolumeSource += (ptrdiff_t)(XSource + NxSource * (YSource + NySource * ZSource));
	VolumeDestination += (ptrdiff_t)(XDestination + NxDestination * (YDestination
		+ NyDestination * ZDestination));
	for (Z = -NzCopy; (Z < 0L); Z++) {
		P = VolumeSource;
		Q = VolumeDestination;
		for (Y = -NyCopy; (Y < 0L); Y++) {
			for (X = -NxCopy; (X < 0L); X++) {
				*Q++ = ConvertDoubleToFloat(*P++);
			}
			P += IgnoredSource;
			Q += IgnoredDestination;
		}
		VolumeSource += NxNySource;
		VolumeDestination += NxNyDestination;
	}
/**/DEBUG_WRITE_LEAVING(CopyDoubleToFloat, "Done")
	return(Status);
} /* end CopyDoubleToFloat */


/*--------------------------------------------------------------------------*/
extern int		CopyDoubleToDouble
				(
					double	*VolumeSource,		/* double input data */
					long	NxSource,			/* width of the input */
					long	NySource,			/* height of the input */
					long	NzSource,			/* depth of the input */
					long	XSource,			/* x coordinate to get from */
					long	YSource,			/* y coordinate to get from */
					long	ZSource,			/* z coordinate to get from */
					double	*VolumeDestination,	/* double output data */
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

/* copies a sub-volume from a double input volume into a portion of a double output volume */
/* input and output may not share their memory space */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin CopyDoubleToDouble */

	double	*P, *Q;
	const ptrdiff_t
			NxNySource = (ptrdiff_t)(NxSource * NySource),
			NxNyDestination = (ptrdiff_t)(NxDestination * NyDestination);
	long	Y, Z;
	int		Status = !ERROR;

#if defined(CODEWARRIOR) && !defined(DEBUG)
#pragma unused(NzSource, NzDestination)
#endif
/**/DEBUG_CHECK_NULL_POINTER(CopyDoubleToDouble, VolumeSource, Status,
/**/	"No input")
/**/DEBUG_CHECK_RANGE_LONG(CopyDoubleToDouble, NxSource, 1L, LONG_MAX, Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyDoubleToDouble, NySource, 1L, LONG_MAX, Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyDoubleToDouble, NzSource, 1L, LONG_MAX, Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyDoubleToDouble, XSource, 0L, NxSource - 1L, Status,
/**/	"Invalid X coordinate")
/**/DEBUG_CHECK_RANGE_LONG(CopyDoubleToDouble, YSource, 0L, NySource - 1L, Status,
/**/	"Invalid Y coordinate")
/**/DEBUG_CHECK_RANGE_LONG(CopyDoubleToDouble, ZSource, 0L, NzSource - 1L, Status,
/**/	"Invalid Z coordinate")
/**/DEBUG_CHECK_NULL_POINTER(CopyDoubleToDouble, VolumeDestination, Status,
/**/	"No output")
/**/DEBUG_CHECK_RANGE_LONG(CopyDoubleToDouble, NxDestination, 1L, LONG_MAX, Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyDoubleToDouble, NyDestination, 1L, LONG_MAX, Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyDoubleToDouble, NzDestination, 1L, LONG_MAX, Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyDoubleToDouble, XDestination, 0L, NxDestination - 1L, Status,
/**/	"Invalid X coordinate")
/**/DEBUG_CHECK_RANGE_LONG(CopyDoubleToDouble, YDestination, 0L, NyDestination - 1L, Status,
/**/	"Invalid Y coordinate")
/**/DEBUG_CHECK_RANGE_LONG(CopyDoubleToDouble, ZDestination, 0L, NzDestination - 1L, Status,
/**/	"Invalid Z coordinate")
/**/DEBUG_CHECK_RANGE_LONG(CopyDoubleToDouble, NxCopy, 1L, NxSource - XSource, Status,
/**/	"Invalid width (negative or excessive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyDoubleToDouble, NyCopy, 1L, NySource - YSource, Status,
/**/	"Invalid height (negative or excessive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyDoubleToDouble, NzCopy, 1L, NzSource - ZSource, Status,
/**/	"Invalid depth (negative or excessive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyDoubleToDouble, NxCopy, 1L, NxDestination - XDestination, Status,
/**/	"Invalid width (negative or excessive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyDoubleToDouble, NyCopy, 1L, NyDestination - YDestination, Status,
/**/	"Invalid height (negative or excessive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyDoubleToDouble, NzCopy, 1L, NzDestination - ZDestination, Status,
/**/	"Invalid depth (negative or excessive)")
/**/DEBUG_RETURN_ON_ERROR(CopyDoubleToDouble, Status)
/**/DEBUG_WRITE_ENTERING(CopyDoubleToDouble,
/**/	"About to copy a double subvolume into a double subvolume")
#ifdef DEBUG
/**/if (((VolumeSource + (ptrdiff_t)(NxSource * NySource * NzSource))
/**/	> (VolumeDestination + (ptrdiff_t)(XDestination + NxDestination
/**/	* (YDestination + NyDestination * ZDestination))))
/**/	&& ((VolumeDestination + (ptrdiff_t)(NxDestination * NyDestination
/**/	* NzDestination)) > (VolumeSource + (ptrdiff_t)(XSource + NxSource
/**/	* (YSource + NySource * ZSource))))) {
/**/	WRITE_WARNING(CopyDoubleToDouble, "Data overlap: memcpy may fail")
/**/}
#endif

	VolumeSource += (ptrdiff_t)(XSource + NxSource * (YSource + NySource * ZSource));
	VolumeDestination += (ptrdiff_t)(XDestination + NxDestination * (YDestination
		+ NyDestination * ZDestination));
	for (Z = -NzCopy; (Z < 0L); Z++) {
		P = VolumeSource;
		Q = VolumeDestination;
		for (Y = -NyCopy; (Y < 0L); Y++) {
			Q = (double *)memcpy(Q, P, (size_t)(NxCopy * (long)sizeof(double)));
			P += (ptrdiff_t)NxSource;
			Q += (ptrdiff_t)NxDestination;
		}
		VolumeSource += NxNySource;
		VolumeDestination += NxNyDestination;
	}
/**/DEBUG_WRITE_LEAVING(CopyDoubleToDouble, "Done")
	return(Status);
} /* end CopyDoubleToDouble */


/*--------------------------------------------------------------------------*/
extern int		CopyFloatToDouble
				(
					float	*VolumeSource,		/* float input data */
					long	NxSource,			/* width of the input */
					long	NySource,			/* height of the input */
					long	NzSource,			/* depth of the input */
					long	XSource,			/* x coordinate to get from */
					long	YSource,			/* y coordinate to get from */
					long	ZSource,			/* z coordinate to get from */
					double	*VolumeDestination,	/* double output data */
					long	NxDestination,		/* width of the output */
					long	NyDestination,		/* height of the output */
					long	NzDestination,		/* depth of the output */
					long	XDestination,		/* x coordinate to put into */
					long	YDestination,		/* y coordinate to put into */
					long	ZDestination,		/* z coordinate to put into */
					long	NxCopy,				/* width of block to copy */
					long	NyCopy,				/* height of block to copy */
					long	NzCopy				/* depth of the block to copy */
				)

/* copies a sub-volume from a float input volume into a portion of a double output volume */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin CopyFloatToDouble */

	double	*Q;
	float	*P;
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
/**/DEBUG_CHECK_NULL_POINTER(CopyFloatToDouble, VolumeSource, Status,
/**/	"No input")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToDouble, NxSource, 1L, LONG_MAX, Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToDouble, NySource, 1L, LONG_MAX, Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToDouble, NzSource, 1L, LONG_MAX, Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToDouble, XSource, 0L, NxSource - 1L, Status,
/**/	"Invalid X coordinate")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToDouble, YSource, 0L, NySource - 1L, Status,
/**/	"Invalid Y coordinate")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToDouble, ZSource, 0L, NzSource - 1L, Status,
/**/	"Invalid Z coordinate")
/**/DEBUG_CHECK_NULL_POINTER(CopyFloatToDouble, VolumeDestination, Status,
/**/	"No output")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToDouble, NxDestination, 1L, LONG_MAX, Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToDouble, NyDestination, 1L, LONG_MAX, Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToDouble, NzDestination, 1L, LONG_MAX, Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToDouble, XDestination, 0L, NxDestination - 1L, Status,
/**/	"Invalid X coordinate")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToDouble, YDestination, 0L, NyDestination - 1L, Status,
/**/	"Invalid Y coordinate")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToDouble, ZDestination, 0L, NzDestination - 1L, Status,
/**/	"Invalid Z coordinate")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToDouble, NxCopy, 1L, NxSource - XSource, Status,
/**/	"Invalid width (negative or excessive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToDouble, NyCopy, 1L, NySource - YSource, Status,
/**/	"Invalid height (negative or excessive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToDouble, NzCopy, 1L, NzSource - ZSource, Status,
/**/	"Invalid depth (negative or excessive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToDouble, NxCopy, 1L, NxDestination - XDestination, Status,
/**/	"Invalid width (negative or excessive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToDouble, NyCopy, 1L, NyDestination - YDestination, Status,
/**/	"Invalid height (negative or excessive)")
/**/DEBUG_CHECK_RANGE_LONG(CopyFloatToDouble, NzCopy, 1L, NzDestination - ZDestination, Status,
/**/	"Invalid depth (negative or excessive)")
/**/DEBUG_RETURN_ON_ERROR(CopyFloatToDouble, Status)
/**/DEBUG_WRITE_ENTERING(CopyFloatToDouble,
/**/	"About to copy a float subvolume into a double subvolume")

	VolumeSource += (ptrdiff_t)(XSource + NxSource * (YSource + NySource * ZSource));
	VolumeDestination += (ptrdiff_t)(XDestination + NxDestination * (YDestination
		+ NyDestination * ZDestination));
	for (Z = -NzCopy; (Z < 0L); Z++) {
		P = VolumeSource;
		Q = VolumeDestination;
		for (Y = -NyCopy; (Y < 0L); Y++) {
			for (X = -NxCopy; (X < 0L); X++) {
				*Q++ = (double)*P++;
			}
			P += IgnoredSource;
			Q += IgnoredDestination;
		}
		VolumeSource += NxNySource;
		VolumeDestination += NxNyDestination;
	}
/**/DEBUG_WRITE_LEAVING(CopyFloatToDouble, "Done")
	return(Status);
} /* end CopyFloatToDouble */

/*--------------------------------------------------------------------------*/
extern int		FreeVolumeDouble
				(
					double	**Volume			/* 3D double array */
				)

/* Frees a 3D array of double */
/* success: return(!ERROR) and set *Volume to NULL; failure: return(ERROR); */

{ /* begin FreeVolumeDouble */

	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(FreeVolumeDouble, *Volume, Status,
/**/	"Nothing to free")
/**/DEBUG_RETURN_ON_ERROR(FreeVolumeDouble, Status)
/**/DEBUG_WRITE_ENTERING(FreeVolumeDouble,
/**/	"About to free a 3D array of double")

	free(*Volume);
	*Volume = (double *)NULL;
/**/DEBUG_WRITE_LEAVING(FreeVolumeDouble, "Done")
	return(Status);
} /* end FreeVolumeDouble */

/*--------------------------------------------------------------------------*/
extern int		GetxDoubleToDouble
				(
					double	*VolumeSource,		/* double input data */
					long	NxSource,			/* width of the input */
					long	NySource,			/* height of the input */
					long	NzSource,			/* depth of the input */
					long	XSource,			/* x coordinate to get from */
					long	YSource,			/* y coordinate to get from */
					long	ZSource,			/* z coordinate to get from */
					double	RowDestination[],	/* double output data */
					long	NxCopy				/* length of the output */
				)

/* copies a horizontal line from a double input volume into a 1D double output array */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin GetxDoubleToDouble */

	int		Status = !ERROR;

#if defined(CODEWARRIOR) && !defined(DEBUG)
#pragma unused(NzSource)
#endif
/**/DEBUG_CHECK_NULL_POINTER(GetxDoubleToDouble, VolumeSource, Status,
/**/	"No input volume")
/**/DEBUG_CHECK_RANGE_LONG(GetxDoubleToDouble, NxSource, 1L, LONG_MAX, Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(GetxDoubleToDouble, NySource, 1L, LONG_MAX, Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(GetxDoubleToDouble, NzSource, 1L, LONG_MAX, Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(GetxDoubleToDouble, XSource, 0L, NxSource - 1L, Status,
/**/	"Invalid X coordinate")
/**/DEBUG_CHECK_RANGE_LONG(GetxDoubleToDouble, YSource, 0L, NySource - 1L, Status,
/**/	"Invalid Y coordinate")
/**/DEBUG_CHECK_RANGE_LONG(GetxDoubleToDouble, ZSource, 0L, NzSource - 1L, Status,
/**/	"Invalid Z coordinate")
/**/DEBUG_CHECK_NULL_POINTER(GetxDoubleToDouble, RowDestination, Status,
/**/	"No output line")
/**/DEBUG_CHECK_RANGE_LONG(GetxDoubleToDouble, NxCopy, 1L, NxSource - XSource, Status,
/**/	"Invalid width (negative or excessive)")
/**/DEBUG_RETURN_ON_ERROR(GetxDoubleToDouble, Status)
/**/DEBUG_WRITE_ENTERING(GetxDoubleToDouble,
/**/	"About to get a (double)row from a (double)volume")

	VolumeSource += (ptrdiff_t)(XSource + NxSource * (YSource + NySource * ZSource));
	RowDestination = (double *)memcpy(RowDestination, VolumeSource, (size_t)(NxCopy
		* (long)sizeof(double)));
/**/DEBUG_WRITE_LEAVING(GetxDoubleToDouble, "Done")
	return(Status);
} /* end GetxDoubleToDouble */

/*--------------------------------------------------------------------------*/
extern int		GetyDoubleToDouble
				(
					double	*VolumeSource,		/* double input data */
					long	NxSource,			/* width of the input */
					long	NySource,			/* height of the input */
					long	NzSource,			/* depth of the input */
					long	XSource,			/* x coordinate to get from */
					long	YSource,			/* y coordinate to get from */
					long	ZSource,			/* z coordinate to get from */
					double	ColumnDestination[],/* double output data */
					long	NyCopy				/* length of the output */
				)

/* copies a vertical line from a double input volume into a 1D double output array */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin GetyDoubleToDouble */

	double	*End;
	int		Status = !ERROR;

#if defined(CODEWARRIOR) && !defined(DEBUG)
#pragma unused(NzSource)
#endif
/**/DEBUG_CHECK_NULL_POINTER(GetyDoubleToDouble, VolumeSource, Status,
/**/	"No input volume")
/**/DEBUG_CHECK_RANGE_LONG(GetyDoubleToDouble, NxSource, 1L, LONG_MAX, Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(GetyDoubleToDouble, NySource, 1L, LONG_MAX, Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(GetyDoubleToDouble, NzSource, 1L, LONG_MAX, Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(GetyDoubleToDouble, XSource, 0L, NxSource - 1L, Status,
/**/	"Invalid X coordinate")
/**/DEBUG_CHECK_RANGE_LONG(GetyDoubleToDouble, YSource, 0L, NySource - 1L, Status,
/**/	"Invalid Y coordinate")
/**/DEBUG_CHECK_RANGE_LONG(GetyDoubleToDouble, ZSource, 0L, NzSource - 1L, Status,
/**/	"Invalid Z coordinate")
/**/DEBUG_CHECK_NULL_POINTER(GetyDoubleToDouble, ColumnDestination, Status,
/**/	"No output line")
/**/DEBUG_CHECK_RANGE_LONG(GetyDoubleToDouble, NyCopy, 1L, NySource - YSource, Status,
/**/	"Invalid height (negative or excessive)")
/**/DEBUG_RETURN_ON_ERROR(GetyDoubleToDouble, Status)
/**/DEBUG_WRITE_ENTERING(GetyDoubleToDouble,
/**/	"About to get a (double)column from a (double)volume")

	VolumeSource += (ptrdiff_t)(XSource + NxSource * (YSource + NySource * ZSource));
	End = ColumnDestination + (ptrdiff_t)NyCopy;
	while (ColumnDestination < End) {
		*ColumnDestination++ = *VolumeSource;
		VolumeSource += NxSource;
	}
/**/DEBUG_WRITE_LEAVING(GetyDoubleToDouble, "Done")
	return(Status);
} /* end GetyDoubleToDouble */

/*--------------------------------------------------------------------------*/
extern int		GetzDoubleToDouble
				(
					double	*VolumeSource,		/* double input data */
					long	NxSource,			/* width of the input */
					long	NySource,			/* height of the input */
					long	NzSource,			/* depth of the input */
					long	XSource,			/* x coordinate to get from */
					long	YSource,			/* y coordinate to get from */
					long	ZSource,			/* z coordinate to get from */
					double	PillarDestination[],/* double output data */
					long	NzCopy				/* length of the output */
				)

/* copies a line perpendicular to the screen from a double input volume into a 1D double output array */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin GetzDoubleToDouble */

	double	*End;
	const ptrdiff_t
			NxNySource = (ptrdiff_t)(NxSource * NySource);
	int		Status = !ERROR;

#if defined(CODEWARRIOR) && !defined(DEBUG)
#pragma unused(NzSource)
#endif
/**/DEBUG_CHECK_NULL_POINTER(GetzDoubleToDouble, VolumeSource, Status,
/**/	"No input volume")
/**/DEBUG_CHECK_RANGE_LONG(GetzDoubleToDouble, NxSource, 1L, LONG_MAX, Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(GetzDoubleToDouble, NySource, 1L, LONG_MAX, Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(GetzDoubleToDouble, NzSource, 1L, LONG_MAX, Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(GetzDoubleToDouble, XSource, 0L, NxSource - 1L, Status,
/**/	"Invalid X coordinate")
/**/DEBUG_CHECK_RANGE_LONG(GetzDoubleToDouble, YSource, 0L, NySource - 1L, Status,
/**/	"Invalid Y coordinate")
/**/DEBUG_CHECK_RANGE_LONG(GetzDoubleToDouble, ZSource, 0L, NzSource - 1L, Status,
/**/	"Invalid Z coordinate")
/**/DEBUG_CHECK_NULL_POINTER(GetzDoubleToDouble, PillarDestination, Status,
/**/	"No output line")
/**/DEBUG_CHECK_RANGE_LONG(GetzDoubleToDouble, NzCopy, 1L, NzSource - ZSource, Status,
/**/	"Invalid depth (negative or excessive)")
/**/DEBUG_RETURN_ON_ERROR(GetzDoubleToDouble, Status)
/**/DEBUG_WRITE_ENTERING(GetzDoubleToDouble,
/**/	"About to get a (double)pillar from a (double)volume")

	VolumeSource += (ptrdiff_t)(XSource + NxSource * (YSource + NySource * ZSource));
	End = PillarDestination + (ptrdiff_t)NzCopy;
	while (PillarDestination < End) {
		*PillarDestination++ = *VolumeSource;
		VolumeSource += NxNySource;
	}
/**/DEBUG_WRITE_LEAVING(GetzDoubleToDouble, "Done")
	return(Status);
} /* end GetzDoubleToDouble */

/*--------------------------------------------------------------------------*/
extern int		PutxDoubleToDouble
				(
					double	*VolumeDestination,	/* double output data */
					long	NxDestination,		/* width of the output */
					long	NyDestination,		/* height of the output */
					long	NzDestination,		/* depth of the output */
					long	XDestination,		/* x coordinate to put into */
					long	YDestination,		/* y coordinate to put into */
					long	ZDestination,		/* z coordinate to put into */
					double	RowSource[],		/* double input data */
					long	NxCopy				/* length of the input */
				)

/* copies a 1D double input array into a horizontal line of a double output volume */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin PutxDoubleToDouble */

	int		Status = !ERROR;

#if defined(CODEWARRIOR) && !defined(DEBUG)
#pragma unused(NzDestination)
#endif
/**/DEBUG_CHECK_NULL_POINTER(PutxDoubleToDouble, VolumeDestination, Status,
/**/	"No output volume")
/**/DEBUG_CHECK_RANGE_LONG(PutxDoubleToDouble, NxDestination, 1L, LONG_MAX, Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(PutxDoubleToDouble, NyDestination, 1L, LONG_MAX, Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(PutxDoubleToDouble, NzDestination, 1L, LONG_MAX, Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(PutxDoubleToDouble, XDestination, 0L, NxDestination - 1L, Status,
/**/	"Invalid X coordinate")
/**/DEBUG_CHECK_RANGE_LONG(PutxDoubleToDouble, YDestination, 0L, NyDestination - 1L, Status,
/**/	"Invalid Y coordinate")
/**/DEBUG_CHECK_RANGE_LONG(PutxDoubleToDouble, ZDestination, 0L, NzDestination - 1L, Status,
/**/	"Invalid Z coordinate")
/**/DEBUG_CHECK_NULL_POINTER(PutxDoubleToDouble, RowSource, Status,
/**/	"No output line")
/**/DEBUG_CHECK_RANGE_LONG(PutxDoubleToDouble, NxCopy, 1L, NxDestination - XDestination, Status,
/**/	"Invalid width (negative or excessive)")
/**/DEBUG_RETURN_ON_ERROR(PutxDoubleToDouble, Status)
/**/DEBUG_WRITE_ENTERING(PutxDoubleToDouble,
/**/	"About to put a (double)row into a (double)volume")

	VolumeDestination += (ptrdiff_t)(XDestination + NxDestination * (YDestination
		+ NyDestination * ZDestination));
	VolumeDestination = (double *)memcpy(VolumeDestination, RowSource, (size_t)(NxCopy
		* (long)sizeof(double)));
/**/DEBUG_WRITE_LEAVING(PutxDoubleToDouble, "Done")
	return(Status);
} /* end PutxDoubleToDouble */

/*--------------------------------------------------------------------------*/
extern int		PutyDoubleToDouble
				(
					double	*VolumeDestination,	/* double output data */
					long	NxDestination,		/* width of the output */
					long	NyDestination,		/* height of the output */
					long	NzDestination,		/* depth of the output */
					long	XDestination,		/* x coordinate to put into */
					long	YDestination,		/* y coordinate to put into */
					long	ZDestination,		/* z coordinate to put into */
					double	ColumnSource[],		/* double input data */
					long	NyCopy				/* length of the input */
				)

/* copies a 1D double input array into a vertical line of a double output volume */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin PutyDoubleToDouble */

	double	*End;
	int		Status = !ERROR;

#if defined(CODEWARRIOR) && !defined(DEBUG)
#pragma unused(NzDestination)
#endif
/**/DEBUG_CHECK_NULL_POINTER(PutyDoubleToDouble, VolumeDestination, Status,
/**/	"No output volume")
/**/DEBUG_CHECK_RANGE_LONG(PutyDoubleToDouble, NxDestination, 1L, LONG_MAX, Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(PutyDoubleToDouble, NyDestination, 1L, LONG_MAX, Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(PutyDoubleToDouble, NzDestination, 1L, LONG_MAX, Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(PutyDoubleToDouble, XDestination, 0L, NxDestination - 1L, Status,
/**/	"Invalid X coordinate")
/**/DEBUG_CHECK_RANGE_LONG(PutyDoubleToDouble, YDestination, 0L, NyDestination - 1L, Status,
/**/	"Invalid Y coordinate")
/**/DEBUG_CHECK_RANGE_LONG(PutyDoubleToDouble, ZDestination, 0L, NzDestination - 1L, Status,
/**/	"Invalid Z coordinate")
/**/DEBUG_CHECK_NULL_POINTER(PutyDoubleToDouble, ColumnSource, Status,
/**/	"No output line")
/**/DEBUG_CHECK_RANGE_LONG(PutyDoubleToDouble, NyCopy, 1L, NyDestination - YDestination, Status,
/**/	"Invalid height (negative or excessive)")
/**/DEBUG_RETURN_ON_ERROR(PutyDoubleToDouble, Status)
/**/DEBUG_WRITE_ENTERING(PutyDoubleToDouble,
/**/	"About to put a (double)column into a (double)volume")

	VolumeDestination += (ptrdiff_t)(XDestination + NxDestination * (YDestination
		+ NyDestination * ZDestination));
	End = ColumnSource + (ptrdiff_t)NyCopy;
	while (ColumnSource < End) {
		*VolumeDestination = *ColumnSource++;
		VolumeDestination += NxDestination;
	}
/**/DEBUG_WRITE_LEAVING(PutyDoubleToDouble, "Done")
	return(Status);
} /* end PutyDoubleToDouble */

/*--------------------------------------------------------------------------*/
extern int		PutzDoubleToDouble
				(
					double	*VolumeDestination,	/* double output data */
					long	NxDestination,		/* width of the output */
					long	NyDestination,		/* height of the output */
					long	NzDestination,		/* depth of the output */
					long	XDestination,		/* x coordinate to put into */
					long	YDestination,		/* y coordinate to put into */
					long	ZDestination,		/* z coordinate to put into */
					double	PillarSource[],		/* double input data */
					long	NzCopy				/* length of the input */
				)

/* copies a 1D double input array into a line perpendicular to the screen of a double output volume */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin PutzDoubleToDouble */

	double	*End;
	const ptrdiff_t
			NxNyDestination = (ptrdiff_t)(NxDestination * NyDestination);
	int		Status = !ERROR;

#if defined(CODEWARRIOR) && !defined(DEBUG)
#pragma unused(NzDestination)
#endif
/**/DEBUG_CHECK_NULL_POINTER(PutzDoubleToDouble, VolumeDestination, Status,
/**/	"No output volume")
/**/DEBUG_CHECK_RANGE_LONG(PutzDoubleToDouble, NxDestination, 1L, LONG_MAX, Status,
/**/	"Invalid width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(PutzDoubleToDouble, NyDestination, 1L, LONG_MAX, Status,
/**/	"Invalid height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(PutzDoubleToDouble, NzDestination, 1L, LONG_MAX, Status,
/**/	"Invalid depth (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(PutzDoubleToDouble, XDestination, 0L, NxDestination - 1L, Status,
/**/	"Invalid X coordinate")
/**/DEBUG_CHECK_RANGE_LONG(PutzDoubleToDouble, YDestination, 0L, NyDestination - 1L, Status,
/**/	"Invalid Y coordinate")
/**/DEBUG_CHECK_RANGE_LONG(PutzDoubleToDouble, ZDestination, 0L, NzDestination - 1L, Status,
/**/	"Invalid Z coordinate")
/**/DEBUG_CHECK_NULL_POINTER(PutzDoubleToDouble, PillarSource, Status,
/**/	"No output line")
/**/DEBUG_CHECK_RANGE_LONG(PutzDoubleToDouble, NzCopy, 1L, NzDestination - ZDestination, Status,
/**/	"Invalid depth (negative or excessive)")
/**/DEBUG_RETURN_ON_ERROR(PutzDoubleToDouble, Status)
/**/DEBUG_WRITE_ENTERING(PutzDoubleToDouble,
/**/	"About to put a (double)pillar into a (double)volume")

	VolumeDestination += (ptrdiff_t)(XDestination + NxDestination * (YDestination
		+ NyDestination * ZDestination));
	End = PillarSource + (ptrdiff_t)NzCopy;
	while (PillarSource < End) {
		*VolumeDestination = *PillarSource++;
		VolumeDestination += NxNyDestination;
	}
/**/DEBUG_WRITE_LEAVING(PutzDoubleToDouble, "Done")
	return(Status);
} /* end PutzDoubleToDouble */

