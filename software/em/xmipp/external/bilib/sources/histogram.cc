/*****************************************************************************
 *	System includes
 ****************************************************************************/
#include	<math.h>
#include	<stddef.h>
#include	<stdio.h>
#include	<stdlib.h>

/*****************************************************************************
 *	Toolbox defines
 ****************************************************************************/
#include	"configs.h"
#include	"debug.h"
#include	"error.h"

/*****************************************************************************
 *	Toolbox includes
 ****************************************************************************/
#include	"histogram.h"
#include	"messagedisplay.h"

/*****************************************************************************
 *	Conditional includes
 ****************************************************************************/
#ifdef		DEBUG
#include	<float.h>
#include	<limits.h>
#include	<string.h>
#endif

/*****************************************************************************
 *	Local struct
 ****************************************************************************/
/*--------------------------------------------------------------------------*/
struct TOrderedTree
{
	float		Value;
	long		Occurence;
	struct TOrderedTree
				*LowerBranch, *HigherBranch, *Trunk;
};

/*****************************************************************************
 *	Declaration of static procedures
 ****************************************************************************/
/*--------------------------------------------------------------------------*/
static void		FlattenAndFreeTree
				(
					struct TOrderedTree
							**Root,				/* root of the ordered tree to flatten */
					double	Frequency[],		/* ordinates */
					float	Value[]				/* abscissa */
				);

/*--------------------------------------------------------------------------*/
static void		FreeTree
				(
					struct TOrderedTree
							**Root				/* root of the ordered tree to free */
				);

/*****************************************************************************
 *	Definition of static procedures
 ****************************************************************************/
/*--------------------------------------------------------------------------*/
static void		FlattenAndFreeTree
				(
					struct TOrderedTree
							**Root,				/* root of the ordered tree to free */
					double	Frequency[],		/* ordinates */
					float	Value[]				/* abscissa */
				)

{ /* begin FlattenAndFreeTree */

	struct TOrderedTree
			*u, *v;

/**/DEBUG_WRITE_ENTERING(FlattenAndFreeTree,
/**/	"About to execute FlattenAndFreeTree")
	u = *Root;
	while (u != (struct TOrderedTree *)NULL) {
		if (u->LowerBranch != (struct TOrderedTree *)NULL) {
			v = u->LowerBranch;
			u->LowerBranch = (struct TOrderedTree *)NULL;
			u = v;
			continue;
		}
		if (u->HigherBranch != (struct TOrderedTree *)NULL) {
			*Value++ = u->Value;
			*Frequency++ = (double)(u->Occurence);
			u->Occurence = -1L;
			v = u->HigherBranch;
			u->HigherBranch = (struct TOrderedTree *)NULL;
			u = v;
			continue;
		}
		if (u->Occurence != -1L) {
			*Value++ = u->Value;
			*Frequency++ = (double)(u->Occurence);
		}
		v = u;
		u = u->Trunk;
		free(v);
	}
	*Root = (struct TOrderedTree *)NULL;
/**/DEBUG_WRITE_LEAVING(FlattenAndFreeTree, "Done")
} /* end FlattenAndFreeTree */

/*--------------------------------------------------------------------------*/
static void		FreeTree
				(
					struct TOrderedTree
							**Root				/* root of the ordered tree to free */
				)

{ /* begin FreeTree */

	struct TOrderedTree
			*u, *v;

/**/DEBUG_WRITE_ENTERING(FreeTree,
/**/	"About to execute FreeTree")
	u = *Root;
	while (u != (struct TOrderedTree *)NULL) {
		if (u->LowerBranch != (struct TOrderedTree *)NULL) {
			v = u->LowerBranch;
			u->LowerBranch = (struct TOrderedTree *)NULL;
			u = v;
			continue;
		}
		if (u->HigherBranch != (struct TOrderedTree *)NULL) {
			v = u->HigherBranch;
			u->HigherBranch = (struct TOrderedTree *)NULL;
			u = v;
			continue;
		}
		v = u;
		u = u->Trunk;
		free(v);
	}
	*Root = (struct TOrderedTree *)NULL;
/**/DEBUG_WRITE_LEAVING(FreeTree, "Done")
} /* end FreeTree */

/*****************************************************************************
 *	Definition of extern procedures
 ****************************************************************************/
/*--------------------------------------------------------------------------*/
extern int		HistogramBuild
				(
					float	*VolumeSource,		/* data to process */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					double	Frequency[],		/* output vector of ordinates */
					float	Value[],			/* output vector of abscissa */
					int		*Status				/* error management */
				)

/* computation of the frequencies of occurence of data values */
/* VolumeSource is a (float)volume of size (Nx x Ny x Nz) */
/* Value[] is a (float)array of length previously determined by HistogramGetSize */
/* the returned content of Value[] is sorted in strict ascending order */
/* Frequency[] is a (double)array of length previously determined by HistogramGetSize */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin HistogramBuild */

	struct TOrderedTree
			*Root;
	struct TOrderedTree
			*u, *v;
	float	*p;
	double	Normalization;
	long	HistogramLength;
	long	i, Nxyz = Nx * Ny * Nz;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(HistogramBuild, VolumeSource, *Status,
/**/	"Missing input VolumeSource")
/**/DEBUG_CHECK_NULL_POINTER(HistogramBuild, Value, *Status,
/**/	"Missing output Value")
/**/DEBUG_CHECK_NULL_POINTER(HistogramBuild, Frequency, *Status,
/**/	"Missing output Frequency")
/**/DEBUG_CHECK_RANGE_LONG(HistogramBuild, Nx, 1L, LONG_MAX, *Status,
/**/	"Invalid volume width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(HistogramBuild, Ny, 1L, LONG_MAX, *Status,
/**/	"Invalid volume height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(HistogramBuild, Nz, 1L, LONG_MAX, *Status,
/**/	"Invalid volume depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(HistogramBuild, *Status)
/**/DEBUG_WRITE_ENTERING(HistogramBuild,
/**/	"About to execute HistogramBuild")

	p = VolumeSource + (ptrdiff_t)(Nxyz - 1L);
	Root = (struct TOrderedTree *)malloc(sizeof(struct TOrderedTree));
	if (Root == (struct TOrderedTree *)NULL) {
		*Status = ERROR;
		WRITE_ERROR(HistogramBuild, "Unable to perform allocation of Root")
/**/	DEBUG_WRITE_LEAVING(HistogramBuild, "Done")
		return(*Status);
	}
	Root->Value = *p;
	HistogramLength = 1L;
	Root->Occurence = 1L;
	Root->LowerBranch = (struct TOrderedTree *)NULL;
	Root->HigherBranch = (struct TOrderedTree *)NULL;
	Root->Trunk = (struct TOrderedTree *)NULL;
	while (VolumeSource <= --p) {
		u = Root;
		do {
			if (*p == u->Value) {
				(u->Occurence)++;
				u = Root;
			}
			else {
				if (*p < u->Value) {
					v = u->LowerBranch;
					if (v == (struct TOrderedTree *)NULL) {
						v = (struct TOrderedTree *)malloc(sizeof(struct TOrderedTree));
						if (v == (struct TOrderedTree *)NULL) {
							*Status = ERROR;
							FreeTree(&Root);
							WRITE_ERROR(HistogramBuild, "Unable to perform allocation of a low leaf")
/**/						DEBUG_WRITE_LEAVING(HistogramBuild, "Done")
							return(*Status);
						}
						v->Value = *p;
						HistogramLength++;
						v->Occurence = 1L;
						v->LowerBranch = (struct TOrderedTree *)NULL;
						v->HigherBranch = (struct TOrderedTree *)NULL;
						v->Trunk = u;
						u->LowerBranch = v;
						u = Root;
					}
					else {
						u = v;
					}
				}
				else {
					v = u->HigherBranch;
					if (v == (struct TOrderedTree *)NULL) {
						v = (struct TOrderedTree *)malloc(sizeof(struct TOrderedTree));
						if (v == (struct TOrderedTree *)NULL) {
							*Status = ERROR;
							FreeTree(&Root);
							WRITE_ERROR(HistogramBuild, "Unable to perform allocation of a high leaf")
/**/						DEBUG_WRITE_LEAVING(HistogramBuild, "Done")
							return(*Status);
						}
						v->Value = *p;
						HistogramLength++;
						v->Occurence = 1L;
						v->LowerBranch = (struct TOrderedTree *)NULL;
						v->HigherBranch = (struct TOrderedTree *)NULL;
						v->Trunk = u;
						u->HigherBranch = v;
						u = Root;
					}
					else {
						u = v;
					}
				}
			}
		} while (u != Root);
	}
	FlattenAndFreeTree(&Root, Frequency, Value);
	Normalization = 1.0 / (double)Nxyz;
	for (i = 0L; (i < HistogramLength); i++) {
		Frequency[i] *= Normalization;
	}
/**/DEBUG_WRITE_LEAVING(HistogramBuild, "Done")
	return(*Status);
} /* end HistogramBuild */

/*--------------------------------------------------------------------------*/
extern int		HistogramEqualize
				(
					double	Frequency[],		/* histogram ordinates */
					float	Value[],			/* histogram abscissa */
					float	EqualizedValue[],	/* output vector of abscissa */
					long	HistogramLength,	/* length of the histogram */
					long	*NumberOfClasses,	/* number of classes, desired -> actual */
					double	Tolerance			/* admissible relative error */
				)

/* construction of the lookup table: Value[k] <-> EqualizedValue[k] */
/* EqualizedValue[] satisfies */
/*		sum(k in K(n0)) Frequency[k] ~=~ sum(k in K(n1)) Frequency[k] */
/*		for all n0, n1 in [0L, NumberOfClasses - 1L] */
/* where ~=~ means "is about equal to", */
/* and where K(n) is a domain such that */
/*		EqualizedValue[k0] = (sum(k1 in K(n)) Frequency[k1] * Value[k1]) */
/*			/ (sum(k2 in K(n)) Frequency[k2]) */
/*		for all k0 in K(n) */
/* under the constraint */
/*		DistinctElements(QuantizedValue[]) == NumberOfClasses */
/* Frequency[] is a (double)array of length HistogramLength */
/* the content of Frequency[] must be strictly positive */
/* the content of Frequency[] must have unit sum */
/* Value[] is a (float)array of length HistogramLength */
/* the content of Value[] must be sorted in strictly ascending order */
/* EqualizedValue[] is a returned (float)array of length HistogramLength */
/* on input, NumberOfClasses indicates the desired number of classes */
/* on output, NumberOfClasses returns the effective number of classes */
/* NumberOfClasses is no greater than (1.0 / max(Frequency[])), and never increases */
/* it may happen that the only solution that satisfies all constraints is undesirable */
/*		e.g.,	Frequency[] = {0.9, 0.1}; */
/*				Value[] = {10.0F, 90.0F}; */
/*				NumberOfClasses = 2L; (desired) */
/*		results in */
/*				QuantizedValues[] = {18.0F, 18.0F}; */
/*				NumberOfClasses = 1L; (actual) */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin HistogramEqualize */

	double	Mean, Norm, Threshold;
	long	i, j, k;
	int		LessClasses;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(HistogramEqualize, Frequency, Status,
/**/	"Missing input Frequency")
/**/DEBUG_CHECK_NULL_POINTER(HistogramEqualize, Value, Status,
/**/	"Missing input Value")
/**/DEBUG_CHECK_NULL_POINTER(HistogramEqualize, EqualizedValue, Status,
/**/	"Missing output EqualizedValue")
/**/DEBUG_CHECK_NULL_POINTER(HistogramEqualize, NumberOfClasses, Status,
/**/	"Missing in-place NumberOfClasses")
/**/DEBUG_CHECK_RANGE_LONG(HistogramEqualize, HistogramLength, 1L, LONG_MAX, Status,
/**/	"Invalid HistogramLength (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(HistogramEqualize, *NumberOfClasses, 1L, LONG_MAX, Status,
/**/	"Invalid NumberOfClasses (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_DOUBLE(HistogramEqualize, Tolerance, 0.0, DBL_MAX, Status,
/**/	"Invalid Tolerance (should be positive)")
/**/DEBUG_RETURN_ON_ERROR(HistogramEqualize, Status)
/**/DEBUG_WRITE_ENTERING(HistogramEqualize,
/**/	"About to execute HistogramEqualize")

	if (HistogramLength < *NumberOfClasses) {
		WRITE_WARNING(HistogramEqualize, "Reducing the number of classes")
		*NumberOfClasses = HistogramLength;
	}
	Threshold = 1.0 / (double)*NumberOfClasses;
	Norm = 0.0;
	for (i = 0L; (i < HistogramLength); i++) {
		if (Frequency[i] <= 0.0) {
			Status = ERROR;
			WRITE_ERROR(HistogramEqualize, "Invalid Frequency[i] (should be strictly positive)")
/**/		DEBUG_WRITE_LEAVING(HistogramEqualize, "Done")
			return(Status);
		}
		if (Threshold < Frequency[i]) {
			WRITE_WARNING(HistogramEqualize, "Reducing the number of classes")
			*NumberOfClasses = (long)floor(1.0 / Frequency[i]);
			Threshold = 1.0 / (double)*NumberOfClasses;
		}
		Norm += Frequency[i];
	}
	if (Tolerance < ((Norm - 1.0) * (Norm - 1.0))) {
		Status = ERROR;
		WRITE_ERROR(HistogramEqualize, "Invalid Frequency[] (should have unit sum)")
/**/	DEBUG_WRITE_LEAVING(HistogramEqualize, "Done")
		return(Status);
	}
	for (i = 1L; (i < HistogramLength); i++) {
		if (Value[i] <= Value[i - 1L]) {
			Status = ERROR;
			WRITE_ERROR(HistogramEqualize, "Invalid Value[] (should be strictly increasing)")
/**/		DEBUG_WRITE_LEAVING(HistogramEqualize, "Done")
			return(Status);
		}
	}
	do {
		LessClasses = FALSE;
		EqualizedValue[0] = (float)Frequency[0];
		for (i = 1L; (i < HistogramLength); i++) {
			EqualizedValue[i] = (float)((double)EqualizedValue[i - 1L] + Frequency[i]);
		}
		k = 0L;
		for (i = 1L; (i < *NumberOfClasses); i++) {
			Mean = 0.0;
			Norm = 0.0;
			Threshold = (float)((double)i / (double)*NumberOfClasses);
			j = k;
			while (EqualizedValue[j] <= Threshold) {
				Mean += Frequency[j] * (double)Value[j];
				Norm += Frequency[j++];
			}
			if (Norm < Tolerance) {
				WRITE_WARNING(HistogramEqualize, "Reducing the number of classes")
				(*NumberOfClasses)--;
				LessClasses = TRUE;
				break;
			}
			Mean /= Norm;
			while (EqualizedValue[k] <= Threshold) {
				EqualizedValue[k++] = (float)Mean;
			}
		}
		if (LessClasses) {
			continue;
		}
		Mean = 0.0;
		Norm = 0.0;
		j = k;
		while (j < HistogramLength) {
			Mean += Frequency[j] * (double)Value[j];
			Norm += Frequency[j++];
		}
		if (Norm < Tolerance) {
			WRITE_WARNING(HistogramEqualize, "Reducing the number of classes")
			(*NumberOfClasses)--;
			LessClasses = TRUE;
			continue;
		}
		Mean /= Norm;
		while (k < HistogramLength) {
			EqualizedValue[k++] = (float)Mean;
		}
	} while (LessClasses);
/**/DEBUG_WRITE_LEAVING(HistogramEqualize, "Done")
	return(Status);
} /* end HistogramEqualize */

/*--------------------------------------------------------------------------*/
extern int		HistogramGetSize
				(
					float	*VolumeSource,		/* data to process */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					long	*HistogramLength,	/* output length of the histogram */
					int		*Status				/* error management */
				)

/* determination of the number of differing data values in a volume */
/* VolumeSource is a (float)volume of size (Nx x Ny x Nz) */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin HistogramGetSize */

	struct TOrderedTree
			*Root;
	struct TOrderedTree
			*u, *v;
	float	*p;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(HistogramGetSize, VolumeSource, *Status,
/**/	"Missing input VolumeSource")
/**/DEBUG_CHECK_NULL_POINTER(HistogramGetSize, HistogramLength, *Status,
/**/	"Missing output HistogramLength")
/**/DEBUG_CHECK_RANGE_LONG(HistogramGetSize, Nx, 1L, LONG_MAX, *Status,
/**/	"Invalid volume width (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(HistogramGetSize, Ny, 1L, LONG_MAX, *Status,
/**/	"Invalid volume height (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(HistogramGetSize, Nz, 1L, LONG_MAX, *Status,
/**/	"Invalid volume depth (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(HistogramGetSize, *Status)
/**/DEBUG_WRITE_ENTERING(HistogramGetSize,
/**/	"About to execute HistogramGetSize")

	p = VolumeSource + (ptrdiff_t)(Nx * Ny * Nz - 1L);
	Root = (struct TOrderedTree *)malloc(sizeof(struct TOrderedTree));
	if (Root == (struct TOrderedTree *)NULL) {
		*Status = ERROR;
		WRITE_ERROR(HistogramGetSize, "Unable to perform allocation of Root")
/**/	DEBUG_WRITE_LEAVING(HistogramGetSize, "Done")
		return(*Status);
	}
	Root->Value = *p;
	*HistogramLength = 1L;
	Root->LowerBranch = (struct TOrderedTree *)NULL;
	Root->HigherBranch = (struct TOrderedTree *)NULL;
	Root->Trunk = (struct TOrderedTree *)NULL;
	while (VolumeSource <= --p) {
		u = Root;
		do {
			if (*p == u->Value) {
				u = Root;
			}
			else {
				if (*p < u->Value) {
					v = u->LowerBranch;
					if (v == (struct TOrderedTree *)NULL) {
						v = (struct TOrderedTree *)malloc(sizeof(struct TOrderedTree));
						if (v == (struct TOrderedTree *)NULL) {
							*Status = ERROR;
							FreeTree(&Root);
							WRITE_ERROR(HistogramGetSize, "Unable to perform allocation of a low leaf")
/**/						DEBUG_WRITE_LEAVING(HistogramGetSize, "Done")
							return(*Status);
						}
						v->Value = *p;
						(*HistogramLength)++;
						v->LowerBranch = (struct TOrderedTree *)NULL;
						v->HigherBranch = (struct TOrderedTree *)NULL;
						v->Trunk = u;
						u->LowerBranch = v;
						u = Root;
					}
					else {
						u = v;
					}
				}
				else {
					v = u->HigherBranch;
					if (v == (struct TOrderedTree *)NULL) {
						v = (struct TOrderedTree *)malloc(sizeof(struct TOrderedTree));
						if (v == (struct TOrderedTree *)NULL) {
							*Status = ERROR;
							FreeTree(&Root);
							WRITE_ERROR(HistogramGetSize, "Unable to perform allocation of a high leaf")
/**/						DEBUG_WRITE_LEAVING(HistogramGetSize, "Done")
							return(*Status);
						}
						v->Value = *p;
						(*HistogramLength)++;
						v->LowerBranch = (struct TOrderedTree *)NULL;
						v->HigherBranch = (struct TOrderedTree *)NULL;
						v->Trunk = u;
						u->HigherBranch = v;
						u = Root;
					}
					else {
						u = v;
					}
				}
			}
		} while (u != Root);
	}
	FreeTree(&Root);
/**/DEBUG_WRITE_LEAVING(HistogramGetSize, "Done")
	return(*Status);
} /* end HistogramGetSize */

/*--------------------------------------------------------------------------*/
extern int		HistogramKMeans
				(
					double	Frequency[],		/* histogram ordinates */
					float	Value[],			/* histogram abscissa */
					float	QuantizedValue[],	/* output vector of abscissa */
					long	HistogramLength,	/* length of the histogram */
					long	*NumberOfClasses,	/* number of classes, desired -> actual */
					double	Tolerance,			/* admissible relative error */
					int		*Status				/* error management */
				)

/* construction of the lookup table */
/*		Value[k] <-> QuantizedValue[k] */
/* minimization of */
/*		sum(k) Frequency[k] * (Value[k] - QuantizedValue[k])^2 */
/* under the constraint */
/*		DistinctElements(QuantizedValue[]) == NumberOfClasses */
/* Frequency[] is a (double)array of length HistogramLength */
/* the content of Frequency[] must be strictly positive */
/* the content of Frequency[] must have unit sum */
/* Value[] is a (float)array of length HistogramLength */
/* the content of Value[] must be sorted in strictly ascending order */
/* QuantizedValue[] is a returned (float)array of length HistogramLength */
/* on input, NumberOfClasses indicates the desired number of classes */
/* on output, NumberOfClasses returns the effective number of classes */
/* NumberOfClasses never increases */
/* important cases that go undetected (unfortunately): */
/* 1. convergence to a non-global optimum */
/*		e.g.,	Frequency[] = {0.25, 0.25, 0.25, 0.25}; */
/*				Value[] = {1.0F, 2.0F, 10.0F, 20.0F}; */
/*				NumberOfClasses = 3L; */
/*		results in the local optimum */
/*				QuantizedValues[] = {1.0F, 2.0F, 15.0F, 15.0F}; */
/*		instead of the true global optimum */
/*				QuantizedValues[] = {1.5F, 1.5F, 10.0F, 20.0F}; */
/* 2. there may be more than one global optimum */
/*		e.g.,	Frequency[] = {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0}; */
/*				Value[] = {-1.0F, 0.0F, 1.0F}; */
/*				NumberOfClasses = 2L; */
/*		results in the global optimum */
/*				QuantizedValues[] = {-0.5F, -0.5F, 1.0F}; */
/*		the other global optimum is ignored */
/*				QuantizedValues[] = {-1.0F, 0.5F, 0.5F}; */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin HistogramKMeans */

	float	*Threshold, *ClassRepresentative;
	double	Mean, Norm;
	float	Old;
	long	i, j, k;
	int		Modified, LessClasses;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(HistogramKMeans, Frequency, *Status,
/**/	"Missing input Frequency")
/**/DEBUG_CHECK_NULL_POINTER(HistogramKMeans, Value, *Status,
/**/	"Missing input Value")
/**/DEBUG_CHECK_NULL_POINTER(HistogramKMeans, QuantizedValue, *Status,
/**/	"Missing output QuantizedValue")
/**/DEBUG_CHECK_NULL_POINTER(HistogramKMeans, NumberOfClasses, *Status,
/**/	"Missing in-place NumberOfClasses")
/**/DEBUG_CHECK_RANGE_LONG(HistogramKMeans, HistogramLength, 1L, LONG_MAX, *Status,
/**/	"Invalid HistogramLength (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(HistogramKMeans, *NumberOfClasses, 1L, LONG_MAX, *Status,
/**/	"Invalid NumberOfClasses (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_DOUBLE(HistogramKMeans, Tolerance, 0.0, DBL_MAX, *Status,
/**/	"Invalid Tolerance (should be positive)")
/**/DEBUG_RETURN_ON_ERROR(HistogramKMeans, *Status)
/**/DEBUG_WRITE_ENTERING(HistogramKMeans,
/**/	"About to execute HistogramKMeans")

	if (HistogramLength < *NumberOfClasses) {
		WRITE_WARNING(HistogramKMeans, "Reducing the number of classes")
		*NumberOfClasses = HistogramLength;
	}
	Norm = 0.0;
	for (i = 0L; (i < HistogramLength); i++) {
		if (Frequency[i] <= 0.0) {
			*Status = ERROR;
			WRITE_ERROR(HistogramKMeans, "Invalid Frequency[i] (should be strictly positive)")
/**/		DEBUG_WRITE_LEAVING(HistogramKMeans, "Done")
			return(*Status);
		}
		Norm += Frequency[i];
	}
	if (Tolerance < ((Norm - 1.0) * (Norm - 1.0))) {
		*Status = ERROR;
		WRITE_ERROR(HistogramKMeans, "Invalid Frequency[] (should have unit sum)")
/**/	DEBUG_WRITE_LEAVING(HistogramKMeans, "Done")
		return(*Status);
	}
	for (i = 1L; (i < HistogramLength); i++) {
		if (Value[i] <= Value[i - 1L]) {
			*Status = ERROR;
			WRITE_ERROR(HistogramKMeans, "Invalid Value[] (should be strictly increasing)")
/**/		DEBUG_WRITE_LEAVING(HistogramKMeans, "Done")
			return(*Status);
		}
	}
	if (*NumberOfClasses == 1L) {
		Mean = 0.0;
		for (k = 0L; (k < HistogramLength); k++) {
			Mean += Frequency[k] * (double)Value[k];
		}
		for (k = 0L; (k < HistogramLength); k++) {
			QuantizedValue[k] = (float)Mean;
		}
/**/	DEBUG_WRITE_LEAVING(HistogramKMeans, "Done")
		return(*Status);
	}
	Threshold = (float *)malloc((size_t)((*NumberOfClasses - 1L) * (long)sizeof(float)));
	if (Threshold == (float *)NULL) {
		*Status = ERROR;
		WRITE_ERROR(HistogramKMeans, "Unable to allocate Threshold")
/**/	DEBUG_WRITE_LEAVING(HistogramKMeans, "Done")
		return(*Status);
	}
	ClassRepresentative = (float *)malloc((size_t)(*NumberOfClasses * (long)sizeof(float)));
	if (ClassRepresentative == (float *)NULL) {
		*Status = ERROR;
		free(Threshold);
		WRITE_ERROR(HistogramKMeans, "Unable to allocate ClassRepresentative")
/**/	DEBUG_WRITE_LEAVING(HistogramKMeans, "Done")
		return(*Status);
	}
	do {
		LessClasses = FALSE;
		if (*NumberOfClasses == 1L) {
			Mean = 0.0;
			for (k = 0L; (k < HistogramLength); k++) {
				Mean += Frequency[k] * (double)Value[k];
			}
			for (k = 0L; (k < HistogramLength); k++) {
				QuantizedValue[k] = (float)Mean;
			}
			free(ClassRepresentative);
			free(Threshold);
/**/		DEBUG_WRITE_LEAVING(HistogramKMeans, "Done")
			return(*Status);
		}
		j = (HistogramLength - 1L) / (*NumberOfClasses - 1L);
		for (i = 0L; (i < *NumberOfClasses); i++) {
			ClassRepresentative[i] = Value[i * j];
		}
		for (i = 0L; (i < (*NumberOfClasses - 1L)); i++) {
			Threshold[i] = (float)(0.5 * ((double)ClassRepresentative[i]
				+ (double)ClassRepresentative[i + 1L]));
		}
		do {
			Modified = FALSE;
			k = 0L;
			for (i = 0L; (i < (*NumberOfClasses - 1L)); i++) {
				Mean = 0.0;
				Norm = 0.0;
				while (Value[k] < Threshold[i]) {
					Mean += Frequency[k] * (double)Value[k];
					Norm += Frequency[k++];
				}
				if (Norm < Tolerance) {
					WRITE_WARNING(HistogramKMeans, "Reducing the number of classes")
					(*NumberOfClasses)--;
					LessClasses = TRUE;
					break;
				}
				Old = ClassRepresentative[i];
				ClassRepresentative[i] = (float)(Mean / Norm);
				Modified = (Modified || (Tolerance < ((Old
					- ClassRepresentative[i]) * (Old - ClassRepresentative[i]))));
			}
			if (LessClasses) {
				break;
			}
			Mean = 0.0;
			Norm = 0.0;
			while (k < HistogramLength) {
				Mean += Frequency[k] * (double)Value[k];
				Norm += Frequency[k++];
			}
			if (Norm < Tolerance) {
				WRITE_WARNING(HistogramKMeans, "Reducing the number of classes")
				(*NumberOfClasses)--;
				LessClasses = TRUE;
				break;
			}
			Old = ClassRepresentative[*NumberOfClasses - 1L];
			ClassRepresentative[*NumberOfClasses - 1L] = (float)(Mean / Norm);
			Modified = (Modified || (Tolerance < ((Old
				- ClassRepresentative[*NumberOfClasses - 1L]) * (Old
				- ClassRepresentative[*NumberOfClasses - 1L]))));
			for (i = 0L; (i < (*NumberOfClasses - 1L)); i++) {
				Old = Threshold[i];
				Threshold[i] = (float)(0.5 * ((double)ClassRepresentative[i]
					+ (double)ClassRepresentative[i + 1L]));
				Modified = (Modified || (Tolerance < ((Old - Threshold[i])
					* (Old - Threshold[i]))));
			}
		} while (Modified);
		if (LessClasses) {
			continue;
		}
		k = 0L;
		for (i = 0L; (i < (*NumberOfClasses - 1L)); i++) {
			while (Value[k] < Threshold[i]) {
				QuantizedValue[k++] = ClassRepresentative[i];
			}
		}
		while (k < HistogramLength) {
			QuantizedValue[k++] = ClassRepresentative[*NumberOfClasses - 1L];
		}
	} while (LessClasses);
	free(ClassRepresentative);
	free(Threshold);
/**/DEBUG_WRITE_LEAVING(HistogramKMeans, "Done")
	return(*Status);
} /* end HistogramKMeans */

