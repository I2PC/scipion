/*****************************************************************************
 *	System includes
 ****************************************************************************/
#include	<float.h>
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
#include	"linearalgebra.h"
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
/*--------------------------------------------------------------------------*/
static double	Pythag
				(
					double	a,					/* 1st side */
					double	b					/* 2nd side */
				);

/*****************************************************************************
 *	Definition of static procedures
 ****************************************************************************/
/*--------------------------------------------------------------------------*/
static double	Pythag
				(
					double	a,					/* 1st side */
					double	b					/* 2nd side */
				)

/* compute sqrt(a^2+b^2) in a stable numerical form */

{ /* begin Pythag */

	double			absa, absb;

/**/DEBUG_WRITE_ENTERING(Pythag,
/**/	"About to apply the Pythagorean theorem")

	absa = fabs(a);
	absb = fabs(b);
	if (absb < absa) {
/**/	DEBUG_WRITE_LEAVING(Pythag, "Done")
		return(absa * sqrt(1.0 + absb * absb / (absa * absa)));
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Pythag, "Done")
		return((absb == 0.0) ? (0.0) : (absb * sqrt(1.0 + absa * absa / (absb * absb))));
	}
} /* end Pythag */

/*****************************************************************************
 *	Definition of extern procedures
 ****************************************************************************/
/*--------------------------------------------------------------------------*/
extern int		AllocateMatrix
				(
					double	**Matrix,			/* double output pointer */
					long	Lines,				/* height of the matrix */
					long	Columns,			/* width of the matrix */
					int		*Status				/* error management */
				)

/* allocate a (double)matrix of size (Lines x Columns) */
/* success: return(!ERROR); failure: return(ERROR) and set *Matrix to NULL */

{ /* begin AllocateMatrix */

	*Status = !ERROR;

/**/DEBUG_CHECK_RANGE_LONG(AllocateMatrix, Lines, 1L, LONG_MAX, *Status,
/**/	"Invalid Lines (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(AllocateMatrix, Columns, 1L, LONG_MAX, *Status,
/**/	"Invalid Columns (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(AllocateMatrix, *Status)
/**/DEBUG_WRITE_ENTERING(AllocateMatrix,
/**/	"About to allocate a matrix")
#ifdef DEBUG
/**/if (*Matrix != (double *)NULL) {
/**/	WRITE_WARNING(AllocateMatrix, "Matrix may have been previously allocated")
/**/}
#endif

	*Matrix = (double *)malloc((size_t)(Lines * Columns * (long)sizeof(double)));
	if (*Matrix == (double *)NULL) {
		*Status = ERROR;
		WRITE_ERROR(AllocateMatrix, "Unable to perform allocation")
	}
/**/DEBUG_WRITE_LEAVING(AllocateMatrix, "Done")
	return(*Status);
} /* end AllocateMatrix */

/*--------------------------------------------------------------------------*/
extern int		AllocateVector
				(
					double	*(Vector[]),		/* double output pointer */
					long	Lines,				/* height of the vector */
					int		*Status				/* error management */
				)

/* allocate a (double)vector of size (Lines x 1) */
/* success: return(!ERROR); failure: return(ERROR) and set *Vector to NULL */

{ /* begin AllocateVector */

	*Status = !ERROR;

/**/DEBUG_CHECK_RANGE_LONG(AllocateVector, Lines, 1L, LONG_MAX, *Status,
/**/	"Invalid Lines (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(AllocateVector, *Status)
/**/DEBUG_WRITE_ENTERING(AllocateVector,
/**/	"About to allocate a vector")
#ifdef DEBUG
/**/if (*Vector != (double *)NULL) {
/**/	WRITE_WARNING(AllocateVector, "Vector may have been previously allocated")
/**/}
#endif

	*Vector = (double *)malloc((size_t)(Lines * (long)sizeof(double)));
	if (*Vector == (double *)NULL) {
		*Status = ERROR;
		WRITE_ERROR(AllocateVector, "Unable to perform allocation")
	}
/**/DEBUG_WRITE_LEAVING(AllocateVector, "Done")
	return(*Status);
} /* end AllocateVector */

/*--------------------------------------------------------------------------*/
extern int		FreeMatrix
				(
					double	**Matrix			/* (double)matrix array */
				)

/* free a (double)matrix */
/* success: return(!ERROR) and set *Matrix to NULL; failure: return(ERROR); */

{ /* begin FreeMatrix */

	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(FreeMatrix, *Matrix, Status,
/**/	"Nothing to free")
/**/DEBUG_RETURN_ON_ERROR(FreeMatrix, Status)
/**/DEBUG_WRITE_ENTERING(FreeMatrix,
/**/	"About to free a matrix")

	free(*Matrix);
	*Matrix = (double *)NULL;
/**/DEBUG_WRITE_LEAVING(FreeMatrix, "Done")
	return(Status);
} /* end FreeMatrix */

/*--------------------------------------------------------------------------*/
extern int		FreeVector
				(
					double	*(Vector[])			/* (double)vector array */
				)

/* free a (double)vector */
/* success: return(!ERROR) and set *Vector to NULL; failure: return(ERROR); */

{ /* begin FreeVector */

	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(FreeVector, *Vector, Status,
/**/	"Nothing to free")
/**/DEBUG_RETURN_ON_ERROR(FreeVector, Status)
/**/DEBUG_WRITE_ENTERING(FreeVector,
/**/	"About to free a vector")

	free(*Vector);
	*Vector = (double *)NULL;
/**/DEBUG_WRITE_LEAVING(FreeVector, "Done")
	return(Status);
} /* end FreeVector */

/*--------------------------------------------------------------------------*/
extern int		FrobeniusNorm
				(
					double	*A,					/* input matrix */
					double	*Norm,				/* resulting scalar norm */
					long	Lines,				/* height of the matrix A */
					long	Columns				/* width of the matrix A */
				)

/* compute the Frobenius norm of the matrix A */
/* the size of A is (Lines x Columns) */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin FrobeniusNorm */

	double	*p, *q;
	double	Sum;
	long	i, j;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(FrobeniusNorm, A, Status,
/**/	"Empty input major matrix")
/**/DEBUG_CHECK_NULL_POINTER(FrobeniusNorm, Norm, Status,
/**/	"Empty output")
/**/DEBUG_CHECK_RANGE_LONG(FrobeniusNorm, Lines, 1L, LONG_MAX, Status,
/**/	"Invalid Lines (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(FrobeniusNorm, Columns, 1L, LONG_MAX, Status,
/**/	"Invalid Columns (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(FrobeniusNorm, Status)
/**/DEBUG_WRITE_ENTERING(FrobeniusNorm,
/**/	"About to compute the Frobenius norm of a matrix")

	Sum = 0.0;
	p = A;
	q = A;
	for (i = Lines; (0L < i); i--) {
		for (j = Columns; (0L < j); j--) {
			Sum += *p++ * *q++;
		}
	}
	*Norm = sqrt(Sum);
/**/DEBUG_WRITE_LEAVING(FrobeniusNorm, "Done")
	return(Status);
} /* end FrobeniusNorm */

/*--------------------------------------------------------------------------*/
extern int		GetIdentitySquareMatrix
				(
					double	*A,					/* output identity matrix */
					long	Size				/* height or width of the square matrix */
				)

/* fill a matrix with identity */
/* the size of the output matrix A is (Size x Size) */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin GetIdentitySquareMatrix */

	long	i, j;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(GetIdentitySquareMatrix, A, Status,
/**/	"Missing output")
/**/DEBUG_CHECK_RANGE_LONG(GetIdentitySquareMatrix, Size, 1L, LONG_MAX, Status,
/**/	"Invalid size (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(GetIdentitySquareMatrix, Status)
/**/DEBUG_WRITE_ENTERING(GetIdentitySquareMatrix,
/**/	"About to build an identity matrix")

	for (i = Size; (0L < i); i--) {
		for (j = Size; (i < j); j--) {
			*A++ = 0.0;
		}
		*A++ = 1.0;
		j--;
		while (0L < j--) {
			*A++ = 0.0;
		}
	}
/**/DEBUG_WRITE_LEAVING(GetIdentitySquareMatrix, "Done")
	return(Status);
} /* end GetIdentitySquareMatrix */

/*--------------------------------------------------------------------------*/
extern int		LinearLeastSquares
				(
					double	*A,					/* matrix */
					long	Lines,				/* height the matrix */
					long	Columns,			/* width of the matrix */
					double	B[],				/* right-hand side vector */
					double	X[],				/* resulting vector */
					double	Tolerance,			/* admissible relative error */
					int		*Status				/* error management */
				)

/* give the least-squares solution to (A * X = B) such that */
/* (A^T * A)^-1 * A^T * B = X is a vector of size (Columns), where A is */
/* a (Lines x Columns) matrix, and where B is a vector of size (Lines) */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin LinearLeastSquares */

	double	*Q, *R;
	double	*p, *q, *x, *y;
	double	s;
	long	i, j;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(LinearLeastSquares, A, *Status,
/**/	"No input matrix")
/**/DEBUG_CHECK_NULL_POINTER(LinearLeastSquares, B, *Status,
/**/	"No input vector")
/**/DEBUG_CHECK_NULL_POINTER(LinearLeastSquares, X, *Status,
/**/	"No output vector")
/**/DEBUG_CHECK_RANGE_LONG(LinearLeastSquares, Lines, 1L, LONG_MAX, *Status,
/**/	"Invalid size (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(LinearLeastSquares, Columns, 1L, LONG_MAX, *Status,
/**/	"Invalid size (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_DOUBLE(LinearLeastSquares, Tolerance, 0.0, DBL_MAX, *Status,
/**/	"Invalid Tolerance (should be positive)")
/**/DEBUG_RETURN_ON_ERROR(LinearLeastSquares, *Status)
/**/DEBUG_WRITE_ENTERING(LinearLeastSquares,
/**/	"About to solve a linear least-squares problem by QR decomposition")

	AllocateMatrix(&Q, Lines, Columns, Status);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(LinearLeastSquares, "Done")
		return(*Status);
	}
	AllocateMatrix(&R, Columns, Columns, Status);
	if (*Status == ERROR) {
		FreeMatrix(&Q);
/**/	DEBUG_WRITE_LEAVING(LinearLeastSquares, "Done")
		return(*Status);
	}
	q = (double *)memcpy(Q, A, (size_t)(Lines * Columns * (long)sizeof(double)));
	QRdecomposition(Q, R, Lines, Columns, Tolerance, Status);
	if (*Status == ERROR) {
		FreeMatrix(&R);
		FreeMatrix(&Q);
/**/	DEBUG_WRITE_LEAVING(LinearLeastSquares, "Done")
		return(*Status);
	}
	x = X;
	for (i = Columns; (0L < i); i--) {
		s = 0.0;
		p = q++;
		y = B;
		for (j = Lines; (0L < j); j--) {
			s += *p * *y++;
			p += (ptrdiff_t)Columns;
		}
		*x++ = s;
	}
	p = R + (ptrdiff_t)(Columns * Columns);
	for (i = Columns; (0L < i); i--) {
		s = *--p;
		p -= (ptrdiff_t)Columns;
		if ((s * s) < Tolerance) {
			*--x = 0.0;
		}
		else {
			*--x /= s;
		}
		q = R + (ptrdiff_t)((i - 2L) * Columns + i - 1L);
		y = X + (ptrdiff_t)(i - 2L);
		for (j = i - 2L; (0L <= j); j--) {
			*y-- -= *q * *x;
			q -= (ptrdiff_t)Columns;
		}
	}
	*Status = FreeMatrix(&R);
	if (*Status == ERROR) {
		FreeMatrix(&Q);
/**/	DEBUG_WRITE_LEAVING(LinearLeastSquares, "Done")
		return(*Status);
	}
	*Status = FreeMatrix(&Q);
/**/DEBUG_WRITE_LEAVING(LinearLeastSquares, "Done")
	return(*Status);
} /* end LinearLeastSquares */

/*--------------------------------------------------------------------------*/
extern int		MatrixAdd
				(
					double	*A,					/* left matrix operand */
					double	*B,					/* right matrix operand */
					double	*X,					/* resulting matrix */
					long	Lines,				/* height of the matrices */
					long	Columns				/* width of the matrices */
				)

/* perform the matrix addition X = A + B */
/* the size of A, B and X is (Lines x Columns) */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin MatrixAdd */

	long	i, j;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(MatrixAdd, A, Status,
/**/	"Empty input left operand")
/**/DEBUG_CHECK_NULL_POINTER(MatrixAdd, B, Status,
/**/	"Empty input right operand")
/**/DEBUG_CHECK_NULL_POINTER(MatrixAdd, X, Status,
/**/	"Missing result")
/**/DEBUG_CHECK_RANGE_LONG(MatrixAdd, Lines, 1L, LONG_MAX, Status,
/**/	"Invalid Lines (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MatrixAdd, Columns, 1L, LONG_MAX, Status,
/**/	"Invalid Columns (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(MatrixAdd, Status)
/**/DEBUG_WRITE_ENTERING(MatrixAdd,
/**/	"About to add two matrices")

	for (i = Lines; (0L < i); i--) {
		for (j = Columns; (0L < j); j--) {
			*X++ = *A++ + *B++;
		}
	}
/**/DEBUG_WRITE_LEAVING(MatrixAdd, "Done")
	return(Status);
} /* end MatrixAdd */

/*--------------------------------------------------------------------------*/
extern int		MatrixConditionNumber
				(
					double	*A,					/* input matrix */
					long	Lines,				/* height of the matrix A */
					long	Columns,			/* width of the matrix A */
					double	*ConditionNumber,	/* resulting scalar condition number */
					double	Tolerance,			/* admissible relative error */
					long	MaxIterations,		/* convergence limit */
					int		*Status				/* error management */
				)

/* compute the condition number of the matrix A */
/* the size of A is (Lines x Columns) */
/* ConditionNumber = -1.0 indicates a singular matrix */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin MatrixConditionNumber */

	double	*w;
	double	*U = (double *)NULL, *W = (double *)NULL, *V = (double *)NULL;
	double	Small, Large;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(MatrixConditionNumber, A, *Status,
/**/	"Empty input")
/**/DEBUG_CHECK_NULL_POINTER(MatrixConditionNumber, ConditionNumber, *Status,
/**/	"Empty output")
/**/DEBUG_CHECK_RANGE_LONG(MatrixConditionNumber, Lines, 1L, LONG_MAX, *Status,
/**/	"Invalid Lines (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MatrixConditionNumber, Columns, 1L, LONG_MAX, *Status,
/**/	"Invalid Columns (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MatrixConditionNumber, MaxIterations, 1L, LONG_MAX, *Status,
/**/	"Invalid iterations (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_DOUBLE(MatrixConditionNumber, Tolerance, 0.0, DBL_MAX, *Status,
/**/	"Invalid Tolerance (should be positive)")
/**/DEBUG_RETURN_ON_ERROR(MatrixConditionNumber, *Status)
/**/DEBUG_WRITE_ENTERING(MatrixConditionNumber,
/**/	"About to compute the condition number of a matrix")

	AllocateMatrix(&U, Lines, Columns, Status);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(MatrixConditionNumber, "Done")
		return(*Status);
	}
	AllocateMatrix(&V, Columns, Columns, Status);
	if (*Status == ERROR) {
		FreeMatrix(&U);
/**/	DEBUG_WRITE_LEAVING(MatrixConditionNumber, "Done")
		return(*Status);
	}
	AllocateVector(&W, Columns, Status);
	if (*Status == ERROR) {
		FreeMatrix(&V);
		FreeMatrix(&U);
/**/	DEBUG_WRITE_LEAVING(MatrixConditionNumber, "Done")
		return(*Status);
	}
	U = (double *)memcpy(U, A, (size_t)(Lines * Columns * (long)sizeof(double)));
	SingularValueDecomposition(U, Lines, Columns, W, V, MaxIterations, Status);
	if (*Status == ERROR) {
		FreeVector(&W);
		FreeMatrix(&V);
		FreeMatrix(&U);
/**/	DEBUG_WRITE_LEAVING(MatrixConditionNumber, "Done")
		return(*Status);
	}
	*Status = FreeMatrix(&V);
	if (*Status == ERROR) {
		FreeVector(&W);
		FreeMatrix(&U);
/**/	DEBUG_WRITE_LEAVING(MatrixConditionNumber, "Done")
		return(*Status);
	}
	*Status = FreeMatrix(&U);
	if (*Status == ERROR) {
		FreeVector(&W);
/**/	DEBUG_WRITE_LEAVING(MatrixConditionNumber, "Done")
		return(*Status);
	}
	w = W;
	Small = *w;
	Large = *w;
	while (0L < --Columns) {
		w++;
		Small = (*w < Small) ? (*w) : (Small);
		Large = (Large < *w) ? (*w) : (Large);
	}
	*ConditionNumber = ((Tolerance * Large) < Small) ? (Large / Small) : (-1.0);
	*Status = FreeVector(&W);
/**/DEBUG_WRITE_LEAVING(MatrixConditionNumber, "Done")
	return(*Status);
} /* end MatrixConditionNumber */

/*--------------------------------------------------------------------------*/
extern int		MatrixGramSchmidtOrthonormalize
				(
					double	*A,					/* input matrix */
					double	*B,					/* resulting matrix */
					long	Lines,				/* height of the matrix A */
					long	Columns,			/* width of the matrix A */
					double	Tolerance			/* admissible relative error */
				)

/* orthonormalize a matrix */
/* the size of A and B is (Lines x Columns) */
/* the orthogonalization proceeds row-wise */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin MatrixGramSchmidtOrthonormalize */

	double	*u, *v, *w, *x;
	double	Sum;
	long	i, j, k;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(MatrixGramSchmidtOrthonormalize, A, Status,
/**/	"Empty input")
/**/DEBUG_CHECK_NULL_POINTER(MatrixGramSchmidtOrthonormalize, B, Status,
/**/	"Empty output")
/**/DEBUG_CHECK_RANGE_LONG(MatrixGramSchmidtOrthonormalize, Lines, 1L, LONG_MAX, Status,
/**/	"Invalid Lines (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MatrixGramSchmidtOrthonormalize, Columns, 1L, LONG_MAX, Status,
/**/	"Invalid Columns (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_DOUBLE(MatrixGramSchmidtOrthonormalize, Tolerance, 0.0, DBL_MAX, Status,
/**/	"Invalid Tolerance (should be positive)")
/**/DEBUG_RETURN_ON_ERROR(MatrixGramSchmidtOrthonormalize, Status)
/**/DEBUG_WRITE_ENTERING(MatrixGramSchmidtOrthonormalize,
/**/	"About to othonormalize a matrix by the Gram-Schmidt procedure")

	u = B;
	i = 0L;
	do {
		v = A;
		Sum = 0.0;
		for (j = Columns; (0L < j); j--) {
			Sum += *v * *v;
			v++;
		}
		if (Tolerance < Sum) {
			Sum = 1.0 / sqrt(Sum);
			for (j = Columns; (0L < j); j--) {
				*u++ = *A++ * Sum;
			}
		}
		else {
			u = (double *)memcpy(u, A, (size_t)(Columns * (long)sizeof(double)));
			u += (ptrdiff_t)Columns;
			A += (ptrdiff_t)Columns;
		}
		i++;
	} while (Sum < Tolerance);
	while (i++ < Lines) {
		u = (double *)memcpy(u, A, (size_t)(Columns * (long)sizeof(double)));
		v = B;
		for (k = i - 1L; (0L < k); k--) {
			w = A;
			x = v;
			Sum = 0.0;
			for (j = Columns; (0L < j); j--) {
				Sum += *w++ * *x++;
			}
			w = u;
			for (j = Columns; (0L < j); j--) {
				*w++ -= *v++ * Sum;
			}
		}
		Sum = 0.0;
		for (j = Columns; (0L < j); j--) {
			Sum += *v * *v;
			v++;
		}
		if (Tolerance < Sum) {
			Sum = 1.0 / sqrt(Sum);
			for (j = Columns; (0L < j); j--) {
				*u++ *= Sum;
			}
		}
		else {
			u += (ptrdiff_t)Columns;
		}
		A += (ptrdiff_t)Columns;
	}
/**/DEBUG_WRITE_LEAVING(MatrixGramSchmidtOrthonormalize, "Done")
	return(Status);
} /* end MatrixGramSchmidtOrthonormalize */

/*--------------------------------------------------------------------------*/
extern int		MatrixMinor
				(
					double	*A,					/* input matrix */
					double	*B,					/* smaller output matrix */
					long	Lines,				/* height of the matrix A */
					long	Columns,			/* width of the matrix A */
					long	i,					/* line to delete */
					long	j					/* column to delete */
				)

/* extract a minor matrix B from the major matrix A */
/* the size of A is (Lines x Columns) */
/* the size of B is ((Lines - 1) x (Columns - 1)) */
/* the line to delete is i (indexing starts from 0) */
/* the column to delete is j (indexing starts from 0) */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin MatrixMinor */

	long	c, l;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(MatrixMinor, A, Status,
/**/	"Empty input major matrix")
/**/DEBUG_CHECK_NULL_POINTER(MatrixMinor, B, Status,
/**/	"Empty output minor matrix")
/**/DEBUG_CHECK_RANGE_LONG(MatrixMinor, Lines, 1L, LONG_MAX, Status,
/**/	"Invalid Lines (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MatrixMinor, Columns, 1L, LONG_MAX, Status,
/**/	"Invalid Columns (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MatrixMinor, i, 0L, Lines - 1L, Status,
/**/	"Invalid line to delete (should be in [0, Lines - 1])")
/**/DEBUG_CHECK_RANGE_LONG(MatrixMinor, j, 0L, Columns - 1L, Status,
/**/	"Invalid column to delete (should be in [0, Columns - 1])")
/**/DEBUG_RETURN_ON_ERROR(MatrixMinor, Status)
/**/DEBUG_WRITE_ENTERING(MatrixMinor,
/**/	"About to extract a minor matrix from a major matrix")

	for (l = 0L; (l < i); l++) {
		for (c = 0L; (c < j); c++) {
			*B++ = *A++;
		}
		A++;
		while (c++ < Columns) {
			*B++ = *A++;
		}
	}
	A += (ptrdiff_t)Columns;
	while (l++ < Lines) {
		for (c = 0L; (c < j); c++) {
			*B++ = *A++;
		}
		A++;
		while (c++ < Columns) {
			*B++ = *A++;
		}
	}
/**/DEBUG_WRITE_LEAVING(MatrixMinor, "Done")
	return(Status);
} /* end MatrixMinor */

/*--------------------------------------------------------------------------*/
extern int		MatrixMultiply
				(
					double	*A,					/* left matrix operand */
					double	*B,					/* right maxtrix operand */
					double	*X,					/* resulting matrix */
					long	Lines,				/* height (left and resulting) */
					long	CommonSize,			/* left width and right height */
					long	Columns				/* width (right and resulting) */
				)

/* perform the matrix multiplication X = A.B */
/* the size of A is (Lines x CommonSize) */
/* the size of B is (CommonSize x Columns) */
/* the size of X is (Lines x Columns) */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin MatrixMultiply */

	double	*a, *b, *p;
	double	Sum;
	long	i, j, k;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(MatrixMultiply, A, Status,
/**/	"Empty input left operand")
/**/DEBUG_CHECK_NULL_POINTER(MatrixMultiply, B, Status,
/**/	"Empty input right operand")
/**/DEBUG_CHECK_NULL_POINTER(MatrixMultiply, X, Status,
/**/	"Missing result")
/**/DEBUG_CHECK_RANGE_LONG(MatrixMultiply, Lines, 1L, LONG_MAX, Status,
/**/	"Invalid Lines (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MatrixMultiply, CommonSize, 1L, LONG_MAX, Status,
/**/	"Invalid CommonSize (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MatrixMultiply, Columns, 1L, LONG_MAX, Status,
/**/	"Invalid Columns (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(MatrixMultiply, Status)
/**/DEBUG_WRITE_ENTERING(MatrixMultiply,
/**/	"About to multiply two matrices together")

	for (i = Lines; (0L < i); i--) {
		p = B;
		for (j = Columns; (0L < j); j--) {
			Sum = 0.0;
			a = A;
			b = p++;
			for (k = CommonSize; (0L < k); k--) {
				Sum += *a++ * *b;
				b += (ptrdiff_t)Columns;
			}
			*X++ = Sum;
		}
		A += (ptrdiff_t)CommonSize;
	}
/**/DEBUG_WRITE_LEAVING(MatrixMultiply, "Done")
	return(Status);
} /* end MatrixMultiply */

/*--------------------------------------------------------------------------*/
extern int		MatrixRank
				(
					double	*A,					/* input matrix */
					long	Lines,				/* height of the matrix A */
					long	Columns,			/* width of the matrix A */
					long	*Rank,				/* resulting scalar condition number */
					double	Tolerance,			/* admissible relative error */
					long	MaxIterations,		/* convergence limit */
					int		*Status				/* error management */
				)

/* compute the rank of the matrix A */
/* the size of A is (Lines x Columns) */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin MatrixRank */

	double	*w;
	double	*U = (double *)NULL, *W = (double *)NULL, *V = (double *)NULL;
	double	Norm;
	long	j;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(MatrixRank, A, *Status,
/**/	"Empty input")
/**/DEBUG_CHECK_NULL_POINTER(MatrixRank, Rank, *Status,
/**/	"Empty output")
/**/DEBUG_CHECK_RANGE_LONG(MatrixRank, Lines, 1L, LONG_MAX, *Status,
/**/	"Invalid Lines (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MatrixRank, Columns, 1L, LONG_MAX, *Status,
/**/	"Invalid Columns (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MatrixRank, MaxIterations, 1L, LONG_MAX, *Status,
/**/	"Invalid iterations (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_DOUBLE(MatrixRank, Tolerance, 0.0, DBL_MAX, *Status,
/**/	"Invalid Tolerance (should be positive)")
/**/DEBUG_RETURN_ON_ERROR(MatrixRank, *Status)
/**/DEBUG_WRITE_ENTERING(MatrixRank,
/**/	"About to compute the rank of a matrix")

	AllocateMatrix(&U, Lines, Columns, Status);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(MatrixRank, "Done")
		return(*Status);
	}
	AllocateMatrix(&V, Columns, Columns, Status);
	if (*Status == ERROR) {
		FreeMatrix(&U);
/**/	DEBUG_WRITE_LEAVING(MatrixRank, "Done")
		return(*Status);
	}
	AllocateVector(&W, Columns, Status);
	if (*Status == ERROR) {
		FreeMatrix(&V);
		FreeMatrix(&U);
/**/	DEBUG_WRITE_LEAVING(MatrixRank, "Done")
		return(*Status);
	}
	U = (double *)memcpy(U, A, (size_t)(Lines * Columns * (long)sizeof(double)));
	SingularValueDecomposition(U, Lines, Columns, W, V, MaxIterations, Status);
	if (*Status == ERROR) {
		FreeVector(&W);
		FreeMatrix(&V);
		FreeMatrix(&U);
/**/	DEBUG_WRITE_LEAVING(MatrixRank, "Done")
		return(*Status);
	}
	*Status = FreeMatrix(&V);
	if (*Status == ERROR) {
		FreeVector(&W);
		FreeMatrix(&U);
/**/	DEBUG_WRITE_LEAVING(MatrixRank, "Done")
		return(*Status);
	}
	*Status = FreeMatrix(&U);
	if (*Status == ERROR) {
		FreeVector(&W);
/**/	DEBUG_WRITE_LEAVING(MatrixRank, "Done")
		return(*Status);
	}
	w = W;
	Norm = 0.0;
	for (j = Columns; (0L < j); j--) {
		Norm += *w * *w;
		w++;
	}
	Norm = sqrt(Norm);
	*Rank = 0L;
	if (Tolerance < Norm) {
		w = W;
		Norm = 1.0 / Norm;
		for (j = Columns; (0L < j); j--) {
			*Rank += ((*w++ * Norm) < Tolerance) ? (0L) : (1L);
		}
	}
	*Status = FreeVector(&W);
/**/DEBUG_WRITE_LEAVING(MatrixRank, "Done")
	return(*Status);
} /* end MatrixRank */

/*--------------------------------------------------------------------------*/
extern int		MatrixSubtract
				(
					double	*A,					/* left matrix operand */
					double	*B,					/* right matrix operand */
					double	*X,					/* resulting matrix */
					long	Lines,				/* height of the matrices */
					long	Columns				/* width of the matrices */
				)

/* perform the matrix subtraction X = A - B */
/* the size of A, B and X is (Lines x Columns) */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin MatrixSubtract */

	long	i, j;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(MatrixSubtract, A, Status,
/**/	"Empty input left operand")
/**/DEBUG_CHECK_NULL_POINTER(MatrixSubtract, B, Status,
/**/	"Empty input right operand")
/**/DEBUG_CHECK_NULL_POINTER(MatrixSubtract, X, Status,
/**/	"Missing result")
/**/DEBUG_CHECK_RANGE_LONG(MatrixSubtract, Lines, 1L, LONG_MAX, Status,
/**/	"Invalid Lines (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MatrixSubtract, Columns, 1L, LONG_MAX, Status,
/**/	"Invalid Columns (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(MatrixSubtract, Status)
/**/DEBUG_WRITE_ENTERING(MatrixSubtract,
/**/	"About to subtract two matrices")

	for (i = Lines; (0L < i); i--) {
		for (j = Columns; (0L < j); j--) {
			*X++ = *A++ - *B++;
		}
	}
/**/DEBUG_WRITE_LEAVING(MatrixSubtract, "Done")
	return(Status);
} /* end MatrixSubtract */

/*--------------------------------------------------------------------------*/
extern int		MatrixTimesVector
				(
					double	*A,					/* left matrix operand */
					double	*B,					/* right vector operand */
					double	*X,					/* resulting vector */
					long	Lines,				/* matrix height */
					long	Columns				/* matrix width and vector height */
				)

/* perform the matrix/vector multiplication X = A.B */
/* the size of A is (Lines x Columns) */
/* the size of B is (Columns x 1) */
/* the size of X is (Lines x 1) */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin MatrixTimesVector */

	double	*p;
	double	Sum;
	long	i, j;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(MatrixTimesVector, A, Status,
/**/	"Empty input left operand")
/**/DEBUG_CHECK_NULL_POINTER(MatrixTimesVector, B, Status,
/**/	"Empty input right operand")
/**/DEBUG_CHECK_NULL_POINTER(MatrixTimesVector, X, Status,
/**/	"Missing result")
/**/DEBUG_CHECK_RANGE_LONG(MatrixTimesVector, Lines, 1L, LONG_MAX, Status,
/**/	"Invalid Lines (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MatrixTimesVector, Columns, 1L, LONG_MAX, Status,
/**/	"Invalid Columns (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(MatrixTimesVector, Status)
/**/DEBUG_WRITE_ENTERING(MatrixTimesVector,
/**/	"About to multiply two matrices together")

	for (i = Lines; (0L < i); i--) {
		Sum = 0.0;
		p = B;
		for (j = Columns; (0L < j); j--) {
			Sum += *A++ * *p++;
		}
		*X++ = Sum;
	}
/**/DEBUG_WRITE_LEAVING(MatrixTimesVector, "Done")
	return(Status);
} /* end MatrixTimesVector */

/*--------------------------------------------------------------------------*/
extern int		MatrixTranspose
				(
					double	*A,					/* input matrix */
					double	*At,				/* transposed output matrix */
					long	Lines,				/* height input, width output */
					long	Columns				/* width input, height output */
				)

/* perform the transposition of a matrix */
/* the size of the input matrix A is (Lines x Columns) */
/* the size of the output matrix At is (Columns x Lines) */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin MatrixTranspose */

	double	*p;
	long	i, j;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(MatrixTranspose, A, Status,
/**/	"No input")
/**/DEBUG_CHECK_NULL_POINTER(MatrixTranspose, At, Status,
/**/	"Missing output")
/**/DEBUG_CHECK_RANGE_LONG(MatrixTranspose, Lines, 1L, LONG_MAX, Status,
/**/	"Invalid Lines (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(MatrixTranspose, Columns, 1L, LONG_MAX, Status,
/**/	"Invalid Columns (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(MatrixTranspose, Status)
/**/DEBUG_WRITE_ENTERING(MatrixTranspose,
/**/	"About to transpose a matrix")

	for (i = Lines; (0L < i); i--) {
		p = At++;
		for (j = Columns; (0L < j); j--) {
			*p = *A++;
			p += (ptrdiff_t)Lines;
		}
	}
/**/DEBUG_WRITE_LEAVING(MatrixTranspose, "Done")
	return(Status);
} /* end MatrixTranspose */

/*--------------------------------------------------------------------------*/
int		multiply_3Matrices(	double	*A,			 /* left matrix operand */
							double	*B,			 /* middle maxtrix operand */
							double	*C,			 /* right maxtrix operand */
							double	*X,			 /* resulting matrix */
							long	Lines,		 /* height (left and resulting) */
							long	CommonSizeH, /* left width and middle hight */
							long	CommonSizeW, /* right hight and middle width */
							long	Columns		 /* width (right and resulting) */ ) {

double	*Help;

Help = (double *)malloc((size_t)(Columns * CommonSizeH) * sizeof(double));
if (Help == (double *)NULL){ 
	WRITE_ERROR( multiply_3Matrices, "ERROR - Not enough memory for Help");
	return(ERROR);
}


if (MatrixMultiply(B, C, Help, CommonSizeH, CommonSizeW, Columns) == ERROR){ 
	WRITE_ERROR( multiply_3Matrices, "Error returned by MatrixMultiply");
	return(ERROR);
}
if (MatrixMultiply(A, Help, X, Lines, CommonSizeH, Columns) == ERROR){
	WRITE_ERROR( multiply_3Matrices, "Error returned by MatrixMultiply"); 
	return(ERROR);
}

free(Help);

return(!ERROR);
}

/*--------------------------------------------------------------------------*/
extern int		QRdecomposition
				(
					double	*Q,					/* in-place matrix */
					double	*R,					/* output matrix */
					long	Lines,				/* height the matrix */
					long	Columns,			/* width of the matrix */
					double	Tolerance,			/* admissible relative error */
					int		*Status				/* error management */
				)

/* decompose the (Lines x Columns) input matrix Q into an orthonormal */
/* output matrix Q of same size (Lines x Columns) and an upper-diagonal */
/* square matrix R of size (Columns x Columns), such that the matrix */
/* product (Q * R) gives the input matrix, and such that the matrix */
/* product (Q^T * Q) gives the identity. */
/* requirement: Columns <= Lines */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin QRdecomposition */

	double	*A;
	double	*a, *b, *c, *p, *q, *u, *v;
	double	s;
	long	i, j, k;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(QRdecomposition, Q, *Status,
/**/	"No input matrix Q")
/**/DEBUG_CHECK_NULL_POINTER(QRdecomposition, R, *Status,
/**/	"No output matrix R")
/**/DEBUG_CHECK_RANGE_LONG(QRdecomposition, Lines, 1L, LONG_MAX, *Status,
/**/	"Invalid size (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(QRdecomposition, Columns, 1L, Lines, *Status,
/**/	"Invalid size (should be [0 < Columns <= Lines])")
/**/DEBUG_CHECK_RANGE_DOUBLE(QRdecomposition, Tolerance, 0.0, DBL_MAX, *Status,
/**/	"Invalid Tolerance (should be positive)")
/**/DEBUG_RETURN_ON_ERROR(QRdecomposition, *Status)
/**/DEBUG_WRITE_ENTERING(QRdecomposition,
/**/	"About to a QR decomposition")

	AllocateMatrix(&A, Lines, Columns, Status);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(QRdecomposition, "Done")
		return(*Status);
	}
	a = (double *)memcpy(A, Q, (size_t)(Lines * Columns * (long)sizeof(double)));
	u = Q;
	for (j = 0L; (j < Columns); a++, u++, j++) {
		p = Q;
		for (k = j; (0 < k); k--) {
			s = 0.0;
			b = a;
			q = p++;
			for (i = Lines; (0 < i); i--) {
				s += *b * *q;
				b += (ptrdiff_t)Columns;
				q += (ptrdiff_t)Columns;
			}
			v = u + (ptrdiff_t)(Lines * Columns);
			for (i = Lines; (0 < i); i--) {
				q -= (ptrdiff_t)Columns;
				v -= (ptrdiff_t)Columns;
				*v -= s * *q;
			}
		}
		s = 0.0;
		v = u;
		for (i = Lines; (0 < i); i--) {
			s += *v * *v;
			v += (ptrdiff_t)Columns;
		}
		if ((s * s) < Tolerance) {
			s = 0.0;
		}
		else {
			s = 1.0 / sqrt((double)s);
		}
		for (i = Lines; (0 < i); i--) {
			v -= (ptrdiff_t)Columns;
			*v *= s;
		}
	}
	a = A;
	for (i = 0L; (i < Columns); Q++, i++) {
		for (j = 0; (j < i); j++) {
			*R++ = 0.0;
		}
		b = a++;
		for (j = i; (j < Columns); j++) {
			s = 0.0;
			c = b++;
			q = Q;
			for (k = 0L; (k < Lines); k++) {
				s += *q * *c;
				c += (ptrdiff_t)Columns;
				q += (ptrdiff_t)Columns;
			}
			*R++ = s;
		}
	}
	*Status = FreeMatrix(&A);
/**/DEBUG_WRITE_LEAVING(QRdecomposition, "Done")
	return(*Status);
} /* end QRdecomposition */

/*--------------------------------------------------------------------------*/
extern int		SingularValueBackSubstitution
				(
					double	*U,					/* input matrix */
					double	W[],				/* vector of singular values */
					double	*V,					/* untransposed orthogonal matrix */
					long	Lines,				/* height the matrix */
					long	Columns,			/* width of the matrix */
					double	B[],				/* input vector */
					double	X[],				/* returned solution */
					int		*Status				/* error management */
				)

/* solve (U.W.Transpose(V)).X == B in terms of X */
/* {U, W, V} are given by SingularValueDecomposition */
/* by convention, set w[i,j]=0 to get (1/w[i,j])=0 */
/* the size of the input matrix U is (Lines x Columns) */
/* the size of the vector (1/W) of singular values is (Columns) */
/* the size of the untransposed orthogonal matrix V is (Columns x Columns) */
/* the size of the input vector B is (Columns) */
/* the size of the output vector X is (Columns) */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin SingularValueBackSubstitution */

	double			*Workspace = (double *)NULL;
	double			*p, *q, *r;
	double			s;
	long			i, j;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(SingularValueBackSubstitution, U, *Status,
/**/	"No input matrix")
/**/DEBUG_CHECK_NULL_POINTER(SingularValueBackSubstitution, W, *Status,
/**/	"No input singular values")
/**/DEBUG_CHECK_NULL_POINTER(SingularValueBackSubstitution, V, *Status,
/**/	"No input untransposed orthogonal matrix")
/**/DEBUG_CHECK_NULL_POINTER(SingularValueBackSubstitution, B, *Status,
/**/	"No input vector")
/**/DEBUG_CHECK_NULL_POINTER(SingularValueBackSubstitution, X, *Status,
/**/	"No output vector")
/**/DEBUG_CHECK_RANGE_LONG(SingularValueBackSubstitution, Lines, 1L, LONG_MAX, *Status,
/**/	"Invalid size (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(SingularValueBackSubstitution, Columns, 1L, LONG_MAX, *Status,
/**/	"Invalid size (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(SingularValueBackSubstitution, *Status)
/**/DEBUG_WRITE_ENTERING(SingularValueBackSubstitution,
/**/	"About to perform SVD back substitution")

	AllocateVector(&Workspace, Columns, Status);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(SingularValueBackSubstitution, "Done")
		return(*Status);
	}
	q = Workspace;
	for (j = Columns; (0L < j); W++, j--) {
		s = 0.0;
		r = U++;
		if (*W != 0.0) {
			p = B;
			for (i = Lines; (0L < i); i--) {
				s += *r * *p++;
				r += (ptrdiff_t)Columns;
			}
			s /= *W;
		}
		*q++ = s;
	}
	for (i = Columns; (0L < i); i--) {
		s = 0.0;
		q = Workspace;
		for (j = Columns; (0L < j); j--) {
			s += *V++ * *q++;
		}
		*X++ = s;
	}
	*Status = FreeVector(&Workspace);
/**/DEBUG_WRITE_LEAVING(SingularValueBackSubstitution, "Done")
	return(*Status);
} /* end SingularValueBackSubstitution */

/*--------------------------------------------------------------------------*/
extern int		SingularValueDecomposition
				(
					double	*U,					/* input matrix (in-place processing) */
					long	Lines,				/* height the matrix */
					long	Columns,			/* width of the matrix */
					double	W[],				/* resulting vector of singular values */
					double	*V,					/* untransposed orthogonal matrix */
					long	MaxIterations,		/* convergence limit */
					int		*Status				/* error management */
				)

/* decomposition of the input matrix into U.W.Transpose(V) */
/* the size of the input matrix U is (Lines x Columns) */
/* the singular values are returned in the vector W of size (Columns) */
/* the size of the returned untransposed orthogonal matrix V is (Columns x Columns) */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin SingularValueDecomposition */

	double	*rv1 = (double *)NULL;
	double	Norm, Scale;
	double	c, f, g, h, s;
	double	x, y, z;
	long	i, its, j, jj, k, l = 0L, nm = 0L;
	int		Flag;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(SingularValueDecomposition, U, *Status,
/**/	"No input matrix")
/**/DEBUG_CHECK_NULL_POINTER(SingularValueDecomposition, W, *Status,
/**/	"No output singular values")
/**/DEBUG_CHECK_NULL_POINTER(SingularValueDecomposition, V, *Status,
/**/	"No output matrix")
/**/DEBUG_CHECK_RANGE_LONG(SingularValueDecomposition, Lines, 1L, LONG_MAX, *Status,
/**/	"Invalid size (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(SingularValueDecomposition, Columns, 1L, LONG_MAX, *Status,
/**/	"Invalid size (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(SingularValueDecomposition, MaxIterations, 1L, LONG_MAX, *Status,
/**/	"Invalid iterations (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(SingularValueDecomposition, *Status)
/**/DEBUG_WRITE_ENTERING(SingularValueDecomposition,
/**/	"About to decompose a matrix into singular values")

	AllocateMatrix(&rv1, Columns, Columns, Status);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(SingularValueDecomposition, "Done")
		return(*Status);
	}
	g = Scale = Norm = 0.0;
	for (i = 0L; (i < Columns); i++) {
		l = i + 1L;
		rv1[i] = Scale * g;
		g = s = Scale = 0.0;
		if (i < Lines) {
			for (k = i; (k < Lines); k++) {
				Scale += fabs(U[k * Columns + i]);
			}
			if (Scale != 0.0) {
				for (k = i; (k < Lines); k++) {
					U[k * Columns + i] /= Scale;
					s += U[k * Columns + i] * U[k * Columns + i];
				}
				f = U[i * Columns + i];
				g = (0.0 <= f) ? (-sqrt(s)) : (sqrt(s));
				h = f * g - s;
				U[i * Columns + i] = f - g;
				for (j = l; (j < Columns); j++) {
					for (s = 0.0, k = i; (k < Lines); k++) {
						s += U[k * Columns + i] * U[k * Columns + j];
					}
					f = s / h;
					for (k = i; (k < Lines); k++) {
						U[k * Columns + j] += f * U[k * Columns + i];
					}
				}
				for (k = i; (k < Lines); k++) {
					U[k * Columns + i] *= Scale;
				}
			}
		}
		W[i] = Scale * g;
		g = s = Scale = 0.0;
		if ((i < Lines) && (i != (Columns - 1L))) {
			for (k = l; (k < Columns); k++) {
				Scale += fabs(U[i * Columns + k]);
			}
			if (Scale != 0.0) {
				for (k = l; (k < Columns); k++) {
					U[i * Columns + k] /= Scale;
					s += U[i * Columns + k] * U[i * Columns + k];
				}
				f = U[i * Columns + l];
				g = (0.0 <= f) ? (-sqrt(s)) : (sqrt(s));
				h = f * g - s;
				U[i * Columns + l] = f - g;
				for (k = l; (k < Columns); k++) {
					rv1[k] = U[i * Columns + k] / h;
				}
				for (j = l; (j < Lines); j++) {
					for (s = 0.0, k = l; (k < Columns); k++) {
						s += U[j * Columns + k] * U[i * Columns + k];
					}
					for (k = l; (k < Columns); k++) {
						U[j * Columns + k] += s * rv1[k];
					}
				}
				for (k = l; (k < Columns); k++) {
					U[i * Columns + k] *= Scale;
				}
			}
		}
		Norm = ((fabs(W[i]) + fabs(rv1[i])) < Norm) ? (Norm) : (fabs(W[i]) + fabs(rv1[i]));
	}
	for (i = Columns - 1L; (0L <= i); i--) {
		if (i < (Columns - 1L)) {
			if (g != 0.0) {
				for (j = l; (j < Columns); j++) {
					V[j * Columns + i] = U[i * Columns + j] / (U[i * Columns + l] * g);
				}
				for (j = l; (j < Columns); j++) {
					for (s = 0.0, k = l; (k < Columns); k++) {
						s += U[i * Columns + k] * V[k * Columns + j];
					}
					for (k = l; (k < Columns); k++) {
						if (s != 0.0) {
							V[k * Columns + j] += s * V[k * Columns + i];
						}
					}
				}
			}
			for (j = l; (j < Columns); j++) {
				V[i * Columns + j] = V[j * Columns + i] = 0.0;
			}
		}
		V[i * Columns + i] = 1.0;
		g = rv1[i];
		l = i;
	}
	for (i = (Lines < Columns) ? (Lines - 1L) : (Columns - 1L); (0L <= i); i--) {
		l = i + 1L;
		g = W[i];
		for (j = l; (j < Columns); j++) {
			U[i * Columns + j] = 0.0;
		}
		if (g != 0.0) {
			g = 1.0 / g;
			for (j = l; (j < Columns); j++) {
				for (s = 0.0, k = l; (k < Lines); k++) {
					s += U[k * Columns + i] * U[k * Columns + j];
				}
				f = s * g / U[i * Columns + i];
				for (k = i; (k < Lines); k++) {
					if (f != 0.0) {
						U[k * Columns + j] += f * U[k * Columns + i];
					}
				}
			}
			for (j = i; (j < Lines); j++) {
				U[j * Columns + i] *= g;
			}
		}
		else {
			for (j = i; (j < Lines); j++) {
				U[j * Columns + i] = 0.0;
			}
		}
		U[i * Columns + i] += 1.0;
	}
	for (k = Columns - 1L; (0L <= k); k--) {
		for (its = 1L; (its <= MaxIterations); its++) {
			Flag = TRUE;
			for (l = k; (0L <= l); l--) {
				nm = l - 1L;
				if ((fabs(rv1[l]) + Norm) == Norm) {
					Flag = FALSE;
					break;
				}
				if ((fabs(W[nm]) + Norm) == Norm) {
					break;
				}
			}
			if (Flag) {
				c = 0.0;
				s = 1.0;
				for (i = l; (i <= k); i++) {
					f = s * rv1[i];
					rv1[i] *= c;
					if ((fabs(f) + Norm) == Norm) {
						break;
					}
					g = W[i];
					h = Pythag(f, g);
					W[i] = h;
					h = 1.0 / h;
					c = g * h;
					s = -f * h;
					for (j = 0L; (j < Lines); j++) {
						y = U[j * Columns + nm];
						z = U[j * Columns + i];
						U[j * Columns + nm] = y * c + z * s;
						U[j * Columns + i] = z * c - y * s;
					}
				}
			}
			z = W[k];
			if (l == k) {
				if (z < 0.0) {
					W[k] = -z;
					for (j = 0L; (j < Columns); j++) {
						V[j * Columns + k] = -V[j * Columns + k];
					}
				}
				break;
			}
			if (its == MaxIterations) {
				*Status = ERROR;
				FreeMatrix(&rv1);
				WRITE_ERROR(SingularValueDecomposition, "No convergence")
/**/			DEBUG_WRITE_LEAVING(SingularValueDecomposition, "Done")
				return(*Status);
			}
			x = W[l];
			nm = k - 1L;
			y = W[nm];
			g = rv1[nm];
			h = rv1[k];
			f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
			g = Pythag(f, 1.0);
			f = ((x - z) * (x + z) + h * ((y / (f + ((0.0 <= f) ? (fabs(g))
				: (-fabs(g))))) - h)) / x;
			c = s = 1.0;
			for (j = l; (j <= nm); j++) {
				i = j + 1L;
				g = rv1[i];
				y = W[i];
				h = s * g;
				g = c * g;
				z = Pythag(f, h);
				rv1[j] = z;
				c = f / z;
				s = h / z;
				f = x * c + g * s;
				g = g * c - x * s;
				h = y * s;
				y *= c;
				for (jj = 0L; (jj < Columns); jj++) {
					x = V[jj * Columns + j];
					z = V[jj * Columns + i];
					V[jj * Columns + j] = x * c + z * s;
					V[jj * Columns + i] = z * c - x * s;
				}
				z = Pythag(f, h);
				W[j] = z;
				if (z != 0.0) {
					z = 1.0 / z;
					c = f * z;
					s = h * z;
				}
				f = c * g + s * y;
				x = c * y - s * g;
				for (jj = 0L; (jj < Lines); jj++) {
					y = U[jj * Columns + j];
					z = U[jj * Columns + i];
					U[jj * Columns + j] = y * c + z * s;
					U[jj * Columns + i] = z * c - y * s;
				}
			}
			rv1[l] = 0.0;
			rv1[k] = f;
			W[k] = x;
		}
	}
	*Status = FreeMatrix(&rv1);
/**/DEBUG_WRITE_LEAVING(SingularValueDecomposition, "Done")
	return(*Status);
} /* end SingularValueDecomposition */

/*--------------------------------------------------------------------------*/
extern int		SquareMatrixDeterminant
				(
					double	*A,					/* input square matrix */
					long	Size,				/* height or width of the square matrix */
					double	*Determinant,		/* resulting scalar determinant */
					int		*Status				/* error management */
				)

/* compute the determinant of a square matrix */
/* the size of the input matrix A is (Size x Size) */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin SquareMatrixDeterminant */

	double	*B = (double *)NULL;
	double	Det;
	double	a0, a1, a2, a3, a4, a5;
	long	i;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(SquareMatrixDeterminant, A, *Status,
/**/	"Missing input")
/**/DEBUG_CHECK_NULL_POINTER(SquareMatrixDeterminant, Determinant, *Status,
/**/	"Missing output")
/**/DEBUG_CHECK_RANGE_LONG(SquareMatrixDeterminant, Size, 1L, LONG_MAX, *Status,
/**/	"Invalid size (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(SquareMatrixDeterminant, *Status)
/**/DEBUG_WRITE_ENTERING(SquareMatrixDeterminant,
/**/	"About to compute a determinant")

	switch (Size) {
		case 1L:
			*Determinant = *A;
			break;
		case 2L:
			*Determinant = A[0] * A[3] - A[1] * A[2];
			break;
		case 3L:
			*Determinant = A[0] * A[4] * A[8] + A[1] * A[5] * A[6] + A[2] * A[3] * A[7]
				- A[0] * A[5] * A[7] - A[1] * A[3] * A[8] - A[2] * A[4] * A[6];
			break;
		case 4L:
			a0 = A[10] * A[15] - A[11] * A[14];
			a1 = A[11] * A[13] - A[9] * A[15];
			a2 = A[9] * A[14] - A[10] * A[13];
			a3 = A[8] * A[15] - A[11] * A[12];
			a4 = A[10] * A[12] - A[8] * A[14];
			a5 = A[8] * A[13] - A[9] * A[12];
			*Determinant = A[0] * A[5] * a0 + A[0] * A[6] * a1 + A[0] * A[7] * a2
				+ A[1] * A[6] * a3 + A[1] * A[7] * a4 + A[2] * A[7] * a5
				- A[1] * A[4] * a0 - A[2] * A[4] * a1 - A[2] * A[5] * a3
				- A[3] * A[4] * a2 - A[3] * A[5] * a4 - A[3] * A[6] * a5;
			break;
		default:
			*Determinant = 0.0;
			AllocateMatrix(&B, Size - 1L, Size - 1L, Status);
			if (*Status == ERROR) {
/**/			DEBUG_WRITE_LEAVING(SquareMatrixDeterminant, "Done")
				return(*Status);
			}
			for (i = 0L; (i < Size); i++) {
				*Status = MatrixMinor(A, B, Size, Size, 0L, i);
				if (*Status == ERROR) {
					FreeMatrix(&B);
/**/				DEBUG_WRITE_LEAVING(SquareMatrixDeterminant, "Done")
					return(*Status);
				}
				SquareMatrixDeterminant(B, Size - 1L, &Det, Status);
				if (*Status == ERROR) {
					FreeMatrix(&B);
/**/				DEBUG_WRITE_LEAVING(SquareMatrixDeterminant, "Done")
					return(*Status);
				}
				if ((i & 1L) != 0L) {
					*Determinant -= A[i] * Det;
				}
				else {
					*Determinant += A[i] * Det;
				}
			}
			*Status = FreeMatrix(&B);
			break;
	}
/**/DEBUG_WRITE_LEAVING(SquareMatrixDeterminant, "Done")
	return(*Status);
} /* end SquareMatrixDeterminant */

/*--------------------------------------------------------------------------*/
extern int		SquareMatrixInvertGauss
				(
					double	*Direct,			/* input matrix */
					double	*Inverse,			/* resulting matrix */
					long	Size,				/* height or width of the square matrix */
					double	Tolerance,			/* admissible relative error */
					int		*Status				/* error management */
				)

/* compute the inverse of a square matrix such that Direct.Inverse = Identity */
/* the size of the matrices Direct and Inverse is (Size x Size) */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin SquareMatrixInvertGauss */

	double	*Workspace = (double *)NULL;
	double	*p, *q, *r, *s, *t, *u, *v;
	double	AbsMax, Max;
	double	Swap, Gain;
	long	i, j, k;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(SquareMatrixInvertGauss, Direct, *Status,
/**/	"No input matrix")
/**/DEBUG_CHECK_NULL_POINTER(SquareMatrixInvertGauss, Inverse, *Status,
/**/	"No output matrix")
/**/DEBUG_CHECK_RANGE_LONG(SquareMatrixInvertGauss, Size, 1L, LONG_MAX, *Status,
/**/	"Invalid size (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_DOUBLE(SquareMatrixInvertGauss, Tolerance, 0.0, DBL_MAX, *Status,
/**/	"Invalid Tolerance (should be positive)")
/**/DEBUG_RETURN_ON_ERROR(SquareMatrixInvertGauss, *Status)
/**/DEBUG_WRITE_ENTERING(SquareMatrixInvertGauss,
/**/	"About to invert a matrix by Gauss-Jordan")

	/* Workspace allocation */
	AllocateMatrix(&Workspace, Size, Size, Status);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(SquareMatrixInvertGauss, "Done")
		return(*Status);
	}
	Workspace = (double *)memcpy(Workspace, Direct,
		(size_t)(Size * Size * (long)sizeof(double)));
	/* Null matrix */
	q = Inverse;
	for (i = Size; (0L < i); i--) {
		for (j = Size; (0L < j); j--) {
			*q++ = 0.0;
		}
	}
	/* Normalize */
	p = Workspace;
	q = Inverse;
	r = p;
	for (i = 0L; (i < Size); i++) {
		Max = *p++;
		AbsMax = fabs(Max);
		for (j = 1L; (j < Size); p++, j++) {
			if (AbsMax < fabs(*p)) {
				Max = *p;
				AbsMax = fabs(Max);
			}
		}
		if (AbsMax < Tolerance) {
			*Status = ERROR;
			FreeMatrix(&Workspace);
			WRITE_ERROR(SquareMatrixInvertGauss, "Singular matrix")
/**/		DEBUG_WRITE_LEAVING(SquareMatrixInvertGauss, "Done")
			return(*Status);
		}
		Max = 1.0 / Max;
		q[i] = Max;
		for (j = Size; (0L < j); j--) {
			*r++ *= Max;
		}
		q += (ptrdiff_t)Size;
	}
	p = Workspace;
	u = Workspace;
	v = Inverse;
	for (j = 0L; (j < Size); j++) {
		/* Select pivot */
		q = p++;
		r = q;
		Max = *r;
		AbsMax = fabs(Max);
		r += (ptrdiff_t)Size;
		for (i = 1L; (i < (Size - j)); i++) {
			if (AbsMax < fabs(*r)) {
				Max = *r;
				AbsMax = fabs(Max);
				q = r;
			}
			r += (ptrdiff_t)Size;
		}
		if (AbsMax < Tolerance) {
			*Status = ERROR;
			FreeMatrix(&Workspace);
			WRITE_ERROR(SquareMatrixInvertGauss, "Singular matrix")
/**/		DEBUG_WRITE_LEAVING(SquareMatrixInvertGauss, "Done")
			return(*Status);
		}
		Max = 1.0 / Max;
		/* Permute lines */
		q -= (ptrdiff_t)j;
		if (q != u) {
			r = u;
			t = v;
			s = t + (q - r);
			for (i = 0L; (i < Size); i++) {
				Swap = *r;
				*r++ = *q;
				*q++ = Swap;
				Swap = *t;
				*t++ = *s;
				*s++ = Swap;
			}
		}
		/* Subtract line */
		q = Workspace + (ptrdiff_t)j;
		r = Inverse;
		for (i = 0L; (i < j); i++) {
			Gain = *q * Max;
			s = u + (ptrdiff_t)j;
			t = v;
			for (k = 0L; (k < j); k++) {
				*r++ -= Gain * *t++;
			}
			while (k++ < Size) {
				*q++ -= Gain * *s++;
				*r++ -= Gain * *t++;
			}
			q += (ptrdiff_t)j;
		}
		i++;
		q += (ptrdiff_t)Size;
		r += (ptrdiff_t)Size;
		while (i++ < Size) {
			Gain = *q * Max;
			s = u + (ptrdiff_t)j;
			t = v;
			for (k = 0L; (k < j); k++) {
				*r++ -= Gain * *t++;
			}
			while (k++ < Size) {
				*q++ -= Gain * *s++;
				*r++ -= Gain * *t++;
			}
			q += (ptrdiff_t)j;
		}
		p += (ptrdiff_t)Size;
		u += (ptrdiff_t)Size;
		v += (ptrdiff_t)Size;
	}
	/* Re-normalize */
	p = Workspace;
	q = Inverse;
	for (i = Size; (0L < i); i--) {
		Gain = 1.0 / *p++;
		for (j = Size; (0L < j); j--) {
			*q++ *= Gain;
		}
		p += (ptrdiff_t)Size;
	}
	*Status = FreeMatrix(&Workspace);
/**/DEBUG_WRITE_LEAVING(SquareMatrixInvertGauss, "Done")
	return(*Status);
} /* end SquareMatrixInvertGauss */

/*--------------------------------------------------------------------------*/
extern int		SquareMatrixSquareRoot
				(
					double	*A,					/* input matrix */
					double	*B,					/* output matrix and initial solution */
					long	Size,				/* height or width of the square matrix */
					double	Tolerance,			/* admissible relative error */
					long	MaxIterations,		/* convergence limit */
					int		*Status				/* error management */
				)

/* compute one out of the 2^Size square roots B of a matrix A such that B.B = A */
/* the size of the matrices A and B is (Size x Size) */
/* B must provide an initial solution */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin SquareMatrixSquareRoot */

	double	*x, *y, *z;
	double	*R = (double *)NULL, *X = (double *)NULL, *Y = (double *)NULL, *Z = (double *)NULL;
	double	Err, Eps, OldEps, Norm;
	long	Iteration;
	long	i, j;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(SquareMatrixSquareRoot, A, *Status,
/**/	"Empty input")
/**/DEBUG_CHECK_NULL_POINTER(SquareMatrixSquareRoot, B, *Status,
/**/	"Empty output")
/**/DEBUG_CHECK_RANGE_LONG(SquareMatrixSquareRoot, Size, 1L, LONG_MAX, *Status,
/**/	"Invalid size (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(SquareMatrixSquareRoot, MaxIterations, 1L, LONG_MAX, *Status,
/**/	"Invalid iterations (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_DOUBLE(SquareMatrixSquareRoot, Tolerance, 0.0, DBL_MAX, *Status,
/**/	"Invalid Tolerance (should be positive)")
/**/DEBUG_RETURN_ON_ERROR(SquareMatrixSquareRoot, *Status)
/**/DEBUG_WRITE_ENTERING(SquareMatrixSquareRoot,
/**/	"About to compute one square root of a matrix")

	*Status = FrobeniusNorm(A, &Norm, Size, Size);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(SquareMatrixSquareRoot, "Done")
		return(*Status);
	}
	if (Norm < Tolerance) {
		for (i = Size; (0L < i); i--) {
			for (j = Size; (0L < j); j--) {
				*B++ = 0.0;
			}
		}
/**/	DEBUG_WRITE_LEAVING(SquareMatrixSquareRoot, "Done")
		return(*Status);
	}
	AllocateMatrix(&R, Size, Size, Status);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(SquareMatrixSquareRoot, "Done")
		return(*Status);
	}
	AllocateMatrix(&X, Size, Size, Status);
	if (*Status == ERROR) {
		FreeMatrix(&R);
/**/	DEBUG_WRITE_LEAVING(SquareMatrixSquareRoot, "Done")
		return(*Status);
	}
	AllocateMatrix(&Y, Size, Size, Status);
	if (*Status == ERROR) {
		FreeMatrix(&X);
		FreeMatrix(&R);
/**/	DEBUG_WRITE_LEAVING(SquareMatrixSquareRoot, "Done")
		return(*Status);
	}
	AllocateMatrix(&Z, Size, Size, Status);
	if (*Status == ERROR) {
		FreeMatrix(&Y);
		FreeMatrix(&X);
		FreeMatrix(&R);
/**/	DEBUG_WRITE_LEAVING(SquareMatrixSquareRoot, "Done")
		return(*Status);
	}
	Iteration = 0L;
	Z = (double *)memcpy(Z, B, (size_t)(Size * Size * (long)sizeof(double)));
	OldEps = DBL_MAX;
	do {
		SquareMatrixInvertGauss(Z, R, Size, Tolerance, Status);
		if (*Status == ERROR) {
			FreeMatrix(&Z);
			FreeMatrix(&Y);
			FreeMatrix(&X);
			FreeMatrix(&R);
/**/		DEBUG_WRITE_LEAVING(SquareMatrixSquareRoot, "Done")
			return(*Status);
		}
		*Status = MatrixMultiply(R, A, X, Size, Size, Size);
		if (*Status == ERROR) {
			FreeMatrix(&Z);
			FreeMatrix(&Y);
			FreeMatrix(&X);
			FreeMatrix(&R);
/**/		DEBUG_WRITE_LEAVING(SquareMatrixSquareRoot, "Done")
			return(*Status);
		}
		*Status = MatrixMultiply(A, R, Y, Size, Size, Size);
		if (*Status == ERROR) {
			FreeMatrix(&Z);
			FreeMatrix(&Y);
			FreeMatrix(&X);
			FreeMatrix(&R);
/**/		DEBUG_WRITE_LEAVING(SquareMatrixSquareRoot, "Done")
			return(*Status);
		}
		x = X;
		y = Y;
		z = Z;
		for (i = Size; (0L < i); i--) {
			for (j = Size; (0L < j); j--) {
				*z += *z + *x++ + *y++;
				*z++ *= 0.25;
			}
		}
		*Status = MatrixMultiply(Z, Z, X, Size, Size, Size);
		if (*Status == ERROR) {
			FreeMatrix(&Z);
			FreeMatrix(&Y);
			FreeMatrix(&X);
			FreeMatrix(&R);
/**/		DEBUG_WRITE_LEAVING(SquareMatrixSquareRoot, "Done")
			return(*Status);
		}
		z = A;
		x = X;
		Eps = 0.0;
		for (i = Size; (0L < i); i--) {
			for (j = Size; (0L < j); j--) {
				Err = *z++ - *x++;
				Eps += Err * Err;
			}
		}
		Eps = sqrt(Eps);
		Eps /= Norm;
		if (Eps <= OldEps) {
			B = (double *)memcpy(B, Z, (size_t)(Size * Size * (long)sizeof(double)));
			OldEps = Eps;
		}
		Iteration++;
	} while ((Tolerance < Eps) && (Iteration < MaxIterations));
	*Status = FreeMatrix(&Z);
	if (*Status == ERROR) {
		FreeMatrix(&Y);
		FreeMatrix(&X);
		FreeMatrix(&R);
/**/	DEBUG_WRITE_LEAVING(SquareMatrixSquareRoot, "Done")
		return(*Status);
	}
	*Status = FreeMatrix(&Y);
	if (*Status == ERROR) {
		FreeMatrix(&X);
		FreeMatrix(&R);
/**/	DEBUG_WRITE_LEAVING(SquareMatrixSquareRoot, "Done")
		return(*Status);
	}
	*Status = FreeMatrix(&X);
	if (*Status == ERROR) {
		FreeMatrix(&R);
/**/	DEBUG_WRITE_LEAVING(SquareMatrixSquareRoot, "Done")
		return(*Status);
	}
	*Status = FreeMatrix(&R);
/**/DEBUG_WRITE_LEAVING(SquareMatrixSquareRoot, "Done")
	return(*Status);
} /* end SquareMatrixSquareRoot */

/*--------------------------------------------------------------------------*/
extern int		SquareMatrixTrace
				(
					double	*A,					/* input matrix */
					double	*Trace,				/* resulting scalar trace */
					long	Size				/* height or width of the square matrix */
				)

/* perform the vector cross product X = A x B */
/* the size of A is (Size x Size) */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin SquareMatrixTrace */

	double	Sum;
	long	i;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(SquareMatrixTrace, A, Status,
/**/	"Empty input")
/**/DEBUG_CHECK_NULL_POINTER(SquareMatrixTrace, Trace, Status,
/**/	"Empty output")
/**/DEBUG_CHECK_RANGE_LONG(SquareMatrixTrace, Size, 1L, LONG_MAX, Status,
/**/	"Invalid size (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(SquareMatrixTrace, Status)
/**/DEBUG_WRITE_ENTERING(SquareMatrixTrace,
/**/	"About to compute the trace of a matrix")

	Sum = 0.0;
	for (i = -Size; (i < 0L); i++) {
		Sum += *A++;
		A += (ptrdiff_t)Size;
	}
	*Trace = Sum;
/**/DEBUG_WRITE_LEAVING(SquareMatrixTrace, "Done")
	return(Status);
} /* end SquareMatrixTrace */

/*--------------------------------------------------------------------------*/
extern int		Vector3DcrossProduct
				(
					double	A[],				/* left vector operand */
					double	B[],				/* right vector operand */
					double	X[]					/* resulting vector */
				)

/* perform the vector cross product X = A x B */
/* the size of A, B and X is (3 x 1) */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin Vector3DcrossProduct */

	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(Vector3DcrossProduct, A, Status,
/**/	"Empty input left operand")
/**/DEBUG_CHECK_NULL_POINTER(Vector3DcrossProduct, B, Status,
/**/	"Empty input right operand")
/**/DEBUG_CHECK_NULL_POINTER(Vector3DcrossProduct, X, Status,
/**/	"Missing result")
/**/DEBUG_RETURN_ON_ERROR(Vector3DcrossProduct, Status)
/**/DEBUG_WRITE_ENTERING(Vector3DcrossProduct,
/**/	"About to compute a cross product")

	*X++ = A[1] * B[2] - A[2] * B[1];
	*X++ = A[2] * B[0] - A[0] * B[2];
	*X = A[0] * B[1] - A[1] * B[0];
/**/DEBUG_WRITE_LEAVING(Vector3DcrossProduct, "Done")
	return(Status);
} /* end Vector3DcrossProduct */

/*--------------------------------------------------------------------------*/
extern int		VectorAdd
				(
					double	A[],				/* left vector operand */
					double	B[],				/* right vector operand */
					double	X[],				/* resulting vector */
					long	Lines				/* height of the vectors */
				)

/* perform the vector addition X = A + B */
/* the vectors A, B and X have (Lines) elements */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin VectorAdd */

	long	i;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(VectorAdd, A, Status,
/**/	"Empty input left operand")
/**/DEBUG_CHECK_NULL_POINTER(VectorAdd, B, Status,
/**/	"Empty input right operand")
/**/DEBUG_CHECK_NULL_POINTER(VectorAdd, X, Status,
/**/	"Missing result")
/**/DEBUG_CHECK_RANGE_LONG(VectorAdd, Lines, 1L, LONG_MAX, Status,
/**/	"Invalid Lines (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(VectorAdd, Status)
/**/DEBUG_WRITE_ENTERING(VectorAdd,
/**/	"About to add two vectors")

	for (i = Lines; (0L < i); i--) {
		*X++ = *A++ + *B++;
	}
/**/DEBUG_WRITE_LEAVING(VectorAdd, "Done")
	return(Status);
} /* end VectorAdd */

/*--------------------------------------------------------------------------*/
extern int		VectorAngle
				(
					double	A[],				/* 1st input vector */
					double	B[],				/* 2nd input vector */
					double	*Angle,				/* resulting scalar angle in rad */
					long	Lines,				/* height of the vectors */
					double	Tolerance			/* admissible relative error */
				)

/* compute the angle between two vectors */
/* the size of A and B is (Lines x 1) */
/* the angle unit is radian */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin VectorAngle */

	double	NormA, NormB, Sum;
	long	i;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(VectorAngle, A, Status,
/**/	"Empty input A")
/**/DEBUG_CHECK_NULL_POINTER(VectorAngle, B, Status,
/**/	"Empty input B")
/**/DEBUG_CHECK_NULL_POINTER(VectorAngle, Angle, Status,
/**/	"Missing result")
/**/DEBUG_CHECK_RANGE_LONG(VectorAngle, Lines, 1L, LONG_MAX, Status,
/**/	"Invalid Lines (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_DOUBLE(VectorAngle, Tolerance, 0.0, DBL_MAX, Status,
/**/	"Invalid Tolerance (should be positive)")
/**/DEBUG_RETURN_ON_ERROR(VectorAngle, Status)
/**/DEBUG_WRITE_ENTERING(VectorAngle,
/**/	"About to compute the angle between {A,B}")

	NormA = 0.0;
	NormB = 0.0;
	Sum = 0.0;
	for (i = Lines; (0L < i); i--) {
		NormA += *A * *A;
		NormB += *B * *B;
		Sum += *A++ * *B++;
	}
	NormA = sqrt(NormA);
	NormB = sqrt(NormB);
	NormA *= NormB;
	if (NormA < Tolerance) {
		WRITE_WARNING(VectorAngle, "Undefined angle")
/**/	DEBUG_WRITE_LEAVING(VectorAngle, "Done")
		return(Status);
	}
	Sum /= NormA;
	if ((Sum < -1.0) || (1.0 < Sum)) {
		Status = ERROR;
		WRITE_ERROR(VectorAngle, "Unexpected numeric error")
/**/	DEBUG_WRITE_LEAVING(VectorAngle, "Done")
		return(Status);
	}
	*Angle = acos(Sum);
/**/DEBUG_WRITE_LEAVING(VectorAngle, "Done")
	return(Status);
} /* end VectorAngle */

/*--------------------------------------------------------------------------*/
extern int		VectorHomogenize
				(
					double	A[],				/* input vector */
					double	X[],				/* output homogenous vector */
					long	Lines,				/* height of the vectors */
					double	Tolerance			/* admissible relative error */
				)

/* convert a vector from non-homogenous to homogenous coordinates */
/* satisfy (Vout = a * Vin) such that (Vout[last] = 1.0) */
/* the size of A and X is (Lines x 1) */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin VectorHomogenize */

	double	Scale;
	long	i;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(VectorHomogenize, A, Status,
/**/	"Empty input A")
/**/DEBUG_CHECK_NULL_POINTER(VectorHomogenize, X, Status,
/**/	"Empty output X")
/**/DEBUG_CHECK_RANGE_LONG(VectorHomogenize, Lines, 1L, LONG_MAX, Status,
/**/	"Invalid Lines (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_DOUBLE(VectorHomogenize, Tolerance, 0.0, DBL_MAX, Status,
/**/	"Invalid Tolerance (should be positive)")
/**/DEBUG_RETURN_ON_ERROR(VectorHomogenize, Status)
/**/DEBUG_WRITE_ENTERING(VectorHomogenize,
/**/	"About to convert a vector into homogenous coordinates")

	Scale = A[Lines - 1L];
	if (fabs(Scale) < Tolerance) {
		Status = ERROR;
		WRITE_ERROR(VectorHomogenize,
			"Unable to perform conversion in homogenous coordinates")
/**/	DEBUG_WRITE_LEAVING(VectorHomogenize, "Done")
		return(Status);
	}
	Scale = 1.0 / Scale;
	for (i = Lines - 1L; (0L < i); i--) {
		*X++ = Scale * *A++;
	}
	*X = 1.0;
/**/DEBUG_WRITE_LEAVING(VectorHomogenize, "Done")
	return(Status);
} /* end VectorHomogenize */

/*--------------------------------------------------------------------------*/
extern int		VectorNorm
				(
					double	A[],				/* input vector */
					double	*Norm,				/* resulting scalar norm */
					long	Lines				/* height of the vector */
				)

/* compute the norm of a vector */
/* the size of A is (Lines x 1) */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin VectorNorm */

	double	Sum;
	long	i;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(VectorNorm, A, Status,
/**/	"Empty input A")
/**/DEBUG_CHECK_NULL_POINTER(VectorNorm, Norm, Status,
/**/	"Missing result")
/**/DEBUG_CHECK_RANGE_LONG(VectorNorm, Lines, 1L, LONG_MAX, Status,
/**/	"Invalid Lines (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(VectorNorm, Status)
/**/DEBUG_WRITE_ENTERING(VectorNorm,
/**/	"About to compute the norm of a vector")

	Sum = 0.0;
	for (i = Lines; (0L < i); i--) {
		Sum += *A * *A;
		A++;
	}
	*Norm = sqrt(Sum);
/**/DEBUG_WRITE_LEAVING(VectorNorm, "Done")
	return(Status);
} /* end VectorNorm */

/*--------------------------------------------------------------------------*/
extern int		VectorNormalize
				(
					double	A[],				/* input vector */
					double	X[],				/* resulting normalized vector */
					long	Lines				/* height of the vector */
				)

/* normalize of a vector to unit length */
/* the size of A and X is (Lines x 1) */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin VectorNormalize */

	double	*p;
	double	Sum;
	long	i;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(VectorNormalize, A, Status,
/**/	"Empty input A")
/**/DEBUG_CHECK_NULL_POINTER(VectorNormalize, X, Status,
/**/	"Missing result")
/**/DEBUG_CHECK_RANGE_LONG(VectorNormalize, Lines, 1L, LONG_MAX, Status,
/**/	"Invalid Lines (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(VectorNormalize, Status)
/**/DEBUG_WRITE_ENTERING(VectorNormalize,
/**/	"About to compute the norm of a vector")

	p = A;
	Sum = 0.0;
	for (i = Lines; (0L < i); i--) {
		Sum += *p * *p;
		p++;
	}
	Sum = sqrt(Sum);
	if (Sum == 0.0) {
		WRITE_WARNING(VectorNormalize, "Degenerate vector")
	}
	else {
		Sum = 1.0 / Sum;
	}
	for (i = Lines; (0L < i); i--) {
		*X++ = Sum * *A++;
	}
/**/DEBUG_WRITE_LEAVING(VectorNormalize, "Done")
	return(Status);
} /* end VectorNormalize */

/*--------------------------------------------------------------------------*/
extern int		VectorScalarProduct
				(
					double	A[],				/* left vector operand */
					double	B[],				/* right vector operand */
					double	*X,					/* resulting scalar product */
					long	Lines				/* height of the vector */
				)

/* perform the vector scalar product X = Transpose(A).B */
/* the size of A and B is (Lines x 1) */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin VectorScalarProduct */

	double	Sum;
	long	i;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(VectorScalarProduct, A, Status,
/**/	"Empty input left operand")
/**/DEBUG_CHECK_NULL_POINTER(VectorScalarProduct, B, Status,
/**/	"Empty input right operand")
/**/DEBUG_CHECK_NULL_POINTER(VectorScalarProduct, X, Status,
/**/	"Missing result")
/**/DEBUG_CHECK_RANGE_LONG(VectorScalarProduct, Lines, 1L, LONG_MAX, Status,
/**/	"Invalid Lines (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(VectorScalarProduct, Status)
/**/DEBUG_WRITE_ENTERING(VectorScalarProduct,
/**/	"About to compute the scalar product of two vectors")

	Sum = 0.0;
	for (i = Lines; (0L < i); i--) {
		Sum += *A++ * *B++;
	}
	*X = Sum;
/**/DEBUG_WRITE_LEAVING(VectorScalarProduct, "Done")
	return(Status);
} /* end VectorScalarProduct */

/*--------------------------------------------------------------------------*/
extern int		VectorScale
				(
					double	A[],				/* input vector */
					double	X[],				/* output vector */
					double	Scale,				/* scale factor */
					long	Lines				/* height of the vector */
				)

/* scale a vector  */
/* satisfy (Vout = Scale * Vin) */
/* the size of A and X is (Lines x 1) */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin VectorScale */

	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(VectorScale, A, Status,
/**/	"Empty input A")
/**/DEBUG_CHECK_NULL_POINTER(VectorScale, X, Status,
/**/	"Empty output X")
/**/DEBUG_CHECK_RANGE_LONG(VectorScale, Lines, 1L, LONG_MAX, Status,
/**/	"Invalid Lines (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(VectorScale, Status)
/**/DEBUG_WRITE_ENTERING(VectorScale,
/**/	"About to scale a vector")

	while (Lines-- > 0L) {
		*X++ = Scale * *A++;
	}
/**/DEBUG_WRITE_LEAVING(VectorScale, "Done")
	return(Status);
} /* end VectorScale */

/*--------------------------------------------------------------------------*/
extern int		VectorSubtract
				(
					double	A[],				/* left vector operand */
					double	B[],				/* right vector operand */
					double	X[],				/* resulting vector */
					long	Lines				/* height of the vectors */
				)

/* perform the vector subtraction X = A - B */
/* the vectors A, B and X have (Lines) elements */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin VectorSubtract */

	long	i;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(VectorSubtract, A, Status,
/**/	"Empty input left operand")
/**/DEBUG_CHECK_NULL_POINTER(VectorSubtract, B, Status,
/**/	"Empty input right operand")
/**/DEBUG_CHECK_NULL_POINTER(VectorSubtract, X, Status,
/**/	"Missing result")
/**/DEBUG_CHECK_RANGE_LONG(VectorSubtract, Lines, 1L, LONG_MAX, Status,
/**/	"Invalid Lines (should be strictly positive)")
/**/DEBUG_RETURN_ON_ERROR(VectorSubtract, Status)
/**/DEBUG_WRITE_ENTERING(VectorSubtract,
/**/	"About to subtract two vectors")

	for (i = Lines; (0L < i); i--) {
		*X++ = *A++ - *B++;
	}
/**/DEBUG_WRITE_LEAVING(VectorSubtract, "Done")
	return(Status);
} /* end VectorSubtract */

