/*****************************************************************************
 *	System includes
 ****************************************************************************/
#include	<math.h>
#include	<stdio.h>

/*****************************************************************************
 *	Toolbox defines
 ****************************************************************************/
#include	"configs.h"
#include	"debug.h"
#include	"error.h"

/*****************************************************************************
 *	Toolbox includes
 ****************************************************************************/
#include	"geometry.h"
#include	"linearalgebra.h"
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
extern int		AffineToRotationMatrix
				(
					double	*A,					/* input 4x4 homogenous affine matrix */
					double	*R,					/* output 4x4 homogenous rotation matrix */
					double	*xRotation,			/* output rotation around x-axis */
					double	*yRotation,			/* output rotation around y-axis */
					double	*zRotation,			/* output rotation around z-axis */
					int		FirstOctant,		/* constrains the rotation */
					double	Tolerance,			/* admissible relative error */
					long	MaxIterations,		/* convergence limit */
					int		*Status				/* error management */
				)

/* approximate a (4 x 4) homogenous affine matrix by a (4 x 4) rotation matrix */
/* the rotation matrix has the form R = Rx.Ry.Rz */
/* A and R are homogenous: (A, R) = {{*, *, *, 0}, {*, *, *, 0}, {*, *, *, 0}, {0, 0, 0, 1}} */
/* the returned rotation angles are given in the unit of radian */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin AffineToRotationMatrix */

	double	*p, *q;
	double	U[3*3], W[3], V[3*3];
	double	Angles[32][4];
	double	InvDeterminant, SecTheta;
	double	Phi, Theta, Psi;
	double	Best;
	long	i, j ,k;
	long	iBest;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(AffineToRotationMatrix, A, *Status,
/**/	"Missing input")
/**/DEBUG_CHECK_NULL_POINTER(AffineToRotationMatrix, R, *Status,
/**/	"Missing output")
/**/DEBUG_CHECK_NULL_POINTER(AffineToRotationMatrix, xRotation, *Status,
/**/	"Missing xRotation")
/**/DEBUG_CHECK_NULL_POINTER(AffineToRotationMatrix, yRotation, *Status,
/**/	"Missing yRotation")
/**/DEBUG_CHECK_NULL_POINTER(AffineToRotationMatrix, zRotation, *Status,
/**/	"Missing zRotation")
/**/DEBUG_CHECK_RANGE_INT(AffineToRotationMatrix, FirstOctant, FALSE, TRUE, *Status,
/**/	"Invalid FirstOctant (should be in [FALSE, TRUE])")
/**/DEBUG_CHECK_RANGE_LONG(AffineToRotationMatrix, MaxIterations, 1L, LONG_MAX, *Status,
/**/	"Invalid iterations (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_DOUBLE(AffineToRotationMatrix, Tolerance, 0.0, DBL_MAX, *Status,
/**/	"Invalid Tolerance (should be positive)")
/**/DEBUG_RETURN_ON_ERROR(AffineToRotationMatrix, *Status)
/**/DEBUG_WRITE_ENTERING(AffineToRotationMatrix,
/**/	"About to approximate an affine matrix by a rotation")

	if ((A[12] != 0.0) || (A[13] != 0.0) || (A[14] != 0.0) || (A[15] != 1.0)) {
		*Status = ERROR;
		WRITE_ERROR(AffineToRotationMatrix, "Non-homogenous matrix")
/**/	DEBUG_WRITE_LEAVING(AffineToRotationMatrix, "Done")
		return(*Status);
	}
	p = A;
	q = U;
	for (i = -3L; (i < 0L); i++) {
		for (j = -3L; (j < 0L); j++) {
			*q++ = *p++;
		}
		p++;
	}
	SingularValueDecomposition(U, 3L, 3L, W, R, MaxIterations, Status);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(AffineToRotationMatrix, "Done")
		return(*Status);
	}
	*Status = MatrixTranspose(R, V, 3L, 3L);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(AffineToRotationMatrix, "Done")
		return(*Status);
	}
	*Status = MatrixMultiply(U, V, R, 3L, 3L, 3L);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(AffineToRotationMatrix, "Done")
		return(*Status);
	}
	SquareMatrixDeterminant(R, 3L, &InvDeterminant, Status);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(AffineToRotationMatrix, "Done")
		return(*Status);
	}
	if (fabs(InvDeterminant) <= Tolerance) {
		*Status = ERROR;
		WRITE_ERROR(AffineToRotationMatrix, "Degenerate Matrix")
/**/	DEBUG_WRITE_LEAVING(AffineToRotationMatrix, "Done")
		return(*Status);
	}
	if (InvDeterminant < 0.0) {
		WRITE_ERROR(AffineToRotationMatrix, "Left-handed matrix")
		*Status = ERROR;
/**/	DEBUG_WRITE_LEAVING(AffineToRotationMatrix, "Done")
		return(*Status);
	}
	InvDeterminant = 1.0 / InvDeterminant;
	Theta = asin(InvDeterminant * R[2]);
	for (k = 0L; (k < 4L); k++) {
		Angles[k][1] = Theta;
		Angles[k + 8L][1] = Theta;
		Angles[k + 16L][1] = Theta;
		Angles[k + 24L][1] = Theta;
	}
	Theta = PI - Theta;
	for (k = 4L; (k < 8L); k++) {
		Angles[k][1] = Theta;
		Angles[k + 8L][1] = Theta;
		Angles[k + 16L][1] = Theta;
		Angles[k + 24L][1] = Theta;
	}
	SecTheta = cos(Theta);
	if (fabs(SecTheta) <= Tolerance) {
		WRITE_WARNING(AffineToRotationMatrix, "Gimbals lock")
		if (0.0 < R[2]) {
			*yRotation = 0.5 * PI;
			*xRotation = 0.5 * atan2(R[7] + R[3], R[4] - R[6]);
			*zRotation = *xRotation;
		}
		else {
			*yRotation = -0.5 * PI;
			*xRotation = 0.5 * atan2(R[7] - R[3], R[4] + R[6]);
			*zRotation = -*xRotation;
		}
	}
	else {
		SecTheta = -1.0 / SecTheta;
		Phi = asin(-InvDeterminant * R[5] * SecTheta);
		for (k = 0L; (k < 8L); k++) {
			Angles[k][0] = Phi;
		}
		Phi = PI - Phi;
		for (k = 8L; (k < 16L); k++) {
			Angles[k][0] = Phi;
		}
		Phi = asin(InvDeterminant * R[5] * SecTheta);
		for (k = 16L; (k < 24L); k++) {
			Angles[k][0] = Phi;
		}
		Phi = PI - Phi;
		for (k = 24L; (k < 32L); k++) {
			Angles[k][0] = Phi;
		}
		Psi = asin(-InvDeterminant * R[1] * SecTheta);
		for (k = 0L; (k < 32L); k += 4L) {
			Angles[k][2] = Psi;
		}
		Psi = PI - Psi;
		for (k = 1L; (k < 32L); k += 4L) {
			Angles[k][2] = Psi;
		}
		Psi = asin(InvDeterminant * R[1] * SecTheta);
		for (k = 2L; (k < 32L); k += 4L) {
			Angles[k][2] = Psi;
		}
		Psi = PI - Psi;
		for (k = 3L; (k < 32L); k += 4L) {
			Angles[k][2] = Psi;
		}
		for (k = 0L; (k < 32L); k++) {
			*Status = GetRotationMatrix(R, Angles[k][0], Angles[k][1], Angles[k][2]);
			if (*Status == ERROR) {
/**/			DEBUG_WRITE_LEAVING(AffineToRotationMatrix, "Done")
				return(*Status);
			}
			Angles[k][3] = 0.0;
			p = A;
			q = R;
			for (i = -3L; (i < 0L); i++) {
				for (j = -3L; (j < 0L); j++) {
					Angles[k][3] += (*p - *q) * (*p - *q);
					p++;
					q++;
				}
				p++;
				q++;
			}
		}
		Best = Angles[0][3];
		iBest = 0L;
		if (!FirstOctant)
			for (k = 1L; (k < 32L); k++)
				if (Angles[k][3] < Best) {
					Best = Angles[k][3];
					iBest = k;
				}
		*xRotation = Angles[iBest][0];
		*yRotation = Angles[iBest][1];
		*zRotation = Angles[iBest][2];
		*Status = GetRotationMatrix(R, Angles[iBest][0], Angles[iBest][1], Angles[iBest][2]);
		if (*Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(AffineToRotationMatrix, "Done")
			return(*Status);
		}
		R[3] = A[3];
		R[7] = A[7];
		R[11] = A[11];
	}
	*xRotation = (*xRotation < 0.0) ? (PI + fmod(*xRotation, PI)) : (fmod(*xRotation, PI));
	*xRotation = (*xRotation <= (0.5 * PI)) ? (*xRotation) : (*xRotation - PI);
	*yRotation = (*yRotation < 0.0) ? (PI + fmod(*yRotation, PI)) : (fmod(*yRotation, PI));
	*yRotation = (*yRotation <= (0.5 * PI)) ? (*yRotation) : (*yRotation - PI);
	*zRotation = (*zRotation < 0.0) ? (PI + fmod(*zRotation, PI)) : (fmod(*zRotation, PI));
	*zRotation = (*zRotation <= (0.5 * PI)) ? (*zRotation) : (*zRotation - PI);
/**/DEBUG_WRITE_LEAVING(AffineToRotationMatrix, "Done")
	return(*Status);
} /* end AffineToRotationMatrix */

/*--------------------------------------------------------------------------*/
extern int		GetRotationMatrix
				(
					double	*R,					/* output 4x4 homogenous matrix */
					double	xRotation,			/* rotation around x-axis */
					double	yRotation,			/* rotation around y-axis */
					double	zRotation			/* rotation around z-axis */
				)

/* fill a rotation matrix R such that R = Rx.Ry.Rz */
/* the size of the output matrix R is (4 x 4) */
/* R is homogenous: R = {{*, *, *, 0}, {*, *, *, 0}, {*, *, *, 0}, {0, 0, 0, 1}} */
/* the rotation angles are given in the unit of radian */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin GetRotationMatrix */

	double	Rx[4][4], Ry[4][4], Rz[4][4];
	double	Rxy[4][4];
	double	Sum;
	long	i, j, k;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(GetRotationMatrix, R, Status,
/**/	"Missing output")
/**/DEBUG_RETURN_ON_ERROR(GetRotationMatrix, Status)
/**/DEBUG_WRITE_ENTERING(GetRotationMatrix,
/**/	"About to build a rotation matrix")

	for (i = 0L; (i < 4L); i++) {
		for (j = 0L; (j < 4L); j++) {
			Rx[i][j] = 0.0;
			Ry[i][j] = 0.0;
			Rz[i][j] = 0.0;
			Rxy[i][j] = 0.0;
		}
		Rx[i][i] = 1.0;
		Ry[i][i] = 1.0;
		Rz[i][i] = 1.0;
		Rxy[i][i] = 1.0;
	}
	Rx[1][1] = cos(xRotation);
	Rx[1][2] = -sin(xRotation);
	Rx[2][1] = sin(xRotation);
	Rx[2][2] = cos(xRotation);
	Ry[0][0] = cos(yRotation);
	Ry[0][2] = sin(yRotation);
	Ry[2][0] = -sin(yRotation);
	Ry[2][2] = cos(yRotation);
	Rz[0][0] = cos(zRotation);
	Rz[0][1] = -sin(zRotation);
	Rz[1][0] = sin(zRotation);
	Rz[1][1] = cos(zRotation);
	for (i = 0L; (i < 3L); i++) {
		for (j = 0L; (j < 3L); j++) {
			Sum = 0.0;
			for (k = 0L; (k < 3L); k++) {
				Sum += Rx[i][k] * Ry[k][j];
			}
			Rxy[i][j] = Sum;
		}
	}
	for (i = 0; (i < 4L); i++) {
		for (j = 0; (j < 4L); j++) {
			Sum = 0.0;
			for (k = 0; (k < 4L); k++) {
				Sum += Rxy[i][k] * Rz[k][j];
			}
			*R++ = Sum;
		}
	}
/**/DEBUG_WRITE_LEAVING(GetRotationMatrix, "Done")
	return(Status);
} /* end GetRotationMatrix */

/*--------------------------------------------------------------------------*/
extern int		GetScalingMatrix
				(
					double	*S,					/* output 4x4 homogenous matrix */
					double	xScale,				/* scaling along x-axis */
					double	yScale,				/* scaling along y-axis */
					double	zScale				/* scaling along z-axis */
				)

/* fill a scaling matrix T */
/* the size of the output matrix T is (4 x 4) */
/* T is homogenous: T = {{sx, 0, 0, 0}, {0, sy, 0, 0}, {0, 0, sz, 0}, {0, 0, 0, 1}} */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin GetScalingMatrix */

	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(GetScalingMatrix, S, Status,
/**/	"Missing output")
/**/DEBUG_RETURN_ON_ERROR(GetScalingMatrix, Status)
/**/DEBUG_WRITE_ENTERING(GetScalingMatrix,
/**/	"About to build a scaling matrix")

	Status = GetIdentitySquareMatrix(S, 4L);
	if (Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(GetScalingMatrix, "Done")
		return(Status);
	}
	S[0] = xScale;
	S[5] = yScale;
	S[10] = zScale;
/**/DEBUG_WRITE_LEAVING(GetScalingMatrix, "Done")
	return(Status);
} /* end GetScalingMatrix */

/*--------------------------------------------------------------------------*/
extern int		GetTranslationMatrix
				(
					double	*T,					/* output 4x4 homogenous matrix */
					double	xTranslation,		/* translation along x-axis */
					double	yTranslation,		/* translation along y-axis */
					double	zTranslation		/* translation along z-axis */
				)

/* fill a translation matrix T */
/* the size of the output matrix T is (4 x 4) */
/* T is homogenous: T = {{1, 0, 0, dx}, {0, 1, 0, dy}, {0, 0, 1, dz}, {0, 0, 0, 1}} */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin GetTranslationMatrix */

	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(GetTranslationMatrix, T, Status,
/**/	"Missing output")
/**/DEBUG_RETURN_ON_ERROR(GetTranslationMatrix, Status)
/**/DEBUG_WRITE_ENTERING(GetTranslationMatrix,
/**/	"About to build a translation matrix")

	Status = GetIdentitySquareMatrix(T, 4L);
	if (Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(GetTranslationMatrix, "Done")
		return(Status);
	}
	T[3] = xTranslation;
	T[7] = yTranslation;
	T[11] = zTranslation;
/**/DEBUG_WRITE_LEAVING(GetTranslationMatrix, "Done")
	return(Status);
} /* end GetTranslationMatrix */

/*--------------------------------------------------------------------------*/
extern int		LineToLineIntersection
				(
					double	P[],				/* 1st point on line1 */
					double	Q[],				/* 2nd point on line1 */
					double	A[],				/* 1st point on line2 */
					double	B[],				/* 2nd point on line2 */
					double	X[],				/* output intersection */
					long	Lines,				/* height of the vectors */
					double	Tolerance,			/* admissible relative error */
					long	MaxIterations,		/* convergence limit */
					int		*Intersection,		/* resulting validity of intersection */
					int		*Status				/* error management */
				)

/* compute the intersection of the line (P,Q) with the line (A,B) */
/* all vectors have (Lines) elements */
/* return Intersection = FALSE if there is no intersection */
/* return Intersection = TRUE if the intersection is valid */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin LineToLineIntersection */

	double	*u, *v, *w, *x;
	double	U[4], V[4], W[2];
	double	Y[2], Z[2];
	double	AB, PQ, PA;
	double	ABAB, ABPQ, PQPQ;
	double	ABPA, PQPA;
	long	i;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(LineToLineIntersection, P, *Status,
/**/	"Empty Line1 1st point")
/**/DEBUG_CHECK_NULL_POINTER(LineToLineIntersection, Q, *Status,
/**/	"Empty Line1 2nd point")
/**/DEBUG_CHECK_NULL_POINTER(LineToLineIntersection, A, *Status,
/**/	"Empty Line2 1st point")
/**/DEBUG_CHECK_NULL_POINTER(LineToLineIntersection, B, *Status,
/**/	"Empty Line2 2nd point")
/**/DEBUG_CHECK_NULL_POINTER(LineToLineIntersection, X, *Status,
/**/	"Empty output")
/**/DEBUG_CHECK_NULL_POINTER(LineToLineIntersection, Intersection, *Status,
/**/	"Missing result")
/**/DEBUG_CHECK_RANGE_LONG(LineToLineIntersection, Lines, 1L, LONG_MAX, *Status,
/**/	"Invalid Lines (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(LineToLineIntersection, MaxIterations, 1L, LONG_MAX, *Status,
/**/	"Invalid iterations (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_DOUBLE(LineToLineIntersection, Tolerance, 0.0, DBL_MAX, *Status,
/**/	"Invalid Tolerance (should be positive)")
/**/DEBUG_RETURN_ON_ERROR(LineToLineIntersection, *Status)
/**/DEBUG_WRITE_ENTERING(LineToLineIntersection,
/**/	"About to determine a (line/line) intersection")

	u = A;
	v = B;
	w = P;
	x = Q;
	ABAB = 0.0;
	ABPQ = 0.0;
	PQPQ = 0.0;
	ABPA = 0.0;
	PQPA = 0.0;
	for (i = Lines; (0L < i); i--) {
		AB = *v++ - *u;
		PQ = *x++ - *w;
		PA = *u++ - *w++;
		ABAB += AB * AB;
		ABPQ -= AB * PQ;
		PQPQ += PQ * PQ;
		ABPA -= AB * PA;
		PQPA += PQ * PA;
	}
	U[0] = ABAB;
	U[1] = ABPQ;
	U[2] = ABPQ;
	U[3] = PQPQ;
	SingularValueDecomposition(U, 2L, 2L, W, V, MaxIterations, Status);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(LineToLineIntersection, "Done")
		return(*Status);
	}
	W[0] = (W[0] < Tolerance) ? (0.0) : (W[0]);
	W[1] = (W[1] < Tolerance) ? (0.0) : (W[1]);
	Y[0] = ABPA;
	Y[1] = PQPA;
	SingularValueBackSubstitution(U, W, V, 2L, 2L, Y, Z, Status);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(LineToLineIntersection, "Done")
		return(*Status);
	}
	ABPQ = 0.0;
	for (i = Lines; (0L < i); i--) {
		AB = *B++ - *A;
		PQ = *Q++ - *P;
		AB = *A++ + Z[0] * AB;
		PQ = *P++ + Z[1] * PQ;
		PA = AB - PQ;
		ABPQ += PA * PA;
		*X++ = 0.5 * (AB + PQ);
	}
	*Intersection = (sqrt(ABPQ) < Tolerance);
/**/DEBUG_WRITE_LEAVING(LineToLineIntersection, "Done")
	return(*Status);
} /* end LineToLineIntersection */

/*--------------------------------------------------------------------------*/
extern int		LineToLineIntersection3D
				(
					double	P[],				/* point on line1 */
					double	q[],				/* direction of line1 */
					double	A[],				/* point on line2 */
					double	b[],				/* direction of line2 */
					double	X[],				/* output intersection */
					double	Tolerance,			/* admissible relative error */
					int		*Intersection		/* resulting validity of intersection */
				)

/* compute the intersection of the line (P,q) with the line (A,b) */
/* all points have 3 elements */
/* all vectors have 3 elements */
/* a line is described by a point and a vector: X = P + t1 * q, X = A + t2 * b */
/* where (t1,t2) are free scalars */
/* the vectors (q, b) must have a normalized unit length */
/* return Intersection = FALSE if there is no intersection */
/* return Intersection = TRUE if the intersection is valid */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin LineToLineIntersection3D */

	double	Y[3], Z[3];
	double	d, qb, PAq, PAb;
	double	tq, tb;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(LineToLineIntersection3D, P, Status,
/**/	"Empty Line1 point")
/**/DEBUG_CHECK_NULL_POINTER(LineToLineIntersection3D, q, Status,
/**/	"Empty Line1 direction")
/**/DEBUG_CHECK_NULL_POINTER(LineToLineIntersection3D, A, Status,
/**/	"Empty Line2 point")
/**/DEBUG_CHECK_NULL_POINTER(LineToLineIntersection3D, b, Status,
/**/	"Empty Line2 direction")
/**/DEBUG_CHECK_NULL_POINTER(LineToLineIntersection3D, X, Status,
/**/	"Empty output")
/**/DEBUG_CHECK_NULL_POINTER(LineToLineIntersection3D, Intersection, Status,
/**/	"Missing result")
/**/DEBUG_CHECK_RANGE_DOUBLE(LineToLineIntersection3D, Tolerance, 0.0, DBL_MAX, Status,
/**/	"Invalid Tolerance (should be positive)")
/**/DEBUG_RETURN_ON_ERROR(LineToLineIntersection3D, Status)
/**/DEBUG_WRITE_ENTERING(LineToLineIntersection3D,
/**/	"About to determine a (line/line) intersection")
#ifdef DEBUG
/**/tq = sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]);
/**/if (Tolerance < fabs(tq - 1.0)) {
/**/	Status = ERROR;
/**/	WRITE_ERROR(LineToLineIntersection3D,
/**/		"Invalid vector q (should have a unit length)")
/**/	DEBUG_WRITE_LEAVING(LineToLineIntersection3D, "Done")
/**/	return(Status);
/**/}
/**/tb = sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
/**/if (Tolerance < fabs(tb - 1.0)) {
/**/	Status = ERROR;
/**/	WRITE_ERROR(LineToLineIntersection3D,
/**/		"Invalid vector b (should have a unit length)")
/**/	DEBUG_WRITE_LEAVING(LineToLineIntersection3D, "Done")
/**/	return(Status);
/**/}
#endif

	X[0] = A[0] - P[0];
	X[1] = A[1] - P[1];
	X[2] = A[2] - P[2];
	PAq = X[0] * q[0] + X[1] * q[1] + X[2] * q[2];
	PAb = X[0] * b[0] + X[1] * b[1] + X[2] * b[2];
	qb = q[0] * b[0] + q[1] * b[1] + q[2] * b[2];
	d = 1.0 - qb * qb;
	if (Tolerance < fabs(d)) {
		tq = (PAq * qb - PAb) / d;
		tb = (PAq - qb * PAb) / d;
		Y[0] = P[0] + tq * q[0];
		Y[1] = P[1] + tq * q[1];
		Y[2] = P[2] + tq * q[2];
		Z[0] = A[0] + tb * b[0];
		Z[1] = A[1] + tb * b[1];
		Z[2] = A[2] + tb * b[2];
		X[0] = Z[0] - Y[0];
		X[1] = Z[1] - Y[1];
		X[2] = Z[2] - Y[2];
		*Intersection = (sqrt(X[0] * X[0] + X[1] * X[1] + X[2] * X[2]) <= Tolerance);
		X[0] = 0.5 * (Y[0] + Z[0]);
		X[1] = 0.5 * (Y[1] + Z[1]);
		X[2] = 0.5 * (Y[2] + Z[2]);
	}
	else {
		d = sqrt(X[0] * X[0] + X[1] * X[1] + X[2] * X[2]);
		if (d <= Tolerance) {
			*Intersection = TRUE;
		}
		else {
			d = 1.0 / d;
			Y[0] = (X[1] * q[2] - X[2] * q[1]) * d;
			Y[1] = (X[2] * q[0] - X[0] * q[2]) * d;
			Y[2] = (X[0] * q[1] - X[1] * q[0]) * d;
			Z[0] = (X[1] * b[2] - X[2] * b[1]) * d;
			Z[1] = (X[2] * b[0] - X[0] * b[2]) * d;
			Z[2] = (X[0] * b[1] - X[1] * b[0]) * d;
			*Intersection = (sqrt(Y[0] * Y[0] + Y[1] * Y[1] + Y[2] * Y[2]
				+ Z[0] * Z[0] + Z[1] * Z[1] + Z[2] * Z[2]) <= Tolerance);
		}
		X[0] = 0.5 * (P[0] + A[0]);
		X[1] = 0.5 * (P[1] + A[1]);
		X[2] = 0.5 * (P[2] + A[2]);
	}
/**/DEBUG_WRITE_LEAVING(LineToLineIntersection3D, "Done")
	return(Status);
} /* end LineToLineIntersection3D */

/*--------------------------------------------------------------------------*/
extern int		LineToPlaneIntersection
				(
					double	P[],				/* 1st point on line */
					double	Q[],				/* 2nd point on line */
					double	A[],				/* 1st point on plane */
					double	B[],				/* 2nd point on plane */
					double	C[],				/* 3rd point on plane */
					double	X[],				/* output intersection */
					long	Lines,				/* height of the vectors */
					double	Tolerance,			/* admissible relative error */
					long	MaxIterations,		/* convergence limit */
					int		*Intersection,		/* resulting validity of intersection */
					int		*Status				/* error management */
				)

/* compute the intersection of the line (P,Q) with the plane (A,B,C) */
/* all vectors have (Lines) elements */
/* return Intersection = FALSE if there is no intersection */
/* return Intersection = TRUE if the intersection is valid */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin LineToPlaneIntersection */

	double	*u, *v, *w, *x, *y;
	double	U[16], V[16], W[4];
	double	Y[4], Z[4];
	double	AB, BC, CA, PQ, PA;
	double	ABPQ, ABCA, ABBC, ABAB;
	double	BCPQ, BCCA, BCBC;
	double	CAPQ, CACA;
	double	PQPQ;
	double	ABPA, BCPA, CAPA, PQPA;
	long	i;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(LineToPlaneIntersection, P, *Status,
/**/	"Empty Line 1st point")
/**/DEBUG_CHECK_NULL_POINTER(LineToPlaneIntersection, Q, *Status,
/**/	"Empty Line 2nd point")
/**/DEBUG_CHECK_NULL_POINTER(LineToPlaneIntersection, A, *Status,
/**/	"Empty Plane 1st point")
/**/DEBUG_CHECK_NULL_POINTER(LineToPlaneIntersection, B, *Status,
/**/	"Empty Plane 2nd point")
/**/DEBUG_CHECK_NULL_POINTER(LineToPlaneIntersection, C, *Status,
/**/	"Empty Plane 3rd point")
/**/DEBUG_CHECK_NULL_POINTER(LineToPlaneIntersection, X, *Status,
/**/	"Empty output")
/**/DEBUG_CHECK_NULL_POINTER(LineToPlaneIntersection, Intersection, *Status,
/**/	"Missing result")
/**/DEBUG_CHECK_RANGE_LONG(LineToPlaneIntersection, Lines, 1L, LONG_MAX, *Status,
/**/	"Invalid Lines (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(LineToPlaneIntersection, MaxIterations, 1L, LONG_MAX, *Status,
/**/	"Invalid iterations (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_DOUBLE(LineToPlaneIntersection, Tolerance, 0.0, DBL_MAX, *Status,
/**/	"Invalid Tolerance (should be positive)")
/**/DEBUG_RETURN_ON_ERROR(LineToPlaneIntersection, *Status)
/**/DEBUG_WRITE_ENTERING(LineToPlaneIntersection,
/**/	"About to determine a (line/plane) intersection")

	u = A;
	v = B;
	w = C;
	x = P;
	y = Q;
	ABAB = 0.0;
	ABBC = 0.0;
	ABCA = 0.0;
	ABPQ = 0.0;
	BCBC = 0.0;
	BCCA = 0.0;
	BCPQ = 0.0;
	CACA = 0.0;
	CAPQ = 0.0;
	PQPQ = 0.0;
	ABPA = 0.0;
	BCPA = 0.0;
	CAPA = 0.0;
	PQPA = 0.0;
	for (i = Lines; (0L < i); i--) {
		AB = *v - *u;
		BC = *w - *v++;
		CA = *u - *w++;
		PQ = *y++ - *x;
		PA = *u++ - *x++;
		ABAB += AB * AB;
		ABBC += AB * BC;
		ABCA += AB * CA;
		ABPQ -= AB * PQ;
		BCBC += BC * BC;
		BCCA += BC * CA;
		BCPQ -= BC * PQ;
		CACA += CA * CA;
		CAPQ -= CA * PQ;
		PQPQ += PQ * PQ;
		ABPA -= AB * PA;
		BCPA -= BC * PA;
		CAPA -= CA * PA;
		PQPA += PQ * PA;
	}
	U[0] = ABAB;
	U[1] = ABBC;
	U[2] = ABCA;
	U[3] = ABPQ;
	U[4] = ABBC;
	U[5] = BCBC;
	U[6] = BCCA;
	U[7] = BCPQ;
	U[8] = ABCA;
	U[9] = BCCA;
	U[10] = CACA;
	U[11] = CAPQ;
	U[12] = ABPQ;
	U[13] = BCPQ;
	U[14] = CAPQ;
	U[15] = PQPQ;
	SingularValueDecomposition(U, 4L, 4L, W, V, MaxIterations, Status);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(LineToPlaneIntersection, "Done")
		return(*Status);
	}
	W[0] = (W[0] < Tolerance) ? (0.0) : (W[0]);
	W[1] = (W[1] < Tolerance) ? (0.0) : (W[1]);
	W[2] = (W[2] < Tolerance) ? (0.0) : (W[2]);
	W[3] = (W[3] < Tolerance) ? (0.0) : (W[3]);
	Y[0] = ABPA;
	Y[1] = BCPA;
	Y[2] = CAPA;
	Y[3] = PQPA;
	SingularValueBackSubstitution(U, W, V, 4L, 4L, Y, Z, Status);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(LineToPlaneIntersection, "Done")
		return(*Status);
	}
	ABPQ = 0.0;
	for (i = Lines; (0L < i); i--) {
		AB = *B - *A;
		BC = *C - *B++;
		CA = *A - *C++;
		PQ = *Q++ - *P;
		AB = *A++ + Z[0] * AB + Z[1] * BC + Z[2] * CA;
		PQ = *P++ + Z[3] * PQ;
		PA = AB - PQ;
		ABPQ += PA * PA;
		*X++ = 0.5 * (AB + PQ);
	}
	*Intersection = (sqrt(ABPQ) <= Tolerance);
/**/DEBUG_WRITE_LEAVING(LineToPlaneIntersection, "Done")
	return(*Status);
} /* end LineToPlaneIntersection */

/*--------------------------------------------------------------------------*/
extern int		LineToPlaneIntersection3D
				(
					double	P[],				/* point on line */
					double	q[],				/* direction of line */
					double	A[],				/* point on plane */
					double	b[],				/* normal of plane */
					double	X[],				/* output intersection */
					double	Tolerance,			/* admissible relative error */
					int		*Intersection		/* resulting validity of intersection */
				)

/* compute the intersection of the line (P,q) with the plane (A,b) */
/* all points have 3 elements */
/* all vectors have 3 elements */
/* a line is described by a point P and a vector q: X = P + t * q, with t a free scalar */
/* a plane is described by a point A and a normal b: <AX, b> = 0 */
/* the vectors (q, b) must have a normalized unit length */
/* return Intersection = FALE if there is no intersection */
/* return Intersection = TRUE if the intersection is valid */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin LineToPlaneIntersection3D */

	double	PAb, qb, t;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(LineToPlaneIntersection3D, P, Status,
/**/	"Empty Line point")
/**/DEBUG_CHECK_NULL_POINTER(LineToPlaneIntersection3D, q, Status,
/**/	"Empty Line direction")
/**/DEBUG_CHECK_NULL_POINTER(LineToPlaneIntersection3D, A, Status,
/**/	"Empty Plane point")
/**/DEBUG_CHECK_NULL_POINTER(LineToPlaneIntersection3D, b, Status,
/**/	"Empty Plane normal")
/**/DEBUG_CHECK_NULL_POINTER(LineToPlaneIntersection3D, X, Status,
/**/	"Empty output")
/**/DEBUG_CHECK_NULL_POINTER(LineToPlaneIntersection3D, Intersection, Status,
/**/	"Missing result")
/**/DEBUG_CHECK_RANGE_DOUBLE(LineToPlaneIntersection3D, Tolerance, 0.0, DBL_MAX, Status,
/**/	"Invalid Tolerance (should be positive)")
/**/DEBUG_RETURN_ON_ERROR(LineToPlaneIntersection3D, Status)
/**/DEBUG_WRITE_ENTERING(LineToPlaneIntersection3D,
/**/	"About to determine a (line/plane) intersection")
#ifdef DEBUG
/**/t = sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]);
/**/if (Tolerance < fabs(t - 1.0)) {
/**/	Status = ERROR;
/**/	WRITE_ERROR(LineToPlaneIntersection3D,
/**/		"Invalid vector q (should have a unit length)")
/**/	DEBUG_WRITE_LEAVING(LineToPlaneIntersection3D, "Done")
/**/	return(Status);
/**/}
/**/t = sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
/**/if (Tolerance < fabs(t - 1.0)) {
/**/	Status = ERROR;
/**/	WRITE_ERROR(LineToPlaneIntersection3D,
/**/		"Invalid vector b (should have a unit length)")
/**/	DEBUG_WRITE_LEAVING(LineToPlaneIntersection3D, "Done")
/**/	return(Status);
/**/}
#endif

	X[0] = A[0] - P[0];
	X[1] = A[1] - P[1];
	X[2] = A[2] - P[2];
	PAb = X[0] * b[0] + X[1] * b[1] + X[2] * b[2];
	qb = q[0] * b[0] + q[1] * b[1] + q[2] * b[2];
	if (fabs(qb) <= Tolerance) {
		*Intersection = (fabs(PAb) <= Tolerance);
		X[0] = A[0];
		X[1] = A[1];
		X[2] = A[2];
	}
	else {
		*Intersection = TRUE;
		t = PAb / qb;
		X[0] = P[0] + t * q[0];
		X[1] = P[1] + t * q[1];
		X[2] = P[2] + t * q[2];
	}
/**/DEBUG_WRITE_LEAVING(LineToPlaneIntersection3D, "Done")
	return(Status);
} /* end LineToPlaneIntersection3D */

/*--------------------------------------------------------------------------*/
extern int		PointsToLine3D
				(
					double	P[],				/* 1st point on input line */
					double	Q[],				/* 2nd point on input line */
					double	A[],				/* point on output line */
					double	b[],				/* direction of output line */
					double	Tolerance			/* admissible relative error */
				)

/* transform the points (P,Q) into the line (A,b) */
/* all points have 3 elements */
/* all vectors have 3 elements */
/* a line is described by a point A and a vector b: X = A + t * b, with t a free scalar */
/* the vector b has a normalized unit length */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin PointsToLine3D */

	double	t;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(PointsToLine3D, P, Status,
/**/	"Empty Line 1st point")
/**/DEBUG_CHECK_NULL_POINTER(PointsToLine3D, Q, Status,
/**/	"Empty Line 2nd point")
/**/DEBUG_CHECK_NULL_POINTER(PointsToLine3D, A, Status,
/**/	"Empty output line point")
/**/DEBUG_CHECK_NULL_POINTER(PointsToLine3D, b, Status,
/**/	"Empty output line direction")
/**/DEBUG_CHECK_RANGE_DOUBLE(PointsToLine3D, Tolerance, 0.0, DBL_MAX, Status,
/**/	"Invalid Tolerance (should be positive)")
/**/DEBUG_RETURN_ON_ERROR(PointsToLine3D, Status)
/**/DEBUG_WRITE_ENTERING(PointsToLine3D,
/**/	"About to transform 2 points into a line")

	A[0] = 0.5 * (P[0] + Q[0]);
	A[1] = 0.5 * (P[1] + Q[1]);
	A[2] = 0.5 * (P[2] + Q[2]);
	b[0] = Q[0] - P[0];
	b[1] = Q[1] - P[1];
	b[2] = Q[2] - P[2];
	t = sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
	if (Tolerance < t) {
		b[0] /= t;
		b[1] /= t;
		b[2] /= t;
	}
	else {
		Status = ERROR;
	}
/**/DEBUG_WRITE_LEAVING(PointsToLine3D, "Done")
	return(Status);
} /* end PointsToLine3D */

/*--------------------------------------------------------------------------*/
extern int		PointsToPlane3D
				(
					double	P[],				/* 1st point on input line */
					double	Q[],				/* 2nd point on input line */
					double	R[],				/* 3rd point on input line */
					double	A[],				/* point on output plane */
					double	b[],				/* normal of output plane */
					double	Tolerance			/* admissible relative error */
				)

/* transform the points (P,Q,R) into the plane (A,b) */
/* all points have 3 elements */
/* all vectors have 3 elements */
/* a plane is described by a point A and a normal b: <AX, b> = 0 */
/* the vectors b has a normalized unit length */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin PointsToPlane3D */

	double	u[3], v[3], w[3];
	double	t;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(PointsToPlane3D, P, Status,
/**/	"Empty Plane 1st point")
/**/DEBUG_CHECK_NULL_POINTER(PointsToPlane3D, Q, Status,
/**/	"Empty Plane 2nd point")
/**/DEBUG_CHECK_NULL_POINTER(PointsToPlane3D, R, Status,
/**/	"Empty Plane 3rd point")
/**/DEBUG_CHECK_NULL_POINTER(PointsToPlane3D, A, Status,
/**/	"Empty output plane point")
/**/DEBUG_CHECK_NULL_POINTER(PointsToPlane3D, b, Status,
/**/	"Empty output plane normal")
/**/DEBUG_CHECK_RANGE_DOUBLE(PointsToPlane3D, Tolerance, 0.0, DBL_MAX, Status,
/**/	"Invalid Tolerance (should be positive)")
/**/DEBUG_RETURN_ON_ERROR(PointsToPlane3D, Status)
/**/DEBUG_WRITE_ENTERING(PointsToPlane3D,
/**/	"About to transform 2 points into a line")

	A[0] = (1.0 / 3.0) * (P[0] + Q[0] + R[0]);
	A[1] = (1.0 / 3.0) * (P[1] + Q[1] + R[1]);
	A[2] = (1.0 / 3.0) * (P[2] + Q[2] + R[2]);
	u[0] = Q[0] - P[0];
	u[1] = Q[1] - P[1];
	u[2] = Q[2] - P[2];
	v[0] = R[0] - Q[0];
	v[1] = R[1] - Q[1];
	v[2] = R[2] - Q[2];
	w[0] = P[0] - R[0];
	w[1] = P[1] - R[1];
	w[2] = P[2] - R[2];
	b[0] = u[1] * v[2] + v[1] * w[2] + w[1] * u[2] - u[2] * v[1] - v[2] * w[1] - w[2] * u[1];
	b[1] = u[2] * v[0] + v[2] * w[0] + w[2] * u[0] - u[0] * v[2] - v[0] * w[2] - w[0] * u[2];
	b[2] = u[0] * v[1] + v[0] * w[1] + w[0] * u[1] - u[1] * v[0] - v[1] * w[0] - w[1] * u[0];
	t = sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
	if (Tolerance < t) {
		b[0] /= t;
		b[1] /= t;
		b[2] /= t;
	}
	else {
		Status = ERROR;
	}
/**/DEBUG_WRITE_LEAVING(PointsToPlane3D, "Done")
	return(Status);
} /* end PointsToPlane3D */

/*--------------------------------------------------------------------------*/
extern int		PointToLineMembership
				(
					double	P[],				/* point to test */
					double	A[],				/* 1st point on line */
					double	B[],				/* 2nd point on line */
					long	Lines,				/* height of the vectors */
					double	Tolerance,			/* admissible relative error */
					int		*Membership,		/* resulting type of intersection */
					int		*Status				/* error management */
				)

/* test the membership of the point P to the line (A,B) */
/* all vectors have (Lines) elements */
/* return Membership = FALSE if P doesn't belong to (A,B) */
/* return Membership = TRUE if P does belong to (A,B) */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin PointToLineMembership */

	double	*u;
	double	*U = (double *)NULL;
	double	Diff, Norm;
	long	i;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(PointToLineMembership, P, *Status,
/**/	"Empty point")
/**/DEBUG_CHECK_NULL_POINTER(PointToLineMembership, A, *Status,
/**/	"Empty line 1st point")
/**/DEBUG_CHECK_NULL_POINTER(PointToLineMembership, B, *Status,
/**/	"Empty line 2nd point")
/**/DEBUG_CHECK_NULL_POINTER(PointToLineMembership, Membership, *Status,
/**/	"Missing result")
/**/DEBUG_CHECK_RANGE_LONG(PointToLineMembership, Lines, 1L, LONG_MAX, *Status,
/**/	"Invalid Lines (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_DOUBLE(PointToLineMembership, Tolerance, 0.0, DBL_MAX, *Status,
/**/	"Invalid Tolerance (should be positive)")
/**/DEBUG_RETURN_ON_ERROR(PointToLineMembership, *Status)
/**/DEBUG_WRITE_ENTERING(PointToLineMembership,
/**/	"About to determine the membership of a point to a line")

	AllocateVector(&U, Lines, Status);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(PointToLineMembership, "Done")
		return(*Status);
	}
	*Status = ProjectPointToLine(P, A, B, U, Lines, Tolerance);
	if (*Status == ERROR) {
		FreeVector(&U);
/**/	DEBUG_WRITE_LEAVING(PointToLineMembership, "Done")
		return(*Status);
	}
	u = U;
	Norm = 0.0;
	for (i = Lines; (0L < i); i--) {
		Diff = *P++ - *u++;
		Norm += Diff * Diff;
	}
	*Membership = (sqrt(Norm) <= Tolerance);
	*Status = FreeVector(&U);
/**/DEBUG_WRITE_LEAVING(PointToLineMembership, "Done")
	return(*Status);
} /* end PointToLineMembership */

/*--------------------------------------------------------------------------*/
extern int		PointToLineMembership3D
				(
					double	P[],				/* point to test */
					double	A[],				/* point on line */
					double	b[],				/* direction of line */
					double	Tolerance,			/* admissible relative error */
					int		*Membership			/* resulting type of intersection */
				)

/* test the membership of the point P to the line (A,b) */
/* all points have 3 elements */
/* all vectors have 3 elements */
/* a line is described by a point A and a vector b: X = A + t * b, with t a free scalar */
/* the vector b must have a normalized unit length */
/* return Membership = FALSE if P doesn't belong to (A,b) */
/* return Membership = TRUE if P does belong to (A,b) */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin PointToLineMembership3D */

	double	v[3], w[3];
	double	t;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(PointToLineMembership3D, P, Status,
/**/	"Empty point")
/**/DEBUG_CHECK_NULL_POINTER(PointToLineMembership3D, A, Status,
/**/	"Empty line point")
/**/DEBUG_CHECK_NULL_POINTER(PointToLineMembership3D, b, Status,
/**/	"Empty line direction")
/**/DEBUG_CHECK_NULL_POINTER(PointToLineMembership3D, Membership, Status,
/**/	"Missing result")
/**/DEBUG_CHECK_RANGE_DOUBLE(PointToLineMembership3D, Tolerance, 0.0, DBL_MAX, Status,
/**/	"Invalid Tolerance (should be positive)")
/**/DEBUG_RETURN_ON_ERROR(PointToLineMembership3D, Status)
/**/DEBUG_WRITE_ENTERING(PointToLineMembership3D,
/**/	"About to determine the membership of a point to a line")
#ifdef DEBUG
/**/t = sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
/**/if (Tolerance < fabs(t - 1.0)) {
/**/	Status = ERROR;
/**/	WRITE_ERROR(PointToLineMembership3D,
/**/		"Invalid vector b (should have a unit length)")
/**/	DEBUG_WRITE_LEAVING(PointToLineMembership3D, "Done")
/**/	return(Status);
/**/}
#endif

	v[0] = A[0] - P[0];
	v[1] = A[1] - P[1];
	v[2] = A[2] - P[2];
	t = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
	if (Tolerance < t) {
		v[0] /= t;
		v[1] /= t;
		v[2] /= t;
		w[0] = v[1] * b[2] - v[2] * b[1];
		w[1] = v[2] * b[0] - v[0] * b[2];
		w[2] = v[0] * b[1] - v[1] * b[0];
		t = sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
		*Membership = (t <= Tolerance);
	}
	else {
		*Membership = TRUE;
	}
/**/DEBUG_WRITE_LEAVING(PointToLineMembership3D, "Done")
	return(Status);
} /* end PointToLineMembership3D */

/*--------------------------------------------------------------------------*/
extern int		PointToPlaneMembership
				(
					double	P[],				/* point to test */
					double	A[],				/* 1st point on plane */
					double	B[],				/* 2nd point on plane */
					double	C[],				/* 2nd point on plane */
					long	Lines,				/* height of the vectors */
					double	Tolerance,			/* admissible relative error */
					long	MaxIterations,		/* convergence limit */
					int		*Membership,		/* resulting type of intersection */
					int		*Status				/* error management */
				)

/* test the membership of the point P to the plane (A,B,C) */
/* all vectors have (Lines) elements */
/* return Membership = FALSE if P doesn't belong to (A,B,C) */
/* return Membership = TRUE if P does belong to (A,B,C) */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin PointToPlaneMembership */

	double	*u;
	double	*U = (double *)NULL;
	double	Diff, Norm;
	long	i;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(PointToPlaneMembership, P, *Status,
/**/	"Empty point")
/**/DEBUG_CHECK_NULL_POINTER(PointToPlaneMembership, A, *Status,
/**/	"Empty plane 1st point")
/**/DEBUG_CHECK_NULL_POINTER(PointToPlaneMembership, B, *Status,
/**/	"Empty plane 2nd point")
/**/DEBUG_CHECK_NULL_POINTER(PointToPlaneMembership, C, *Status,
/**/	"Empty plane 3rd point")
/**/DEBUG_CHECK_NULL_POINTER(PointToPlaneMembership, Membership, *Status,
/**/	"Missing result")
/**/DEBUG_CHECK_RANGE_LONG(PointToPlaneMembership, Lines, 1L, LONG_MAX, *Status,
/**/	"Invalid Lines (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(PointToPlaneMembership, MaxIterations, 1L, LONG_MAX, *Status,
/**/	"Invalid iterations (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_DOUBLE(PointToPlaneMembership, Tolerance, 0.0, DBL_MAX, *Status,
/**/	"Invalid Tolerance (should be positive)")
/**/DEBUG_RETURN_ON_ERROR(PointToPlaneMembership, *Status)
/**/DEBUG_WRITE_ENTERING(PointToPlaneMembership,
/**/	"About to determine the membership of a point to a plane")

	AllocateVector(&U, Lines, Status);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(PointToPlaneMembership, "Done")
		return(*Status);
	}
	ProjectPointToPlane(P, A, B, C, U, Lines, Tolerance, MaxIterations, Status);
	if (*Status == ERROR) {
		FreeVector(&U);
/**/	DEBUG_WRITE_LEAVING(PointToPlaneMembership, "Done")
		return(*Status);
	}
	u = U;
	Norm = 0.0;
	for (i = Lines; (0L < i); i--) {
		Diff = *P++ - *u++;
		Norm += Diff * Diff;
	}
	*Membership = (sqrt(Norm) <= Tolerance);
	*Status = FreeVector(&U);
/**/DEBUG_WRITE_LEAVING(PointToPlaneMembership, "Done")
	return(*Status);
} /* end PointToPlaneMembership */

/*--------------------------------------------------------------------------*/
extern int		PointToPlaneMembership3D
				(
					double	P[],				/* point to test */
					double	A[],				/* point on plane */
					double	b[],				/* normal of plane */
					double	Tolerance,			/* admissible relative error */
					int		*Membership			/* resulting type of intersection */
				)

/* test the membership of the point P to the plane (A,b) */
/* all points have 3 elements */
/* all vectors have 3 elements */
/* a plane is described by a point A and a normal b: <AX, b> = 0 */
/* the vectors b must have a normalized unit length */
/* return Membership = FALSE if P doesn't belong to (A,b) */
/* return Membership = TRUE if P does belong to (A,b) */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin PointToPlaneMembership3D */

	double	u[3];
	double	t;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(PointToPlaneMembership3D, P, Status,
/**/	"Empty point")
/**/DEBUG_CHECK_NULL_POINTER(PointToPlaneMembership3D, A, Status,
/**/	"Empty plane point")
/**/DEBUG_CHECK_NULL_POINTER(PointToPlaneMembership3D, b, Status,
/**/	"Empty plane normal")
/**/DEBUG_CHECK_NULL_POINTER(PointToPlaneMembership3D, Membership, Status,
/**/	"Missing result")
/**/DEBUG_CHECK_RANGE_DOUBLE(PointToPlaneMembership3D, Tolerance, 0.0, DBL_MAX, Status,
/**/	"Invalid Tolerance (should be positive)")
/**/DEBUG_RETURN_ON_ERROR(PointToPlaneMembership3D, Status)
/**/DEBUG_WRITE_ENTERING(PointToPlaneMembership3D,
/**/	"About to determine the membership of a point to a plane")
#ifdef DEBUG
/**/t = sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
/**/if (Tolerance < fabs(t - 1.0)) {
/**/	Status = ERROR;
/**/	WRITE_ERROR(PointToPlaneMembership3D,
/**/		"Invalid vector b (should have a unit length)")
/**/	DEBUG_WRITE_LEAVING(PointToPlaneMembership3D, "Done")
/**/	return(Status);
/**/}
#endif

	u[0] = A[0] - P[0];
	u[1] = A[1] - P[1];
	u[2] = A[2] - P[2];
	t = u[0] * b[0] + u[1] * b[1] + u[2] * b[2];
	*Membership = (fabs(t) <= Tolerance);
/**/DEBUG_WRITE_LEAVING(PointToPlaneMembership3D, "Done")
	return(Status);
} /* end PointToPlaneMembership3D */

/*--------------------------------------------------------------------------*/
extern int		ProjectPointToLine
				(
					double	P[],				/* point to project */
					double	A[],				/* 1st point on line */
					double	B[],				/* 2nd point on line */
					double	X[],				/* resulting projection */
					long	Lines,				/* height of the vectors */
					double	Tolerance			/* admissible relative error */
				)

/* project the point P onto the line (A,B) */
/* all vectors have (Lines) elements */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin ProjectPointToLine */

	double	*u, *v, *w;
	double	ABAB, APAB;
	double	AB, AP;
	long	i;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(ProjectPointToLine, P, Status,
/**/	"Empty Point")
/**/DEBUG_CHECK_NULL_POINTER(ProjectPointToLine, A, Status,
/**/	"Empty Plane 1st point")
/**/DEBUG_CHECK_NULL_POINTER(ProjectPointToLine, B, Status,
/**/	"Empty Plane 2nd point")
/**/DEBUG_CHECK_NULL_POINTER(ProjectPointToLine, X, Status,
/**/	"Empty output")
/**/DEBUG_CHECK_RANGE_LONG(ProjectPointToLine, Lines, 1L, LONG_MAX, Status,
/**/	"Invalid Lines (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_DOUBLE(ProjectPointToLine, Tolerance, 0.0, DBL_MAX, Status,
/**/	"Invalid Tolerance (should be positive)")
/**/DEBUG_RETURN_ON_ERROR(ProjectPointToLine, Status)
/**/DEBUG_WRITE_ENTERING(ProjectPointToLine,
/**/	"About to project a point onto a plane")

	u = P;
	v = A;
	w = B;
	APAB = 0.0;
	ABAB = 0.0;
	for (i = Lines; (0L < i); i--) {
		AP = *u++ - *v;
		AB = *w++ - *v++;
		APAB += AB * AP;
		ABAB += AB * AB;
	}
	if (ABAB < Tolerance) {
		for (i = Lines; (0L < i); i--) {
			*X++ = 0.5 * (*A++ + *B++);
		}
	}
	else {
		APAB /= ABAB;
		for (i = Lines; (0L < i); i--) {
			*X++ = *A + APAB * (*B++ - *A);
			A++;
		}
	}
/**/DEBUG_WRITE_LEAVING(ProjectPointToLine, "Done")
	return(Status);
} /* end ProjectPointToLine */

/*--------------------------------------------------------------------------*/
extern int		ProjectPointToLine3D
				(
					double	P[],				/* point to project */
					double	A[],				/* point on line */
					double	b[],				/* direction of line */
					double	X[],				/* resulting projection */
					double	Tolerance			/* admissible relative error */
				)

/* project the point P onto the line (A,b) */
/* all points have 3 elements */
/* all vectors have 3 elements */
/* a line is described by a point A and a vector b: X = A + t * b, with t a free scalar */
/* the vector b must have a normalized unit length */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin ProjectPointToLine3D */

	double	t;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(ProjectPointToLine3D, P, Status,
/**/	"Empty Point")
/**/DEBUG_CHECK_NULL_POINTER(ProjectPointToLine3D, A, Status,
/**/	"Empty Plane point")
/**/DEBUG_CHECK_NULL_POINTER(ProjectPointToLine3D, b, Status,
/**/	"Empty Plane normal")
/**/DEBUG_CHECK_NULL_POINTER(ProjectPointToLine3D, X, Status,
/**/	"Empty output")
/**/DEBUG_CHECK_RANGE_DOUBLE(ProjectPointToLine3D, Tolerance, 0.0, DBL_MAX, Status,
/**/	"Invalid Tolerance (should be positive)")
/**/DEBUG_RETURN_ON_ERROR(ProjectPointToLine3D, Status)
/**/DEBUG_WRITE_ENTERING(ProjectPointToLine3D,
/**/	"About to project a point onto a plane")
#ifdef DEBUG
/**/t = sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
/**/if (Tolerance < fabs(t - 1.0)) {
/**/	Status = ERROR;
/**/	WRITE_ERROR(ProjectPointToLine3D,
/**/		"Invalid vector b (should have a unit length)")
/**/	DEBUG_WRITE_LEAVING(ProjectPointToLine3D, "Done")
/**/	return(Status);
/**/}
#endif

	X[0] = A[0] - P[0];
	X[1] = A[1] - P[1];
	X[2] = A[2] - P[2];
	t = X[0] * b[0] + X[1] * b[1] + X[2] * b[2];
	X[0] = P[0] + t * b[0];
	X[1] = P[1] + t * b[1];
	X[2] = P[2] + t * b[2];
/**/DEBUG_WRITE_LEAVING(ProjectPointToLine3D, "Done")
	return(Status);
} /* end ProjectPointToLine3D */

/*--------------------------------------------------------------------------*/
extern int		ProjectPointToPlane
				(
					double	P[],				/* point to project */
					double	A[],				/* 1st point on plane */
					double	B[],				/* 2nd point on plane */
					double	C[],				/* 3rd point on plane */
					double	X[],				/* resulting projection */
					long	Lines,				/* height of the vectors */
					double	Tolerance,			/* admissible relative error */
					long	MaxIterations,		/* convergence limit */
					int		*Status				/* error management */
				)

/* project the point P onto the plane (A,B,C) */
/* all vectors have (Lines) elements */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin ProjectPointToPlane */

	double	*u, *v, *w, *x;
	double	U[9], V[9], W[3], Y[3], Z[3];
	double	APAB, ABAB;
	double	APBC, ABBC, BCBC;
	double	APCA, ABCA, BCCA, CACA;
	double	AP, AB, BC, CA;
	long	i;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(ProjectPointToPlane, P, *Status,
/**/	"Empty Point")
/**/DEBUG_CHECK_NULL_POINTER(ProjectPointToPlane, A, *Status,
/**/	"Empty Plane 1st point")
/**/DEBUG_CHECK_NULL_POINTER(ProjectPointToPlane, B, *Status,
/**/	"Empty Plane 2nd point")
/**/DEBUG_CHECK_NULL_POINTER(ProjectPointToPlane, X, *Status,
/**/	"Empty output")
/**/DEBUG_CHECK_RANGE_LONG(ProjectPointToPlane, Lines, 1L, LONG_MAX, *Status,
/**/	"Invalid Lines (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(ProjectPointToPlane, MaxIterations, 1L, LONG_MAX, *Status,
/**/	"Invalid iterations (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_DOUBLE(ProjectPointToPlane, Tolerance, 0.0, DBL_MAX, *Status,
/**/	"Invalid Tolerance (should be positive)")
/**/DEBUG_RETURN_ON_ERROR(ProjectPointToPlane, *Status)
/**/DEBUG_WRITE_ENTERING(ProjectPointToPlane,
/**/	"About to project a point onto a plane")

	u = P;
	v = A;
	w = B;
	x = C;
	APAB = 0.0;
	APBC = 0.0;
	APCA = 0.0;
	ABAB = 0.0;
	ABBC = 0.0;
	BCBC = 0.0;
	ABCA = 0.0;
	BCCA = 0.0;
	CACA = 0.0;
	for (i = Lines; (0L < i); i--) {
		AP = *u++ - *v;
		AB = *w - *v;
		BC = *x - *w++;
		CA = *v++ - *x++;
		APAB += AP * AB;
		APBC += AP * BC;
		APCA += AP * CA;
		ABAB += AB * AB;
		ABBC += AB * BC;
		BCBC += BC * BC;
		ABCA += AB * CA;
		BCCA += BC * CA;
		CACA += CA * CA;
	}
	U[0] = ABAB;
	U[1] = ABBC;
	U[2] = ABCA;
	U[3] = ABBC;
	U[4] = BCBC;
	U[5] = BCCA;
	U[6] = ABCA;
	U[7] = BCCA;
	U[8] = CACA;
	SingularValueDecomposition(U, 3L, 3L, W, V, MaxIterations, Status);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(ProjectPointToPlane, "Done")
		return(*Status);
	}
	W[0] = (W[0] < Tolerance) ? (0.0) : (W[0]);
	W[1] = (W[1] < Tolerance) ? (0.0) : (W[1]);
	W[2] = (W[2] < Tolerance) ? (0.0) : (W[2]);
	Y[0] = APAB;
	Y[1] = APBC;
	Y[2] = APCA;
	SingularValueBackSubstitution(U, W, V, 3L, 3L, Y, Z, Status);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(ProjectPointToPlane, "Done")
		return(*Status);
	}
	for (i = Lines; (0L < i); i--) {
		*X++ = *A + Z[0] * (*B - *A) + Z[1] * (*C - *B) + Z[2] * (*A - *C);
		A++;
		B++;
		C++;
	}
/**/DEBUG_WRITE_LEAVING(ProjectPointToPlane, "Done")
	return(*Status);
} /* end ProjectPointToPlane */

/*--------------------------------------------------------------------------*/
extern int		ProjectPointToPlane3D
				(
					double	P[],				/* point to project */
					double	A[],				/* point on plane */
					double	b[],				/* normal of plane */
					double	X[],				/* resulting projection */
					double	Tolerance			/* admissible relative error */
				)

/* project the point P onto the plane (A,b) */
/* all points have 3 elements */
/* all vectors have 3 elements */
/* a plane is described by a point A and a normal b: <AX, b> = 0 */
/* the vectors b must have a normalized unit length */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin ProjectPointToPlane3D */

	double	t;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(ProjectPointToPlane3D, P, Status,
/**/	"Empty Point")
/**/DEBUG_CHECK_NULL_POINTER(ProjectPointToPlane3D, A, Status,
/**/	"Empty Plane point")
/**/DEBUG_CHECK_NULL_POINTER(ProjectPointToPlane3D, b, Status,
/**/	"Empty Plane normal")
/**/DEBUG_CHECK_NULL_POINTER(ProjectPointToPlane3D, X, Status,
/**/	"Empty output")
/**/DEBUG_CHECK_RANGE_DOUBLE(ProjectPointToPlane3D, Tolerance, 0.0, DBL_MAX, Status,
/**/	"Invalid Tolerance (should be positive)")
/**/DEBUG_RETURN_ON_ERROR(ProjectPointToPlane3D, Status)
/**/DEBUG_WRITE_ENTERING(ProjectPointToPlane3D,
/**/	"About to project a point onto a plane")
#ifdef DEBUG
/**/t = sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
/**/if (Tolerance < fabs(t - 1.0)) {
/**/	Status = ERROR;
/**/	WRITE_ERROR(ProjectPointToPlane3D,
/**/		"Invalid vector b (should have a unit length)")
/**/	DEBUG_WRITE_LEAVING(ProjectPointToPlane3D, "Done")
/**/	return(Status);
/**/}
#endif

	X[0] = A[0] - P[0];
	X[1] = A[1] - P[1];
	X[2] = A[2] - P[2];
	t = X[0] * b[0] + X[1] * b[1] + X[2] * b[2];
	X[0] = A[0] - t * b[0];
	X[1] = A[1] - t * b[1];
	X[2] = A[2] - t * b[2];
/**/DEBUG_WRITE_LEAVING(ProjectPointToPlane3D, "Done")
	return(Status);
} /* end ProjectPointToPlane3D */

/*--------------------------------------------------------------------------*/
extern int		TestColinearVector
				(
					double	U[],				/* 1st vector */
					double	V[],				/* 2nd vector */
					long	Lines,				/* height of the vector */
					double	Tolerance,			/* admissible relative error */
					int		*Colinear			/* resulting test of colinearity */
				)

/* determine whether the vectors U and V are colinear */
/* the vectors U and V have (Lines) elements */
/* return Colinear = -1 if at least one of the vectors is degenerate */
/* return Colinear = 0 if the vectors are not colinear */
/* return Colinear = 1 if the vectors are colinear */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin TestColinearVector */

	double	*u, *v;
	double	NormU, NormV;
	double	Diff, Norm;
	long	i;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(TestColinearVector, U, Status,
/**/	"Empty input left operand")
/**/DEBUG_CHECK_NULL_POINTER(TestColinearVector, V, Status,
/**/	"Empty input right operand")
/**/DEBUG_CHECK_NULL_POINTER(TestColinearVector, Colinear, Status,
/**/	"Missing result")
/**/DEBUG_CHECK_RANGE_LONG(TestColinearVector, Lines, 1L, LONG_MAX, Status,
/**/	"Invalid Lines (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_DOUBLE(TestColinearVector, Tolerance, 0.0, DBL_MAX, Status,
/**/	"Invalid Tolerance (should be positive)")
/**/DEBUG_RETURN_ON_ERROR(TestColinearVector, Status)
/**/DEBUG_WRITE_ENTERING(TestColinearVector,
/**/	"About to test for colinearity")

	u = U;
	v = V;
	NormU = 0.0;
	NormV = 0.0;
	for (i = Lines; (0L < i); i--) {
		NormU += *u * *u;
		NormV += *v * *v;
		u++;
		v++;
	}
	NormU = sqrt(NormU);
	NormV = sqrt(NormV);
	if (NormU <= Tolerance) {
		*Colinear = -1;
/**/	DEBUG_WRITE_LEAVING(TestColinearVector, "Done")
		return(Status);
	}
	if (NormV <= Tolerance) {
		*Colinear = -1;
/**/	DEBUG_WRITE_LEAVING(TestColinearVector, "Done")
		return(Status);
	}
	NormU = 1.0 / NormU;
	NormV = 1.0 / NormV;
	Norm = 0.0;
	for (i = Lines; (0L < i); i--) {
		Diff = *U++ * NormU - *V++ * NormV;
		Norm += Diff * Diff;
	}
	*Colinear = (sqrt(Norm) < Tolerance) ? (1) : (0);
/**/DEBUG_WRITE_LEAVING(TestColinearVector, "Done")
	return(Status);
} /* end TestColinearVector */

/*--------------------------------------------------------------------------*/
extern int		TestCoplanarVector
				(
					double	U[],				/* 1st vector */
					double	V[],				/* 2nd vector */
					double	W[],				/* 3rd vector */
					long	Lines,				/* height of the vectors */
					double	Tolerance,			/* admissible relative error */
					long	MaxIterations,		/* convergence limit */
					int		*Coplanar,			/* resulting test of coplanarity */
					int		*Status				/* error management */
				)

/* determine whether the vectors (U,V,W) are coplanar */
/* the vectors (U,V,W) have (Lines) elements */
/* return Coplanar = -1 in degenerate cases */
/* return Coplanar = 0 if the vectors are not coplanar */
/* return Coplanar = 1 if the vectors are coplanar */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin TestCoplanarVector */

	double	*u, *v, *w, *x;
	double	*M = (double *)NULL;
	double	X[9], Y[3];
	long	i;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(TestCoplanarVector, U, *Status,
/**/	"Empty 1st vector")
/**/DEBUG_CHECK_NULL_POINTER(TestCoplanarVector, V, *Status,
/**/	"Empty 2nd vector")
/**/DEBUG_CHECK_NULL_POINTER(TestCoplanarVector, W, *Status,
/**/	"Empty 3rd vector")
/**/DEBUG_CHECK_NULL_POINTER(TestCoplanarVector, Coplanar, *Status,
/**/	"Missing result")
/**/DEBUG_CHECK_RANGE_LONG(TestCoplanarVector, Lines, 1L, LONG_MAX, *Status,
/**/	"Invalid Lines (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(TestCoplanarVector, MaxIterations, 1L, LONG_MAX, *Status,
/**/	"Invalid iterations (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_DOUBLE(TestCoplanarVector, Tolerance, 0.0, DBL_MAX, *Status,
/**/	"Invalid Tolerance (should be positive)")
/**/DEBUG_RETURN_ON_ERROR(TestCoplanarVector, *Status)
/**/DEBUG_WRITE_ENTERING(TestCoplanarVector,
/**/	"About to determine a (line/line) intersection")

	*Status = TestColinearVector(U, V, Lines, Tolerance, Coplanar);
	switch (*Coplanar) {
		case 1:
			*Coplanar = -1;
		case -1:
/**/		DEBUG_WRITE_LEAVING(TestCoplanarVector, "Done")
			return(*Status);
		case 0:
			break;
		default:
			*Status = ERROR;
			WRITE_ERROR(TestCoplanarVector, "Unexpected internal error")
/**/		DEBUG_WRITE_LEAVING(TestCoplanarVector, "Done")
			return(*Status);
	}
	*Status = TestColinearVector(V, W, Lines, Tolerance, Coplanar);
	switch (*Coplanar) {
		case 1:
			*Coplanar = -1;
		case -1:
/**/		DEBUG_WRITE_LEAVING(TestCoplanarVector, "Done")
			return(*Status);
		case 0:
			break;
		default:
			*Status = ERROR;
			WRITE_ERROR(TestCoplanarVector, "Unexpected internal error")
/**/		DEBUG_WRITE_LEAVING(TestCoplanarVector, "Done")
			return(*Status);
	}
	*Status = TestColinearVector(W, U, Lines, Tolerance, Coplanar);
	switch (*Coplanar) {
		case 1:
			*Coplanar = -1;
		case -1:
/**/		DEBUG_WRITE_LEAVING(TestCoplanarVector, "Done")
			return(*Status);
		case 0:
			break;
		default:
			*Status = ERROR;
			WRITE_ERROR(TestCoplanarVector, "Unexpected internal error")
/**/		DEBUG_WRITE_LEAVING(TestCoplanarVector, "Done")
			return(*Status);
	}
	AllocateMatrix(&M, Lines, 3L, Status);
	if (*Status == ERROR) {
/**/	DEBUG_WRITE_LEAVING(TestCoplanarVector, "Done")
		return(*Status);
	}
	u = U;
	v = V;
	w = W;
	x = M;
	for (i = Lines; (0L < i); i--) {
		*x++ = *u++;
		*x++ = *v++;
		*x++ = *w++;
	}
	SingularValueDecomposition(M, Lines, 3L, Y, X, MaxIterations, Status);
	if (*Status == ERROR) {
		FreeMatrix(&M);
/**/	DEBUG_WRITE_LEAVING(TestCoplanarVector, "Done")
		return(*Status);
	}
	if (((Y[0] <= Tolerance) && (Tolerance < Y[1]) && (Tolerance < Y[2]))
		|| ((Tolerance < Y[0]) && (Y[1] <= Tolerance) && (Tolerance < Y[2]))
		|| ((Tolerance < Y[0]) && (Tolerance < Y[1]) && (Y[2] <= Tolerance))) {
		*Coplanar = 1;
	}
	else if (((Tolerance < Y[0]) && (Tolerance < Y[1]) && (Tolerance < Y[2]))
		|| ((Tolerance < Y[0]) && (Tolerance < Y[1]) && (Tolerance < Y[2]))
		|| ((Tolerance < Y[0]) && (Tolerance < Y[1]) && (Tolerance < Y[2]))) {
		*Coplanar = 0;
	}
	else {
		WRITE_WARNING(TestCoplanarVector, "Unexpected condition")
		*Coplanar = -1;
	}
	*Status = FreeMatrix(&M);
/**/DEBUG_WRITE_LEAVING(TestCoplanarVector, "Done")
	return(*Status);
} /* end TestCoplanarVector */

