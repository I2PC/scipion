/**@name Geometry */
//@{
/*--------------------------------------------------------------------------*/
/** Affine -> Rotation matrix.
    Approximate a (4 x 4) homogenous affine matrix by a (4 x 4) rotation
    matrix. The rotation matrix has the form R = Rx.Ry.Rz.
    
    A and R are homogenous: (A, R) = {{*, *, *, 0}, {*, *, *, 0},
       {*, *, *, 0}, {0, 0, 0, 1}}
    
    The returned rotation angles are given in the unit of radian.
    
    success: return(!ERROR); failure: return(ERROR)
*/
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
				);

/*--------------------------------------------------------------------------*/
/** Produce Euler rotation matrix (X,Y,Z).
    Fill a rotation matrix R such that R = Rx.Ry.Rz.
    The size of the output matrix R is (4 x 4). R is homogenous:
    R = {{*, *, *, 0}, {*, *, *, 0}, {*, *, *, 0}, {0, 0, 0, 1}}.
    The rotation angles are given in the unit of radian.
    
    success: return(!ERROR); failure: return(ERROR)
*/
extern int		GetRotationMatrix
				(
					double	*R,					/* output 4x4 homogenous matrix */
					double	xRotation,			/* rotation around x-axis */
					double	yRotation,			/* rotation around y-axis */
					double	zRotation			/* rotation around z-axis */
				);

/*--------------------------------------------------------------------------*/
/** Produce scaling matrix.
    Fill a scaling matrix T. The size of the output matrix T is (4 x 4).
    T is homogenous:
    T = {{sx, 0, 0, 0}, {0, sy, 0, 0}, {0, 0, sz, 0}, {0, 0, 0, 1}}
    
    success: return(!ERROR); failure: return(ERROR)
*/
extern int		GetScalingMatrix
				(
					double	*S,					/* output 4x4 homogenous matrix */
					double	xScale,				/* scaling along x-axis */
					double	yScale,				/* scaling along y-axis */
					double	zScale				/* scaling along z-axis */
				);

/*--------------------------------------------------------------------------*/
/** Produce a translation matrix.
    Fill a translation matrix T. The size of the output matrix T is (4 x 4)
    T is homogenous:
    T = {{1, 0, 0, dx}, {0, 1, 0, dy}, {0, 0, 1, dz}, {0, 0, 0, 1}}
    
    success: return(!ERROR); failure: return(ERROR)
*/
extern int		GetTranslationMatrix
				(
					double	*T,					/* output 4x4 homogenous matrix */
					double	xTranslation,		/* translation along x-axis */
					double	yTranslation,		/* translation along y-axis */
					double	zTranslation		/* translation along z-axis */
				);

/*--------------------------------------------------------------------------*/
/** Do two lines intersect (n-D)?.
    Compute the intersection of the line (P,Q) with the line (A,B).
    All vectors have (Lines) elements. I.e., the dimension of the space
    is (Lines).
    
    return Intersection = FALSE if there is no intersection.
    return Intersection = TRUE if the intersection is valid.
    
    success: return(!ERROR); failure: return(ERROR)
*/
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
				);

/*--------------------------------------------------------------------------*/
/** Do two lines intersect (3-D)?
    Compute the intersection of the line (P,q) with the line (A,b).
    All points have 3 elements. All vectors have 3 elements.
    
    A line is described by a point and a vector: X = P + t1 * q,
    X = A + t2 * b where (t1,t2) are free scalars.
    
    The vectors (q, b) must have a normalized unit length.
    
    return Intersection = FALSE if there is no intersection.
    return Intersection = TRUE if the intersection is valid.
    
    success: return(!ERROR); failure: return(ERROR)
*/
extern int		LineToLineIntersection3D
				(
					double	P[],				/* point on line1 */
					double	q[],				/* direction of line1 */
					double	A[],				/* point on line2 */
					double	b[],				/* direction of line2 */
					double	X[],				/* output intersection */
					double	Tolerance,			/* admissible relative error */
					int		*Intersection		/* resulting validity of intersection */
				);

/*--------------------------------------------------------------------------*/
/** Does a line intersect a plane (n-D)?
    Compute the intersection of the line (P,Q) with the plane (A,B,C).
    All vectors have (Lines) elements. I.e., the dimension of the space
    is (Lines).
    
    return Intersection = FALSE if there is no intersection.
    return Intersection = TRUE if the intersection is valid.
    
    success: return(!ERROR); failure: return(ERROR)
*/
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
				);

/*--------------------------------------------------------------------------*/
/** Does a line intersect a plane (3-D)?
    Compute the intersection of the line (P,q) with the plane (A,b).
    All points have 3 elements. All vectors have 3 elements.
    
    A line is described by a point P and a vector q: X = P + t * q,
    with t a free scalar. A plane is described by a point A and a normal b:
    <AX, b> = 0. The vectors (q, b) must have a normalized unit length.
    
    return Intersection = FALE if there is no intersection.
    return Intersection = TRUE if the intersection is valid.
    
    success: return(!ERROR); failure: return(ERROR)
*/
extern int		LineToPlaneIntersection3D
				(
					double	P[],				/* point on line */
					double	q[],				/* direction of line */
					double	A[],				/* point on plane */
					double	b[],				/* normal of plane */
					double	X[],				/* output intersection */
					double	Tolerance,			/* admissible relative error */
					int		*Intersection		/* resulting validity of intersection */
				);

/*--------------------------------------------------------------------------*/
/** What is the line that passes through two points (3-D)?.
    Transform the points (P,Q) into the line (A,b).
    All points have 3 elements. All vectors have 3 elements.
    
    A line is described by a point A and a vector b: X = A + t * b,
    with t a free scalar. The vector b has a normalized unit length.
    
    success: return(!ERROR); failure: return(ERROR)
*/
extern int		PointsToLine3D
				(
					double	P[],				/* 1st point on input line */
					double	Q[],				/* 2nd point on input line */
					double	A[],				/* point on output line */
					double	b[],				/* direction of output line */
					double	Tolerance			/* admissible relative error */
				);

/*--------------------------------------------------------------------------*/
/** What is the plane passing through 3 points (3-D)?
    Transform the points (P,Q,R) into the plane (A,b).
    All points have 3 elements. All vectors have 3 elements.
    
    A plane is described by a point A and a normal b: <AX, b> = 0.
    The vectors b has a normalized unit length.
    
    success: return(!ERROR); failure: return(ERROR)
*/
extern int		PointsToPlane3D
				(
					double	P[],				/* 1st point on input line */
					double	Q[],				/* 2nd point on input line */
					double	R[],				/* 3rd point on input line */
					double	A[],				/* point on output plane */
					double	b[],				/* normal of output plane */
					double	Tolerance			/* admissible relative error */
				);

/*--------------------------------------------------------------------------*/
/** Does a point belong to a line (n-D)?.
    Test the membership of the point P to the line (A,B).
    All vectors have (Lines) elements.
    
    return Membership = FALSE if P doesn't belong to (A,B).
    return Membership = TRUE if P does belong to (A,B).
    
    success: return(!ERROR); failure: return(ERROR)
*/
extern int		PointToLineMembership
				(
					double	P[],				/* point to test */
					double	A[],				/* 1st point on line */
					double	B[],				/* 2nd point on line */
					long	Lines,				/* height of the vectors */
					double	Tolerance,			/* admissible relative error */
					int		*Membership,		/* resulting type of intersection */
					int		*Status				/* error management */
				);

/*--------------------------------------------------------------------------*/
/** Does a point belong to a line (3-D)?
    Test the membership of the point P to the line (A,b).
    All points have 3 elements. All vectors have 3 elements.
    
    A line is described by a point A and a vector b: X = A + t * b,
    with t a free scalar. The vector b must have a normalized unit length.
    
    return Membership = FALSE if P doesn't belong to (A,b).
    return Membership = TRUE if P does belong to (A,b).
    
    success: return(!ERROR); failure: return(ERROR)
*/
extern int		PointToLineMembership3D
				(
					double	P[],				/* point to test */
					double	A[],				/* point on line */
					double	b[],				/* direction of line */
					double	Tolerance,			/* admissible relative error */
					int		*Membership			/* resulting type of intersection */
				);

/*--------------------------------------------------------------------------*/
/** Does a point belong to a plane (n-D)?.
    Test the membership of the point P to the plane (A,B,C).
    All vectors have (Lines) elements.
    
    return Membership = FALSE if P doesn't belong to (A,B,C).
    return Membership = TRUE if P does belong to (A,B,C)
    
    success: return(!ERROR); failure: return(ERROR)
*/
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
				);

/*--------------------------------------------------------------------------*/
/** Does a point belong to a plane (3-D)?
    Test the membership of the point P to the plane (A,b).
    All points have 3 elements. All vectors have 3 elements.
    
    A plane is described by a point A and a normal b: <AX, b> = 0.
    The vectors b must have a normalized unit length.
    
    return Membership = FALSE if P doesn't belong to (A,b).
    return Membership = TRUE if P does belong to (A,b)
    
    success: return(!ERROR); failure: return(ERROR)
*/
extern int		PointToPlaneMembership3D
				(
					double	P[],				/* point to test */
					double	A[],				/* point on plane */
					double	b[],				/* normal of plane */
					double	Tolerance,			/* admissible relative error */
					int		*Membership			/* resulting type of intersection */
				);

/*--------------------------------------------------------------------------*/
/** Project a point onto a line (n-D).
    Project the point P onto the line (A,B).
    All vectors have (Lines) elements.
    
    success: return(!ERROR); failure: return(ERROR)
*/
extern int		ProjectPointToLine
				(
					double	P[],				/* point to project */
					double	A[],				/* 1st point on line */
					double	B[],				/* 2nd point on line */
					double	X[],				/* resulting projection */
					long	Lines,				/* height of the vectors */
					double	Tolerance			/* admissible relative error */
				);

/*--------------------------------------------------------------------------*/
/** Project a point onto a line (3-D).
    Project the point P onto the line (A,b).
    All points have 3 elements. All vectors have 3 elements.
    
    A line is described by a point A and a vector b: X = A + t * b,
    with t a free scalar. The vector b must have a normalized unit length.
    
    success: return(!ERROR); failure: return(ERROR)
*/
extern int		ProjectPointToLine3D
				(
					double	P[],				/* point to project */
					double	A[],				/* point on line */
					double	b[],				/* direction of line */
					double	X[],				/* resulting projection */
					double	Tolerance			/* admissible relative error */
				);

/*--------------------------------------------------------------------------*/
/** Project a point onto a plane (n-D).
    Project the point P onto the plane (A,B,C).
    All vectors have (Lines) elements.
    
    success: return(!ERROR); failure: return(ERROR)
*/
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
				);

/*--------------------------------------------------------------------------*/
/** Project a point onto a plane (3-D).
    Project the point P onto the plane (A,b).
    All points have 3 elements. All vectors have 3 elements.
    
    A plane is described by a point A and a normal b: <AX, b> = 0.
    The vectors b must have a normalized unit length.
    
    success: return(!ERROR); failure: return(ERROR)
*/

extern int		ProjectPointToPlane3D
				(
					double	P[],				/* point to project */
					double	A[],				/* point on plane */
					double	b[],				/* normal of plane */
					double	X[],				/* resulting projection */
					double	Tolerance			/* admissible relative error */
				);

/*--------------------------------------------------------------------------*/
/** Are two vectors colinear (n-D)?.
    Determine whether the vectors U and V are colinear.
    The vectors U and V have (Lines) elements.
    
    return Colinear = -1 if at least one of the vectors is degenerate.
    return Colinear = 0 if the vectors are not colinear.
    return Colinear = 1 if the vectors are colinear.
    
    success: return(!ERROR); failure: return(ERROR)
*/
extern int		TestColinearVector
				(
					double	U[],				/* 1st vector */
					double	V[],				/* 2nd vector */
					long	Lines,				/* height of the vector */
					double	Tolerance,			/* admissible relative error */
					int		*Colinear			/* resulting test of colinearity */
				);

/*--------------------------------------------------------------------------*/
/** Are three vectors coplanar (n-D)?
    Determine whether the vectors (U,V,W) are coplanar.
    The vectors (U,V,W) have (Lines) elements.
    
    return Coplanar = -1 in degenerate cases.
    return Coplanar = 0 if the vectors are not coplanar.
    return Coplanar = 1 if the vectors are coplanar.
    
    success: return(!ERROR); failure: return(ERROR)
*/
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
				);
//@}