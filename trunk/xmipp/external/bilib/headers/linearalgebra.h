#ifndef _LINEARALGEBRA
#define _LINEARALGEBRA
/**@defgroup BilibLinearAlgebra Linear algebra
   @ingroup BilibLibrary */
//@{
/*--------------------------------------------------------------------------*/
extern int  AllocateMatrix
(
    double **Matrix,   /* double output pointer */
    long Lines,    /* width of the matrix */
    long Columns,   /* height of the matrix */
    int  *Status    /* error management */
);

/*--------------------------------------------------------------------------*/
extern int  AllocateVector
(
    double *(Vector[]),  /* double output pointer */
    long Lines,    /* height of the vector */
    int  *Status    /* error management */
);

/*--------------------------------------------------------------------------*/
extern int  FreeMatrix
(
    double **Matrix   /* (double)matrix array */
);

/*--------------------------------------------------------------------------*/
extern int  FreeVector
(
    double *(Vector[])   /* (double)vector array */
);

/*--------------------------------------------------------------------------*/
/** Compute the Frobenius norm of the matrix A.
    The size of A is (Lines x Columns).
    success: return(!ERROR); failure: return(ERROR) */
extern int  FrobeniusNorm
(
    double *A,     /* input matrix */
    double *Norm,    /* resulting scalar norm */
    long Lines,    /* height of the matrix A */
    long Columns    /* width of the matrix A */
);

/*--------------------------------------------------------------------------*/
extern int  GetIdentitySquareMatrix
(
    double *A,     /* output identity matrix */
    long Size    /* height or width of the square matrix */
);

/*--------------------------------------------------------------------------*/
extern int  LinearLeastSquares
(
    double *A,     /* matrix */
    long Lines,    /* height the matrix */
    long Columns,   /* width of the matrix */
    double B[],    /* right-hand side vector */
    double X[],    /* resulting vector */
    double Tolerance,   /* admissible relative error */
    int  *Status    /* error management */
);

/*--------------------------------------------------------------------------*/
extern int  MatrixAdd
(
    double *A,     /* left matrix operand */
    double *B,     /* right matrix operand */
    double *X,     /* resulting matrix */
    long Lines,    /* height of the matrices */
    long Columns    /* width of the matrices */
);

/*--------------------------------------------------------------------------*/
/** Compute the condition number of the matrix A.
    The size of A is (Lines x Columns).
    ConditionNumber = -1.0 indicates a singular matrix.

    success: return(!ERROR); failure: return(ERROR) */
extern int  MatrixConditionNumber
(
    double *A,     /* input matrix */
    long Lines,    /* height of the matrix A */
    long Columns,   /* width of the matrix A */
    double *ConditionNumber, /* resulting scalar condition number */
    double Tolerance,   /* admissible relative error */
    long MaxIterations,  /* convergence limit */
    int  *Status    /* error management */
);

/*--------------------------------------------------------------------------*/
/** Orthonormalize a matrix.
    The size of A and B is (Lines x Columns).
    The orthogonalization proceeds row-wise.

    success: return(!ERROR); failure: return(ERROR) */
extern int  MatrixGramSchmidtOrthonormalize
(
    double *A,     /* input matrix */
    double *B,     /* resulting matrix */
    long Lines,    /* height of the matrix A */
    long Columns,   /* width of the matrix A */
    double Tolerance   /* admissible relative error */
);

/*--------------------------------------------------------------------------*/
/** Extract a minor matrix B from the major matrix A.
    The size of A is (Lines x Columns).
    The size of B is ((Lines - 1) x (Columns - 1)).
    The line to delete is i (indexing starts from 0).
    The column to delete is j (indexing starts from 0).

    success: return(!ERROR); failure: return(ERROR) */
extern int  MatrixMinor
(
    double *A,     /* input matrix */
    double *B,     /* smaller output matrix */
    long Lines,    /* height of the matrix A */
    long Columns,   /* width of the matrix A */
    long i,     /* line to delete */
    long j     /* column to delete */
);

/*--------------------------------------------------------------------------*/
extern int  MatrixMultiply
(
    double *A,     /* left matrix operand */
    double *B,     /* right matrix operand */
    double *X,     /* resulting matrix */
    long Lines,    /* height (left and resulting) */
    long CommonSize,   /* left width and right height */
    long Columns    /* width (right and resulting) */
);

/*--------------------------------------------------------------------------*/
/** Compute the rank of the matrix A.
    The size of A is (Lines x Columns).

    success: return(!ERROR); failure: return(ERROR) */
extern int  MatrixRank
(
    double *A,     /* input square matrix */
    long Lines,    /* height of the matrix A */
    long Columns,   /* width of the matrix A */
    long *Rank,    /* resulting scalar condition number */
    double Tolerance,   /* admissible relative error */
    long MaxIterations,  /* convergence limit */
    int  *Status    /* error management */
);

/*--------------------------------------------------------------------------*/
extern int  MatrixSubtract
(
    double *A,     /* left matrix operand */
    double *B,     /* right matrix operand */
    double *X,     /* resulting matrix */
    long Lines,    /* height of the matrices */
    long Columns    /* width of the matrices */
);

/*--------------------------------------------------------------------------*/
extern int  MatrixTimesVector
(
    double *A,     /* left matrix operand */
    double *B,     /* right vector operand */
    double *X,     /* resulting vector */
    long Lines,    /* matrix height */
    long Columns    /* matrix width and vector height */
);

/*--------------------------------------------------------------------------*/
extern int  MatrixTranspose
(
    double *A,     /* input matrix */
    double *At,    /* transposed output matrix */
    long Lines,    /* height input, width output */
    long Columns    /* width input, height output */
);

/*--------------------------------------------------------------------------*/
/** Multiply 3 Matrices.
    X=A*B*C; A is of size Lines x CommonSizeH,
             B is of size CommonSizeH x CommonSizeW,
             C is of size CommonSizeW x Columns.
             Finally, X is of size Lines x Columns.
*/
int  multiply_3Matrices(double *A,
                        double *B,
                        double *C,
                        double *X,
                        long Lines,
                        long CommonSizeH,
                        long CommonSizeW,
                        long Columns
                       );

int multiply_4Matrices( double  *A,
                        double  *B,
                        double  *C,
                        double  *D,
                        double  *X,
                        long    Lines,
                        long    CommonSizeH1,
                        long    CommonSizeW1,
                        long    CommonSizeW2,
                        long    Columns);

int multiply_5Matrices( double  *A,
                        double  *B,
                        double  *C,
                        double  *D,
                        double  *E,
                        double  *X,
                        long    Lines,
                        long    CommonSizeH1,
                        long    CommonSizeW1,
                        long    CommonSizeW2,
                        long    CommonSizeH2,
                        long    Columns );

/*--------------------------------------------------------------------------*/
/** Decompose the (Lines x Columns) input matrix Q into an orthonormal.
    Output matrix Q of same size (Lines x Columns) and an upper-diagonal.
    Square matrix R of size (Columns x Columns), such that the matrix.
    Product (Q * R) gives the input matrix, and such that the matrix.
    Product (Q^T * Q) gives the identity.

    requirement: Columns <= Lines

    success: return(!ERROR); failure: return(ERROR) */
extern int  QRdecomposition
(
    double *Q,     /* in-place matrix */
    double *R,     /* output matrix */
    long Lines,    /* height the matrix */
    long Columns,   /* width of the matrix */
    double Tolerance,   /* admissible relative error */
    int  *Status    /* error management */
);

/*--------------------------------------------------------------------------*/
extern int  SingularValueBackSubstitution
(
    double *U,     /* input matrix (in-place processing) */
    double W[],    /* vector of singular values */
    double *V,     /* untransposed orthogonal matrix */
    long Lines,    /* height the matrix */
    long Columns,   /* width of the matrix */
    double B[],    /* input vector */
    double X[],    /* returned solution */
    int  *Status    /* error management */
);

/*--------------------------------------------------------------------------*/
#define SVDMAXITER 1000
extern int  SingularValueDecomposition
(
    double *U,     /* input matrix (in-place processing) */
    long Lines,    /* height the matrix */
    long Columns,   /* width of the matrix */
    double W[],    /* resulting vector of singular values */
    double *V,     /* untransposed orthogonal matrix */
    long MaxIterations,  /* convergence limit */
    int  *Status    /* error management */
);

/*--------------------------------------------------------------------------*/
extern int  SquareMatrixDeterminant
(
    double *A,     /* input square matrix */
    long Size,    /* height or width of the square matrix */
    double *Determinant,  /* resulting scalar determinant */
    int  *Status    /* error management */
);

/*--------------------------------------------------------------------------*/
extern int  SquareMatrixInvertGauss
(
    double *Direct,   /* input matrix */
    double *Inverse,   /* resulting matrix */
    long Size,    /* height or width of the square matrix */
    double Tolerance,   /* admissible relative error */
    int  *Status    /* error management */
);

/*--------------------------------------------------------------------------*/
/** Compute one out of the 2^Size square roots B of a matrix A such that
    B.B = A.
    The size of the matrices A and B is (Size x Size).
    B must provide an initial solution.

    success: return(!ERROR); failure: return(ERROR) */
extern int  SquareMatrixSquareRoot
(
    double *A,     /* input matrix */
    double *B,     /* output matrix and initial solution */
    long Size,    /* height or width of the square matrix */
    double Tolerance,   /* admissible relative error */
    long MaxIterations,  /* convergence limit */
    int  *Status    /* error management */
);

/*--------------------------------------------------------------------------*/
/** Perform the vector cross product X = A x B.
    The size of A is (Size x Size)

    success: return(!ERROR); failure: return(ERROR) */
extern int  SquareMatrixTrace
(
    double *A,     /* input matrix */
    double *Trace,    /* resulting scalar trace */
    long Size    /* height or width of the square matrix */
);

/*--------------------------------------------------------------------------*/
extern int  Vector3DcrossProduct
(
    double A[],    /* left vector operand */
    double B[],    /* right vector operand */
    double X[]     /* resulting vector */
);

/*--------------------------------------------------------------------------*/
extern int  VectorAdd
(
    double A[],    /* left vector operand */
    double B[],    /* right vector operand */
    double X[],    /* resulting vector */
    long Lines    /* height of the vectors */
);

/*--------------------------------------------------------------------------*/
/** Compute the angle between two vectors *(n-D).
    The size of A and B is (Lines x 1). The angle unit is radian.

    success: return(!ERROR); failure: return(ERROR) */
extern int  VectorAngle
(
    double A[],    /* 1st input vector */
    double B[],    /* 2nd input vector */
    double *Angle,    /* resulting scalar angle in rad */
    long Lines,    /* height of the vectors */
    double Tolerance   /* admissible relative error */
);

/*--------------------------------------------------------------------------*/
/** Convert a vector from non-homogenous to homogenous coordinates.
    Satisfy (Vout = a * Vin) such that (Vout[last] = 1.0).
    The size of A and X is (Lines x 1).

    success: return(!ERROR); failure: return(ERROR) */
extern int  VectorHomogenize
(
    double A[],    /* input vector */
    double X[],    /* output homogenous vector */
    long Lines,    /* height of the vectors */
    double Tolerance   /* admissible relative error */
);

/*--------------------------------------------------------------------------*/
extern int  VectorNorm
(
    double A[],    /* input vector */
    double *Norm,    /* resulting scalar norm */
    long Lines    /* height of the vector */
);

/*--------------------------------------------------------------------------*/
extern int  VectorNormalize
(
    double A[],    /* input vector */
    double X[],    /* resulting normalized vector */
    long Lines    /* height of the vector */
);

/*--------------------------------------------------------------------------*/
extern int  VectorScalarProduct
(
    double A[],    /* left vector operand */
    double B[],    /* right vector operand */
    double *X,     /* resulting scalar product */
    long Lines    /* height of the vector */
);

/*--------------------------------------------------------------------------*/
extern int  VectorScale
(
    double A[],    /* input vector */
    double X[],    /* output vector */
    double Scale,    /* scale factor */
    long Lines    /* height of the vector */
);

/*--------------------------------------------------------------------------*/
extern int  VectorSubtract
(
    double A[],    /* left vector operand */
    double B[],    /* right vector operand */
    double X[],    /* resulting vector */
    long Lines    /* height of the vectors */
);
//@}
#endif
