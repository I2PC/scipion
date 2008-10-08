/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/

#ifndef MATRIX2D_H
#define MATRIX2D_H

#include <iostream>
#include <string>
#include <complex>

#include <external/bilib/headers/pyramidtools.h>
#include "matrix1d.h"
#include "multidimensional_array.h"
#include "error.h"

#ifndef SWIG
template<typename T>
void multiplyElements(const Matrix2D<T>& op1, const Matrix2D<T>& op2,
    Matrix2D<T>& result);

template<typename T>
void multiplyMatrix(const Matrix2D<T>& op1, const Matrix2D<T>& op2,
	Matrix2D<T>& result);

template<typename T>
void solve(const Matrix2D<T>& A, const Matrix1D<T>& b, Matrix1D<T>& result);

template<typename T>
void solve(const Matrix2D<T>& A, const Matrix2D<T>& b, Matrix2D<T>& result);

template<typename T>
void solveBySVD(const Matrix2D<T>& A,
                const Matrix1D<T>& b,
                Matrix1D< double >& result,
                double tolerance);

template<typename T>
void ludcmp(const Matrix2D<T>& A, Matrix2D<T>& LU, Matrix1D< int >& indx, T& d);

template<typename T>
void lubksb(const Matrix2D<T>& LU, Matrix1D< int >& indx, Matrix1D<T>& b);

template<typename T>
void svdcmp(const Matrix2D< T >& a,
            Matrix2D< double >& u,
            Matrix1D< double >& w,
            Matrix2D< double >& v);

void svbksb(Matrix2D< double >& u,
            Matrix1D< double >& w,
            Matrix2D< double >& v,
            Matrix1D< double >& b,
            Matrix1D< double >& x);

template<typename T>
void applyGeometry(Matrix2D<T>& m2,
                Matrix2D< double > A,
                const Matrix2D<T>& m1,
                bool inv,
                bool wrap,
                T outside = (T) 0);

template<typename T>
void applyGeometryBSpline(Matrix2D<T>& m2,
                        Matrix2D< double > A,
                        const Matrix2D<T>& m1,
                        int Splinedegree,
                        bool inv,
                        bool wrap,
                        T outside = (T) 0);

Matrix2D< double > rotation2DMatrix(double ang);

Matrix2D< double > translation2DMatrix(const Matrix1D< double > v);

int bestPrecision(float F, int _width);
#endif

/// @defgroup Matrices Matrices
/// @ingroup MultidimensionalArrays

/** @defgroup MatricesSpeedUp Speed up macros
 * @ingroup Matrices
 *
 * This macros are defined to allow high speed in critical parts of your
 * program. They shouldn't be used systematically as usually there is no
 * checking on the correctness of the operation you are performing. Speed comes
 * from three facts: first, they are macros and no function call is performed
 * (although most of the critical functions are inline functions), there is no
 * checking on the correctness of the operation (it could be wrong and you are
 * not warned of it), and destination vectors are not returned saving time in
 * the copy constructor and in the creation/destruction of temporary vectors.
 */

/** @defgroup MatricesSizeShape Size and shape
 * @ingroup MatricesSpeedUp
 *
 * Although they are not defined here you can also use STARTINGX and FINISHINGX
 * (defined for Matrix1D)
 */

/** TRUE if both arrays have the same shape
 * @ingroup MatricesSizeShape
 *
 * Two arrays have the same shape if they have the same size and the same
 * starting point. Be aware that this is a macro which simplifies to a boolean.
 */
#define SAME_SHAPE2D(v1, v2) \
    (XSIZE(v1) == XSIZE(v2) && \
     YSIZE(v1) == YSIZE(v2) && \
     STARTINGX(v1) == STARTINGX(v2) && \
     STARTINGY(v1) == STARTINGY(v2))

/** For all elements in the array
 * @ingroup MatricesSizeShape
 *
 * This macro is used to generate loops for the matrix in an easy way. It
 * defines internal indexes 'i' and 'j' which ranges the matrix using its
 * mathematical definition (ie, logical access).
 *
 * @code
 * FOR_ALL_ELEMENTS_IN_MATRIX2D(m)
 * {
 *     std::cout << m(i, j) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_ELEMENTS_IN_MATRIX2D(m) \
    for (int i=STARTINGY(m); i<=FINISHINGY(m); i++) \
        for (int j=STARTINGX(m); j<=FINISHINGX(m); j++)

/** For all elements in the array between corners
 * @ingroup MatricesSizeShape
 *
 * This macro is used to generate loops for a matrix in an easy manner. It needs
 * an externally defined Matrix1D< double > r(2). Then YY(r) and XX(r) range
 * from (int) YY(corner1) to (int)YY(corner2), (int) XX(corner1) to (int)
 * XX(corner2) (included limits) respectively. Notice that corner1 and corner2
 * need only be Matrix1D.
 *
 * @code
 * Matrix1D< double > corner1(2), corner2(2);
 * Matrix1D< int > r(2);
 * XX(corner1) = -1;
 * XX(corner2) = 1;
 * YY(corner1) = -2;
 * YY(corner2) = 2;
 *
 * FOR_ALL_ELEMENTS_IN_MATRIX2D_BETWEEN(corner1, corner2)
 * {
 *     std::cout << v(r) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_ELEMENTS_IN_MATRIX2D_BETWEEN(corner1, corner2) \
    for (YY(r)=YY((corner1)); YY(r)<=YY((corner2)); YY(r)++) \
        for (XX(r)=XX((corner1)); XX(r)<=XX((corner2)); XX(r)++)

/** For all elements in common
 * @ingroup MatricesSizeShape
 *
 * This macro is used to generate loops for all the elements logically in common
 * between two images in an easy manner. Then i and j (locally defined) range
 * from MAX(STARTINGY(V1), STARTINGY(V2)) to MIN(FINISHINGY(V1),
 * FINISHINGY(V2)), MAX(STARTINGX(V1), STARTINGX(V2)) to MIN(FINISHINGX(V1),
 * FINISHINGX(V2)) (included limits) respectively. You need to define
 * SPEED_UP_temps.
 *
 * @code
 * Matrix2D< double > m1(10, 10), m2(20, 20);
 * m1.setXmippOrigin();
 * m2.setXmippOrigin();
 *
 * FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX2D(m1, m2)
 * {
 *     ...
 * }
 * @endcode
 */
#define FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX2D(m1, m2) \
    ispduptmp2 = XMIPP_MAX(STARTINGY(m1), STARTINGY(m2)); \
    ispduptmp3 = XMIPP_MIN(FINISHINGY(m1), FINISHINGY(m2)); \
    ispduptmp4 = XMIPP_MAX(STARTINGX(m1), STARTINGX(m2)); \
    ispduptmp5 = XMIPP_MIN(FINISHINGX(m1), FINISHINGX(m2)); \
    for (int i=ispduptmp2; i<=ispduptmp3; i++) \
        for (int j=ispduptmp4; j<=ispduptmp5; j++)

/** For all elements in the array, accessed physically
 * @ingroup MatricesSizeShape
 *
 * This macro is used to generate loops for the matrix in an easy way using
 * physical indexes. It defines internal indexes 'i' and 'j' which ranges the
 * matrix using its physical definition.
 *
 * @code
 * FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(m)
 * {
 *     std::cout << DIRECT_MAT_ELEM(m, i, j) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(m) \
    for (int i=0; i<YSIZE(m); i++) \
        for (int j=0; j<XSIZE(m); j++)

/// @defgroup MatricesMemory Memory access
/// @ingroup MatricesSpeedUp

/** @defgroup MatricesArithmetic Arithmethic operations
 * @ingroup MatricesSpeedUp
 *
 * The vectors and matrices involved in these macros should be created with the
 * correct size before entering in them. These macros allow a fast operation on
 * R2 and R3 vectors, and small size matrices. These macros need some temporary
 * variables. You must "load" them by calling the macro SPEED_UP_temps at the
 * beginning of your function
 */

/** Matrix (3x3) by vector (3x1) (a=M*b)
 * @ingroup MatricesArithmetic
 *
 * You must "load" the temporary variables, and create the result vector with
 * the appropiate size. You can reuse the vector b to store the results (that
 * is, M3x3_BY_V3x1(b, M, b);, is allowed).
 *
 * @code
 * double example
 * {
 *     SPEED_UP_temps;
 *
 *     Matrix1D< double > a(3), b(3);
 *     Matrix2D< double > M(3, 3);
 *
 *     M.init_random(0, 1);
 *     b.init_random(0, 1);
 *     M3x3_BY_V3x1(a, M, b);
 *
 *     return a.sum();
 * }
 * @endcode
 */
#define M3x3_BY_V3x1(a, M, b) { \
        spduptmp0 = dMij(M, 0, 0) * XX(b) + dMij(M, 0, 1) * YY(b) + dMij(M, 0, 2) \
                    * ZZ(b); \
        spduptmp1 = dMij(M, 1, 0) * XX(b) + dMij(M, 1, 1) * YY(b) + dMij(M, 1, 2) \
                    * ZZ(b); \
        spduptmp2 = dMij(M, 2, 0) * XX(b) + dMij(M, 2, 1) * YY(b) + dMij(M, 2, 2) \
                    * ZZ(b); \
        XX(a) = spduptmp0; YY(a) = spduptmp1; ZZ(a) = spduptmp2; }


/** Matrix (3x3) by Matrix (3x3) (A=B*C)
 * @ingroup MatricesArithmetic
 *
 * You must "load" the temporary variables, and create the result vector with
 * the appropiate size. You can reuse any of the multiplicands to store the
 * results (that is, M3x3_BY_M3x3(A, A, B);, is allowed).
 */
#define M3x3_BY_M3x3(A, B, C) { \
        spduptmp0 = dMij(B, 0, 0) * dMij(C, 0, 0) + dMij(B, 0, 1) * dMij(C, 1, 0) \
                    + dMij(B, 0, 2) * dMij(C, 2, 0); \
        spduptmp1 = dMij(B, 0, 0) * dMij(C, 0, 1) + dMij(B, 0, 1) * dMij(C, 1, 1) \
                    + dMij(B, 0, 2) * dMij(C, 2, 1); \
        spduptmp2 = dMij(B, 0, 0) * dMij(C, 0, 2) + dMij(B, 0, 1) * dMij(C, 1, 2) \
                    + dMij(B, 0, 2) * dMij(C, 2, 2); \
        spduptmp3 = dMij(B, 1, 0) * dMij(C, 0, 0) + dMij(B, 1, 1) * dMij(C, 1, 0) \
                    + dMij(B, 1, 2) * dMij(C, 2, 0); \
        spduptmp4 = dMij(B, 1, 0) * dMij(C, 0, 1) + dMij(B, 1, 1) * dMij(C, 1, 1) \
                    + dMij(B, 1, 2) * dMij(C, 2, 1); \
        spduptmp5 = dMij(B, 1, 0) * dMij(C, 0, 2) + dMij(B, 1, 1) * dMij(C, 1, 2) \
                    + dMij(B, 1, 2) * dMij(C, 2, 2); \
        spduptmp6 = dMij(B, 2, 0) * dMij(C, 0, 0) + dMij(B, 2, 1) * dMij(C, 1, 0) \
                    + dMij(B, 2, 2) * dMij(C, 2, 0); \
        spduptmp7 = dMij(B, 2, 0) * dMij(C, 0, 1) + dMij(B, 2, 1) * dMij(C, 1, 1) \
                    + dMij(B, 2, 2) * dMij(C, 2, 1); \
        spduptmp8 = dMij(B, 2, 0) * dMij(C, 0, 2) + dMij(B, 2, 1) * dMij(C, 1, 2) \
                    + dMij(B, 2, 2) * dMij(C, 2, 2); \
        dMij(A, 0, 0) = spduptmp0; \
        dMij(A, 0, 1) = spduptmp1; \
        dMij(A, 0, 2) = spduptmp2; \
        dMij(A, 1, 0) = spduptmp3; \
        dMij(A, 1, 1) = spduptmp4; \
        dMij(A, 1, 2) = spduptmp5; \
        dMij(A, 2, 0) = spduptmp6; \
        dMij(A, 2, 1) = spduptmp7; \
        dMij(A, 2, 2) = spduptmp8; }

/** Matrix (2x2) by vector (2x1) (a=M*b)
 * @ingroup MatricesArithmetic
 *
 * You must "load" the temporary variables, and create the result vector with
 * the appropiate size. You can reuse the vector b to store the results (that
 * is, M2x2_BY_V2x1(b, M, b);, is allowed).
 *
 * @code
 * double example
 * {
 *     SPEED_UP_temps;
 *
 *     Matrix1D< double > a(2), b(2);
 *     Matrix2D< double > M(2, 2);
 *
 *     M.init_random(0, 1);
 *     b.init_random(0, 1);
 *
 *     M2x2_BY_V2x1(a, M, b);
 *
 *     return a.sum();
 * }
 * @endcode
 */
#define M2x2_BY_V2x1(a, M, b) { \
        spduptmp0 = dMij(M, 0, 0) * XX(b) + dMij(M, 0, 1) * YY(b); \
        spduptmp1 = dMij(M, 1, 0) * XX(b) + dMij(M, 1, 1) * YY(b); \
        XX(a) = spduptmp0; \
        YY(a) = spduptmp1; }

/** Matrix (2x2) by constant (M2=M1*k)
 * @ingroup MatricesArithmetic
 *
 * You must create the result matrix with the appropiate size. You can reuse
 * the matrix M1 to store the results (that is, M2x2_BY_CT(M, M, k);, is
 * allowed).
 */
#define M2x2_BY_CT(M2, M1, k) { \
        dMij(M2, 0, 0) = dMij(M1, 0, 0) * k; \
        dMij(M2, 0, 1) = dMij(M1, 0, 1) * k; \
        dMij(M2, 1, 0) = dMij(M1, 1, 0) * k; \
        dMij(M2, 1, 1) = dMij(M1, 1, 1) * k; }

/** Matrix (3x3) by constant (M2=M1*k)
 * @ingroup MatricesArithmetic
 *
 * You must create the result matrix with the appropiate size. You can reuse the
 * matrix M1 to store the results (that is, M2x2_BY_CT(M, M, k);, is allowed).
 */
#define M3x3_BY_CT(M2, M1, k) { \
        dMij(M2, 0, 0) = dMij(M1, 0, 0) * k; \
        dMij(M2, 0, 1) = dMij(M1, 0, 1) * k; \
        dMij(M2, 0, 2) = dMij(M1, 0, 2) * k; \
        dMij(M2, 1, 0) = dMij(M1, 1, 0) * k; \
        dMij(M2, 1, 1) = dMij(M1, 1, 1) * k; \
        dMij(M2, 1, 2) = dMij(M1, 1, 2) * k; \
        dMij(M2, 2, 0) = dMij(M1, 2, 0) * k; \
        dMij(M2, 2, 1) = dMij(M1, 2, 1) * k; \
        dMij(M2, 2, 2) = dMij(M1, 2, 2) * k; }

/** Inverse of a matrix (2x2)
 * @ingroup MatricesArithmetic
 *
 * Input and output matrix cannot be the same one. The output is supposed to be
 * already resized.
 */
#define M2x2_INV(Ainv, A) { \
        spduptmp0 = 1.0 / (dMij(A, 0, 0) * dMij(A, 1, 1) - dMij(A, 0, 1) \
                           * dMij(A, 1, 0)); \
        dMij(Ainv, 0, 0) = dMij(A, 1, 1); \
        dMij(Ainv, 0, 1) = -dMij(A, 0, 1); \
        dMij(Ainv, 1, 0) = -dMij(A, 1, 0); \
        dMij(Ainv, 1, 1) =  dMij(A, 0, 0); \
        M2x2_BY_CT(Ainv, Ainv, spduptmp0); }
//@}

/// Template class for Xmipp matrices
/// @ingroup Matrices
template<typename T>
class Matrix2D: public MultidimArray<T>
{
public:
    /// Matrix multiplication in a algebraic way
    friend void multiplyMatrix<>(const Matrix2D<T>& op1, const Matrix2D<T>& op2, Matrix2D<T>& result);

    /** The following functions are directly an interface to the Numerical
     * Recipes ones with the same name, but with a better interface LU
     * decomposition
     */
    friend void ludcmp<>(const Matrix2D<T>& A, Matrix2D<T>& LU, Matrix1D< int >& indx, T& d);

    /// Solving an equation system based on LU. Remember to free LU outside
    friend void lubksb<>(const Matrix2D<T>& LU, Matrix1D< int >& indx, Matrix1D<T>& b);

    /// @defgroup MatricesConstructors Constructors
    /// @ingroup Matrices

    /** Empty constructor
     * @ingroup MatricesConstructors
     */
    Matrix2D(): MultidimArray<T>()
    {
    }

    /** Dimension constructor
     * @ingroup MatricesConstructors
     *
     * The dimension constructor creates a matrix with memory associated (but
     * not assigned to anything, could be full of garbage) origin=0, size=the
     * given one. e careful that first number is the Y dimension (number of
     * rows), and the second the X dimension (number of columns).
     *
     * @code
     * Matrix2D< double > v1(6, 3);
     * @endcode
     */
    Matrix2D(int Ydim, int Xdim): MultidimArray<T>()
    {
    	resize(Ydim,Xdim);
    }

    /** Copy constructor
     * @ingroup MatricesConstructors
     *
     * The created matrix is a perfect copy of the input matrix but with a
     * different memory assignment.
     *
     * @code
     * Matrix2D< double > m2(m1);
     * @endcode
     */
    Matrix2D(const Matrix2D<T>& v)
    {
        *this = v;
    }

    /** Clear.
     * @ingroup MatricesConstructors
     */
     void clear()
     {
        MultidimArray<T>::clear();
     }

    /// @defgroup MatricesInitialization Initialization
    /// @ingroup Matrices

    /** Alias a multidimarray.
     * @ingroup MatricesInitialization
     *
     * Treat the multidimarray as if it were a matrix. The data is not copied
     * into new memory, but a pointer to the multidimarray is copied.
     * You should not make any operation on this matrix such that the
     * memory locations are changed
     */
     void alias(const MultidimArray<T> &m)
     {
         copyShape(m);
         this->data=m.data;
         this->destroyData=false;
     }

    /** Identity matrix of current size
     * @ingroup MatricesInitialization
     *
     * If actually the matrix is not squared then an identity matrix is
     * generated of size (Xdim x Xdim).
     *
     * @code
     * m.initIdentity();
     * @endcode
     */
    void initIdentity()
    {
        initIdentity(XSIZE(*this));
    }

    /** Identity matrix of a given size
     * @ingroup MatricesInitialization
     *
     * A (dim x dim) identity matrix is generated.
     *
     * @code
     * m.initIdentity(3);
     * @endcode
     */
    void initIdentity(int dim)
    {
        initIdentity(dim, dim);
    }

    /** Identity matrix of a given size
     * @ingroup MatricesInitialization
     *
     * A (dimX x dimY) identity matrix is generated. That is, any element i, j
     * of the matrix such that i = j is equal to 1.
     *
     * @code
     * m.initIdentity(2, 3);
     * @endcode
     */
    void initIdentity(int Ydim, int Xdim)
    {
        if (Xdim == 0 || Ydim == 0)
        {
            clear();
            return;
        }

        resize(Ydim, Xdim);
        FOR_ALL_ELEMENTS_IN_MATRIX2D(*this)
	        DIRECT_MAT_ELEM(*this, i, j) = (T)(i == j);
    }

    /** Zero initialisation with a new dimension
     * @ingroup MatricesInitialization
     *
     * Be careful to the size order (Ydim, Xdim).
     *
     * @code
     * v1.initZeros(6, 3);
     * @endcode
     */
    void initZeros(int Ydim, int Xdim)
    {
        resize(Ydim, Xdim);
        initConstant((T) 0);
    }

    /** Zero initialisation with current dimension
     * @ingroup MatricesInitialization
     */
    void initZeros()
    {
        MultidimArray<T>::initZeros();
    }

    /** Zero initialisation with current dimension
     * @ingroup MatricesInitialization
     */
    template <typename T1>
    void initZeros(const Matrix2D<T1> &m)
    {
        MultidimArray<T>::initZeros(m);
    }

    /** Makes a matrix from a vector
     * @ingroup MatricesInitialization
     *
     * The origin of the matrix is set such that it has one of the index origins
     * (X or Y) to the same value as the vector, and the other set to 0
     * according to the shape.
     *
     * @code
     * Matrix2D< double > m = fromVector(v);
     * @endcode
     */
    void fromVector(const Matrix1D<T>& op1)
    {
        // Null vector => Null matrix
        if (XSIZE(op1) == 0)
        {
            clear();
            return;
        }

        // Look at shape and copy values
        if (op1.isRow())
        {
            resize(1, XSIZE(op1));

            for (int j = 0; j < XSIZE(op1); j++)
                DIRECT_MAT_ELEM(*this, 0, j) = DIRECT_VEC_ELEM(op1, j);

            STARTINGX(*this) = STARTINGX(op1);
            STARTINGY(*this) = 0;
        }
        else
        {
            resize(XSIZE(op1), 1);

            for (int i = 0; i < XSIZE(op1); i++)
                DIRECT_MAT_ELEM(*this, i, 0) = DIRECT_VEC_ELEM(op1, i);

            STARTINGX(*this) = 0;
            STARTINGY(*this) = STARTINGX(op1);
        }
    }

    /** Makes a vector from a matrix
     * @ingroup MatricesInitialization
     *
     * An exception is thrown if the matrix is not a single row or a single
     * column. The origin of the vector is set according to the one of the
     * matrix.
     *
     * @code
     * Matrix1D< double > v;
     * m.toVector(v);
     * @endcode
     */
    void toVector(Matrix1D<T>& op1) const
    {
        // Null matrix => Null vector
        if (XSIZE(*this) == 0 || YSIZE(*this) == 0)
        {
            op1.clear();
            return;
        }

        // If matrix is not a vector, produce an error
        if (XSIZE(*this) != 1 && (YSIZE(*this) != 1))
            REPORT_ERROR(1102,
                         "To_vector: Matrix cannot be converted to vector");

        // Look at shape and copy values
        if (YSIZE(*this) == 1)
        {
            // Row vector
            op1.resize(XSIZE(*this));

            for (int j = 0; j < XSIZE(*this); j++)
                DIRECT_VEC_ELEM(op1, j) = DIRECT_MAT_ELEM(*this, 0, j);

            op1.setRow();
            STARTINGX(op1) = STARTINGX(*this);
        }
        else
        {
            // Column vector
            op1.resize(YSIZE(*this));

            for (int i = 0; i < YSIZE(*this); i++)
                DIRECT_VEC_ELEM(op1, i) = DIRECT_MAT_ELEM(*this, i, 0);

            op1.setCol();
            STARTINGX(op1) = STARTINGY(*this);
        }
    }

    /** @defgroup MatricesSize Size and shape
     * @ingroup Matrices
     *
     * The shape of a matrix is defined by its origin and its size. The size is
     * clear, and the origin is the logical position of the first real position
     * of the array. For instance, if we have a matrix of dimension (5,3) =
     * (Ydim, Xdim) and origin (-2,-1), this means that the array is
     * representing the logical positions
     *
     * @code
     * [(-2,-1) (-2,0) (-2,1)
     *  (-1,-1) (-1,0) (-1,1)
     *  ( 0,-1) ( 0,0) ( 0,1)
     *  ( 1,-1) ( 1,0) ( 1,1)
     *  ( 2,-1) ( 2,0) ( 2,1)]
     * @endcode
     *
     * we could access to any of these positions (Ex: v(-2,1)=3;) and actually
     * any try to access to a position related to 5 (Ex: v(4,1)=3;), although
     * it physically exists, is not logically correct and hence it will throw an
     * exception. The startingX and finishingX positions for this sample vector
     * are -1 and 1 respectively, while for Y are -2 and 2. The "for" iterations
     * through the matrix should include these two values if you want to cover
     * the whole matrix.
     *
     * @code
     * for (int i=STARTINGY(m); i<=FINISHINGY(m); i++)
     *     for (int j=STARTINGX(m); j<=FINISHINGX(m); j++)
     *         MAT_ELEM(m, i, j) += 1;
     * @endcode
     */

    /** Resize to a given size
     * @ingroup MatricesSize
     *
     * This function resize the actual array to the given size. The origin is
     * not modified. If the actual array is larger than the pattern then the
     * values outside the new size are lost, if it is smaller then 0's are
     * added. An exception is thrown if there is no memory.
     *
     * @code
     * m1.resize(3, 2);
     * @endcode
     */
    void resize(int Ydim, int Xdim)
    {
    	MultidimArray<T>::resize(1, Ydim, Xdim);
    }

    /** Resize taking the shape from another matrix
     * @ingroup MatricesSize
     */
    template <typename T1>
    void resize(const Matrix2D<T1> &M)
    {
    	MultidimArray<T>::resize(M);
    }

    /** Produce an array suitable for working with Numerical Recipes
     * @ingroup MatricesSize
     *
     * This function must be used only as a preparation for routines which need
     * that the first physical index is 1 and not 0 as it usually is in C. New
     * memory is needed to hold the new double pointer array.
     */
    T** adaptForNumericalRecipes() const
    {
        T** m = NULL;
        ask_Tmatrix(m, 1, YSIZE(*this), 1, XSIZE(*this));

        FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(*this)
            m[i+1][j+1] = DIRECT_MAT_ELEM(*this, i, j);

        return m;
    }

    /** Produce a pointer suitable for working with Numerical Recipes (2)
     * @ingroup MatricesSize
     *
     * This function meets the same goal as the one before, however this one
     * work with 2D arrays as a single pointer. The first element of the array
     * is pointed by result[1*Xdim+1], and in general result[i*Xdim+j]
     */
    T* adaptForNumericalRecipes2() const
    {
        return MULTIDIM_ARRAY(*this) - 1 - XSIZE(*this);
    }

    /** Load from numerical recipes result.
     * @ingroup MatricesSize
     */
    void loadFromNumericalRecipes(T** m, int Ydim, int Xdim)
    {
        resize(Ydim, Xdim);

        for (int i = 1; i <= Ydim; i++)
            for (int j = 1; j <= Xdim; j++)
                (*this)(i - 1, j - 1) = m[i][j];
    }

    /** Kill an array produced for numerical recipes
     * @ingroup MatricesSize
     *
     * The allocated memory is freed.
     */
    void killAdaptationForNumericalRecipes(T** m) const
    {
        free_Tmatrix(m, 1, YSIZE(*this), 1, XSIZE(*this));
    }

    /** Kill an array produced for numerical recipes, 2.
     * @ingroup MatricesSize
     *
     * Nothing needs to be done in fact.
     */
    void killAdaptationForNumericalRecipes2(T** m) const
        {}

    /** Outside
     * @ingroup MatricesSize
     *
     * TRUE if the logical index given is outside the definition region of this
     * array.
     */
    bool outside(int i, int j) const
    {
        return (j < STARTINGX(*this) || j > FINISHINGX(*this) ||
                i < STARTINGY(*this) || i > FINISHINGY(*this));
    }

    /** Outside
     * @ingroup MatricesSize
     *
     * TRUE if the logical index given is outside the definition region of this
     * array.
     */
    bool outside(const Matrix1D<double> &r) const
    {
        if (XSIZE(r) < 2)
            REPORT_ERROR(1, "Outside: index vector has not got enough components");
    
        return (XX(r) < STARTINGX(*this) || XX(r) > FINISHINGX(*this) ||
                YY(r) < STARTINGY(*this) || YY(r) > FINISHINGY(*this));
    }

    /** IsCorner
     * @ingroup MatricesSize
     *
     * TRUE if the logical index given is a corner of the definition region of this
     * array.
     */
    bool isCorner(const Matrix1D< double >& v) const
    {
        if (XSIZE(v) < 2)
            REPORT_ERROR(1, "isCorner: index vector has got not enough components");

        return ((XX(v) == STARTINGX(*this)  && YY(v) == STARTINGY(*this))  ||
                (XX(v) == STARTINGX(*this)  && YY(v) == FINISHINGY(*this)) ||
                (XX(v) == FINISHINGX(*this) && YY(v) == STARTINGY(*this))  ||
                (XX(v) == FINISHINGX(*this) && YY(v) == FINISHINGY(*this)));
    }

    /** Move origin to
     * @ingroup MatricesSize
     *
     * This function adjust logical indexes such that the Xmipp origin of the
     * array moves to the specified position. For instance, an array whose x
     * indexes go from -1 to 1, if we move the origin to 4, then the x indexes
     * go from 3 to 5. This is very useful for convolution operations where you
     * only need to move the logical starting of the array.
     */
    void moveOriginTo(int i, int j)
    {
        STARTINGY(*this) = i + FIRST_XMIPP_INDEX(YSIZE(*this));
        STARTINGX(*this) = j + FIRST_XMIPP_INDEX(XSIZE(*this));
    }

    /** Same shape
     * @ingroup MatricesSize
     *
     * Returns true if this object has got the same shape (origin and size) than
     * the argument
     */
    bool sameShape(const Matrix2D<T>& op) const
    {
        return SAME_SHAPE2D(*this, op);
    }

    /** Print shape shape
     * @ingroup MatricesSize
    */
    void printShape(std::ostream& out=std::cout) const
    {
	out << "Size(Y,X): " << YSIZE(*this) << "x" << XSIZE(*this)
            << " i=[" << STARTINGY(*this) << ".." << FINISHINGY(*this) << "]"
            << " j=[" << STARTINGX(*this) << ".." << FINISHINGX(*this) << "]";
    }

    /** @defgroup MatricesMemory Memory access
     * @ingroup Matrices
     *
     * This functions allows you to access to the matrix elements.
     */

    /** Matrix element access via index
     * @ingroup MatricesMemory
     *
     * Returns the value of a matrix logical position. In our example we could
     * access from v(-2,-1) to v(2,1). The elements can be used either by value
     * or by reference. An exception is thrown if the index is outside the
     * logical range. The first argument is the Y position and the second the X
     * position.
     *
     * @code
     * m(-2, 1) = 1;
     * val = m(-2, 1);
     * @endcode
     */
    T& operator()(int i, int j) const
    {
        if (i < STARTINGY(*this) || i > FINISHINGY(*this))
            REPORT_ERROR(1103, static_cast< std::string >
	       ("Matrix subscript (i) out of range i=")+integerToString(i));

        if (j < STARTINGX(*this) || j > FINISHINGX(*this))
            REPORT_ERROR(1103, static_cast< std::string >
	    	("Matrix subscript (j) out of range j=")+integerToString(j));

        return MAT_ELEM(*this, i, j);
    }

    /** Get the pixel at (i, j)
     * @ingroup MatricesMemory
     *
     * Logical access
     */
    T getPixel(int i, int j) const
    {
        return (*this)(i, j);
    }

    /** Set the pixel at (i,j)
     * @ingroup MatricesMemory
     *
     * Logical access
     */
    void setPixel(int i, int j, T val)
    {
        (*this)(i, j) = val;
    }

    /** Matrix element access via double vector
     * @ingroup MatricesMemory
     *
     * Returns the value of a matrix logical position, but this time the element
     * position is determined by a R2 vector. The elements can be used either by
     * value or by reference. An exception is thrown if the index is outside the
     * logical range. Pay attention in the following example that we are
     * accessing the same element as in the previous function but, now we have
     * to give first the X position instead of the Y one because we are building
     * first a vector of the form (x,y).
     *
     * @code
     * m(vectorR2(1, -2)) = 1;
     * val = m(vectorR2(1, -2));
     * @endcode
     */
    T& operator()(const Matrix1D< double >& v) const
    {
        return MAT_ELEM(*this, ROUND(YY(v)), ROUND(XX(v)));
    }

    /** Matrix element access via intger vector
     * @ingroup MatricesMemory
     */
    T& operator()(const Matrix1D< int >& v) const
    {
        return MAT_ELEM(*this, YY(v), XX(v));
    }

    /** Interpolates the value of the 2D matrix M at the point (x,y)
     * @ingroup MatricesMemory
     *
     * Bilinear interpolation. (x,y) are in logical coordinates.
     */
    T interpolatedElement(double x, double y, T outside_value = (T) 0) const
    {
        int x0 = FLOOR(x);
        double fx = x - x0;
        int x1 = x0 + 1;
        int y0 = FLOOR(y);
        double fy = y - y0;
        int y1 = y0 + 1;

        T d00 = outside(y0, x0) ? outside_value : MAT_ELEM(*this, y0, x0);
        T d10 = outside(y1, x0) ? outside_value : MAT_ELEM(*this, y1, x0);
        T d11 = outside(y1, x1) ? outside_value : MAT_ELEM(*this, y1, x1);
        T d01 = outside(y0, x1) ? outside_value : MAT_ELEM(*this, y0, x1);

        double d0 = (T) LIN_INTERP(fx, (double) d00, (double) d01);
        double d1 = (T) LIN_INTERP(fx, (double) d10, (double) d11);
        return (T) LIN_INTERP(fy, d0, d1);
    }

    /** Interpolates the value of the 2D matrix M at the point (x,y) knowing
     * that this image is a set of B-spline coefficients
     * @ingroup MatricesMemory
     *
     * (x,y) are in logical coordinates
     *
     * To interpolate using splines you must first produce the Bspline
     * coefficients. An example to interpolate an image at (0.5,0.5) using
     * splines would be:
     *
     * @code
     * Matrix2D< double > Bspline_coeffs;
     * myImage.produceSplineCoefficients(Bspline_coeffs, 3);
     * interpolated_value = Bspline_coeffs.interpolatedElementBSpline(0.5,
     * 0.5,3);
     * @endcode
     */
    T interpolatedElementBSpline(double x, double y, int SplineDegree = 3) const
    {
        int SplineDegree_1 = SplineDegree - 1;

        // Logical to physical
        y -= STARTINGY(*this);
        x -= STARTINGX(*this);

        int lmax = XSIZE(*this);
        int mmax = YSIZE(*this);
        int l1 = CEIL(x - SplineDegree_1);
        int l2 = l1 + SplineDegree;
        int m1 = CEIL(y - SplineDegree_1);
        int m2 = m1 + SplineDegree;

        double columns = 0.0;
        for (int m = m1; m <= m2; m++)
        {
	    int equivalent_m=m;
	    if      (m<0)             equivalent_m=-m-1;
	    else if (m>=YSIZE(*this)) equivalent_m=2*YSIZE(*this)-m-1;
            int row_m = XSIZE(*this) * equivalent_m;
            double rows = 0.0;
            for (int l = l1; l <= l2; l++)
            {
                double xminusl = x - (double) l;
		int equivalent_l=l;
		if      (l<0)             equivalent_l=-l-1;
		else if (l>=XSIZE(*this)) equivalent_l=2*XSIZE(*this)-l-1;
                double Coeff = (double) DIRECT_MAT_ELEM(*this,equivalent_m,equivalent_l);
                switch (SplineDegree)
                {
                case 2:
                    rows += Coeff * Bspline02(xminusl);
                    break;

                case 3:
                    rows += Coeff * Bspline03(xminusl);
                    break;

                case 4:
                    rows += Coeff * Bspline04(xminusl);
                    break;

                case 5:
                    rows += Coeff * Bspline05(xminusl);
                    break;

                case 6:
                    rows += Coeff * Bspline06(xminusl);
                    break;

                case 7:
                    rows += Coeff * Bspline07(xminusl);
                    break;

                case 8:
                    rows += Coeff * Bspline08(xminusl);
                    break;

                case 9:
                    rows += Coeff * Bspline09(xminusl);
                    break;
                }
            }

            double yminusm = y - (double) m;
            switch (SplineDegree)
            {
            case 2:
                columns += rows * Bspline02(yminusm);
                break;

            case 3:
                columns += rows * Bspline03(yminusm);
                break;

            case 4:
                columns += rows * Bspline04(yminusm);
                break;

            case 5:
                columns += rows * Bspline05(yminusm);
                break;

            case 6:
                columns += rows * Bspline06(yminusm);
                break;

            case 7:
                columns += rows * Bspline07(yminusm);
                break;

            case 8:
                columns += rows * Bspline08(yminusm);
                break;

            case 9:
                columns += rows * Bspline09(yminusm);
                break;
            }
        }
        return (T) columns;
    }

    /** Extracts the profile between two points
     * @ingroup MatricesMemory
     *
     * Given two logical indexes, this function returns samples of the line that
     * joins them. This is done by bilinear interpolation. The number of samples
     * in the line is N.
     */
    void profile(int x0, int y0, int xF, int yF, int N,
                 Matrix1D< double >& profile) const
    {
        profile.initZeros(N);
        double tx_step = (double)(xF - x0) / (N - 1);
        double ty_step = (double)(yF - y0) / (N - 1);
        double tx = x0, ty = y0;

        for (int i = 0; i < N; i++)
        {
            profile(i) = interpolatedElement(tx, ty);
            tx += tx_step;
            ty += ty_step;
        }
    }

    /** Get row
     * @ingroup MatricesMemory
     *
     * This function returns a row vector corresponding to the choosen row
     * inside matrix, the numbering of the rows is also logical not physical.
     *
     * @code
     * std::vector< double > v;
     * m.getRow(-2, v);
     * @endcode
     */
    void getRow(int i, Matrix1D<T>& v) const
    {
        if (XSIZE(*this) == 0 || YSIZE(*this) == 0)
        {
            v.clear();
            return;
        }

        if (i < STARTINGY(*this) || i > FINISHINGY(*this))
            REPORT_ERROR(1103,
                         "getRow: Matrix subscript (i) greater than matrix dimension");

        v.resize(XSIZE(*this));
        STARTINGX(v) = STARTINGX(*this);

        for (int j = STARTINGX(*this); j <= FINISHINGX(*this); j++)
            VEC_ELEM(v, j) = MAT_ELEM(*this, i, j);

        v.setRow();
    }

    /** Return row. The same as previous.
     * @ingroup MatricesMemory
      */
    Matrix1D<T> Row(int i) const
    {
        Matrix1D<T> aux;
        getRow(i, aux);

        return aux;
    }

    /** Get Column
     * @ingroup MatricesMemory
     *
     * This function returns a column vector corresponding to the choosen column
     * inside matrix, the numbering of the column is also logical not physical.
     *
     * @code
     * std::vector< double > v;
     * m.getCol(-1, v);
     * @endcode
     */
    void getCol(int j, Matrix1D<T>& v) const
    {
        if (XSIZE(*this) == 0 || YSIZE(*this) == 0)
        {
            v.clear();
            return;
        }

        if (j < STARTINGX(*this) || j > FINISHINGX(*this))
            REPORT_ERROR(1103,
                         "getCol: Matrix subscript (j) greater than matrix dimension");

        v.resize(YSIZE(*this));
        STARTINGX(v)  = STARTINGY(*this);

        for (int i = STARTINGY(*this); i <= FINISHINGY(*this); i++)
            VEC_ELEM(v, i) = MAT_ELEM(*this, i, j);

        v.setCol();
    }

    /** Return Column. The same as previous.
     * @ingroup MatricesMemory
     */
    Matrix1D<T> Col(int i) const
    {
        Matrix1D<T> aux;
        getCol(i, aux);

        return aux;
    }

    /** Set Row
     * @ingroup MatricesMemory
     *
     * This function sets a row vector corresponding to the choosen row inside
     * matrix, the numbering of the rows is also logical not physical.
     *
     * @code
     * m.setRow(-2, m.row(1)); // Copies row 1 in row -2
     * @endcode
     */
    void setRow(int i, const Matrix1D<T>& v)
    {
        if (XSIZE(*this) == 0 || YSIZE(*this) == 0)
            REPORT_ERROR(1, "setRow: Target matrix is empty");

        if (i < STARTINGY(*this) || i > FINISHINGY(*this))
            REPORT_ERROR(1103, "setRow: Matrix subscript (i) out of range");

        if (XSIZE(v) != XSIZE(*this))
            REPORT_ERROR(1102,
                         "setRow: Vector dimension different from matrix one");

        if (!v.isRow())
            REPORT_ERROR(1107, "setRow: Not a row vector in assignment");

        i = i - STARTINGY(*this);
        for (int j = 0; j < XSIZE(*this); j++)
            DIRECT_MAT_ELEM(*this, i, j) = DIRECT_VEC_ELEM(v, j);
    }

    /** Set Column
     * @ingroup MatricesMemory
     *
     * This function sets a column vector corresponding to the choosen column
     * inside matrix, the numbering of the column is also logical not physical.
     *
     * @code
     * m.setCol(-1, (m.row(1)).transpose()); // Copies row 1 in column -1
     * @endcode
     */
    void setCol(int j, const Matrix1D<T>& v)
    {
        if (XSIZE(*this) == 0 || YSIZE(*this) == 0)
            REPORT_ERROR(1, "setCol: Target matrix is empty");

        if (j < STARTINGX(*this) || j > FINISHINGX(*this))
            REPORT_ERROR(1103, "setCol: Matrix subscript (j) out of range");

        if (XSIZE(v) != YSIZE(*this))
            REPORT_ERROR(1102,
                         "setCol: Vector dimension different from matrix one");

        if (!v.isCol())
            REPORT_ERROR(1107, "setCol: Not a column vector in assignment");

        j = j - STARTINGX(*this);
        for (int i = 0; i < YSIZE(*this); i++)
            DIRECT_MAT_ELEM(*this, i, j) = DIRECT_VEC_ELEM(v, i);
    }

    /** Logical to physical index translation
     * @ingroup MatricesMemory
     *
     * This function returns the physical position of a logical one.
     *
     * @code
     * m.toPhysical(i_log, j_log, i_phys, j_phys);
     * @endcode
     */
    void toPhysical(int i_log, int j_log, int& i_phys, int& j_phys) const
    {
        i_phys = i_log - STARTINGY(*this);
        j_phys = j_log - STARTINGX(*this);
    }

    /** Physical to logical index translation
     * @ingroup MatricesMemory
     *
     * This function returns the logical position of a physical one.
     *
     * @code
     * m.toLogical(i_phys, j_phys, i_log, j_log);
     * @endcode
     */
    void toLogical(int i_phys, int j_phys, int &i_log, int& j_log) const
    {
        i_log = i_phys + STARTINGY(*this);
        j_log = j_phys + STARTINGX(*this);
    }

    /// @defgroup MatrixOperators Operators
    /// @ingroup Matrices

    /** Assignment.
     * @ingroup MatrixOperators
     *
     * You can build as complex assignment expressions as you like. Multiple
     * assignment is allowed.
     *
     * @code
     * v1 = v2 + v3;
     * v1 = v2 = v3;
     * @endcode
     */
    Matrix2D<T>& operator=(const Matrix2D<T>& op1)
    {
	if (&op1 != this)
	{
            resize(op1);
            T* ptr=NULL;
	    unsigned long int n;
            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
	    	*ptr=DIRECT_MULTIDIM_ELEM(op1,n);
	}

	return *this;
    }

    /** Unary minus.
     * @ingroup MatrixOperators
     *
     * It is used to build arithmetic expressions. You can make a minus
     * of anything as long as it is correct semantically.
     */
    Matrix2D<T> operator-() const
    {
        Matrix2D<T> tmp(*this);
	T* ptr;
	unsigned long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(tmp,n,ptr)
            *ptr = -(*ptr);
        return tmp;
    }

    /** v3 = v1 + v2.
     * @ingroup MatrixOperators
     */
    Matrix2D<T> operator+(const Matrix2D<T>& op1) const
    {
        Matrix2D<T> tmp;
        arrayByArray(*this, op1, tmp, '+');
        return tmp;
    }

    /** v3 = v1 - v2.
     * @ingroup MatrixOperators
     */
    Matrix2D<T> operator-(const Matrix2D<T>& op1) const
    {
        Matrix2D<T> tmp;
        arrayByArray(*this, op1, tmp, '-');
        return tmp;
    }

    /** v3 = v1 * v2.
     * @ingroup MatrixOperators
     */
    Matrix2D<T> operator*(const Matrix2D<T>& op1) const
    {
        Matrix2D<T> tmp;
        multiplyMatrix(*this, op1, tmp);
        return tmp;
    }

    /** v3 = v1 / v2.
     * @ingroup MatrixOperators
     */
    Matrix2D<T> operator/(const Matrix2D<T>& op1) const
    {
        Matrix2D<T> tmp;
        arrayByArray(*this, op1, tmp, '/');
        return tmp;
    }

    /** v3 += v2.
     * @ingroup MatrixOperators
     */
    void operator+=(const Matrix2D<T>& op1)
    {
        arrayByArray(*this, op1, *this, '+');
    }

    /** v3 -= v2.
     * @ingroup MatrixOperators
     */
    void operator-=(const Matrix2D<T>& op1)
    {
        arrayByArray(*this, op1, *this, '-');
    }

    /** v3 *= v2.
     * @ingroup MatrixOperators
     */
    void operator*=(const Matrix2D<T>& op1)
    {
    	Matrix2D<T> tmp;
        multiplyElements(*this, op1, tmp);
	*this=tmp;
    }

    /** v3 /= v2.
     * @ingroup MatrixOperators
     */
    void operator/=(const Matrix2D<T>& op1)
    {
        arrayByArray(*this, op1, *this, '/');
    }

    /** v3 = v1 + k.
     * @ingroup MatrixOperators
     */
    Matrix2D<T> operator+(T op1) const
    {
        Matrix2D<T> tmp;
        arrayByScalar(*this, op1, tmp, '+');
        return tmp;
    }

    /** v3 = v1 - k.
     * @ingroup MatrixOperators
     */
    Matrix2D<T> operator-(T op1) const
    {
        Matrix2D<T> tmp;
        arrayByScalar(*this, op1, tmp, '-');
        return tmp;
    }

    /** v3 = v1 * k.
     * @ingroup MatrixOperators
     */
    Matrix2D<T> operator*(T op1) const
    {
        Matrix2D<T> tmp;
        arrayByScalar(*this, op1, tmp, '*');
        return tmp;
    }

    /** v3 = v1 / k.
     * @ingroup MatrixOperators
     */
    Matrix2D<T> operator/(T op1) const
    {
        Matrix2D<T> tmp;
        arrayByScalar(*this, op1, tmp, '/');
        return tmp;
    }

    /** v3 += k.
     * @ingroup MatrixOperators
     *
     * This function is not ported to Python.
     */

    void operator+=(const T& op1)
    {
        arrayByScalar(*this, op1, *this, '+');
    }

    /** v3 -= k.
     * @ingroup MatrixOperators
     *
     * This function is not ported to Python.
     */
    void operator-=(const T& op1)
    {
        arrayByScalar(*this, op1, *this, '-');
    }

    /** v3 *= k.
     * @ingroup MatrixOperators
     *
     * This function is not ported to Python.
     */
    void operator*=(const T& op1)
    {
        arrayByScalar(*this, op1, *this, '*');
    }

    /** v3 /= k.
     * @ingroup MatrixOperators
     *
     * This function is not ported to Python.
     */
    void operator/=(const T& op1)
    {
        arrayByScalar(*this, op1, *this, '/');
    }

    /** v3 = k + v2.
     * @ingroup MatrixOperators
     */
    friend Matrix2D<T> operator+(T op1, const Matrix2D<T>& op2)
    {
        Matrix2D<T> tmp;
        scalarByArray(op1, op2, tmp, '+');
        return tmp;
    }

    /** v3 = k - v2.
     * @ingroup MatrixOperators
     */
    friend Matrix2D<T> operator-(T op1, const Matrix2D<T>& op2)
    {
        Matrix2D<T> tmp;
        scalarByArray(op1, op2, tmp, '-');
        return tmp;
    }

    /** v3 = k * v2.
     * @ingroup MatrixOperators
     */
    friend Matrix2D<T> operator*(T op1, const Matrix2D<T>& op2)
    {
        Matrix2D<T> tmp;
        scalarByArray(op1, op2, tmp, '*');
        return tmp;
    }

    /** v3 = k / v2
     * @ingroup MatrixOperators
     */
    friend Matrix2D<T> operator/(T op1, const Matrix2D<T>& op2)
    {
        Matrix2D<T> tmp;
        scalarByArray(op1, op2, tmp, '/');
        return tmp;
    }

    /// @defgroup MatricesUtilities Utilities
    /// @ingroup Matrices

    /** Algebraic transpose of matrix
     * @ingroup MatricesUtilities
     *
     * You can use the transpose in as complex expressions as you like. The
     * origin of the matrix is not changed.
     *
     * @code
     * m2 = m1.transpose();
     * @endcode
     */
    Matrix2D<T> transpose() const
    {
        T aux;
        Matrix2D<T> result(XSIZE(*this), YSIZE(*this));

        FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(result)
            DIRECT_MAT_ELEM(result, i, j) = DIRECT_MAT_ELEM(*this, j, i);

        STARTINGX(result) = STARTINGX(*this);
        STARTINGY(result) = STARTINGY(*this);

        return result;
    }

    /** Determinant of a matrix
     * @ingroup MatricesUtilities
     *
     * An exception is thrown if the matrix is not squared or it is empty.
     *
     * @code
     * double det = m.det();
     * @endcode
     */
    T det() const
    {
        // (see Numerical Recipes, Chapter 2 Section 5)
        if (XSIZE(*this) == 0)
            REPORT_ERROR(1108, "determinant: Matrix is empty");

        if (XSIZE(*this) != YSIZE(*this))
            REPORT_ERROR(1109, "determinant: Matrix is not squared");

        for (int i = 0; i < YSIZE(*this); i++)
        {
            bool all_zeros = true;
            for (int j = 0; j < XSIZE(*this); j++)
                if (ABS(DIRECT_MAT_ELEM(*this, i, j)) > XMIPP_EQUAL_ACCURACY)
                {
                    all_zeros = false;
                    break;
                }

            if (all_zeros)
                return 0;
        }

        // Perform decomposition
        Matrix1D< int > indx;
        T d;
        Matrix2D<T> LU;
        ludcmp(*this, LU, indx, d);

        // Calculate determinant
        for (int i = 0; i < XSIZE(*this); i++)
            d *= (T) LU(i , i);

        return d;
    }

    /** Inverse of a matrix
     * @ingroup MatricesUtilities
     *
     * The matrix is inverted using a SVD decomposition. In fact the
     * pseudoinverse is returned.
     *
     * @code
     * Matrix2D< double > m1_inv;
     * m1.inv(m1_inv);
     * @endcode
     */
    void inv(Matrix2D<T>& result) const
    {
        if (XSIZE(*this) == 0)
            REPORT_ERROR(1108, "Inverse: Matrix is empty");

        // Perform SVD decomposition
        Matrix2D< double > u, v;
        Matrix1D< double > w;
        svdcmp(*this, u, w, v); // *this = U * W * V^t

        double tol = MultidimArray<T>::computeMax() *
            XMIPP_MAX(XSIZE(*this),YSIZE(*this)) * 1e-14;
        result.initZeros(XSIZE(*this), YSIZE(*this));

        // Compute W^-1
        bool invertible = false;
        for (int i = 0; i < XSIZE(w); i++)
        {
            if (ABS(DIRECT_VEC_ELEM(w, i)) > tol)
            {
                DIRECT_VEC_ELEM(w, i) = 1.0 / DIRECT_VEC_ELEM(w, i);
                invertible = true;
            }
            else
                DIRECT_VEC_ELEM(w, i) = 0.0;
        }

        if (!invertible)
            return;

        // Compute V*W^-1
        FOR_ALL_ELEMENTS_IN_MATRIX2D(v)
            DIRECT_MAT_ELEM(v, i, j) *= DIRECT_VEC_ELEM(w, j);

        // Compute Inverse
        for (int i = 0; i < XSIZE(*this); i++)
            for (int j = 0; j < YSIZE(*this); j++)
                for (int k = 0; k < XSIZE(*this); k++)
                    DIRECT_MAT_ELEM(result, i, j) += (T)
                        (DIRECT_MAT_ELEM(v, i, k) * DIRECT_MAT_ELEM(u, j, k));
    }

    /** Inverse of a matrix
     * @ingroup MatricesUtilities
     */
    Matrix2D<T> inv() const
    {
        Matrix2D<T> result;
        inv(result);

        return result;
    }

    /** Solve equation system
     * @ingroup MatricesUtilities
     *
     * The equation system is defined by Ax=b, it is solved for x. Exceptions
     * are thrown if the equation system is not solvable.
     *
     * @code
     * Matrix1D< double > x = m.inv();
     * @endcode
     */
    friend void solve<>(const Matrix2D<T>& A, const Matrix1D<T>& b,
        Matrix1D<T>& result);

    /** Solve equation system.
     * @ingroup MatricesUtilities
     *
     * The same as before but now, b is a matrix and so, x is also a matrix. A
     * must be a square matrix.
     */
    friend void solve<>(const Matrix2D<T>& A, const Matrix2D<T>& b,
        Matrix2D<T>& result);

    /** Put a window to matrix
     * @ingroup MatricesUtilities
     *
     * The matrix is windowed within the two positions given to this function.
     * Indexes always refer to logical indexes. If a position is outside the
     * actual matrix range then the matrix is padded with init_value until the
     * new position is reached. In the following examples suppose that m1 is the
     * following and that the origin is (-1,-1).
     *
     * @code
     *      [1 2 3               [1 2 3 0
     * m1 =  4 5 6    --->  m1 =  4 5 6 0
     *       7 8 9]               7 8 9 0]
     *
     * @endcode
     *
     * @code
     * m1.window(-1, -1, 1, 2);
     * @endcode
     */
    void window(int y0, int x0, int yF, int xF, T init_value = 0)
    {
        Matrix2D<T> result(yF - y0 + 1, xF - x0 + 1);
        STARTINGY(result) = y0;
        STARTINGX(result) = x0;

        FOR_ALL_ELEMENTS_IN_MATRIX2D(result)
        if (j >= STARTINGX(*this) && j <= FINISHINGX(*this) &&
            i >= STARTINGY(*this) && i <= FINISHINGY(*this))
            MAT_ELEM(result, i, j) = MAT_ELEM(*this, i, j);
        else
            MAT_ELEM(result, i, j) = init_value;

        *this = result;
    }

    /** Computes the center of mass of an image.
     * @ingroup MatricesUtilities
     */
    void centerOfMass(Matrix1D< double >& center, void* mask=NULL)
    {
		center.initZeros(2);
		double mass = 0;
		Matrix2D< int>* imask = (Matrix2D< int >*) mask;

		FOR_ALL_ELEMENTS_IN_MATRIX2D(*this)
		{
            	if ((imask == NULL || MAT_ELEM(*imask, i, j)) &&
			MAT_ELEM(*this, i, j) > 0)
            	{
        		XX(center) += j * MAT_ELEM(*this, i, j);
        		YY(center) += i * MAT_ELEM(*this, i, j);
        		mass += MAT_ELEM(*this, i, j);
            	}
		}

		if (mass != 0)
            	center /= mass;
    }

    /** Adjust the range of the array to a given one within a mask.
     * @ingroup MatricesUtilities
     *
     * A linear operation is performed on the values of the array such that
     * after it, the values of the array are comprissed between the two values
     * set. The actual array is modified itself. The linear transformation
	 * is computed within the mask, but it is applied everywhere.
     *
     * @code
     * v.rangeAdjust(0, 1, mask);
     * // The array is now ranging from 0 to 1
     * @endcode
     */
    // This function must be explictly implemented outside
    void rangeAdjust(T minF, T maxF, Matrix2D<int> &mask)
    {
        if (MULTIDIM_SIZE(*this) <= 0)
            return;

        double min0, max0;
		bool first=true;
		FOR_ALL_ELEMENTS_IN_MATRIX2D(*this)
		{
			if (mask(i,j))
			{
				T val=(*this)(i,j);
				if (first)
				{
					min0=max0=(double)val;
					first=false;
				}
				else
				{
					min0=XMIPP_MIN(min0,val);
					max0=XMIPP_MAX(max0,val);
				}
			}
		}

        // If max0==min0, it means that the vector is a constant one, so the
        // only possible transformation is to a fixed minF
        double slope;
        if (max0 != min0)
            slope = static_cast< double >(maxF - minF) /
                    static_cast< double >(max0 - min0);
        else
            slope = 0;

        T* ptr=NULL;
	    unsigned long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
            *ptr = minF + static_cast< T >(slope *
                static_cast< double >(*ptr - min0));
    }

    /** Adjust the range of the array to a given one.
     * @ingroup MatricesUtilities
     *
     * A linear operation is performed on the values of the array such that
     * after it, the values of the array are comprissed between the two values
     * set. The actual array is modified itself
     *
     * @code
     * v.rangeAdjust(0, 1);
     * // The array is now ranging from 0 to 1
     * @endcode
     */
    void rangeAdjust(T minF, T maxF)
    {
        MultidimArray<T>::rangeAdjust(minF,maxF);
    }

    /** Compute statistics.
     * @ingroup MatricesUtilities
     */
    void computeStats(double& avg, double& stddev, T& minval, T& maxval) const
    {
    	MultidimArray<T>::computeStats(avg, stddev, minval, maxval);
    }
    
    /** Compute statistics within region.
     * @ingroup MatricesUtilities
     *
     * The region is specified by two corners.
     */
    void computeStats(double& avg,
                       double& stddev,
                       T& min_val,
                       T& max_val,
                       const Matrix1D< int >& corner1,
                       const Matrix1D< int >& corner2) const
    {
	min_val = max_val = (*this)(corner1);

	Matrix1D< double > r(3);
	double N = 0, sum = 0, sum2 = 0;

	FOR_ALL_ELEMENTS_IN_MATRIX2D_BETWEEN(corner1, corner2)
	{
            sum += (*this)(r);
            sum2 += (*this)(r) * (*this)(r);
            N++;

            if ((*this)(r) < min_val)
        	min_val = (*this)(r);
            else if ((*this)(r) > max_val)
        	max_val = (*this)(r);
	}

	if (N != 0)
	{
            avg = sum / N;
            stddev = sqrt(sum2 / N - avg * avg);
	}
	else
	{
            avg = stddev = 0;
	}
    }

    /** @defgroup MatricesGeometrical Geometrical Transformations
     * @ingroup Matrices
     *
     * In all geometrical transformations a periodic extension of the matrix is
     * supposed, ie, if a pixel goes out on the left, it is entering on the
     * right, ...
     */

    /** Applies a geometrical transformation
     * @ingroup MatricesGeometrical
     *
     * Any geometrical transformation defined by the matrix A (double (3x3)!!
     * ie, in homogeneous R2 coordinates) is applied to the matrix M1. The
     * result is stored in m2 (it cannot be the same as the input matrix). An
     * exception is thrown if the transformation matrix is not 3x3. The result
     * matrix is resized to the same dimensions as M1 if M2 is empty (0x0) at
     * the beginning, if it is not, ie, if M2 has got some size then only those
     * values in the matrix are filled, this is very useful for resizing the
     * matrix, then you manually resize the output matrix to the desired size
     * and then call this routine.
     *
     * The relationship between the output coordinates and the input ones are
     *
     * @code
     * out = A * in
     * (x,y) = A * (x',y')
     * @endcode
     *
     * This function works independently from the logical indexing of each
     * matrix, it sets the logical center and the physical center of the image
     * and work with these 2 coordinate spaces. At the end the original logical
     * indexing of each matrix is kept.
     *
     * The procedure followed goes from coordinates in the output matrix to the
     * ones in the input one, so the inverse of the A matrix is needed. There is
     * a flag telling if the given transformation matrix is already the inverse
     * one or the normal one. If it is the normal one internally the matrix is
     * inversed. If you are to do many "rotations" then some time is spent in
     * inverting the matrix. Normally the matrix is the normal one.
     *
     * There is something else to tell about the geometrical tranformation. The
     * value of the pixel in the output matrix is computed via bilinear
     * interpolation in the input matrix. If any of the pixels participating in
     * the interpolation falls outside the input matrix, then automatically the
     * corresponding output pixel is set to 0, unless that the wrap flag has
     * been set to 1. In this case if the pixel falls out by the right hand then
     * it is "wrapped" and the corresponding pixel in the left hand is used. The
     * same is appliable to top-bottom. Usually wrap mode is off. Wrap mode is
     * interesting for translations but not for rotations, for example.
     *
     * The inverse mode and wrapping mode should be taken by default by the
     * routine, g++ seems to have problems with template functions outside a
     * class with default parameters. So, I'm sorry, you will have to put them
     * always. The usual combination is applyGeometry(..., IS_NOT_INV, DONT_WRAP).
     * Although you can also use the constants IS_INV, or WRAP.
     *
     * m2 cannot be the same matrix as m1.
     *
     * @code
     * Matrix2D< double > A(3, 3);
     * A.initIdentity;
     * applyGeometry(m2, A, m1);
     * @endcode
     */
    friend void applyGeometry<>(Matrix2D<T>& m2, Matrix2D< double > A,
                             const Matrix2D<T>& m1, bool inv, bool wrap,
                             T outside);

    /** Apply geom with B-spline interpolation
     * @ingroup MatricesGeometrical
     */
    friend void applyGeometryBSpline<>(Matrix2D<T>& m2, Matrix2D< double > A,
                                     const Matrix2D<T>& m1, int Splinedegree,
                                     bool inv, bool wrap, T outside);

    /** Self apply geom
     * @ingroup MatricesGeometrical
     *
     * Same as the previous one, but the result is kept in this object
     */
    void selfApplyGeometry(Matrix2D< double > A, bool inv, bool wrap,
                         T outside = (T) 0)
    {
        Matrix2D<T> aux;
        applyGeometry(aux, A, *this, inv, wrap, outside);

        *this = aux;
    }

    /** Self apply geom with Bspline interpolation
     *@ingroup MatricesGeometrical
     */
    void selfApplyGeometryBSpline(Matrix2D< double > A, int SplineDegree,
                                 bool inv, bool wrap, T outside = (T) 0)
    {
        Matrix2D<T> aux;
        applyGeometryBSpline(aux, A, *this, SplineDegree, inv, wrap, outside);

        *this = aux;
    }

#define IS_INV true
#define IS_NOT_INV false
#define DONT_WRAP false
#define WRAP true

    /** Rotate matrix
     * @ingroup MatricesGeometrical
     *
     * The angle must be in degrees. The result cannot be this object.
     *
     * @code
     * m1.rotate(60, m2);
     * @endcode
     */
    void rotate(double ang, Matrix2D<T>& result, bool wrap = DONT_WRAP) const
    {
        Matrix2D< double > temp = rotation2DMatrix(ang);
        applyGeometry(result, temp, *this, IS_NOT_INV, wrap);
    }

    /** Rotate matrix
     * @ingroup MatricesGeometrical
     *
     * Same as the previous one, but the result is kept in this object
     */
    void selfRotate(double ang, bool wrap = DONT_WRAP)
    {
        Matrix2D<T> aux;
        rotate(ang, aux, wrap);
        *this = aux;
    }

    /** Rotate matrix (using Bspline interpolation)
     * @ingroup MatricesGeometrical
     *
     * The angle must be in degrees. The result cannot be this object.
     *
     * @code
     * m1.rotate(60, m2);
     * @endcode
     */
    void rotateBSpline(int Splinedegree,
                        double ang, Matrix2D<T>& result, bool wrap = DONT_WRAP) const
    {
        Matrix2D< double > temp = rotation2DMatrix(ang);
        applyGeometryBSpline(result, temp, *this, Splinedegree, IS_NOT_INV, wrap);
    }

    /** Rotate matrix (using Bspline interpolation).
     * @ingroup MatricesGeometrical
     *
     * Same as the previous one, but the result is kept in this object
     */
    void selfRotateBSpline(int Splinedegree, double ang, bool wrap = DONT_WRAP)
    {
        Matrix2D<T> aux;
        rotateBSpline(Splinedegree, ang, aux, wrap);
        *this = aux;
    }

    /** Translate matrix
     * @ingroup MatricesGeometrical
     *
     * The displacement is given as a R2 vector of the form (shift_X,shift_Y).
     * The result cannot be this object.
     *
     * @code
     * // m1 is shifted 2 pixels down and stored in m2
     * m2 = m1.translate(vectorR2(0, 2));
     * @endcode
     */
    void translate(const Matrix1D< double >& v, Matrix2D<T>& result,
                   bool wrap = WRAP) const
    {
        Matrix2D< double > temp = translation2DMatrix(v);
        applyGeometry(result, temp, *this, IS_NOT_INV, wrap);
    }

    /** Translate matrix
     * @ingroup MatricesGeometrical
     *
     * Same as the previous one, but the result is kept in this object
     */
    void selfTranslate(const Matrix1D< double >& v, bool wrap = WRAP)
    {
        Matrix2D<T> aux;
        translate(v, aux, wrap);
        *this = aux;
    }

    /** Translate matrix (using Bspline interpolation)
     * @ingroup MatricesGeometrical
     *
     * The displacement is given as a R2 vector of the form (shift_X,shift_Y).
     * The result cannot be this object.
     *
     * @code
     * // m1 is shifted 2 pixels down and stored in m2
     * m2 = m1.translate(vectorR2(0, 2));
     * @endcode
     */
    void translateBSpline(int Splinedegree,
                           const Matrix1D< double >& v, Matrix2D<T>& result,
                           bool wrap = WRAP) const
    {
        Matrix2D< double > temp = translation2DMatrix(v);
        applyGeometryBSpline(result, temp, *this, Splinedegree, IS_NOT_INV, wrap);
    }

    /** Translate matrix (using Bspline interpolation)
     * @ingroup MatricesGeometrical
     *
     * Same as the previous one, but the result is kept in this object
     */
    void selfTranslateBSpline(int Splinedegree,
                              const Matrix1D< double >& v, bool wrap = WRAP)
    {
        Matrix2D<T> aux;
        translateBSpline(Splinedegree, v, aux, wrap);
        *this = aux;
    }

    /** Translate center of mass to center
     * @ingroup MatricesGeometrical
     *
     * If the input has very high values, it is better to rescale it to be
     * between 0 and 1.
     */
    void selfTranslateCenterOfMassToCenter(bool wrap = WRAP)
    {
        MultidimArray<T>::setXmippOrigin();
        Matrix1D< double > center;
        centerOfMass(center);
        center *= -1;
        selfTranslate(center, wrap);
    }

    /** Translate center of mass to center (using Bspline interpolation)
     * @ingroup MatricesGeometrical
     *
     * If the input has very high values, it is better to rescale it to be
     * between 0 and 1.
     */
    void selfTranslateCenterOfMassToCenterBSpline(
        int Splinedegree, bool wrap = WRAP)
    {
        MultidimArray<T>::setXmippOrigin();
        Matrix1D< double > center;
        centerOfMass(center);
        center *= -1;
        selfTranslateBSpline(Splinedegree, center, wrap);
    }

    /** Scales to a new size
     * @ingroup MatricesGeometrical
     *
     * The matrix is scaled (resampled) to fill a new size. It is not the same
     * as "window" in this same class. The size can be larger or smaller than
     * the actual one. But the result matrix cannot be this object.
     *
     * @code
     * m1.scaleToSize(128, 128);
     * @endcode
     */
    void scaleToSize(int Ydim, int Xdim, Matrix2D<T>& result) const
    {
        Matrix2D< double > temp(3, 3);
        result.resize(Ydim, Xdim);
        temp.initIdentity();

        DIRECT_MAT_ELEM(temp, 0, 0) = (double) Xdim / (double) XSIZE(*this);
        DIRECT_MAT_ELEM(temp, 1, 1) = (double) Ydim / (double) YSIZE(*this);

        applyGeometry(result, temp, *this, IS_NOT_INV, WRAP);
    }

    /** Scales to a new size
     * @ingroup MatricesGeometrical
     *
     * Same as the previous one, but the result is kept in this object.
     */
    void selfScaleToSize(int Ydim, int Xdim)
    {
        Matrix2D<T> aux;
        scaleToSize(Ydim, Xdim, aux);
        *this = aux;
    }

    /** Scales to a new size (using Bspline interpolation)
     * @ingroup MatricesGeometrical
     *
     * The matrix is scaled (resampled) to fill a new size. It is not the same
     * as "window" in this same class. The size can be larger or smaller than
     * the actual one. But the result matrix cannot be this object.
     *
     * @code
     * m1.scaleToSize(128, 128);
     * @endcode
     */
    void scaleToSizeBSpline(int Splinedegree,
                               int Ydim, int Xdim, Matrix2D<T>& result) const
    {
        Matrix2D< double > temp(3, 3);
        result.resize(Ydim, Xdim);
        temp.initIdentity();

        DIRECT_MAT_ELEM(temp, 0, 0) = (double) Xdim / (double) XSIZE(*this);
        DIRECT_MAT_ELEM(temp, 1, 1) = (double) Ydim / (double) YSIZE(*this);

        applyGeometryBSpline(result, temp, *this, Splinedegree, IS_NOT_INV, WRAP);
    }

    /** Scales to a new size (using Bspline interpolation)
     * @ingroup MatricesGeometrical
     *
     * Same as the previous one, but the result is kept in this object.
     */
    void selfScaleToSizeBSpline(int Splinedegree, int Ydim, int Xdim)
    {
        Matrix2D<T> aux;
        scaleToSizeBSpline(Splinedegree, Ydim, Xdim, aux);
        *this = aux;
    }

    /** Superpixel reduce
     * @ingroup MatricesGeometrical
     *
     * This function reduces the given image averaging in the superpixel of the
     * given size. The last columns and rows are removed if the size of the
     * original image is not an exact multiple of the given superpixel size
     */
    void superpixelReduce(Matrix2D<T>& result, int size = 2) const
    {
        result.initZeros(YSIZE(*this) / size, XSIZE(*this) / size);
        int size2 = size * size;

        FOR_ALL_ELEMENTS_IN_MATRIX2D(result)
        {
            for (int ii = 0; ii < size; ii++)
                for (int jj = 0; jj < size; jj++)
                    DIRECT_MAT_ELEM(result, i, j) +=
                        DIRECT_MAT_ELEM(*this, size * i + ii, size * j + jj);

            DIRECT_MAT_ELEM(result, i, j) /= size2;
        }
    }

    /** Superpixel expand
     * @ingroup MatricesGeometrical
     *
     * This function copies each pixel to a new image as many times as the size
     * of the superpixel.
     */
    void superpixelExpand(Matrix2D<T>& result, int size = 2) const
    {
        result.initZeros(YSIZE(*this) * size, XSIZE(*this) * size);

        FOR_ALL_ELEMENTS_IN_MATRIX2D(*this)
        {
            for (int ii = 0; ii < size; ii++)
                for (int jj = 0; jj < size; jj++)
                    DIRECT_MAT_ELEM(result, size * i + ii, size * j + jj) =
                        DIRECT_MAT_ELEM(*this, i, j);
        }
    }

    /** Reduce the image by 2 using a BSpline pyramid
     * @ingroup MatricesGeometrical
     */
    void pyramidReduce(Matrix2D< double >& result, int levels = 1) const
    {
        Matrix2D< double > aux, aux2;
        MultidimArray<T>::produceSplineCoefficients(aux, 3);

        for (int i = 0; i < levels; i++)
        {
            aux.reduceBSpline(aux2, 3);
            aux = aux2;
        }

        aux2.produceImageFromSplineCoefficients(result, 3);
    }

    /** Expand the image by 2 using a BSpline pyramid
     * @ingroup MatricesGeometrical
     */
    void pyramidExpand(Matrix2D< double >& result, int levels = 1) const
    {
        Matrix2D< double > aux, aux2;
        MultidimArray<T>::produceSplineCoefficients(aux, 3);

        for (int i = 0; i < levels; i++)
        {
            aux.expandBSpline(aux2, 3);
            aux = aux2;
        }

        aux2.produceImageFromSplineCoefficients(result, 3);
    }

    /** Expand a set of B-spline coefficients
     * @ingroup MatricesGeometrical
     *
     * Knowing that this matrix is a set of B-spline coefficients, produce the
     * expanded set of B-spline coefficients using the two-scale relationship.
     */
    void expandBSpline(Matrix2D< double >& expanded, int SplineDegree = 3) const
    {
        double g[200]; // Coefficients of the reduce filter
        long ng; // Number of coefficients of the reduce filter
        double h[200]; // Coefficients of the expansion filter
        long nh; // Number of coefficients of the expansion filter
        short IsCentered; // Equal TRUE if the filter is a centered spline

        // Get the filter
        if (GetPyramidFilter("Centered Spline", SplineDegree,
                             g, &ng, h, &nh, &IsCentered))
            REPORT_ERROR(1, "Unable to load the filter coefficients");

        Matrix2D< double > aux;
        typeCast(*this, aux);
        expanded.resize(2 * YSIZE(aux), 2 * XSIZE(aux));

        Expand_2D(MULTIDIM_ARRAY(aux), XSIZE(aux), YSIZE(aux),
                  MULTIDIM_ARRAY(expanded), h, nh, IsCentered);
    }

    /** Reduce a set of B-spline coefficients
     * @ingroup MatricesGeometrical
     *
     * Knowing that this matrix is a set of B-spline coefficients, produce the
     * reduced set of B-spline coefficients using the two-scale relationship.
     */
    void reduceBSpline(Matrix2D< double >& reduced, int SplineDegree = 3) const
    {
        double g[200]; // Coefficients of the reduce filter
        long ng; // Number of coefficients of the reduce filter
        double h[200]; // Coefficients of the expansion filter
        long nh; // Number of coefficients of the expansion filter
        short IsCentered; // Equal TRUE if the filter is a centered spline

        // Get the filter
        if (GetPyramidFilter("Centered Spline", SplineDegree,
                             g, &ng, h, &nh, &IsCentered))
            REPORT_ERROR(1, "Unable to load the filter coefficients");

        Matrix2D< double>  aux;
        typeCast(*this, aux);

        if (XSIZE(aux) % 2 != 0 && YSIZE(aux) % 2 != 0)
            aux.resize(YSIZE(aux) - 1, XSIZE(aux) - 1);
        else if (YSIZE(aux) % 2 != 0)
            aux.resize(YSIZE(aux) - 1, XSIZE(aux));
        else if (XSIZE(aux) % 2 != 0)
            aux.resize(YSIZE(aux), XSIZE(aux) - 1);

        reduced.resize(YSIZE(aux) / 2, XSIZE(aux) / 2);

        Reduce_2D(MULTIDIM_ARRAY(aux), XSIZE(aux), YSIZE(aux),
                  MULTIDIM_ARRAY(reduced), g, ng, IsCentered);
    }

    /** @defgroup MatricesAlgebraic Algebraic operations
     * @ingroup Matrices
     *
     * NOTICE!!!: the matrix by matrix multiplier operator (*) has been
     * redefined in this class and it represents the true matrix by matrix
     * algebraic multiplication, if you want an element-wise multiplication you
     * have to use the function defined in this section "multiplyElements"
     */

    /** Matrix multiplication element by element
     * @ingroup MatricesAlgebraic
     *
     * @code
     * multiplyElements(m1, m2, mResult);
     * @endcode
     */
    friend void multiplyElements(const Matrix2D<T>& op1, const Matrix2D<T>& op2,
        Matrix2D<T>& result)
    {
        arrayByArray(op1, op2, result, '*');
    }

    /** Matrix by vector multiplication
     * @ingroup MatricesAlgebraic
     *
     * @code
     * v2 = A*v1;
     * @endcode
     */
    Matrix1D<T> operator*(const Matrix1D<T>& op1) const
    {
        Matrix1D<T> result;

        if (XSIZE(*this) != XSIZE(op1))
            REPORT_ERROR(1102, "Not compatible sizes in matrix by vector");

        if (!op1.isCol())
            REPORT_ERROR(1102, "Vector is not a column");

        result.initZeros(YSIZE(*this));

        for (int i = 0; i < YSIZE(*this); i++)
            for (int j = 0; j < XSIZE(op1); j++)
                DIRECT_VEC_ELEM(result, i) += DIRECT_MAT_ELEM(*this, i, j) *
                                              DIRECT_VEC_ELEM(op1, j);

        result.setCol();
        STARTINGX(result) = STARTINGY(*this);

        return result;
    }

    /// @defgroup MatricesTypes Matrix types
    /// @ingroup Matrices

    /** True if the matrix is diagonal
     * @ingroup MatricesTypes
     *
     * @code
     * if (m.isDiagonal())
     *     std::cout << "The matrix is diagonal\n";
     * @endcode
     */
    bool isDiagonal() const
    {
        if (XSIZE(*this) != YSIZE(*this))
            return false;

        FOR_ALL_ELEMENTS_IN_MATRIX2D(*this)
        if (i != j && ABS(DIRECT_MAT_ELEM(*this, i, j)) >
            XMIPP_EQUAL_ACCURACY)
            return false;

        return true;
    }

    /** True if the matrix is scalar
     * @ingroup MatricesTypes
     *
     * A scalar matrix is diagonal and all the values at the diagonal are the
     * same
     */
    bool isScalar() const
    {
        if (!isDiagonal())
            return false;

        for (int i = 1; i < YSIZE(*this); i++)
            if (ABS(DIRECT_MAT_ELEM(*this, i, i) - DIRECT_MAT_ELEM(*this, 0, 0))
                > XMIPP_EQUAL_ACCURACY)
                return false;

        return true;
    }

    /** True if the matrix is identity
     * @ingroup MatricesTypes
     *
     * @code
     * if (m.isIdentity())
     *     std::cout << "The matrix is identity\n";
     * @endcode
     */
    bool isIdentity() const
    {
        return isScalar() &&
               ABS(DIRECT_MAT_ELEM(*this, 0, 0) - (T) 1) < XMIPP_EQUAL_ACCURACY;
    }

    /** Maximum element
     * @ingroup MatricesTypes
     *
     * This function returns the index of the maximum element of an array(i,j).
     * Returns -1 if the array is empty
     */
    void maxIndex(int& imax, int& jmax) const
    {
        if (XSIZE(*this) == 0)
        {
            imax = jmax = -1;
            return;
        }

        imax = jmax = 0;
        T maxval = MAT_ELEM(*this, imax, jmax);

        FOR_ALL_ELEMENTS_IN_MATRIX2D(*this)
            if (MAT_ELEM(*this, i, j) > maxval)
            {
                maxval = MAT_ELEM(*this, i, j);
                imax = i;
                jmax = j;
            }
    }

    /** Minimum element
     * @ingroup MatricesTypes
     *
     * This function returns the index of the minimum element of an array(i,j).
     * Returns -1 if the array is empty
     */
    void minIndex(int& imin, int& jmin) const
    {
        if (XSIZE(*this) == 0)
        {
            imin = jmin = -1;
            return;
        }

        imin = jmin = 0;
        T minval = MAT_ELEM(*this, imin, jmin);

        FOR_ALL_ELEMENTS_IN_MATRIX2D(*this)
            if (MAT_ELEM(*this, i, j) > minval)
            {
                minval = MAT_ELEM(*this, i, j);
                imin = i;
                jmin = j;
            }
    }
};

// Specializations case for complex numbers
template<>
std::complex< double > Matrix2D< std::complex< double > >::interpolatedElement(
    double x, double y, std::complex< double > outside_value) const;

/** @defgroup MatricesRelated Related functions
 * @ingroup Matrices
 *
 * These functions are not methods of Matrix2D
 */

/// @defgroup MatricesGeometry Geometry with matrices
/// @ingroup MatricesRelated

/** Creates a rotational matrix (3x3) for images
 * @ingroup MatricesGeometry
 *
 * The rotation angle is in degrees.
 *
 * @code
 * m = rotation2DMatrix(60);
 * @endcode
 */
Matrix2D< double > rotation2DMatrix(double ang);

/** Creates a translational matrix (3x3) for images
 * @ingroup MatricesGeometry
 *
 * The shift is given as a R2 vector (shift_X, shift_Y). An exception is thrown
 * if the displacement is not a R2 vector.
 *
 * @code
 * // Displacement of 1 pixel to the right
 * m = translation2DMatrix(vectorR2(1, 0));
 * @endcode
 */
Matrix2D< double > translation2DMatrix(const Matrix1D< double > v);

/** Creates a rotational matrix (4x4) for volumes around system axis
 * @ingroup MatricesGeometry
 *
 * The rotation angle is in degrees, and the rotational axis is either 'X', 'Y'
 * or 'Z'. An exception is thrown if the axis given is not one of these.
 *
 * The returned matrices are respectively alpha degrees around Z
 *
 * @code
 * [ cos(A) -sin(A)     0   ]
 * [ sin(A)  cos(A)     0   ]
 * [   0       0        1   ]
 * @endcode
 *
 * alpha degrees around Y
 * @code
 * [ cos(A)    0    -sin(A) ]
 * [   0       1       0    ]
 * [ sin(A)    0     cos(A) ]
 * @endcode
 *
 * alpha degrees around X
 * @code
 * [   1       0       0    ]
 * [   0     cos(A) -sin(A) ]
 * [   0     sin(A)  cos(A) ]
 * @endcode
 *
 * @code
 * m = rotation3DMatrix(60, 'X');
 * @endcode
 */
Matrix2D< double > rotation3DMatrix(double ang, char axis);

/** Creates a rotational matrix (4x4) for volumes around any axis
 * @ingroup MatricesGeometry
 *
 * The rotation angle is in degrees, and the rotational axis is given as a R3
 * vector. An exception is thrown if the axis is not a R3 vector. The axis needs
 * not to be unitary.
 *
 * @code
 * m = rotation3DMatrix(60, vectorR3(1, 1, 1));
 * @endcode
 */
Matrix2D< double > rotation3DMatrix(double ang, const Matrix1D< double >& axis);

/** Matrix which transforms the given axis into Z
 * @ingroup MatricesGeometry
 *
 * A geometrical transformation matrix (4x4) is returned such that the given
 * axis is rotated until it is aligned with the Z axis. This is very useful in
 * order to produce rotational matrices, for instance, around any axis.
 *
 * @code
 * Matrix2D< double > A = alignWithZ(axis);
 * return A.transpose() * rotation3DMatrix(ang, 'Z') * A;
 * @endcode
 *
 * The returned matrix is such that A*axis=Z, where Z and axis are column
 * vectors.
 */
Matrix2D< double > alignWithZ(const Matrix1D< double >& axis);

/** Creates a translational matrix (4x4) for volumes
 * @ingroup MatricesGeometry
 *
 * The shift is given as a R3 vector (shift_X, shift_Y, shift_Z). An exception
 * is thrown if the displacement is not a R3 vector.
 *
 * @code
 * // Displacement of 2 pixels down
 * m = translation3DMatrix(vectorR3(0, 0, 2));
 * @endcode
 */
Matrix2D< double > translation3DMatrix(const Matrix1D< double >& v);

/** Creates a scaling matrix (4x4) for volumes
 * @ingroup MatricesGeometry
 *
 * The scaling factors for the different axis must be given as a vector. So
 * that, XX(sc)=scale for X axis, YY(sc)=...
 */
Matrix2D< double > scale3DMatrix(const Matrix1D< double >& sc);

/// @defgroup MatricesMisc Miscellaneous
/// @ingroup MatricesRelated

/** Matrix equality.
 * @ingroup MatricesMisc */
template<typename T>
bool operator==(const Matrix2D<T>& op1, const Matrix2D<T>& op2)
{
    return op1.equal(op2);
}

/** Matrix inequality.
 * @ingroup MatricesMisc */
template<typename T>
bool operator!=(const Matrix2D<T>& op1, const Matrix2D<T>& op2)
{
    return !(op1==op2);
}

/** Reduce both matrices to a common size
 * @ingroup MatricesMisc
 *
 * Search the range of logical indexes for which both matrices have got valid
 * values, and cut both to that size, the corresponding origin is automatically
 * computed.
 *
 * @code
 * Matrix2D< double > m1(5, 3);
 * m1.startingX() = -2;
 * m1.startingY() = -2;
 *
 * Matrix2D< double > m2(4, 4);
 * m2.startingX() = 0;
 * m2.startingY() = 0;
 * cutToCommonSize(m1, m2);
 * // m1 and m2 range from (0,0)=(y,x) to (2,0)
 * @endcode
 */
template<typename T>
void cutToCommonSize(Matrix2D<T>& m1, Matrix2D<T>& m2)
{
    int y0 = XMIPP_MAX(STARTINGY(m1), STARTINGY(m2));
    int yF = XMIPP_MIN(FINISHINGY(m1), FINISHINGY(m2));
    int x0 = XMIPP_MAX(STARTINGX(m1), STARTINGX(m2));
    int xF = XMIPP_MIN(FINISHINGX(m1), FINISHINGX(m2));

    m1.window(y0, x0, yF, xF);
    m2.window(y0, x0, yF, xF);
}

/** Does a radial average of a matrix, around the pixel where is the origin
 * @ingroup MatricesMisc
 *
 * A vector is returned where:
 * - the first element is the mean of the pixels whose distance to the origin
 * is (0-1),
 * - the second element is the mean of the pixels whose distance to the origin
 * is (1-2)
 * - and so on.
 *
 * A second vector radial_count is returned containing the number of pixels over
 * which each radial average was calculated. if rounding=true,
 * element=round(distance).
 * - so the first element is the mean of the voxels whose distance to the origin
 * is (0.5-1.5),
 * - the second element is the mean of the voxels whose distance to the origin
 * is (1.5-2.5)
 * - and so on.
 */
template<typename T>
void radialAverage(const Matrix2D< T >& m,
                   const Matrix1D< int >& center_of_rot,
                   Matrix1D< T >& radial_mean,
                   Matrix1D< int >& radial_count,
                   const bool& rounding = false)
{
    Matrix1D< double > idx(2);

    /* First determine the maximum distance that one should expect,
       to set the dimension of the radial average vector */
    Matrix1D< int > distances(4);
    double y = STARTINGY(m) - YY(center_of_rot);
    double x = STARTINGX(m) - XX(center_of_rot);

    distances(0) = (int) FLOOR(sqrt(x * x + y * y));
    x = FINISHINGX(m) - XX(center_of_rot);
    y = STARTINGY(m) - YY(center_of_rot);

    distances(1) = (int) FLOOR(sqrt(x * x + y * y));
    x = STARTINGX(m) - XX(center_of_rot);
    y = FINISHINGY(m) - YY(center_of_rot);

    distances(2) = (int) FLOOR(sqrt(x * x + y * y));
    x = FINISHINGX(m) - XX(center_of_rot);
    y = FINISHINGY(m) - YY(center_of_rot);

    distances(3) = (int) FLOOR(sqrt(x * x + y * y));
    int dim = (int) CEIL(distances.computeMax()) + 1;

    if (rounding)
        dim++;

    // Define the vectors
    radial_mean.resize(dim);
    radial_mean.initZeros();
    radial_count.resize(dim);
    radial_count.initZeros();

    /* Perform the radial sum and count pixels that contribute to
       every distance */
    FOR_ALL_ELEMENTS_IN_MATRIX2D(m)
    {
        YY(idx) = i - YY(center_of_rot);
        XX(idx) = j - XX(center_of_rot);

        // Determine distance to the center
        int distance;
        if (rounding)
            distance = (int) ROUND(idx.module());
        else
            distance = (int) FLOOR(idx.module());

        // Sum te value to the pixels with the same distance
        radial_mean(distance) += m(i, j);

        // Count the pixel
        radial_count(distance)++;
    }

    // Perfor the mean
    FOR_ALL_ELEMENTS_IN_MATRIX1D(radial_mean)
    {
        radial_mean(i) /= (T) radial_count(i);
    }
}

/** Solve equation system, nonnegative solution
 * @ingroup MatricesMisc
 *
 * The equation system is defined by Ax=b, it is solved for x. x is forced to be
 * nonnegative. It is designed to cope with large equation systems. This
 * function is borrowed from LAPACK nnls.
 *
 * The norm of the vector Ax-b is returned.
 */
double solveNonNegative(const Matrix2D< double >& A, const Matrix1D< double >& b,
                    Matrix1D< double >& result);

/** Solve equation system, symmetric positive-definite matrix
 * @ingroup MatricesMisc
 *
 * The equation system is defined by Ax=b, it is solved for x. This method can
 * only be applied if A is positive-definite matrix and symmetric. It applies a
 * Cholesky factorization and backsubstitution (see Numerical Recipes).
 */
void solveViaCholesky(const Matrix2D< double >& A,
                        const Matrix1D< double >& b,
                        Matrix1D< double >& result);

/** Evaluate quadratic form
 * @ingroup MatricesMisc
 *
 * Given x, c and H this function returns the value of the quadratic form
 * val=c^t*x+0.5*x^t*H^t*H*x and the gradient of the quadratic form at x
 * grad=c+H*x.
 *
 * Exceptions are thrown if the vectors and matrices do not have consistent
 * dimensions.
 */
void evaluateQuadratic(const Matrix1D< double >& x, const Matrix1D< double >& c,
                    const Matrix2D< double >& H, double& val,
                    Matrix1D< double >& grad);

/** Solves Quadratic programming subproblem
 * @ingroup MatricesMisc
 *
 * @code
 * min 0.5*x'Cx + d'x   subject to:  A*x <= b
 *  x                                Aeq*x=beq
 *                                   bl<=x<=bu
 * @endcode
 */
void quadraticProgramming(const Matrix2D< double >& C, const Matrix1D< double >& d,
              const Matrix2D< double >& A, const Matrix1D< double >& b,
              const Matrix2D< double >& Aeq, const Matrix1D< double >& beq,
              Matrix1D< double >& bl, Matrix1D< double >& bu,
              Matrix1D< double >& x);


/** Solves the least square problem
 * @ingroup MatricesMisc
 *
 * @code
 * min 0.5*(Norm(C*x-d))   subject to:  A*x <= b
 * x                                    Aeq*x=beq
 *                                      bl<=x<=bu
 * @endcode
 */
void leastSquare(const Matrix2D< double >& C, const Matrix1D< double >& d,
            const Matrix2D< double >& A, const Matrix1D< double >& b,
            const Matrix2D< double >& Aeq, const Matrix1D< double >& beq,
            Matrix1D< double >& bl, Matrix1D< double >& bu,
            Matrix1D< double >& x);

/** Solves the regularized least squares problem
 * @ingroup MatricesMisc
 *
 * @code
 * min Norm(A*x-d) + lambda * norm (G*x) 
 * @endcode
 *
 * Give an empty G matrix (NULL matrix) if G is the identity matrix
 * If AtA is not a NULL matrix, then the product AtA is not computed.
 */
void regularizedLeastSquare(const Matrix2D< double >& A,
    const Matrix1D< double >& d, double lambda,
    const Matrix2D< double >& G, Matrix1D< double >& x);

/// Show matrix
template<typename T>
std::ostream& operator<<(std::ostream& ostrm, const Matrix2D<T>& v)
{
    if (XSIZE(v) == 0 || YSIZE(v) == 0)
        ostrm << "NULL matrix\n";
    else
    {
        ostrm << std::endl;
        double max_val = ABS(DIRECT_MAT_ELEM(v, 0, 0));

    	T* ptr;
	unsigned long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(v,n,ptr)
            max_val = XMIPP_MAX(max_val, ABS(*ptr));

        int prec = bestPrecision(max_val, 10);

        for (int i = STARTINGY(v); i <= FINISHINGY(v); i++)
        {
            for (int j = STARTINGX(v); j <= FINISHINGX(v); j++)
            {
                if (typeid(T)==typeid(std::complex<double>))
                    ostrm << MAT_ELEM(v, i, j) << ' ';
                else
                    ostrm << floatToString((double) MAT_ELEM(v, i, j), 10, prec) << ' ';
            }
            ostrm << std::endl;
        }
    }

    return ostrm;
}

// TODO Document
template<typename T>
void solve(const Matrix2D<T>& A, const Matrix1D<T>& b, Matrix1D<T>& result)
{
    if (XSIZE(A) == 0)
        REPORT_ERROR(1108, "Solve: Matrix is empty");

    if (XSIZE(A) != YSIZE(A))
        REPORT_ERROR(1109, "Solve: Matrix is not squared");

    if (XSIZE(A) != XSIZE(b))
        REPORT_ERROR(1102, "Solve: Different sizes of Matrix and Vector");

    if (b.isRow())
        REPORT_ERROR(1107, "Solve: Not correct vector shape");

    // Perform LU decomposition and then solve
    Matrix1D< int > indx;
    T d;
    Matrix2D<T> LU;
    ludcmp(A, LU, indx, d);
    result = b;
    lubksb(LU, indx, result);
}

// TODO Document
template<typename T>
void solveBySVD(const Matrix2D< T >& A, const Matrix1D< T >& b,
                  Matrix1D< double >& result, double tolerance)
{
    if (XSIZE(A) == 0)
        REPORT_ERROR(1108, "Solve: Matrix is empty");

    if (XSIZE(A) != YSIZE(A))
        REPORT_ERROR(1109, "Solve: Matrix is not squared");

    if (XSIZE(A) != XSIZE(b))
        REPORT_ERROR(1102, "Solve: Different sizes of Matrix and Vector");

    if (b.isRow())
        REPORT_ERROR(1107, "Solve: Not correct vector shape");

    // First perform de single value decomposition
    // Xmipp interface that calls to svdcmp of numerical recipes
    Matrix2D< double > u, v;
    Matrix1D< double > w;
    svdcmp(A, u, w, v);

    // Here is checked if eigenvalues of the svd decomposition are acceptable
    // If a value is lower than tolerance, the it's zeroed, as this increases
    // the precision of the routine.
    FOR_ALL_ELEMENTS_IN_MATRIX1D(w)
    if (w(i) < tolerance)
        w(i) = 0;

    // Set size of matrices
    result.resize(XSIZE(b));

    // Xmipp interface that calls to svdksb of numerical recipes
    Matrix1D< double > bd;
    typeCast(b, bd);
    svbksb(u, w, v, bd, result);
}

// TODO Document
template<typename T>
void solve(const Matrix2D<T>& A, const Matrix2D<T>& b, Matrix2D<T>& result)
{
    if (XSIZE(A) == 0)
        REPORT_ERROR(1108, "Solve: Matrix is empty");

    if (XSIZE(A) != YSIZE(A))
        REPORT_ERROR(1109, "Solve: Matrix is not squared");

    if (YSIZE(A) != YSIZE(b))
        REPORT_ERROR(1102, "Solve: Different sizes of A and b");

    // Solve
    result = b;
    Matrix2D<T> Aux = A;
    gaussj(Aux.adaptForNumericalRecipes2(), Aux.ydim,
           result.adaptForNumericalRecipes2(), b.xdim);
}

// TODO Document
template<typename T>
void ludcmp(const Matrix2D<T>& A, Matrix2D<T>& LU, Matrix1D< int >& indx, T& d)
{
    LU = A;
    indx.resize(XSIZE(A));
    ludcmp(LU.adaptForNumericalRecipes2(), XSIZE(A),
           indx.adaptForNumericalRecipes(), &d);
}

// TODO Document
template<typename T>
void lubksb(const Matrix2D<T>& LU, Matrix1D< int >& indx, Matrix1D<T>& b)
{
    lubksb(LU.adaptForNumericalRecipes2(), XSIZE(indx),
           indx.adaptForNumericalRecipes(),
           b.adaptForNumericalRecipes());
}

// TODO Document
#define VIA_BILIB
template<typename T>
void svdcmp(const Matrix2D< T >& a,
            Matrix2D< double >& u,
            Matrix1D< double >& w,
            Matrix2D< double >& v)
{
    // svdcmp only works with double
    typeCast(a, u);

    // Set size of matrices
    w.initZeros(u.colNumber());
    v.initZeros(u.colNumber(), u.colNumber());

    // Call to the numerical recipes routine
#ifdef VIA_NR
    svdcmp(MULTIDIM_ARRAY(u),
           u.rowNumber(), u.colNumber(),
           MULTIDIM_ARRAY(w),
           MULTIDIM_ARRAY(v));
#endif

#ifdef VIA_BILIB
    int status;
    SingularValueDecomposition(MULTIDIM_ARRAY(u),
                               u.rowNumber(), u.colNumber(),
                               MULTIDIM_ARRAY(w),
                               MULTIDIM_ARRAY(v),
                               5000, &status);
#endif
}

#undef VIA_NR
#undef VIA_BILIB

// TODO Document
//#define DEBUG_APPLYGEO
template<typename T>
void applyGeometry(Matrix2D<T>& M2, Matrix2D< double > A, const Matrix2D<T>& M1, bool inv,
                bool wrap, T outside)
{
    int m1, n1, m2, n2;
    double x, y, xp, yp;
    double minxp, minyp, maxxp, maxyp;
    int cen_x, cen_y, cen_xp, cen_yp;
    double wx, wy; // Weights in X,Y directions for bilinear interpolation
    int Xdim, Ydim;

    if ((XSIZE(A) != 3) || (YSIZE(A) != 3))
        REPORT_ERROR(1102, "Apply_geom: geometrical transformation is not 3x3");

    if (A.isIdentity())
    {
        M2 = M1;
        return;
    }
    if (XSIZE(M1) == 0)
    {
        M2.clear();
        return;
    }

    if (!inv)
        A = A.inv();

    // For scalings the output matrix is resized outside to the final
    // size instead of being resized inside the routine with the
    // same size as the input matrix
    if (XSIZE(M2) == 0)
        M2.resize(M1);

    if (outside != 0.)
    {
        // Initialise output matrix with value=outside
        FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(M2)
        {
            DIRECT_MAT_ELEM(M2, i, j) = outside;
        }
    }

    // Find center and limits of image
    cen_y  = (int)(YSIZE(M2) / 2);
    cen_x  = (int)(XSIZE(M2) / 2);
    cen_yp = (int)(YSIZE(M1) / 2);
    cen_xp = (int)(XSIZE(M1) / 2);
    minxp  = -cen_xp;
    minyp  = -cen_yp;
    maxxp  = XSIZE(M1) - cen_xp - 1;
    maxyp  = YSIZE(M1) - cen_yp - 1;
    Xdim   = XSIZE(M1);
    Ydim   = YSIZE(M1);

    // Now we go from the output image to the input image, ie, for any pixel
    // in the output image we calculate which are the corresponding ones in
    // the original image, make an interpolation with them and put this value
    // at the output pixel
    //#define DEBUG_APPLYGEO
#ifdef DEBUG_APPLYGEO
    std::cout << "A\n" << A << std::endl
              << "(cen_x ,cen_y )=(" << cen_x  << "," << cen_y  << ")\n"
              << "(cen_xp,cen_yp)=(" << cen_xp << "," << cen_yp << ")\n"
              << "(min_xp,min_yp)=(" << minxp  << "," << minyp  << ")\n"
              << "(max_xp,max_yp)=(" << maxxp  << "," << maxyp  << ")\n";
#endif
    //#undef DEBUG_APPLYGEO

    for (int i = 0; i < YSIZE(M2); i++)
    {
        // Calculate position of the beginning of the row in the output image
        x = -cen_x;
        y = i - cen_y;

        // For OldXmipp origins with even XSIZE & YSIZE:
        //      x= -cen_x+0.5;
        //      y=i-cen_y+0.5;

        // Calculate this position in the input image according to the
        // geometrical transformation
        // they are related by
        // coords_output(=x,y) = A * coords_input (=xp,yp)
        xp = x * dMij(A, 0, 0) + y * dMij(A, 0, 1) + dMij(A, 0, 2);
        yp = x * dMij(A, 1, 0) + y * dMij(A, 1, 1) + dMij(A, 1, 2);

        for (int j = 0; j < XSIZE(M2); j++)
        {
            bool interp;
            T tmp;

#ifdef DEBUG_APPLYGEO
            std::cout << "Computing (" << i << "," << j << ")\n";
            std::cout << "   (y, x) =(" << y << "," << x << ")\n"
                      << "   before wrapping (y',x')=(" << yp << "," << xp << ") "
                      << std::endl;
#endif
            // If the point is outside the image, apply a periodic extension
            // of the image, what exits by one side enters by the other
            interp = true;
            if (wrap)
            {
                if (xp < minxp - XMIPP_EQUAL_ACCURACY ||
                    xp > maxxp + XMIPP_EQUAL_ACCURACY)
                    xp = realWRAP(xp, minxp - 0.5, maxxp + 0.5);

                if (yp < minyp - XMIPP_EQUAL_ACCURACY ||
                    yp > maxyp + XMIPP_EQUAL_ACCURACY)

                    yp = realWRAP(yp, minyp - 0.5, maxyp + 0.5);
            }
            else
            {
                if (xp < minxp - XMIPP_EQUAL_ACCURACY ||
                    xp > maxxp + XMIPP_EQUAL_ACCURACY)
                    interp = false;

                if (yp < minyp - XMIPP_EQUAL_ACCURACY ||
                    yp > maxyp + XMIPP_EQUAL_ACCURACY)
                    interp = false;
            }

#ifdef DEBUG_APPLYGEO
            std::cout << "   after wrapping (y',x')=(" << yp << "," << xp << ") "
                      << std::endl;
            std::cout << "   Interp = " << interp << std::endl;
            x++;
#endif

            if (interp)
            {
                // Calculate the integer position in input image, be careful
                // that it is not the nearest but the one at the top left corner
                // of the interpolation square. Ie, (0.7,0.7) would give (0,0)
                // Calculate also weights for point m1+1,n1+1
                wx = xp + cen_xp;
                m1 = (int) wx;
                wx = wx - m1;
                m2 = m1 + 1;
                wy = yp + cen_yp;
                n1 = (int) wy;
                wy = wy - n1;
                n2 = n1 + 1;

                // m2 and n2 can be out by 1 so wrap must be check here
                if (wrap)
                {
                    if (m2 >= Xdim)
                        m2 = 0;
                    if (n2 >= Ydim)
                        n2 = 0;
                }

#ifdef DEBUG_APPLYGEO
                std::cout << "   From (" << n1 << "," << m1 << ") and ("
                          << n2 << "," << m2 << ")\n";
                std::cout << "   wx= " << wx << " wy= " << wy << std::endl;
#endif

                // Perform interpolation
                // if wx == 0 means that the rightest point is useless for this
                // interpolation, and even it might not be defined if m1=xdim-1
                // The same can be said for wy.
                tmp  = (T)((1 - wy) * (1 - wx) * dMij(M1, n1, m1));

                if (wx != 0 && m2 < M1.xdim)
                    tmp += (T)((1 - wy) * wx * dMij(M1, n1, m2));

                if (wy != 0 && n2 < M1.ydim)
                {
                    tmp += (T)(wy * (1 - wx) * dMij(M1, n2, m1));

                    if (wx != 0 && m2 < M1.xdim)
                        tmp += (T)(wy * wx * dMij(M1, n2, m2));
                }

                dMij(M2, i, j) = tmp;

#ifdef DEBUG_APPYGEO
                std::cout << "   val= " << tmp << std::endl;
#endif

            }

            // Compute new point inside input image
            xp += dMij(A, 0, 0);
            yp += dMij(A, 1, 0);
        }
    }
}

#undef DEBUG_APPLYGEO

// TODO Document
//#define DEBUG
template<typename T>
void applyGeometryBSpline(Matrix2D<T>& M2, Matrix2D< double > A, const Matrix2D<T>& M1,
                        int Splinedegree, bool inv, bool wrap, T outside)
{
    int m1, n1, m2, n2;
    double x, y, xp, yp;
    double minxp, minyp, maxxp, maxyp;
    int cen_x, cen_y, cen_xp, cen_yp;

    if ((XSIZE(A) != 3) || (YSIZE(A) != 3))
        REPORT_ERROR(1102, "Apply_geom: geometrical transformation is not 3x3");

    if (A.isIdentity())
    {
        M2 = M1;
        return;
    }

    if (XSIZE(M1) == 0)
    {
        M2.clear();
        return;
    }

    if (!inv)
        A = A.inv();

    // For scalings the output matrix is resized outside to the final
    // size instead of being resized inside the routine with the
    // same size as the input matrix
    if (XSIZE(M2) == 0)
        M2.resize(M1);

    // Initialise output matrix with value=outside
    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(M2)
    {
        DIRECT_MAT_ELEM(M2, i, j) = outside;
    }

    // Find center and limits of image
    cen_y  = (int)(YSIZE(M2) / 2);
    cen_x  = (int)(XSIZE(M2) / 2);
    cen_yp = (int)(YSIZE(M1) / 2);
    cen_xp = (int)(XSIZE(M1) / 2);
    minxp  = -cen_xp;
    minyp  = -cen_yp;
    maxxp  = XSIZE(M1) - cen_xp - 1;
    maxyp  = YSIZE(M1) - cen_yp - 1;

    // Build the B-spline coefficients
    Matrix2D< double > Bcoeffs;
    M1.produceSplineCoefficients(Bcoeffs, Splinedegree);
    STARTINGX(Bcoeffs) = (int) minxp;
    STARTINGY(Bcoeffs) = (int) minyp;

    // Now we go from the output image to the input image, ie, for any pixel
    // in the output image we calculate which are the corresponding ones in
    // the original image, make an interpolation with them and put this value
    // at the output pixel
    for (int i = 0; i < YSIZE(M2); i++)
    {
        // Calculate position of the beginning of the row in the output image
        x = -cen_x;
        y = i - cen_y;

        // Calculate this position in the input image according to the
        // geometrical transformation
        // they are related by
        // coords_output(=x,y) = A * coords_input (=xp,yp)
        xp = x * dMij(A, 0, 0) + y * dMij(A, 0, 1) + dMij(A, 0, 2);
        yp = x * dMij(A, 1, 0) + y * dMij(A, 1, 1) + dMij(A, 1, 2);

        for (int j = 0; j < XSIZE(M2); j++)
        {
            bool interp;
            T tmp;

#ifdef DEBUG_APPLYGEO
            std::cout << "Computing (" << i << "," << j << ")\n";
            std::cout << "   (y, x) =(" << y << "," << x << ")\n"
                      << "   before wrapping (y',x')=(" << yp << "," << xp << ") "
                      << std::endl;
#endif

            // If the point is outside the image, apply a periodic extension
            // of the image, what exits by one side enters by the other
            interp = true;
            if (wrap)
            {
                if (xp < minxp - XMIPP_EQUAL_ACCURACY ||
                    xp > maxxp + XMIPP_EQUAL_ACCURACY)
                    xp = realWRAP(xp, minxp - 0.5, maxxp + 0.5);

                if (yp < minyp - XMIPP_EQUAL_ACCURACY ||
                    yp > maxyp + XMIPP_EQUAL_ACCURACY)
                    yp = realWRAP(yp, minyp - 0.5, maxyp + 0.5);
            }
            else
            {
                if (xp < minxp - XMIPP_EQUAL_ACCURACY ||
                    xp > maxxp + XMIPP_EQUAL_ACCURACY)
                    interp = false;

                if (yp < minyp - XMIPP_EQUAL_ACCURACY ||
                    yp > maxyp + XMIPP_EQUAL_ACCURACY)
                    interp = false;
            }

#ifdef DEBUG_APPLYGEO
            std::cout << "   after wrapping (y',x')=(" << yp << "," << xp << ") "
                      << std::endl;
            std::cout << "   Interp = " << interp << std::endl;
            x++;
#endif

            if (interp)
            {
                dMij(M2, i, j) = (T) Bcoeffs.interpolatedElementBSpline(
                                     xp, yp, Splinedegree);

#ifdef DEBUG_APPLYGEO
                std::cout << "   val= " << dMij(M2, i, j) << std::endl;
#endif

            }

            // Compute new point inside input image
            xp += dMij(A, 0, 0);
            yp += dMij(A, 1, 0);
        }
    }
}

#undef DEBUG_APPLYGEO

// TODO Document
template<>
inline void Matrix2D< std::complex< double > >::scaleToSizeBSpline(
    int Splinedegree, int Ydim, int Xdim,
    Matrix2D< std::complex< double > > & result) const
{
    Matrix2D< double > re, im, scre, scim;

    re.resize(YSIZE(*this), XSIZE(*this));
    im.resize(YSIZE(*this), XSIZE(*this));

    Complex2RealImag(MULTIDIM_ARRAY(*this),
                     MULTIDIM_ARRAY(re), MULTIDIM_ARRAY(im),
                     MULTIDIM_SIZE(*this));

    re.scaleToSizeBSpline(Splinedegree, Ydim, Xdim, scre);
    im.scaleToSizeBSpline(Splinedegree, Ydim, Xdim, scim);

    result.resize(Ydim, Xdim);

    RealImag2Complex(MULTIDIM_ARRAY(re), MULTIDIM_ARRAY(im),
                     MULTIDIM_ARRAY(result), MULTIDIM_SIZE(re));
}

// TODO Document
template<>
bool Matrix2D< std::complex< double > >::isDiagonal() const;

// TODO Document
template<>
bool Matrix2D< std::complex<double> >::isScalar() const;

// TODO Document
template<typename T>
void multiplyMatrix(const Matrix2D<T>& op1, const Matrix2D<T>& op2, Matrix2D<T>& result)
{
    if (XSIZE(op1) != YSIZE(op2))
        REPORT_ERROR(1102, "Not compatible sizes in matrix multiplication");

    result.initZeros(YSIZE(op1), XSIZE(op2));
    for (int i = 0; i < YSIZE(op1); i++)
        for (int j = 0; j < XSIZE(op2); j++)
            for (int k = 0; k < XSIZE(op1); k++)
                DIRECT_MAT_ELEM(result, i, j) += DIRECT_MAT_ELEM(op1, i, k) *
                                                 DIRECT_MAT_ELEM(op2, k, j);

    STARTINGY(result) = STARTINGY(op1);
    STARTINGX(result) = STARTINGX(op1);
}

// TODO Document
template<typename T>
Matrix1D<T> Matrix1D<T>::operator*(const Matrix2D<T>& M)
{
    Matrix1D<T> result;

    if (XSIZE(*this) != YSIZE(M))
        REPORT_ERROR(1102, "Not compatible sizes in matrix by vector");

    if (!isRow())
        REPORT_ERROR(1102, "Vector is not a row");

    result.initZeros(XSIZE(M));
    for (int j = 0; j < XSIZE(M); j++)
        for (int i = 0; i < YSIZE(M); i++)
            DIRECT_VEC_ELEM(result, j) += DIRECT_VEC_ELEM(*this, i) *
                                          DIRECT_MAT_ELEM(M, i, j);

    result.setRow();
    STARTINGX(result) = STARTINGX(M);

    return result;
}

#endif
