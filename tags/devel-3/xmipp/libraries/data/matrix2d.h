/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#ifndef MATRIX2D_H
#define MATRIX2D_H

#include <iostream>
#include <string>
#include <complex>

#include <external/bilib/headers/pyramidtools.h>
#include "multidimensional_array.h"
#include "matrix1d.h"
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

Matrix2D< double > rotation2DMatrix(double ang);

Matrix2D< double > translation2DMatrix(const Matrix1D< double > &v);

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

/** Inverse of a matrix (3x3)
 * @ingroup MatricesArithmetic
 *
 * Input and output matrix cannot be the same one. The output is supposed to be
 * already resized.
 */
#define M3x3_INV(Ainv, A) { \
        dMij(Ainv, 0, 0) =   dMij(A, 2, 2)*dMij(A, 1, 1)-dMij(A, 2, 1)*dMij(A, 1, 2); \
        dMij(Ainv, 0, 1) = -(dMij(A, 2, 2)*dMij(A, 0, 1)-dMij(A, 2, 1)*dMij(A, 0, 2)); \
        dMij(Ainv, 0, 2) =   dMij(A, 1, 2)*dMij(A, 0, 1)-dMij(A, 1, 1)*dMij(A, 0, 2); \
        dMij(Ainv, 1, 0) = -(dMij(A, 2, 2)*dMij(A, 1, 0)-dMij(A, 2, 0)*dMij(A, 1, 2)); \
        dMij(Ainv, 1, 1) =   dMij(A, 2, 2)*dMij(A, 0, 0)-dMij(A, 2, 0)*dMij(A, 0, 2); \
        dMij(Ainv, 1, 2) = -(dMij(A, 1, 2)*dMij(A, 0, 0)-dMij(A, 1, 0)*dMij(A, 0, 2)); \
        dMij(Ainv, 2, 0) =   dMij(A, 2, 1)*dMij(A, 1, 0)-dMij(A, 2, 0)*dMij(A, 1, 1); \
        dMij(Ainv, 2, 1) = -(dMij(A, 2, 1)*dMij(A, 0, 0)-dMij(A, 2, 0)*dMij(A, 0, 1)); \
        dMij(Ainv, 2, 2) =   dMij(A, 1, 1)*dMij(A, 0, 0)-dMij(A, 1, 0)*dMij(A, 0, 1); \
        spduptmp0 = 1.0 / (dMij(A, 0, 0)*dMij(Ainv, 0, 0)+dMij(A, 1, 0)*dMij(Ainv, 0, 1)+\
            dMij(A, 2, 0)*dMij(Ainv, 0, 2)); \
        M3x3_BY_CT(Ainv, Ainv, spduptmp0); }
//@}

/// Template class for Xmipp matrices
/// @ingroup Matrices
template<typename T>
class Matrix2D: public MultidimArray<T>
{
public:

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
        for (int i=0; i<YSIZE(*this); i++)
            for (int j=0; j<XSIZE(*this); j++)
                DIRECT_MAT_ELEM(*this, i, j) = (T)(i == j);

        /* I dont know why this does not work...
        FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(*this)
        {
            DIRECT_MAT_ELEM(*this, i, j) = (T)(i == j);
        }
        */
    }

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
        MultidimArray<T>::resize(1, 1, Ydim, Xdim);
    }

    /** Resize taking the shape from another matrix
     * @ingroup MatricesSize
     */
    template <typename T1>
    void resize(const Matrix2D<T1> &M)
    {
        MultidimArray<T>::resize(M);
    }

    /** Resize taking the shape from another matrix which
        is given as a MultidimArray
     * @ingroup MatricesSize
     */
    template <typename T1>
    void resize(const MultidimArray<T1> &M)
    {
        MultidimArray<T>::resize(M);
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

    /** @defgroup MatricesMemory Memory access
     * @ingroup Matrices
     *
     * This functions allows you to access to the matrix elements.
     */

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
            profile(i) = interpolatedElement2D(tx, ty);
            tx += tx_step;
            ty += ty_step;
        }
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


    /** @defgroup MatricesAlgebraic Algebraic operations
     * @ingroup Matrices
     *
     * NOTICE!!!: the matrix by matrix multiplier operator (*) has been
     * redefined in this class and it represents the true matrix by matrix
     * algebraic multiplication, if you want an element-wise multiplication you
     * have to use the function defined in this section "multiplyElements"
     */

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


};


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

/** Creates a rotational matrix (3x3) for images
 * @ingroup MatricesGeometry
 *
 * The rotation angle is in degrees.
 * m must have been already resized to 3x3
 *
 * @code
 *  rotation2DMatrix(60,m);
 * @endcode
 */
void rotation2DMatrix(double ang, Matrix2D< double > &m);

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
Matrix2D< double > translation2DMatrix(const Matrix1D< double > &v);

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
                ostrm << floatToString((double) MAT_ELEM(v, i, j), 10, prec) << ' ';
            }
            ostrm << std::endl;
        }
    }

    return ostrm;
}

// Specialization for complex matrices
std::ostream& operator<<(std::ostream& ostrm,
    const Matrix2D< std::complex<double> >& v);

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
    gaussj(Aux.adaptForNumericalRecipes22D(), Aux.ydim,
           result.adaptForNumericalRecipes22D(), b.xdim);
}

// TODO Document
template<typename T>
void ludcmp(const Matrix2D<T>& A, Matrix2D<T>& LU, Matrix1D< int >& indx, T& d)
{
    LU = A;
    indx.resize(XSIZE(A));
    ludcmp(LU.adaptForNumericalRecipes22D(), XSIZE(A),
           indx.adaptForNumericalRecipes2D(), &d);
}

// TODO Document
template<typename T>
void lubksb(const Matrix2D<T>& LU, Matrix1D< int >& indx, Matrix1D<T>& b)
{
    lubksb(LU.adaptForNumericalRecipes22D(), XSIZE(indx),
           indx.adaptForNumericalRecipes2D(),
           b.adaptForNumericalRecipes2D());
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


#undef DEBUG_APPLYGEO


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
