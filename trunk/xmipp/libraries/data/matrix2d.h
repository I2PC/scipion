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

#ifndef MATRIX2D_H_
#define MATRIX2D_H_

#include <external/bilib/types/tsplinebasis.h>
#include <external/bilib/types/tboundaryconvention.h>
#include <external/bilib/headers/linearalgebra.h>
#include <external/bilib/headers/changebasis.h>
#include <external/bilib/headers/kernel.h>
#include <external/bilib/headers/pyramidtools.h>
#include "matrix1d.h"

/** @defgroup Matrices Matrix2D Matrices
 * @ingroup DataLibrary
 */

/** @defgroup MatricesSpeedUp Matrices speed up macros
 * @ingroup Matrices
 */

/** Array access.
 * @ingroup MatricesSpeedUp
 *
 * This macro gives you access to the array (T)
 */
#define MATRIX2D_ARRAY(m) ((m).mdata)

/** For all elements in the array
 * @ingroup MatricesSpeedUp
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
    for (int i=0; i<(m).mdimy; i++) \
        for (int j=0; j<(m).mdimx; j++)

/** Access to a matrix element
 * @ingroup MatricesSpeedUp
 * v is the array, i and j define the element v_ij.
 *
  * @code
 * MAT_ELEM(m, 0, 0) = 1;
 * val = MAT_ELEM(m, 0, 0);
 * @endcode
 */
#define MAT_ELEM(m,i,j) ((m).mdata[(i)*(m).mdimx+(j)])

/** X dimension of the matrix
 * @ingroup MatricesSpeedUp
 */
#define MAT_XSIZE(m) ((m).mdimx)

/** Y dimension of the matrix
 * @ingroup MatricesSpeedUp
 */
#define MAT_YSIZE(m) ((m).mdimy)

/** Matrix element: Element access
 * @ingroup MatricesSpeedUp
 *
 * This is just a redefinition
 * of the function above
 * @code
 * dMij(m, -2, 1) = 1;
 * val = dMij(m, -2, 1);
 * @endcode
 */
#define dMij(m, i, j)  MAT_ELEM(m, i, j)

/** Matrix (3x3) by vector (3x1) (a=M*b)
 * @ingroup MatricesSpeedUp
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
 * @ingroup MatricesSpeedUp
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
 * @ingroup MatricesSpeedUp
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
 * @ingroup MatricesSpeedUp
 *
 * You must create the result matrix with the appropriate size. You can reuse
 * the matrix M1 to store the results (that is, M2x2_BY_CT(M, M, k);, is
 * allowed).
 */
#define M2x2_BY_CT(M2, M1, k) { \
        dMij(M2, 0, 0) = dMij(M1, 0, 0) * k; \
        dMij(M2, 0, 1) = dMij(M1, 0, 1) * k; \
        dMij(M2, 1, 0) = dMij(M1, 1, 0) * k; \
        dMij(M2, 1, 1) = dMij(M1, 1, 1) * k; }

/** Matrix (3x3) by constant (M2=M1*k)
 * @ingroup MatricesSpeedUp
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
 * @ingroup MatricesSpeedUp
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
 * @ingroup MatricesSpeedUp
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

// Forward declarations
template<typename T>
class Matrix1D;
template<typename T>
class Matrix2D;

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

/** Matrix2D class */
template<typename T>
class Matrix2D
{
public:
    // The array itself
    T* mdata;

    // Destroy data
    bool destroyData;

    // Number of elements in X
    int mdimx;

    // Number of elements in Y
    int mdimy;

    // Total number of elements
    int mdim;

    /// @defgroup MatrixConstructors Constructors
    /// @ingroup Matrices
    /** Empty constructor
     * @ingroup MatrixConstructors
     */
    Matrix2D()
    {
        coreInit();
    }

    /** Dimension constructor
     * @ingroup MatrixConstructors
     */
    Matrix2D(int Ydim, int Xdim)
    {
        coreInit();
        resize(Ydim, Xdim);
    }

    /** Copy constructor
     * @ingroup MatrixConstructors
     */
    Matrix2D(const Matrix2D<T>& v)
    {
        coreInit();
        *this = v;
    }

    /** Destructor.
     * @ingroup MatrixConstructors
     */
    ~Matrix2D()
    {
        coreDeallocate();
    }

    /** Assignment.
     * @ingroup MatrixConstructors
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
            for (int i = 0; i< op1.mdim; i++)
                mdata[i] = op1.mdata[i];
        }

        return *this;
    }

    /// @defgroup MatrixCore Core memory operations for Matrix2D
    /** Clear.
     * @ingroup MatrixCore
     */
    void clear()
    {
        coreDeallocate();
        coreInit();
    }

    /** Core init.
     * @ingroup MatrixCore
     * Initialize everything to 0
     */
    void coreInit()
    {
        mdimx=mdimy=mdim=0;
        mdata=NULL;
        destroyData=true;
    }

    /** Core allocate.
     * @ingroup MatrixCore
     */
    void coreAllocate(int _mdimy, int _mdimx)
    {
        if (_mdimy <= 0 ||_mdimx<=0)
        {
            clear();
            return;
        }

        mdimx=_mdimx;
        mdimy=_mdimy;
        mdim=_mdimx*_mdimy;
        mdata = new T [mdim];
        if (mdata == NULL)
            REPORT_ERROR(1001, "Allocate: No space left");
    }

    /** Core deallocate.
     * @ingroup MatrixCore
     * Free all mdata.
     */
    void coreDeallocate()
    {
        if (mdata != NULL && destroyData)
            delete[] mdata;
        mdata=NULL;
    }

    /// @defgroup MatrixSize Size and shape of Matrix2D
    /** Resize to a given size
     * @ingroup MatrixSize
     */
    void resize(int Ydim, int Xdim)
    {

        if (Xdim == mdimx && Ydim == mdimy)
            return;

        if (Xdim <= 0 || Ydim <= 0)
        {
            clear();
            return;
        }

        T * new_mdata;
        size_t YXdim=Ydim*Xdim;

        try
        {
            new_mdata = new T [YXdim];
        }
        catch (std::bad_alloc &)
        {
            REPORT_ERROR(1001, "Allocate: No space left");
        }

        // Copy needed elements, fill with 0 if necessary
        for (int i = 0; i < Ydim; i++)
            for (int j = 0; j < Xdim; j++)
            {
                T val;
                if (i >= mdimy)
                    val = 0;
                else if (j >= mdimx)
                    val = 0;
                else
                    val = mdata[i*mdimx + j];
                new_mdata[i*Xdim+j] = val;
            }

        // deallocate old vector
        coreDeallocate();

        // assign *this vector to the newly created
        mdata = new_mdata;
        mdimx = Xdim;
        mdimy = Ydim;
        mdim = Xdim * Ydim;
    }

    /** Resize according to a pattern.
     * @ingroup MatrixSize
     *
     * This function resize the actual array to the same size and origin
     * as the input pattern. If the actual array is larger than the pattern
     * then the trailing values are lost, if it is smaller then 0's are
     * added at the end
     *
     * @code
     * v2.resize(v1);
     * // v2 has got now the same structure as v1
     * @endcode
     */
    template<typename T1>
    void resize(const Matrix2D<T1> &v)
    {
        if (mdimx != v.mdimx || mdimy != v.mdimy)
            resize(v.mdimy, v.mdimx);
    }

    /** Extract submatrix and assign to this object.
     * @ingroup MatrixSize
     */
    void submatrix(int i0, int j0, int iF, int jF)
    {
        if (i0<0 || j0<0 || iF>=MAT_YSIZE(*this) || jF>=MAT_XSIZE(*this))
            REPORT_ERROR(1,"Submatrix indexes out of bounds");
        Matrix2D<T> result(iF - i0 + 1, jF - j0 + 1);

        FOR_ALL_ELEMENTS_IN_MATRIX2D(result)
        MAT_ELEM(result, i, j) = MAT_ELEM(*this, i+i0, j+j0);

        *this = result;
    }

    /** Same shape.
     * @ingroup MatrixSize
     *
     * Returns true if this object has got the same shape (origin and size)
     * than the argument
     */
    template <typename T1>
    bool sameShape(const Matrix2D<T1>& op) const
    {
        return ((mdimx == op.mdimx) && (mdimy == op.mdimy));
    }

    /** X dimension
     * @ingroup MatrixSize
     *
     * Returns X dimension
     */
    inline int Xdim() const
    {
        return mdimx;
    }

    /** Y dimension
     * @ingroup MatrixSize
     *
     * Returns Y dimension
     */
    inline int Ydim() const
    {
        return mdimy;
    }

    /// @defgroup MatrixInitialization Initialization of Matrix2D values
    /** Same value in all components.
     * @ingroup MatrixInitialization
     *
     * The constant must be of a type compatible with the array type, ie,
     * you cannot  assign a double to an integer array without a casting.
     * It is not an error if the array is empty, then nothing is done.
     *
     * @code
     * v.initConstant(3.14);
     * @endcode
     */
    void initConstant(T val)
    {
        for (int j = 0; j < mdim; j++)
        {
            mdata[j] = val;
        }
    }

    /** Initialize to zeros with current size.
     * @ingroup MatrixInitialization
     *
     * All values are set to 0. The current size and origin are kept. It is not
     * an error if the array is empty, then nothing is done.
     *
     * @code
     * v.initZeros();
     * @endcode
     */
    void initZeros()
    {
        initConstant(static_cast< T >(0));
    }

    /** Initialize to zeros with a given size.
     * @ingroup MatrixInitialization
     */
    void initZeros(int Ydim, int Xdim)
    {
        resize(Ydim, Xdim);
        initZeros();
    }

    /** Initialize to zeros following a pattern.
      * @ingroup MatrixInitialization
      *
      * All values are set to 0, and the origin and size of the pattern are
      * adopted.
      *
      * @code
      * v2.initZeros(v1);
      * @endcode
      */
    template <typename T1>
    void initZeros(const Matrix2D<T1>& op)
    {
        resize(op);
        initConstant(static_cast< T >(0));
    }

    /** 2D Identity matrix of current size
     * @ingroup MatrixInitialization
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
        initIdentity(MAT_XSIZE(*this));
    }

    /** 2D Identity matrix of a given size
     * @ingroup MatrixInitialization
     *
     * A (dim x dim) identity matrix is generated.
     *
     * @code
     * m.initIdentity(3);
     * @endcode
     */
    void initIdentity(int dim)
    {
        if (dim == 0)
        {
            clear();
            return;
        }

        resize(dim, dim);
        for (int i = 0; i < mdimy; i++)
            for (int j = 0; j < mdimx; j++)
                (*this)(i,j) = (T)(i == j);
    }

    /// @defgroup MatrixOperators Operators for Matrix2D

    /** Matrix element access
     * @ingroup MatrixOperators
     */
    T& operator()(int i, int j) const
    {
        return MAT_ELEM((*this),i,j);
    }
    /** Parenthesis operator for phyton
     * @ingroup MatrixOperators
    */
    void setVal(T val,int y, int x)
    {
        MAT_ELEM((*this),y,x)=val;
    }
    /** Parenthesis operator for phyton
     * @ingroup MatrixOperators
    */
    T getVal( int y, int x)
    {
        return MAT_ELEM((*this),y,x);
    }

    /** v3 = v1 * k.
     * @ingroup MatrixOperators
     */
    Matrix2D<T> operator*(T op1) const
    {
        Matrix2D<T> tmp(*this);
        for (int i=0; i < mdim; i++)
            tmp.mdata[i] = mdata[i] * op1;
        return tmp;
    }

    /** v3 = v1 / k.
     * @ingroup MatrixOperators
     */
    Matrix2D<T> operator/(T op1) const
    {
        Matrix2D<T> tmp(*this);
        for (int i=0; i < mdim; i++)
            tmp.mdata[i] = mdata[i] / op1;
        return tmp;
    }

    /** v3 = k * v2.
     * @ingroup MatrixOperators
     */
    friend Matrix2D<T> operator*(T op1, const Matrix2D<T>& op2)
    {
        Matrix2D<T> tmp(op2);
        for (int i=0; i < op2.mdim; i++)
            tmp.mdata[i] = op1 * op2.mdata[i];
        return tmp;
    }

    /** v3 *= k.
      * @ingroup MatrixOperators
      */
    void operator*=(T op1)
    {
        for (int i=0; i < mdim; i++)
            mdata[i] *= op1;
    }

    /** v3 /= k.
      * @ingroup MatrixOperators
      */
    void operator/=(T op1)
    {
        for (int i=0; i < mdim; i++)
            mdata[i] /= op1;
    }

    /** Matrix by vector multiplication
     * @ingroup MatrixOperators
     *
     * @code
     * v2 = A*v1;
     * @endcode
     */
    Matrix1D<T> operator*(const Matrix1D<T>& op1) const
    {
        Matrix1D<T> result;

        if (mdimx != op1.size())
            REPORT_ERROR(1102, "Not compatible sizes in matrix by vector");

        if (!op1.isCol())
            REPORT_ERROR(1102, "Vector is not a column");

        result.initZeros(mdimy);

        for (int i = 0; i < mdimy; i++)
            for (int j = 0; j < op1.size(); j++)
                result(i) += (*this)(i, j) * op1(j);

        result.setCol();
        return result;
    }

    /** Matrix by Matrix multiplication
     * @ingroup MatrixOperators
     *
     * @code
     * C = A*B;
     * @endcode
     */
    Matrix2D<T> operator*(const Matrix2D<T>& op1) const
    {
        Matrix2D<T> result;
        if (mdimx != op1.mdimy)
            REPORT_ERROR(1102, "Not compatible sizes in matrix multiplication");

        result.initZeros(mdimy, op1.mdimx);
        for (int i = 0; i < mdimy; i++)
            for (int j = 0; j < op1.mdimx; j++)
                for (int k = 0; k < mdimx; k++)
                    result(i, j) += (*this)(i, k) * op1(k, j);
        return result;
    }

    /** Matrix summation
     * @ingroup MatrixOperators
     *
     * @code
     * C = A + B;
     * @endcode
     */
    Matrix2D<T> operator+(const Matrix2D<T>& op1) const
    {
        Matrix2D<T> result;
        if (mdimx != op1.mdimx || mdimy != op1.mdimy)
            REPORT_ERROR(1102, "Not same sizes in matrix summation");

        result.initZeros(mdimy, mdimx);
        for (int i = 0; i < mdimy; i++)
            for (int j = 0; j < mdimx; j++)
                result(i, j) = (*this)(i, j) + op1(i, j);

        return result;
    }

    /** Matrix summation
     * @ingroup MatrixOperators
     *
     * @code
     * A += B;
     * @endcode
     */
    void operator+=(const Matrix2D<T>& op1) const
    {
        if (mdimx != op1.mdimx || mdimy != op1.mdimy)
            REPORT_ERROR(1102, "Not same sizes in matrix summation");

        for (int i = 0; i < mdimy; i++)
            for (int j = 0; j < mdimx; j++)
                MAT_ELEM(*this,i, j) += MAT_ELEM(op1, i, j);
    }

    /** Matrix subtraction
     * @ingroup MatrixOperators
     *
     * @code
     * C = A - B;
     * @endcode
     */
    Matrix2D<T> operator-(const Matrix2D<T>& op1) const
    {
        Matrix2D<T> result;
        if (mdimx != op1.mdimx || mdimy != op1.mdimy)
            REPORT_ERROR(1102, "Not same sizes in matrix summation");

        result.initZeros(mdimy, mdimx);
        for (int i = 0; i < mdimy; i++)
            for (int j = 0; j < mdimx; j++)
                result(i, j) = (*this)(i, j) - op1(i, j);

        return result;
    }

    /** Matrix substraction
     * @ingroup MatrixOperators
     *
     * @code
     * A -= B;
     * @endcode
     */
    void operator-=(const Matrix2D<T>& op1) const
    {
        if (mdimx != op1.mdimx || mdimy != op1.mdimy)
            REPORT_ERROR(1102, "Not same sizes in matrix summation");

        for (int i = 0; i < mdimy; i++)
            for (int j = 0; j < mdimx; j++)
                MAT_ELEM(*this,i, j) -= MAT_ELEM(op1, i, j);
    }

    /** Equality.
     * @ingroup MatrixOperators
     *
     * Returns true if this object has got the same shape (origin and size)
     * than the argument and the same values (within accuracy).
     */
    bool equal(const Matrix2D<T>& op,
               double accuracy = XMIPP_EQUAL_ACCURACY) const
    {
        if (!sameShape(op))
            return false;
        for (int i = 0; i < mdimy; i++)
            for (int j = 0; j < mdimx; j++)
                if (ABS( (*this)(i,j) - op(i,j) ) > accuracy)
                    return false;
        return true;
    }

    /// @defgroup MatrixUtilities Utilities for Matrix2D
    /** Maximum of the values in the array.
      * @ingroup MatrixUtilities
      *
      * The returned value is of the same type as the type of the array.
      */
    T computeMax() const
    {
        if (mdim <= 0)
            return static_cast< T >(0);

        T maxval = mdata[0];
        for (int n = 0; n < mdim; n++)
            if (mdata[n] > maxval)
                maxval = mdata[n];
        return maxval;
    }

    /** Minimum of the values in the array.
       * @ingroup MatrixUtilities
       *
       * The returned value is of the same type as the type of the array.
       */
    T computeMin() const
    {
        if (mdim <= 0)
            return static_cast< T >(0);

        T minval = mdata[0];
        for (int n = 0; n < mdim; n++)
            if (mdata[n] < minval)
                minval = mdata[n];
        return minval;
    }

    /** Produce a 2D array suitable for working with Numerical Recipes
    * @ingroup MatrixUtilities
    *
    * This function must be used only as a preparation for routines which need
    * that the first physical index is 1 and not 0 as it usually is in C. New
    * memory is needed to hold the new double pointer array.
    */
    T** adaptForNumericalRecipes() const
    {
        T** m = NULL;
        ask_Tmatrix(m, 1, mdimy, 1, mdimx);

        for (int i = 0; i < mdimy; i++)
            for (int j = 0; j < mdimx; j++)
                m[i+1][j+1] = mdata[i*mdimx + j];

        return m;
    }

    /** Produce a 1D pointer suitable for working with Numerical Recipes (2)
     * @ingroup MatrixUtilities
     *
     * This function meets the same goal as the one before, however this one
     * work with 2D arrays as a single pointer. The first element of the array
     * is pointed by result[1*Xdim+1], and in general result[i*Xdim+j]
     */
    T* adaptForNumericalRecipes2() const
    {
        return mdata - 1 - mdimx;
    }

    /** Load 2D array from numerical recipes result.
     * @ingroup MatrixUtilities
     */
    void loadFromNumericalRecipes(T** m, int Ydim, int Xdim)
    {
        resize(Ydim, Xdim);

        for (int i = 1; i <= Ydim; i++)
            for (int j = 1; j <= Xdim; j++)
                (*this)(i - 1, j - 1) = m[i][j];
    }

    /** Kill a 2D array produced for numerical recipes
     * @ingroup MatrixUtilities
     *
     * The allocated memory is freed.
     */
    void killAdaptationForNumericalRecipes(T** m) const
    {
        free_Tmatrix(m, 1, mdimy, 1, mdimx);
    }

    /** Kill a 2D array produced for numerical recipes, 2.
     * @ingroup MatrixUtilities
     *
     * Nothing needs to be done.
     */
    void killAdaptationForNumericalRecipes2(T** m) const
        {}

    /** Write this matrix to file
      * @ingroup MatrixUtilities
      */
    void write(const FileName &fn) const
    {
        std::ofstream fhOut;
        fhOut.open(fn.c_str());
        if (!fhOut)
            REPORT_ERROR(1,(std::string)"Cannot open "+fn+" for output");
        fhOut << *this;
        fhOut.close();
    }

    /** Show matrix
      * @ingroup MatrixUtilities
      */
    friend std::ostream& operator<<(std::ostream& ostrm, const Matrix2D<T>& v)
    {
        if (v.Xdim() == 0 || v.Ydim() == 0)
            ostrm << "NULL matrix\n";
        else
        {
            ostrm << std::endl;
            double max_val = v.computeMax();
            int prec = bestPrecision(max_val, 10);

            for (int i = 0; i < v.Ydim(); i++)
            {
                for (int j = 0; j < v.Xdim(); j++)
                {
                    ostrm << floatToString((double) v(i, j), 10, prec) << ' ';
                }
                ostrm << std::endl;
            }
        }

        return ostrm;
    }

    /** Makes a matrix from a vector
     * @ingroup MatrixUtilities
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
        if (op1.size() == 0)
        {
            clear();
            return;
        }

        // Look at shape and copy values
        if (op1.isRow())
        {
            resize(1, op1.size());

            for (int j = 0; j < op1.size(); j++)
                (*this)(0, j) = op1( j);
        }
        else
        {
            resize(op1.size(), 1);

            for (int i = 0; i < op1.size(); i++)
                (*this)(i, 0) = op1( i);
        }
    }

    /** Makes a vector from a matrix
     * @ingroup MatrixUtilities
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
        if (mdimx == 0 || mdimy == 0)
        {
            op1.clear();
            return;
        }

        // If matrix is not a vector, produce an error
        if (!(mdimx == 1 || mdimy == 1))
            REPORT_ERROR(1102,
                         "To_vector: Matrix cannot be converted to vector");

        // Look at shape and copy values
        if (mdimy == 1)
        {
            // Row vector
            op1.resize(mdimx);

            for (int j = 0; j < mdimx; j++)
                op1(j) = (*this)(0, j);

            op1.setRow();
        }
        else
        {
            // Column vector
            op1.resize(mdimy);

            for (int i = 0; i < mdimy; i++)
                op1(i) = (*this)(i, 0);

            op1.setCol();
        }
    }

    /**Copy matrix to stl::vector
     * @ingroup MatrixUtilities
     */
    void copyToVector(std::vector<T> &v)
    {
        v.assign(mdata, mdata+mdim);
    }
    /**Copy stl::vector to matrix
     * @ingroup MatrixUtilities
      */
    void copyFromVector(std::vector<T> &v,int Xdim, int Ydim)
    {
        resize(Ydim, Xdim);
        copy( v.begin(), v.begin()+v.size(), mdata);
    }

    /** Get row
     * @ingroup MatrixUtilities
     *
     * This function returns a row vector corresponding to the choosen
     * row inside the nth 2D matrix, the numbering of the rows is also
     * logical not physical.
     *
     * @code
     * std::vector< double > v;
     * m.getRow(-2, v);
     * @endcode
     */
    void getRow(int i, Matrix1D<T>& v) const
    {
        if (mdimx == 0 || mdimy == 0)
        {
            v.clear();
            return;
        }

        if (i < 0 || i >= mdimy)
            REPORT_ERROR(1103, "getRow: Matrix subscript (i) greater than matrix dimension");

        v.resize(mdimx);
        for (int j = 0; j < mdimx; j++)
            v(j) = (*this)(i, j);

        v.setRow();
    }

    /** Get Column
     * @ingroup MatrixUtilities
     *
     * This function returns a column vector corresponding to the
     * choosen column.
     *
     * @code
     * std::vector< double > v;
     * m.getCol(-1, v);
     * @endcode
     */
    void getCol(int j, Matrix1D<T>& v) const
    {
        if (mdimx == 0 || mdimy == 0)
        {
            v.clear();
            return;
        }

        if (j < 0 || j >= mdimx)
            REPORT_ERROR(1103,"getCol: Matrix subscript (j) greater than matrix dimension");

        v.resize(mdimy);
        for (int i = 0; i < mdimy; i++)
            v(i) = (*this)(i, j);

        v.setCol();
    }

    /** Set Row
     * @ingroup MatrixUtilities
     *
     * This function sets a row vector corresponding to the choosen row in the 2D Matrix
     *
     * @code
     * m.setRow(-2, m.row(1)); // Copies row 1 in row -2
     * @endcode
     */
    void setRow(int i, const Matrix1D<T>& v)
    {
        if (mdimx == 0 || mdimy == 0)
            REPORT_ERROR(1, "setRow: Target matrix is empty");

        if (i < 0 || i >= mdimy)
            REPORT_ERROR(1103, "setRow: Matrix subscript (i) out of range");

        if (v.size() != mdimx)
            REPORT_ERROR(1102,
                         "setRow: Vector dimension different from matrix one");

        if (!v.isRow())
            REPORT_ERROR(1107, "setRow: Not a row vector in assignment");

        for (int j = 0; j < mdimx; j++)
            (*this)(i, j) = v(j);
    }

    /** Set Column
     * @ingroup MatrixUtilities
     *
     * This function sets a column vector corresponding to the choosen column
     * inside matrix.
     *
     * @code
     * m.setCol(0, (m.row(1)).transpose()); // Copies row 1 in column 0
     * @endcode
     */
    void setCol(int j, const Matrix1D<T>& v)
    {
        if (mdimx == 0 || mdimy == 0)
            REPORT_ERROR(1, "setCol: Target matrix is empty");

        if (j < 0 || j>= mdimx)
            REPORT_ERROR(1103, "setCol: Matrix subscript (j) out of range");

        if (v.size() != mdimy)
            REPORT_ERROR(1102,
                         "setCol: Vector dimension different from matrix one");

        if (!v.isCol())
            REPORT_ERROR(1107, "setCol: Not a column vector in assignment");

        for (int i = 0; i < mdimy; i++)
            (*this)(i, j) = v(i);
    }

    /** Determinant of a matrix
     * @ingroup MatrixUtilities
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
        if (mdimx == 0 || mdimy == 0)
            REPORT_ERROR(1108, "determinant: Matrix is empty");

        if (mdimx != mdimy)
            REPORT_ERROR(1109, "determinant: Matrix is not squared");

        for (int i = 0; i < mdimy; i++)
        {
            bool all_zeros = true;
            for (int j = 0; j < mdimx; j++)
                if (ABS(MAT_ELEM((*this),i, j)) > XMIPP_EQUAL_ACCURACY)
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
        for (int i = 0; i < mdimx; i++)
            d *= (T) MAT_ELEM(LU,i , i);

        return d;
    }

    /** Algebraic transpose of a Matrix
     * @ingroup MatrixUtilities
     *
     * You can use the transpose in as complex expressions as you like. The
     * origin of the vector is not changed.
     *
     * @code
     * v2 = v1.transpose();
     * @endcode
     */
    Matrix2D<T> transpose() const
    {
        Matrix2D<T> result(mdimx, mdimy);
        FOR_ALL_ELEMENTS_IN_MATRIX2D(result)
        MAT_ELEM(result,i,j) = MAT_ELEM((*this),j,i);
        return result;
    }

    /** Inverse of a matrix
     * @ingroup MatrixUtilities
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

        if (mdimx == 0 || mdimy == 0)
            REPORT_ERROR(1108, "Inverse: Matrix is empty");

        // Perform SVD decomposition
        Matrix2D< double > u, v;
        Matrix1D< double > w;
        svdcmp(*this, u, w, v); // *this = U * W * V^t

        double tol = computeMax() * XMIPP_MAX(mdimx, mdimy) * 1e-14;
        result.initZeros(mdimx, mdimy);

        // Compute W^-1
        bool invertible = false;
        FOR_ALL_ELEMENTS_IN_MATRIX1D(w)
        {
            if (ABS(VEC_ELEM(w,i)) > tol)
            {
                VEC_ELEM(w,i) = 1.0 / VEC_ELEM(w,i);
                invertible = true;
            }
            else
                VEC_ELEM(w,i) = 0.0;
        }

        if (!invertible)
            return;

        // Compute V*W^-1
        FOR_ALL_ELEMENTS_IN_MATRIX2D(v)
        MAT_ELEM(v,i,j) *= VEC_ELEM(w,j);

        // Compute Inverse
        for (int i = 0; i < mdimx; i++)
            for (int j = 0; j < mdimy; j++)
                for (int k = 0; k < mdimx; k++)
                    MAT_ELEM(result,i,j) += (T) MAT_ELEM(v,i,k) * MAT_ELEM(u,j,k);
    }

    /** Inverse of a matrix
     * @ingroup MatrixUtilities
     */
    Matrix2D<T> inv() const
    {
        Matrix2D<T> result;
        inv(result);

        return result;
    }

    /** True if the matrix is diagonal
     * @ingroup MatrixUtilities
     *
     * @code
     * if (m.isDiagonal())
     *     std::cout << "The matrix is diagonal\n";
     * @endcode
     */
    bool isDiagonal() const
    {
        if (mdimx != mdimy)
            return false;

        for (int i = 0; i < mdimy; i++)
            for (int j = 0; j < mdimx; j++)
                if (i != j && ABS((*this)(i, j)) >
                    XMIPP_EQUAL_ACCURACY)
                    return false;
        return true;
    }

    /** True if the matrix is scalar
     * @ingroup MatrixUtilities
     *
     * A scalar matrix is diagonal and all the values at the diagonal are the
     * same
     */
    bool isScalar() const
    {
        if (!isDiagonal())
            return false;

        for (int i = 1; i < mdimy; i++)
            if (ABS((*this)(i, i) - (*this)(0, 0))
                > XMIPP_EQUAL_ACCURACY)
                return false;

        return true;
    }

    /** True if the matrix is identity
     * @ingroup MatrixUtilities
     *
     * @code
     * if (m.isIdentity())
     *     std::cout << "The matrix is identity\n";
     * @endcode
     */
    bool isIdentity() const
    {
        return isScalar() &&
               ABS((*this)(0, 0) - (T) 1) < XMIPP_EQUAL_ACCURACY;
    }

};

// Implementation of the vector*matrix
// Documented in matrix1D.h
template<typename T>
Matrix1D<T> Matrix1D<T>::operator*(const Matrix2D<T>& M)
{
    Matrix1D<T> result;

    if (VEC_XSIZE(*this) != MAT_YSIZE(M))
        REPORT_ERROR(1102, "Not compatible sizes in matrix by vector");

    if (!isRow())
        REPORT_ERROR(1102, "Vector is not a row");

    result.initZeros(MAT_XSIZE(M));
    for (int j = 0; j < MAT_XSIZE(M); j++)
        for (int i = 0; i < MAT_YSIZE(M); i++)
            VEC_ELEM(result,j) += VEC_ELEM(*this,i) * MAT_ELEM(M,i, j);

    result.setRow();
    return result;
}

/**@defgroup MatrixRelated Related functions
 * @ingroup Matrices
 *
 * These functions are not methods of Matrix2D
 */

/** LU Decomposition
 * @ingroup MatrixRelated
 */
template<typename T>
void ludcmp(const Matrix2D<T>& A, Matrix2D<T>& LU, Matrix1D< int >& indx, T& d)
{
    LU = A;
    indx.resize(A.mdimx);
    ludcmp(LU.adaptForNumericalRecipes2(), A.mdimx,
           indx.adaptForNumericalRecipes(), &d);
}

/** LU Backsubstitution
 * @ingroup MatrixRelated
 */
template<typename T>
void lubksb(const Matrix2D<T>& LU, Matrix1D< int >& indx, Matrix1D<T>& b)
{
    lubksb(LU.adaptForNumericalRecipes2(), indx.size(),
           indx.adaptForNumericalRecipes(),
           b.adaptForNumericalRecipes());
}

/** SVD Backsubstitution
 * @ingroup MatrixRelated
 */
void svbksb(Matrix2D< double >& u,
            Matrix1D< double >& w,
            Matrix2D< double >& v,
            Matrix1D< double >& b,
            Matrix1D< double >& x);

/** SVD Decomposition
 * @ingroup MatrixRelated
 */
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
    w.initZeros(u.mdimx);
    v.initZeros(u.mdimx, u.mdimx);

    // Call to the numerical recipes routine
#ifdef VIA_NR

    svdcmp(u.mdata,
           u.mdimy, u.mdimx,
           w.vdata,
           v.vdata);
#endif

#ifdef VIA_BILIB

    int status;
    SingularValueDecomposition(u.mdata,
                               u.mdimy, u.mdimx,
                               w.vdata,
                               v.mdata,
                               5000, &status);
#endif
}
#undef VIA_NR
#undef VIA_BILIB

/** Conversion from one type to another.
 * @ingroup MatrixRelated
 *
 * If we have an integer array and we need a double one, we can use this
 * function. The conversion is done through a type casting of each element
 * If n >= 0, only the nth volumes will be converted, otherwise all NSIZE volumes
 */
template<typename T1, typename T2>
void typeCast(const Matrix2D<T1>& v1,  Matrix2D<T2>& v2)
{
    if (v1.mdim == 0)
    {
        v2.clear();
        return;
    }

    v2.resize(v1);
    for (int n = 0; n < v1.mdim; n++)
        v2.mdata[n] = static_cast< T2 > (v1.mdata[n]);

}

#endif /* MATRIX2D_H_ */
