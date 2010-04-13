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

#include "matrix1d.h"

/** @defgroup Matrices Matrices speed up macros
 * @ingroup XXX
 *
 */

/** Matrix (3x3) by vector (3x1) (a=M*b)
 * @ingroup Matrices
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
 * @ingroup Matrices
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
 * @ingroup Matrices
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
 * @ingroup Matrices
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
 * @ingroup Matrices
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


template<typename T> class Matrix1D;
template<typename T> class Matrix2D;

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
class Matrix2D
{
public:
    /* The array itself.
    */
    T* data;

    // Destroy data
    bool destroyData;

    // Number of elements in X
    int xdim;

    // Number of elements in Y
    int ydim;

    // Total number of elements
    int xydim;

    /// @defgroup VectorsConstructors Constructors
    /// @ingroup Vectors
    /** Empty constructor
     * @ingroup VectorsConstructors
     *
     * The empty constructor creates a vector with no memory associated,
     * origin=0, size=0, no statistics, ... You can choose between a column
     * vector (by default), or a row one.
     *
     * @code
     * Matrix1D< double > v1;
     * Matrix1D< double > v1(true);
     * // both are examples of empty column vectors
     *
     * Matrix1D< int > v1(false);
     * // empty row vector
     * @endcode
     */
    Matrix2D()
    {
    	coreInit();
    }

    /** Dimension constructor
     * @ingroup VectorsConstructors
     *
     * The dimension constructor creates a vector with memory associated (but
     * not assigned to anything, could be full of garbage) origin=0, size=the
     * given one. You can choose between a column vector (by default), or a row
     * one.
     *
     * @code
     * Matrix1D< double > v1(6);
     * Matrix1D< double > v1(6, 'y');
     * // both are examples of column vectors of dimensions 6
     *
     * Matrix1D< int > v1('n');
     * // empty row vector
     * @endcode
     */
    Matrix2D(int Ydim, int Xdim)
    {
        resize(Ydim, Xdim);
    }

    /** Copy constructor
     * @ingroup VectorsConstructors
     *
     * The created vector is a perfect copy of the input vector but with a
     * different memory assignment.
     *
     * @code
     * Matrix1D< double > v2(v1);
     * @endcode
     */
    Matrix2D(const Matrix2D<T>& v)
    {
        coreInit();
        *this = v;
    }

    /** Destructor.
     * @ingroup MultidimArrayConstructors
     */
     ~Matrix2D()
     {
        coreDeallocate();
     }

    /** Clear.
     * @ingroup MultidimArrayConstructors
     */
     void clear()
     {
        coreDeallocate();
        coreInit();
     }

    /// @defgroup MultidimArrayCore Core memory operations for MultidimArrays
    /// @ingroup MultidimensionalArrays
    /** Core init.
     * @ingroup MultidimArrayCore
     * Initialize everything to 0
     */
    void coreInit()
    {
        xdim=ydim=xydim=0;
        data=NULL;
        destroyData=true;
    }

    /** Core allocate.
     * @ingroup MultidimArrayCore
     */
    void coreAllocate(int _ydim, int _xdim)
    {
        if (_ydim <= 0 ||_xdim<=0)
        {
            clear();
            return;
        }

        xdim=_xdim;
        ydim=_ydim;
        xydim=_xdim*_ydim;
        data = new T [xydim];
        if (data == NULL)
            REPORT_ERROR(1001, "Allocate: No space left");
    }

    /** Core deallocate.
     * @ingroup MultidimArrayCore
     * Free all data.
     */
    void coreDeallocate()
    {
        if (data != NULL && destroyData)
            delete[] data;
        data=NULL;
    }

    /** Resize to a given size
     * @ingroup MultidimSize
     *
     * This function resize the actual array to the given size. The origin is
     * not modified. If the actual array is larger than the pattern then the
     * values outside the new size are lost, if it is smaller then 0's are
     * added. An exception is thrown if there is no memory.
     *
     * @code
     * V1.resize(3, 3, 2);
     * @endcode
     */
    void resize(int Ydim, int Xdim)
    {
        if (Xdim == xdim && Ydim == ydim)
            return;

        if (Xdim <= 0 || Ydim <= 0)
        {
            clear();
            return;
        }


        T * new_data;
        size_t YXdim=Ydim*Xdim;

        try
        {
        	new_data = new T [YXdim];
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
                if (i >= ydim)
                    val = 0;
                else if (j >= xdim)
                    val = 0;
                else
                    val = data[i*xdim + j];
                new_data[i*Xdim+j] = val;
            }

		// deallocate old vector
		coreDeallocate();

		// assign *this vector to the newly created
		data = new_data;
		xdim = Xdim;
		ydim = Ydim;
		xydim = Xdim * Ydim;
    }

    /** Resize according to a pattern.
     * @ingroup MultidimSize
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
        if (xdim != v.xdim || ydim != v.ydim)
            resize(v.ydim, v.xdim);
    }


    /** Same value in all components.
     * @ingroup Initialization
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
    	for (int j = 0; j < xydim; j++)
    	{
    		data[j] = val;
    	}
    }

    /** Initialize to zeros with current size.
     * @ingroup Initialization
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
     * @ingroup Initialization
     */
    void initZeros(int Ydim, int Xdim)
    {
        resize(Ydim, Xdim);
        initZeros();
    }

    /** 2D Identity matrix of current size
     * @ingroup Initialization
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
        initIdentity(xdim, xdim);
    }

    /** 2D Identity matrix of a given size
     * @ingroup Initialization
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

    /** 2D Identity matrix of a given size
     * @ingroup Initialization
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
		for (int i = 0; i < ydim; i++)
			for (int j = 0; j < xdim; j++)
				(*this)(i,j) = (T)(i == j);
    }

    /** Same shape.
     * @ingroup MultidimSize
     *
     * Returns true if this object has got the same shape (origin and size)
     * than the argument
     */
    template <typename T1>
    bool sameShape(const Matrix2D<T1>& op) const
    {
        return ((xdim == op.xdim) && (ydim == op.ydim));
    }

    /// @defgroup VectorsMemory Memory access
    /// @ingroup Vectors

    /** Vector element access
     * @ingroup VectorsMemory
     *
     * Returns the value of a vector logical position. In our example we could
     * access from v(-2) to v(2). The elements can be used either by value or by
     * reference.
     *
     * @code
     * v(-2) = 1;
     * val = v(-2);
     * @endcode
     */
    T& operator()(int i, int j) const
    {
        return data[i*xdim + j];
    }

    /** Assignment.
     * @ingroup VectorOperators
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
            for (int i = 0; i< op1.xydim; i++)
				data[i] = op1.data[i];
        }

        return *this;
    }


    /** Maximum of the values in the array.
      * @ingroup Statistics
      *
      * The returned value is of the same type as the type of the array.
      */
     T computeMax() const
     {
         if (xydim <= 0)
             return static_cast< T >(0);

         T maxval = data[0];
         for (int n = 0; n < xydim; n++)
        	 if (data[n] > maxval)
        		 maxval = data[n];
         return maxval;
     }

     /** Minimum of the values in the array.
        * @ingroup Statistics
        *
        * The returned value is of the same type as the type of the array.
        */
       T computeMin() const
       {
           if (xydim <= 0)
               return static_cast< T >(0);

           T minval = data[0];
           for (int n = 0; n < xydim; n++)
          	 if (data[n] < minval)
          		 minval = data[n];
           return minval;
       }

       /** Produce a 2D array suitable for working with Numerical Recipes
     * @ingroup MatricesSize
     *
     * This function must be used only as a preparation for routines which need
     * that the first physical index is 1 and not 0 as it usually is in C. New
     * memory is needed to hold the new double pointer array.
     */
    T** adaptForNumericalRecipes2D(unsigned long n = 0) const
    {
        T** m = NULL;
        ask_Tmatrix(m, 1, ydim, 1, xdim);

		for (int i = 0; i < ydim; i++)
			for (int j = 0; j < xdim; j++)
				m[i+1][j+1] = data[i*xdim + j];

        return m;
    }

    /** Produce a 1D pointer suitable for working with Numerical Recipes (2)
     * @ingroup MultidimSize
     *
     * This function meets the same goal as the one before, however this one
     * work with 2D arrays as a single pointer. The first element of the array
     * is pointed by result[1*Xdim+1], and in general result[i*Xdim+j]
     */
    T* adaptForNumericalRecipes22D() const
    {
        return data - 1 - xdim;
    }

    /** Load 2D array from numerical recipes result.
     * @ingroup MatricesSize
     */
    void loadFromNumericalRecipes2D(T** m, int Ydim, int Xdim)
    {
        resize(Ydim, Xdim);

        for (int i = 1; i <= Ydim; i++)
            for (int j = 1; j <= Xdim; j++)
                (*this)(i - 1, j - 1) = m[i][j];
    }

    /** Kill a 2D array produced for numerical recipes
     * @ingroup MatricesSize
     *
     * The allocated memory is freed.
     */
    void killAdaptationForNumericalRecipes2D(T** m) const
    {
        free_Tmatrix(m, 1, ydim, 1, xdim);
    }

    /** Kill a 2D array produced for numerical recipes, 2.
     * @ingroup MatricesSize
     *
     * Nothing needs to be done.
     */
    void killAdaptationForNumericalRecipes22D(T** m) const
        {}

    /// @defgroup MatricesUtilities
    /// @ingroup Matrices

    /** True if the matrix is diagonal
     * @ingroup MatricesUtilities
     *
     * @code
     * if (m.isDiagonal())
     *     std::cout << "The matrix is diagonal\n";
     * @endcode
     */
    bool isDiagonal() const
    {
        if (xdim != ydim)
            return false;

        for (int i = 0; i < ydim; i++)
        	for (int j = 0; j < xdim; j++)
        		if (i != j && ABS((*this)(i, j)) >
        		            XMIPP_EQUAL_ACCURACY)
        		            return false;
        return true;
    }

    /** True if the matrix is scalar
     * @ingroup MatricesUtilities
     *
     * A scalar matrix is diagonal and all the values at the diagonal are the
     * same
     */
    bool isScalar() const
    {
        if (!isDiagonal())
            return false;

        for (int i = 1; i < ydim; i++)
            if (ABS((*this)(i, i) - (*this)(0, 0))
                > XMIPP_EQUAL_ACCURACY)
                return false;

        return true;
    }

    /** True if the matrix is identity
     * @ingroup MatricesUtilities
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

    /** Makes a matrix from a vector
     * @ingroup MatricesUtilities
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
        if (op1.xdim == 0)
        {
            clear();
            return;
        }

        // Look at shape and copy values
        if (op1.isRow())
        {
            resize(1, op1.xdim);

            for (int j = 0; j < op1.xdim; j++)
                (*this)(0, j) = op1( j);
        }
        else
        {
            resize(op1.xdim, 1);

            for (int i = 0; i < op1.xdim; i++)
            	(*this)(i, 0) = op1( i);
        }
    }

    /** Makes a vector from a matrix
     * @ingroup MatricesUtilities
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
        if (xdim == 0 || ydim == 0)
        {
            op1.clear();
            return;
        }

        // If matrix is not a vector, produce an error
        if (!(xdim == 1 || ydim == 1))
            REPORT_ERROR(1102,
                         "To_vector: Matrix cannot be converted to vector");

        // Look at shape and copy values
        if (ydim == 1)
        {
            // Row vector
            op1.resize(xdim);

            for (int j = 0; j < xdim; j++)
				op1(j) = (*this)(0, j);

            op1.setRow();
        }
        else
        {
            // Column vector
            op1.resize(ydim);

            for (int i = 0; i < ydim; i++)
				op1(i) = (*this)(i, 0);

            op1.setCol();
		}
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
        if (xdim == 0 || ydim == 0)
            REPORT_ERROR(1108, "determinant: Matrix is empty");

        if (xdim != ydim)
            REPORT_ERROR(1109, "determinant: Matrix is not squared");

        for (int i = 0; i < ydim; i++)
        {
            bool all_zeros = true;
            for (int j = 0; j < xdim; j++)
                if (ABS((*this)(i, j)) > XMIPP_EQUAL_ACCURACY)
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
        for (int i = 0; i < xdim; i++)
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

        if (xdim == 0 || ydim == 0)
            REPORT_ERROR(1108, "Inverse: Matrix is empty");

        // Perform SVD decomposition
        Matrix2D< double > u, v;
        Matrix1D< double > w;
        svdcmp(*this, u, w, v); // *this = U * W * V^t

        //FIXME define omputeMax
        double tol = computeMax() * XMIPP_MAX(xdim, ydim) * 1e-14;
        result.initZeros(xdim, ydim);

        // Compute W^-1
        bool invertible = false;
        for (int i = 0; i < w.xdim; i++)
        {
            if (ABS(w(i)) > tol)
            {
            	w(i) = 1.0 / w(i);
                invertible = true;
            }
            else
            	w(i) = 0.0;
        }

        if (!invertible)
            return;

        // Compute V*W^-1
        for (int i = 0; i < v.ydim; i++)
        	for (int j = 0; j < v.xdim; j++)
        		v(i,j) = w(j);

        // Compute Inverse
        for (int i = 0; i < xdim; i++)
            for (int j = 0; j < ydim; j++)
                for (int k = 0; k < xdim; k++)
					result(i,j) += (T) v(i,k) * u(j,k);
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


}; // class Matrix2D


///FIXME: implement these as real multiplications!!
/**@defgroup MatricesRelated Related functions
 * @ingroup MultidimensionalArrays
 *
 * These functions are not methods of MultidimArray
 */

/* Matrix x Matrix multiplication
 * @ingroup MatricesRelated
 * This is a matrix multiplication, not an element-by-element multiplication!
 *
 * result = A * A;
 */
template<typename T>
void multiplyMatrixbyMatrix(const Matrix2D<T>& op1, const Matrix2D<T>& op2, Matrix2D<T>& result)
{
     if (op1.xdim != op2.ydim)
        REPORT_ERROR(1102, "Not compatible sizes in matrix multiplication");

    result.initZeros(op1.ydim, op2.xdim);
    for (int i = 0; i < op1.ydim; i++)
        for (int j = 0; j < op2.xdim; j++)
            for (int k = 0; k < op1.xdim; k++)
                result(i, j) += op1(i, k) * op2(k, j);

 }

/* Matrix x Vector multiplication
 * @ingroup MatricesRelated
 * This is a matrix multiplication, not an element-by-element multiplication!
 *
 * result = A * b;
 * b should be a column vector, result will be a row vector
 *
 */
template<typename T>
void multiplyMatrixbyVector(const Matrix2D<T>& A, const Matrix1D<T>& b, Matrix1D<T>& result)
{

    if (A.xdim != b.xdim)
        REPORT_ERROR(1102, "Not compatible sizes in matrix by vector");

    if (!b.isCol())
        REPORT_ERROR(1102, "Vector is not a column");

    result.initZeros(A.ydim);

    for (int i = 0; i < A.ydim; i++)
        for (int j = 0; j < b.xdim; j++)
            result(i) += A(i, j) * b(j);

    result.setCol();
    return result;


}

/* Vector x Matrix multiplication
 * @ingroup MatricesRelated
 * This is a matrix multiplication, not an element-by-element multiplication!
 *
 * result = b * A;
 * b should be a row vector, result will be a column vector
 */
template<typename T>
void multiplyVectorbyMatrix(const Matrix1D<T>& b, const Matrix2D<T>& A, Matrix1D<T>& result)
{
    if (b.xdim != A.ydim)
        REPORT_ERROR(1102, "Not compatible sizes in matrix by vector");

    if (!b.isRow())
        REPORT_ERROR(1102, "Vector is not a row");

    result.initZeros(A.xdim);
    for (int j = 0; j < A.xdim; j++)
        for (int i = 0; i < A.ydim; i++)
            result(j) += b(i) * A( i, j);

    result.setRow();
    return result;

}

// TODO Document
template<typename T>
void ludcmp(const Matrix2D<T>& A, Matrix2D<T>& LU, Matrix1D< int >& indx, T& d)
{
    LU = A;
    indx.resize(A.xdim);
    ludcmp(LU.adaptForNumericalRecipes2(), A.xdim,
           indx.adaptForNumericalRecipes(), &d);
}

// TODO Document
template<typename T>
void lubksb(const Matrix2D<T>& LU, Matrix1D< int >& indx, Matrix1D<T>& b)
{
    lubksb(LU.adaptForNumericalRecipes2(), indx.xdim,
           indx.adaptForNumericalRecipes(),
           b.adaptForNumericalRecipes());
}

// TODO Document
void svbksb(Matrix2D< double >& u,
            Matrix1D< double >& w,
            Matrix2D< double >& v,
            Matrix1D< double >& b,
            Matrix1D< double >& x);

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
    w.initZeros(u.xdim);
    v.initZeros(u.xdim, u.xdim);

    // Call to the numerical recipes routine
#ifdef VIA_NR
    svdcmp(u.data,
           u.ydim, u.xdim,
           w.data,
           v.data);
#endif

#ifdef VIA_BILIB
    int status;
    SingularValueDecomposition(u.data,
                               u.ydim, u.xdim,
                               w.data,
                               v.data,
                               5000, &status);
#endif
}

#undef VIA_NR
#undef VIA_BILIB

/** Conversion from one type to another.
 * @ingroup MultidimFunctions
 *
 * If we have an integer array and we need a double one, we can use this
 * function. The conversion is done through a type casting of each element
 * If n >= 0, only the nth volumes will be converted, otherwise all NSIZE volumes
 */
template<typename T1, typename T2>
void typeCast(const Matrix2D<T1>& v1,  Matrix2D<T2>& v2)
{
    if (v1.xydim == 0)
    {
        v2.clear();
        return;
    }

    v2.resize(v1.ydim, v1.xdim);
	for (int n = 0; n < v1.xydim; n++)
		v2.data[n] = static_cast< T2 > (v1.data[n]);

}


#endif /* MATRIX2D_H_ */
