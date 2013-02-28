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

#ifndef MATRIX1D_H_
#define MATRIX1D_H_

#include <stdlib.h>
#include <cmath>
#include "xmipp_funcs.h"
#include "numerical_recipes.h"

extern int bestPrecision(float F, int _width);
extern std::string floatToString(float F, int _width, int _prec);

template<typename T>
class Matrix2D;

/** @defgroup Vectors Matrix1D Vectors
 * @ingroup DataLibrary
 */
//@{
/** @name Vectors speed up macros
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
//@{
/** Array access.
 * This macro gives you access to the array (T)
 */
#define MATRIX1D_ARRAY(v) ((v).vdata)

/** For all elements in the array
 * This macro is used to generate loops for the vector in an easy manner. It
 * defines an internal index 'i' which ranges the vector using its mathematical
 * definition (ie, logical access).
 *
 * @code
 * FOR_ALL_ELEMENTS_IN_MATRIX1D(v)
 * {
 *     std::cout << v(i) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_ELEMENTS_IN_MATRIX1D(v) \
    for (size_t i=0; i<v.vdim; i++)

/** X dimension of the matrix
 */
#define VEC_XSIZE(m) ((m).vdim)

/** Access to X component
 * @code
 * XX(v) = 1;
 * val = XX(v);
 * @endcode
 */
#define XX(v) (v).vdata[0]

/** Access to Y component
 * @code
 * YY(v) = 1;
 * val = YY(v);
 * @endcode
 */
#define YY(v) (v).vdata[1]

/** Access to Z component
 * @code
 * ZZ(v) = 1;
 * val = ZZ(v);
 * @endcode
 */
#define ZZ(v) (v).vdata[2]

/** Creates vector in R2
 * The vector must be created beforehand to the correct size. After this macro
 * the vector is (x, y) in R2.
 *
 * @code
 * MultidimArray< double > v(2);
 * VECTOR_R2(v, 1, 2);
 * @endcode
 */
#define VECTOR_R2(v, x, y) { \
        XX(v) = x; YY(v) = y; }

/** Creates vector in R3
 * The vector must be created beforehand to the correct size. After this macro
 * the vector is (x, y, z) in R3.
 *
 * @code
 * MultidimArray< double > v(3);
 * VECTOR_R2(v, 1, 2, 1);
 * @endcode
 */
#define VECTOR_R3(v, x, y, z) { \
        XX(v) = x; YY(v) = y; ZZ(v) = z;}

/** Adding two R2 vectors (a=b+c)
 * @code
 * MultidimArray< double > a(2), b(2), c(2);
 * ...;
 * V2_PLUS_V2(a, b, c);
 * @endcode
 */
#define V2_PLUS_V2(a, b, c) { \
        XX(a) = XX(b) + XX(c); \
        YY(a) = YY(b) + YY(c); }

/** Substracting two R2 vectors (a=b-c)
 * @code
 * MultidimArray< double > a(2), b(2), c(2);
 * ...;
 * V2_MINUS_V2(a, b, c);
 * @endcode
 */
#define V2_MINUS_V2(a, b, c) { \
        XX(a) = XX(b) - XX(c); \
        YY(a) = YY(b) - YY(c); }

/** Adding/substracting a constant to a R2 vector (a=b-k).
 * @code
 * MultidimArray< double > a(2), b(2);
 * double k;
 * ...;
 * V2_PLUS_CT(a, b, k);
 *
 * MultidimArray< double > a(2), b(2);
 * double k;
 * ...;
 * V2_PLUS_CT(a, b, -k);
 * @endcode
 */
#define V2_PLUS_CT(a, b, k) { \
        XX(a) = XX(b) + (k); \
        YY(a) = YY(b) + (k); }

/** Multiplying/dividing by a constant a R2 vector (a=b*k)
 * @code
 * MultidimArray< double > a(2), b(2);
 * double k;
 * ...;
 * V2_BY_CT(a, b, k);
 *
 * MultidimArray< double > a(2), b(2);
 * double k;
 * ...;
 * V2_BY_CT(a, b, 1/k);
 * @endcode
 */
#define V2_BY_CT(a, b, k) { \
        XX(a) = XX(b) * (k); \
        YY(a) = YY(b) * (k); }

/** Adding two R3 vectors (a=b+c)
 * @code
 * MultidimArray< double > a(3), b(3), c(3);
 * ...;
 * V3_PLUS_V3(a, b, c);
 * @endcode
 */
#define V3_PLUS_V3(a, b, c) { \
        XX(a) = XX(b) + XX(c); \
        YY(a) = YY(b) + YY(c); \
        ZZ(a) = ZZ(b) + ZZ(c); }

/** Substracting two R3 vectors (a=b-c)
 * @code
 * MultidimArray< double > a(3), b(3), c(3);
 * ...;
 * V3_MINUS_V3(a, b, c);
 * @endcode
 */
#define V3_MINUS_V3(a, b, c) { \
        XX(a) = XX(b) - XX(c); \
        YY(a) = YY(b) - YY(c); \
        ZZ(a) = ZZ(b) - ZZ(c); }

/** Adding/substracting a constant to a R3 vector (a=b-k)
 * @code
 * MultidimArray< double > a(3), b(3);
 * double k;
 * ...;
 * V3_PLUS_CT(a, b, k);
 *
 * MultidimArray< double > a(3), b(3);
 * double k;
 * ...;
 * V3_PLUS_CT(a, b, -k);
 * @endcode
 */
#define V3_PLUS_CT(a, b, c) { \
        XX(a) = XX(b) + (c); \
        YY(a) = YY(b) + (c); \
        ZZ(a) = ZZ(b) + (c); }

/** Multiplying/dividing by a constant a R3 vector (a=b*k)
 * @code
 * MultidimArray< double > a(3), b(3);
 * double k;
 * ...;
 * V3_BY_CT(a, b, k);
 *
 * MultidimArray< double > a(3), b(3);
 * double k;
 * ...;
 * V3_BY_CT(a, b, 1/k);
 * @endcode
 */
#define V3_BY_CT(a, b, c) { \
        XX(a) = XX(b) * (c); \
        YY(a) = YY(b) * (c); \
        ZZ(a) = ZZ(b) * (c); }

/** Direct access to vector element
 */
#define VEC_ELEM(v,i) ((v).vdata[(i)])
#define dMi(v, i) ((v).vdata[(i)])

#define VEC_SWAP(v, i, j, aux) {aux = dMi(v, i); dMi(v, i) = dMi(v, j); dMi(v, j) = aux; }

//@}

/** Matrix1D class.*/
template<typename T>
class Matrix1D
{
public:
    /// The array itself
    T* vdata;

    /// Destroy data
    bool destroyData;

    /// Number of elements
    size_t vdim;

    /// <0=column vector (default), 1=row vector
    bool row;

    /// @name Constructors
    //@{
    /** Empty constructor
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
    Matrix1D(bool column = true)
    {
        coreInit();
        row = !column;
    }

    /** Dimension constructor
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
    template<class T1>
    Matrix1D(T1 dim, bool column = true)
    {
        coreInit();
        row = !column;
        initZeros(dim);
    }

    /** Copy constructor
     *
     * The created vector is a perfect copy of the input vector but with a
     * different memory assignment.
     *
     * @code
     * Matrix1D< double > v2(v1);
     * @endcode
     */
    Matrix1D(const Matrix1D<T>& v)
    {
        coreInit();
        *this = v;
        row = v.row;
    }

    /** Destructor.
     */
    ~Matrix1D()
    {
        coreDeallocate();
    }

    /** Assignment.
     *
     * You can build as complex assignment expressions as you like. Multiple
     * assignment is allowed.
     *
     * @code
     * v1 = v2 + v3;
     * v1 = v2 = v3;
     * @endcode
     */
    Matrix1D<T>& operator=(const Matrix1D<T>& op1)
    {
        if (&op1 != this)
        {
            resizeNoCopy(op1);
            memcpy(vdata, op1.vdata, vdim * sizeof(T));
            row = op1.row;
        }

        return *this;
    }

    /** Equals.
     */
    bool operator==(const Matrix1D<T>& op1) const
    {
        if (row != op1.row)
            return false;

        for (size_t i = 0; i < vdim; ++i)
            if (!XMIPP_EQUAL_REAL(vdata[i], op1.vdata[i]))
                return false;
        return true;
    }

    /** Assignment.
     */
    Matrix1D<T>& operator=(const std::vector<T>& op1)
    {
        resizeNoCopy(op1.size());
        memcpy(&vdata[0],&(op1[0]),vdim*sizeof(T));
        row = false;
        return *this;
    }
    //@}

    /// @name Core memory operations for Matrix1D
    //@{
    /** Clear.
     */
    void clear()
    {
        coreDeallocate();
        coreInit();
    }

    /** Core init.
     * Initialize everything to 0
     */
    void coreInit()
    {
        vdim = 0;
        row = false;
        vdata = NULL;
        destroyData = true;
    }

    /** Core allocate.
     */
    inline void coreAllocate(int _vdim)
    {
        if (_vdim <= 0)
        {
            clear();
            return;
        }

        vdim = _vdim;
        vdata = new T[vdim];
        memset(vdata, 0, vdim * sizeof(T));
        if (vdata == NULL
           )
            REPORT_ERROR(ERR_MEM_NOTENOUGH, "Allocate: No space left");
    }

    /** Core deallocate.
     * Free all vdata.
     */
    inline void coreDeallocate()
    {
        if (vdata != NULL && destroyData)
            delete[] vdata;
        vdata = NULL;
    }
    //@}

    ///@name Size and shape of Matrix1D
    //@{
    /** Resize to a given size
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
    inline void resize(size_t Xdim, bool copy = true)
    {
        if (Xdim == vdim)
            return;

        if (Xdim <= 0)
        {
            clear();
            return;
        }

        T * new_vdata;
        try
        {
            new_vdata = new T[(size_t) Xdim]; //TODO: Kino: Valgrind points here there is something wrong, but I could'nt resolve
            memset(new_vdata, 0, Xdim * sizeof(T));
        }
        catch (std::bad_alloc &)
        {
            REPORT_ERROR(ERR_MEM_NOTENOUGH, "Allocate: No space left");
        }

        // Copy needed elements, fill with 0 if necessary
        if (copy)
        {
        	T zero=0; // Useful for complexes
            for (size_t j = 0; j < Xdim; j++)
            {
                T *val=NULL;
                if (j >= vdim)
                    val = &zero;
                else
                    val = &vdata[j];
                new_vdata[j] = *val;
            }
        }

        // deallocate old vector
        coreDeallocate();

        // assign *this vector to the newly created
        vdata = new_vdata;
        vdim = Xdim;

    }

    /** Resize a single 1D image with no copy
     */
    inline void resizeNoCopy(int Xdim)
    {
        resize(Xdim, false);
    }

    /** Resize according to a pattern.
     *
     * This function resize the actual array to the same size
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
    void resize(const Matrix1D<T1> &v)
    {
        if (vdim != v.vdim)
            resize(v.vdim);
    }

    /** Resize a single 1D image with no copy
     */
    template<typename T1>
    void resizeNoCopy(const Matrix1D<T1> &v)
    {
        if (vdim != v.vdim)
            resize(v.vdim, false);
    }

    /** Same shape.
     *
     * Returns true if this object has got the same shape (origin and size)
     * than the argument
     */
    template<typename T1>
    bool sameShape(const Matrix1D<T1>& op) const
    {
        return (vdim == op.vdim);
    }

    /** Returns the size of this vector
     *
     * @code
     * int nn = a.size();
     * @endcode
     */
    inline size_t size() const
    {
        return vdim;
    }

    /** True if vector is a row.
     *
     * @code
     * if (v.isRow())
     *     std::cout << "v is a row vector\n";
     * @endcode
     */
    int isRow() const
    {
        return row;
    }

    /** True if vector is a column
     *
     * @code
     * if (v.isCol())
     *     std::cout << "v is a column vector\n";
     * @endcode
     */
    int isCol() const
    {
        return !row;
    }

    /** Forces the vector to be a row vector
     *
     * @code
     * v.setRow();
     * @endcode
     */
    void setRow()
    {
        row = true;
    }

    /** Forces the vector to be a column vector
     *
     * @code
     * v.setCol();
     * @endcode
     */
    void setCol()
    {
        row = false;
    }
    //@}

    /// @name Initialization of Matrix1D values
    //@{
    /** Same value in all components.
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
        for (size_t j = 0; j < vdim; j++)
            vdata[j] = val;
    }

    /** Initialize to constant with a given size.
     */
    void initConstant(size_t Xdim, T val)
    {
        if (vdim != Xdim)
            resizeNoCopy(Xdim);
        initConstant(val);
    }

    /** Enumerate starting from 0 to size
     */
    void enumerate()
    {
      for (size_t j = 0; j < vdim; j++)
          vdata[j] = j;
    }

    /** Initialize to zeros with current size.
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
        memset(vdata, 0, vdim * sizeof(T));
    }

    /** Initialize to zeros with a given size.
     */
    void initZeros(size_t Xdim)
    {
        if (vdim != Xdim)
            resizeNoCopy(Xdim);
        memset(vdata, 0, vdim * sizeof(T));
    }

    /** Initialize to zeros following a pattern.
     *
     * All values are set to 0, and the origin and size of the pattern are
     * adopted.
     *
     * @code
     * v2.initZeros(v1);
     * @endcode
     */
    template<typename T1>
    void initZeros(const Matrix1D<T1>& op)
    {
        if (vdim != op.vdim)
            resize(op);
        memset(vdata, 0, vdim * sizeof(T));
    }
    //@}

    /** Check if any element of the array is NaN.
     */
    bool isAnyNaN()
    {
        for (size_t j = 0; j < vdim; j++)
            if (ISNAN(vdata[j]))
                return true;
        return false;
    }

    /// @name Matrix1D operators
    //@{
    /** v3 = v1 * k.
     */
    Matrix1D<T> operator*(T op1) const
    {
        Matrix1D<T> tmp(*this);
        T *ptr1 = &VEC_ELEM(tmp,0);
        const T *ptr2 = &VEC_ELEM(*this,0);
        size_t iBlockMax = vdim / 4;
        for (size_t i = 0; i < iBlockMax; i++)
        {
            (*ptr1++) = (*ptr2++) * op1;
            (*ptr1++) = (*ptr2++) * op1;
            (*ptr1++) = (*ptr2++) * op1;
            (*ptr1++) = (*ptr2++) * op1;
        }
        for (size_t i = iBlockMax * 4; i < vdim; ++i)
            (*ptr1++) = (*ptr2++) * op1;
        return tmp;
    }

    /** v3 = v1 / k.
     */
    Matrix1D<T> operator/(T op1) const
    {
        Matrix1D<T> tmp(*this);
        T iop1 = 1 / op1;
        T *ptr1 = &VEC_ELEM(tmp,0);
        const T *ptr2 = &VEC_ELEM(*this,0);
        size_t iBlockMax = vdim / 4;
        for (size_t i = 0; i < iBlockMax; i++)
        {
            (*ptr1++) = (*ptr2++) * iop1;
            (*ptr1++) = (*ptr2++) * iop1;
            (*ptr1++) = (*ptr2++) * iop1;
            (*ptr1++) = (*ptr2++) * iop1;
        }
        for (size_t i = iBlockMax * 4; i < vdim; ++i)
            (*ptr1++) = (*ptr2++) * iop1;
        return tmp;
    }

    /** v3 = v1 + k.
     */
    Matrix1D<T> operator+(T op1) const
    {
        Matrix1D<T> tmp(*this);
        T *ptr1 = &VEC_ELEM(tmp,0);
        const T *ptr2 = &VEC_ELEM(*this,0);
        size_t iBlockMax = vdim / 4;
        for (size_t i = 0; i < iBlockMax; i++)
        {
            (*ptr1++) = (*ptr2++) + op1;
            (*ptr1++) = (*ptr2++) + op1;
            (*ptr1++) = (*ptr2++) + op1;
            (*ptr1++) = (*ptr2++) + op1;
        }
        for (size_t i = iBlockMax * 4; i < vdim; ++i)
            (*ptr1++) = (*ptr2++) + op1;
        return tmp;
    }

    /** v3 = v1 - k.
     */
    Matrix1D<T> operator-(T op1) const
    {
        Matrix1D<T> tmp(*this);
        T *ptr1 = &VEC_ELEM(tmp,0);
        const T *ptr2 = &VEC_ELEM(*this,0);
        size_t iBlockMax = vdim / 4;
        for (size_t i = 0; i < iBlockMax; i++)
        {
            (*ptr1++) = (*ptr2++) - op1;
            (*ptr1++) = (*ptr2++) - op1;
            (*ptr1++) = (*ptr2++) - op1;
            (*ptr1++) = (*ptr2++) - op1;
        }
        for (size_t i = iBlockMax * 4; i < vdim; ++i)
            (*ptr1++) = (*ptr2++) - op1;
        return tmp;
    }

    /** v3 = k * v2.
     */
    friend Matrix1D<T> operator*(T op1, const Matrix1D<T>& op2)
    {
        Matrix1D<T> tmp(op2);
        T *ptr1 = &VEC_ELEM(tmp,0);
        const T *ptr2 = &VEC_ELEM(op2,0);
        size_t iBlockMax = op2.vdim / 4;
        for (size_t i = 0; i < iBlockMax; i++)
        {
            (*ptr1++) = (*ptr2++) * op1;
            (*ptr1++) = (*ptr2++) * op1;
            (*ptr1++) = (*ptr2++) * op1;
            (*ptr1++) = (*ptr2++) * op1;
        }
        for (size_t i = iBlockMax * 4; i < op2.vdim; ++i)
            (*ptr1++) = (*ptr2++) * op1;
        return tmp;
    }

    /** v3 = k / v2.
     */
    friend Matrix1D<T> operator/(T op1, const Matrix1D<T>& op2)
    {
        Matrix1D<T> tmp(op2);
        T *ptr1 = &VEC_ELEM(tmp,0);
        const T *ptr2 = &VEC_ELEM(op2,0);
        size_t iBlockMax = op2.vdim / 4;
        for (size_t i = 0; i < iBlockMax; i++)
        {
            (*ptr1++) = op1 / (*ptr2++);
            (*ptr1++) = op1 / (*ptr2++);
            (*ptr1++) = op1 / (*ptr2++);
            (*ptr1++) = op1 / (*ptr2++);
        }
        for (size_t i = iBlockMax * 4; i < op2.vdim; ++i)
            (*ptr1++) = op1 / (*ptr2++);
        return tmp;
    }

    /** v3 = k + v2.
     */
    friend Matrix1D<T> operator+(T op1, const Matrix1D<T>& op2)
    {
        Matrix1D<T> tmp(op2);
        T *ptr1 = &VEC_ELEM(tmp,0);
        const T *ptr2 = &VEC_ELEM(op2,0);
        size_t iBlockMax = op2.vdim / 4;
        for (size_t i = 0; i < iBlockMax; i++)
        {
            (*ptr1++) = (*ptr2++) + op1;
            (*ptr1++) = (*ptr2++) + op1;
            (*ptr1++) = (*ptr2++) + op1;
            (*ptr1++) = (*ptr2++) + op1;
        }
        for (size_t i = iBlockMax * 4; i < op2.vdim; ++i)
            (*ptr1++) = (*ptr2++) + op1;
        return tmp;
    }

    /** v3 = k - v2.
     */
    friend Matrix1D<T> operator-(T op1, const Matrix1D<T>& op2)
    {
        Matrix1D<T> tmp(op2);
        T *ptr1 = &VEC_ELEM(tmp,0);
        const T *ptr2 = &VEC_ELEM(op2,0);
        size_t iBlockMax = op2.vdim / 4;
        for (size_t i = 0; i < iBlockMax; i++)
        {
            (*ptr1++) = op1 - (*ptr2++);
            (*ptr1++) = op1 - (*ptr2++);
            (*ptr1++) = op1 - (*ptr2++);
            (*ptr1++) = op1 - (*ptr2++);
        }
        for (size_t i = iBlockMax * 4; i < op2.vdim; ++i)
            (*ptr1++) = op1 - (*ptr2++);
        return tmp;
    }

    /** Vector summation
     *
     * @code
     * A += B;
     * @endcode
     */
    void operator+=(const Matrix1D<T>& op1) const
    {
        if (vdim != op1.vdim)
            REPORT_ERROR(ERR_MATRIX_SIZE, "Not same sizes in vector summation");

        T *ptr1 = &VEC_ELEM(*this,0);
        const T *ptr2 = &VEC_ELEM(op1,0);
        size_t iBlockMax = vdim / 4;
        for (size_t i = 0; i < iBlockMax; i++)
        {
            (*ptr1++) += (*ptr2++);
            (*ptr1++) += (*ptr2++);
            (*ptr1++) += (*ptr2++);
            (*ptr1++) += (*ptr2++);
        }
        for (size_t i = iBlockMax * 4; i < vdim; ++i)
            (*ptr1++) += (*ptr2++);
    }

    /** Vector substraction
     *
     * @code
     * A -= B;
     * @endcode
     */
    void operator-=(const Matrix1D<T>& op1) const
    {
        if (vdim != op1.vdim)
            REPORT_ERROR(ERR_MATRIX_SIZE, "Not same sizes in vector summation");

        T *ptr1 = &VEC_ELEM(*this,0);
        const T *ptr2 = &VEC_ELEM(op1,0);
        size_t iBlockMax = vdim / 4;
        for (size_t i = 0; i < iBlockMax; i++)
        {
            (*ptr1++) -= (*ptr2++);
            (*ptr1++) -= (*ptr2++);
            (*ptr1++) -= (*ptr2++);
            (*ptr1++) -= (*ptr2++);
        }
        for (size_t i = iBlockMax * 4; i < vdim; ++i)
            (*ptr1++) -= (*ptr2++);
    }

    /** v3 *= k.
     */
    void operator*=(T op1)
    {
        T *ptr1 = &VEC_ELEM(*this,0);
        size_t iBlockMax = vdim / 4;
        for (size_t i = 0; i < iBlockMax; i++)
        {
            (*ptr1++) *= op1;
            (*ptr1++) *= op1;
            (*ptr1++) *= op1;
            (*ptr1++) *= op1;
        }
        for (size_t i = iBlockMax * 4; i < vdim; ++i)
            (*ptr1++) *= op1;
    }

    /** v3 /= k.
     */
    void operator/=(T op1)
    {
        T iop1 = 1 / op1;
        T * ptr1 = &VEC_ELEM(*this,0);
        size_t iBlockMax = vdim / 4;
        for (size_t i = 0; i < iBlockMax; i++)
        {
            (*ptr1++) *= iop1;
            (*ptr1++) *= iop1;
            (*ptr1++) *= iop1;
            (*ptr1++) *= iop1;
        }
        for (size_t i = iBlockMax * 4; i < vdim; ++i)
            (*ptr1++) *= iop1;
    }

    /** v3 += k.
     */
    void operator+=(T op1)
    {
        T *ptr1 = &VEC_ELEM(*this,0);
        size_t iBlockMax = vdim / 4;
        for (size_t i = 0; i < iBlockMax; i++)
        {
            (*ptr1++) += op1;
            (*ptr1++) += op1;
            (*ptr1++) += op1;
            (*ptr1++) += op1;
        }
        for (size_t i = iBlockMax * 4; i < vdim; ++i)
            (*ptr1++) += op1;
    }

    /** v3 -= k.
     */
    void operator-=(T op1)
    {
        T *ptr1 = &VEC_ELEM(*this,0);
        size_t iBlockMax = vdim / 4;
        for (size_t i = 0; i < iBlockMax; i++)
        {
            (*ptr1++) -= op1;
            (*ptr1++) -= op1;
            (*ptr1++) -= op1;
            (*ptr1++) -= op1;
        }
        for (size_t i = iBlockMax * 4; i < vdim; ++i)
            (*ptr1++) -= op1;
    }

    /** v3 = v1 * v2.
     */
    Matrix1D<T> operator*(const Matrix1D<T>& op1) const
    {
        Matrix1D<T> tmp(op1);
        T *ptr1 = &VEC_ELEM(tmp,0);
        const T *ptr2 = &VEC_ELEM(*this,0);
        const T *ptr3 = &VEC_ELEM(op1,0);
        size_t iBlockMax = vdim / 4;
        for (size_t i = 0; i < iBlockMax; i++)
        {
            (*ptr1++) = (*ptr2++) * (*ptr3++);
            (*ptr1++) = (*ptr2++) * (*ptr3++);
            (*ptr1++) = (*ptr2++) * (*ptr3++);
            (*ptr1++) = (*ptr2++) * (*ptr3++);
        }
        for (size_t i = iBlockMax * 4; i < vdim; ++i)
            (*ptr1++) = (*ptr2++) * (*ptr3++);
        return tmp;
    }

    /** v3 = v1 / v2.
     */
    Matrix1D<T> operator/(const Matrix1D<T>& op1) const
    {
        Matrix1D<T> tmp(op1);
        T *ptr1 = &VEC_ELEM(tmp,0);
        const T *ptr2 = &VEC_ELEM(*this,0);
        const T *ptr3 = &VEC_ELEM(op1,0);
        size_t iBlockMax = vdim / 4;
        for (size_t i = 0; i < iBlockMax; i++)
        {
            (*ptr1++) = (*ptr2++) / (*ptr3++);
            (*ptr1++) = (*ptr2++) / (*ptr3++);
            (*ptr1++) = (*ptr2++) / (*ptr3++);
            (*ptr1++) = (*ptr2++) / (*ptr3++);
        }
        for (size_t i = iBlockMax * 4; i < vdim; ++i)
            (*ptr1++) = (*ptr2++) / (*ptr3++);
        return tmp;
    }
    /** v3 = v1 + v2.
     */
    Matrix1D<T> operator+(const Matrix1D<T>& op1) const
    {
        Matrix1D<T> tmp(op1);
        T *ptr1 = &VEC_ELEM(tmp,0);
        const T *ptr2 = &VEC_ELEM(*this,0);
        const T *ptr3 = &VEC_ELEM(op1,0);
        size_t iBlockMax = vdim / 4;
        for (size_t i = 0; i < iBlockMax; i++)
        {
            (*ptr1++) = (*ptr2++) + (*ptr3++);
            (*ptr1++) = (*ptr2++) + (*ptr3++);
            (*ptr1++) = (*ptr2++) + (*ptr3++);
            (*ptr1++) = (*ptr2++) + (*ptr3++);
        }
        for (size_t i = iBlockMax * 4; i < vdim; ++i)
            (*ptr1++) = (*ptr2++) + (*ptr3++);
        return tmp;
    }

    /** v3 = v1 - v2.
     */
    Matrix1D<T> operator-(const Matrix1D<T>& op1) const
    {
        Matrix1D<T> tmp(op1);
        T *ptr1 = &VEC_ELEM(tmp,0);
        const T *ptr2 = &VEC_ELEM(*this,0);
        const T *ptr3 = &VEC_ELEM(op1,0);
        size_t iBlockMax = vdim / 4;
        for (size_t i = 0; i < iBlockMax; i++)
        {
            (*ptr1++) = (*ptr2++) - (*ptr3++);
            (*ptr1++) = (*ptr2++) - (*ptr3++);
            (*ptr1++) = (*ptr2++) - (*ptr3++);
            (*ptr1++) = (*ptr2++) - (*ptr3++);
        }
        for (size_t i = iBlockMax * 4; i < vdim; ++i)
            (*ptr1++) = (*ptr2++) - (*ptr3++);
        return tmp;
    }

    /** v3 *= v2.
     */
    void operator*=(const Matrix1D<T>& op1)
    {
        T *ptr1 = &VEC_ELEM(*this,0);
        const T *ptr2 = &VEC_ELEM(op1,0);
        size_t iBlockMax = vdim / 4;
        for (size_t i = 0; i < iBlockMax; i++)
        {
            (*ptr1++) *= (*ptr2++);
            (*ptr1++) *= (*ptr2++);
            (*ptr1++) *= (*ptr2++);
            (*ptr1++) *= (*ptr2++);
        }
        for (size_t i = iBlockMax * 4; i < vdim; ++i)
            (*ptr1++) *= (*ptr2++);
    }

    /** v3 /= v2.
     */
    void operator/=(const Matrix1D<T>& op1)
    {
        T *ptr1 = &VEC_ELEM(*this,0);
        const T *ptr2 = &VEC_ELEM(op1,0);
        size_t iBlockMax = vdim / 4;
        for (size_t i = 0; i < iBlockMax; i++)
        {
            (*ptr1++) /= (*ptr2++);
            (*ptr1++) /= (*ptr2++);
            (*ptr1++) /= (*ptr2++);
            (*ptr1++) /= (*ptr2++);
        }
        for (size_t i = iBlockMax * 4; i < vdim; ++i)
            (*ptr1++) /= (*ptr2++);
    }

    /** v3 += v2.
     */
    void operator+=(const Matrix1D<T>& op1)
    {
        T *ptr1 = &VEC_ELEM(*this,0);
        const T *ptr2 = &VEC_ELEM(op1,0);
        size_t iBlockMax = vdim / 4;
        for (size_t i = 0; i < iBlockMax; i++)
        {
            (*ptr1++) += (*ptr2++);
            (*ptr1++) += (*ptr2++);
            (*ptr1++) += (*ptr2++);
            (*ptr1++) += (*ptr2++);
        }
        for (size_t i = iBlockMax * 4; i < vdim; ++i)
            (*ptr1++) += (*ptr2++);
    }

    /** v3 -= v2.
     */
    void operator-=(const Matrix1D<T>& op1)
    {
        T *ptr1 = &VEC_ELEM(*this,0);
        const T *ptr2 = &VEC_ELEM(op1,0);
        size_t iBlockMax = vdim / 4;
        for (size_t i = 0; i < iBlockMax; i++)
        {
            (*ptr1++) -= (*ptr2++);
            (*ptr1++) -= (*ptr2++);
            (*ptr1++) -= (*ptr2++);
            (*ptr1++) -= (*ptr2++);
        }
        for (size_t i = iBlockMax * 4; i < vdim; ++i)
            (*ptr1++) -= (*ptr2++);
    }

    /** Unary minus.
     *
     * It is used to build arithmetic expressions. You can make a minus
     * of anything as long as it is correct semantically.
     *
     * @code
     * v1 = -v2;
     * v1 = -v2.transpose();
     * @endcode
     */
    Matrix1D<T> operator-() const
    {
        Matrix1D<T> tmp(*this);
        T *ptr1 = &VEC_ELEM(tmp,0);
        const T *ptr2 = &VEC_ELEM(*this,0);
        size_t iBlockMax = vdim / 4;
        for (size_t i = 0; i < iBlockMax; i++)
        {
            (*ptr1++) = -(*ptr2++);
            (*ptr1++) = -(*ptr2++);
            (*ptr1++) = -(*ptr2++);
            (*ptr1++) = -(*ptr2++);
        }
        for (size_t i = iBlockMax * 4; i < vdim; ++i)
            (*ptr1++) = -(*ptr2++);
        return tmp;
    }

    /** Vector by matrix
     *
     * Algebraic vector by matrix multiplication. This function is actually
     * implemented in xmippMatrices2D
     */
    Matrix1D<T> operator*(const Matrix2D<T>& M);

    /** Vector element access
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
    T& operator()(int i) const
    {
        return vdata[i];
    }
    //@}

    /// @name Utilities for Matrix1D
    //@{

    /** Produce a vector suitable for working with Numerical Recipes
     *
     * This function must be used only as a preparation for routines which need
     * that the first physical index is 1 and not 0 as it usually is in C. In
     * fact the vector provided for Numerical recipes is exactly this same one
     * but with the indexes changed.
     *
     * This function is not ported to Python.
     */
    T* adaptForNumericalRecipes() const
    {
        return MATRIX1D_ARRAY(*this) - 1;
    }

    /** Kill an array produced for Numerical Recipes.
     *
     * Nothing needs to be done in fact.
     *
     * This function is not ported to Python.
     */
    void killAdaptationForNumericalRecipes(T* m) const
        {}

    /** CEILING
     *
     * Applies a CEILING (look for the nearest larger integer) to each
     * array element.
     */
    void selfCEIL()
    {
        T *ptr1 = &VEC_ELEM(*this,0);
        size_t iBlockMax = vdim / 4;
        for (size_t i = 0; i < iBlockMax; i++)
        {
            *ptr1 = ceil(*ptr1);
            ++ptr1;
            *ptr1 = ceil(*ptr1);
            ++ptr1;
            *ptr1 = ceil(*ptr1);
            ++ptr1;
            *ptr1 = ceil(*ptr1);
            ++ptr1;
        }
        for (size_t i = iBlockMax * 4; i < vdim; ++i)
        {
            *ptr1 = ceil(*ptr1);
            ++ptr1;
        }
    }

    /** FLOOR
     *
     * Applies a FLOOR (look for the nearest larger integer) to each
     * array element.
     */
    void selfFLOOR()
    {
        T *ptr1 = &VEC_ELEM(*this,0);
        size_t iBlockMax = vdim / 4;
        for (size_t i = 0; i < iBlockMax; i++)
        {
            *ptr1 = floor(*ptr1);
            ++ptr1;
            *ptr1 = floor(*ptr1);
            ++ptr1;
            *ptr1 = floor(*ptr1);
            ++ptr1;
            *ptr1 = floor(*ptr1);
            ++ptr1;
        }
        for (size_t i = iBlockMax * 4; i < vdim; ++i)
        {
            *ptr1 = floor(*ptr1);
            ++ptr1;
        }
    }

    /** ROUND
     *
     * Applies a ROUND (look for the nearest larger integer) to each
     * array element.
     */
    void selfROUND()
    {
        T *ptr1 = &VEC_ELEM(*this,0);
        size_t iBlockMax = vdim / 4;
        for (size_t i = 0; i < iBlockMax; i++)
        {
            *ptr1 = round(*ptr1);
            ++ptr1;
            *ptr1 = round(*ptr1);
            ++ptr1;
            *ptr1 = round(*ptr1);
            ++ptr1;
            *ptr1 = round(*ptr1);
            ++ptr1;
        }
        for (size_t i = iBlockMax * 4; i < vdim; ++i)
        {
            *ptr1 = round(*ptr1);
            ++ptr1;
        }
    }

    /** Sort 1D vector elements
     *
     * Sort in ascending order the vector elements. You can use the "reverse"
     * function to sort in descending order.
     *
     * @code
     * v2 = v1.sort();
     * @endcode
     */
    Matrix1D<T> sort() const
    {
        Matrix1D<T> temp;
        Matrix1D<double> aux;

        if (vdim == 0)
            return temp;

        // Initialise vdata
        typeCast(*this, aux);

        // Sort
        double * aux_array = aux.adaptForNumericalRecipes();
        qcksrt(vdim, aux_array);

        typeCast(aux, temp);
        return temp;
    }

    /** Maximum element */
    T computeMax() const
    {
        if (vdim == 0)
            return 0;

        T maxval = VEC_ELEM(*this,0);
        for (size_t j = 0; j < vdim; ++j)
            if (VEC_ELEM(*this,j) > maxval)
                maxval = VEC_ELEM(*this,j);
        return maxval;
    }

    /** Index for the maximum element.
     *
     * This function returns the index of the maximum element of an matrix1d.
     * Returns -1 if the array is empty
     */
    void maxIndex(int& jmax) const
    {
        if (vdim == 0)
        {
            jmax = -1;
            return;
        }

        jmax = 0;
        T maxval = VEC_ELEM(*this, 0);
        for (size_t j = 0; j < vdim; ++j)
            if (VEC_ELEM(*this,j) > maxval)
            {
                jmax = j;
                maxval = VEC_ELEM(*this,j);
            }
    }

    /** Index for the minimum element.
     *
     * This function returns the index of the minimum element of an matrix1d.
     * Returns -1 if the array is empty
     */
    void minIndex(int& jmin) const
    {
        if (vdim == 0)
        {
            jmin = -1;
            return;
        }

        jmin = 0;
        T minval = VEC_ELEM(*this, 0);
        for (size_t j = 0; j < vdim; ++j)
            if (VEC_ELEM(*this,j) < minval)
            {
                jmin = j;
                minval = VEC_ELEM(*this,j);
            }
    }

    /** Algebraic transpose of vector
     *
     * You can use the transpose in as complex expressions as you like. The
     * origin of the vector is not changed.
     *
     * @code
     * v2 = v1.transpose();
     * @endcode
     */
    Matrix1D<T> transpose() const
    {
        Matrix1D<T> temp(*this);
        temp.selfTranspose();
        return temp;
    }

    /** Algebraic transpose of vector
     *
     * The same as before but the result is stored in this same object.
     */
    inline void selfTranspose()
    {
        row = !row;
    }

    /** Sum of vector values.
     *
     * This function returns the sum of all internal values.
     *
     * @code
     * double sum = m.sum();
     * @endcode
     */
    double sum(bool average = false) const
    {
        double sum = 0;
        const T *ptr1 = &VEC_ELEM(*this,0);
        size_t iBlockMax = vdim / 4;
        for (size_t i = 0; i < iBlockMax; i++)
        {
            sum += (*ptr1++);
            sum += (*ptr1++);
            sum += (*ptr1++);
            sum += (*ptr1++);
        }
        for (size_t i = iBlockMax * 4; i < vdim; ++i)
            sum += (*ptr1++);
        if (average)
            return sum / (double) vdim;
        else
            return sum;
    }

    /** Sum of squared vector values.
     *
     * This function returns the sum of all internal values to the second
     * power.
     *
     * @code
     * double sum2 = m.sum2();
     * @endcode
     */
    double sum2() const
    {
        double sum = 0;
        const T *ptr1 = &VEC_ELEM(*this,0);
        double val;
        size_t iBlockMax = vdim / 4;
        for (size_t i = 0; i < iBlockMax; i++)
        {
            val = *ptr1;
            sum += val * val;
            ++ptr1;
            val = *ptr1;
            sum += val * val;
            ++ptr1;
            val = *ptr1;
            sum += val * val;
            ++ptr1;
            val = *ptr1;
            sum += val * val;
            ++ptr1;
        }
        for (size_t i = iBlockMax * 4; i < vdim; ++i)
        {
            val = *ptr1;
            sum += val * val;
            ++ptr1;
        }
        return sum;
    }

    /** Dot product <*this,op1>.
     */
    double dotProduct(const Matrix1D<T> &op1) const
    {
        double sum = 0;
        const T *ptr1 = &VEC_ELEM(*this,0);
        const T *ptr2 = &VEC_ELEM(op1,0);
        size_t iBlockMax = vdim / 4;
        for (size_t i = 0; i < iBlockMax; i++)
        {
            sum += (*ptr1++) * (*ptr2++);
            sum += (*ptr1++) * (*ptr2++);
            sum += (*ptr1++) * (*ptr2++);
            sum += (*ptr1++) * (*ptr2++);
        }
        for (size_t i = iBlockMax * 4; i < vdim; ++i)
            sum += (*ptr1++) * (*ptr2++);
        return sum;
    }

    /** Module of the vector
     *
     * This module is defined as the square root of the sum of the squared
     * components. Euclidean norm of the vector.
     *
     * @code
     * double mod = v.module();
     * @endcode
     */
    inline double module() const
    {
        return sqrt(sum2());
    }

    /** Angle of the vector
     *
     * Supposing this vector is in R2 this function returns the angle of this
     * vector with X axis, ie, atan2(YY(v), XX(v))
     */
    double angle()
    {
        return atan2((double) YY(*this), (double) XX(*this));
    }

    /** Normalize this vector, store the result here
     */
    void selfNormalize()
    {
        double m = module();
        if (fabs(m) > XMIPP_EQUAL_ACCURACY)
        {
            T im = (T) (1.0 / m);
            *this *= im;
        }
        else
            initZeros();
    }

    /** Reverse vector values, keep in this object.
     */
    void selfReverse()
    {
        for (int j = 0; j <= (int) (vdim - 1) / 2; j++)
        {
            T aux;
            SWAP(vdata[j], vdata[vdim-1-j], aux);
        }
    }

    /** Show using gnuplot
     *
     * This function uses gnuplot to plot this vector. You must supply the
     * xlabel and title.
     */
    void showWithGnuPlot(const std::string& xlabel, const std::string& title)
    {
        FileName fn_tmp;
        fn_tmp.initRandom(10);
        Matrix1D<T>::write(static_cast<std::string>("PPP") + fn_tmp + ".txt");

        std::ofstream fh_gplot;
        fh_gplot.open(
            (static_cast<std::string>("PPP") + fn_tmp + ".gpl").c_str());
        if (!fh_gplot)
            REPORT_ERROR(
                ERR_UNCLASSIFIED,
                static_cast<std::string>("vector::showWithGnuPlot: Cannot open PPP") + fn_tmp + ".gpl for output");
        fh_gplot << "set xlabel \"" + xlabel + "\"\n";
        fh_gplot
        << "plot \"PPP" + fn_tmp + ".txt\" title \"" + title
        + "\" w l\n";
        fh_gplot << "pause 300 \"\"\n";
        fh_gplot.close();
        system(
            (static_cast<std::string>("(gnuplot PPP") + fn_tmp
             + ".gpl; rm PPP" + fn_tmp + ".txt PPP" + fn_tmp
             + ".gpl) &").c_str());
    }

    /** Compute numerical derivative
     *
     * The numerical derivative is of the same size as the input vector.
     * However, the first two and the last two samples are set to 0,
     * because the numerical method is not able to correctly estimate the
     * derivative there.
     */
    void numericalDerivative(Matrix1D<double> &result)
    {
        const double i12 = 1.0 / 12.0;
        result.initZeros(*this);
        for (int i = 2; i <= vdim - 2; i++)
            result(i) = i12
                        * (-(*this)(i + 2) + 8 * (*this)(i + 1) - 8 * (*this)(i - 1)
                           + (*this)(i + 2));
    }

    /** Read from an ASCII file.
     *
     * The array must be previously resized to the correct size.
     */
    void read(const FileName& fn)
    {
        std::ifstream in;
        in.open(fn.c_str(), std::ios::in);

        if (!in)
            REPORT_ERROR(
                ERR_IO_NOTOPEN,
                static_cast< std::string >("MultidimArray::read: File " + fn + " not found"));

        in >> *this;
        in.close();
    }

    /** Input from input stream.
     *
     * Actual size of the array is used to know how many values must be read.
     *
     * @code
     * v.resize(3);
     * std::cin >> v;
     * @endcode
     *
     * This function is not ported to Python.
     */
    // This function must be explictly implemented outside
    friend std::istream& operator>>(std::istream& in, Matrix1D<T>& v)
    {
        for (int i = 0; i < VEC_XSIZE(v); i++)
            in >> VEC_ELEM(v,i);
        return in;
    }

    /** Output to output stream.*/
    friend std::ostream& operator<<(std::ostream& ostrm, const Matrix1D<T>& v)
    {
        if (v.vdim == 0)
            ostrm << "NULL Array\n";
        else
            ostrm << std::endl;

        double max_val = ABS(v.vdata[0]);

        for (size_t j = 0; j < v.vdim; j++)
        {
            max_val = XMIPP_MAX(max_val, v.vdata[j]);
        }

        int prec = bestPrecision(max_val, 10);

        if (v.row)
            for (size_t j = 0; j < v.vdim; j++)
            {
                ostrm << floatToString((double) v.vdata[j], 10, prec)
                << std::endl;
            }
        else
            for (size_t j = 0; j < v.vdim; j++)
            {
                ostrm << floatToString((double) v.vdata[j], 10, prec) << " ";
            }

        return ostrm;
    }

    /** Write to an ASCII file.
     */
    void write(const FileName& fn) const
    {
        std::ofstream out;
        out.open(fn.c_str(), std::ios::out);
        if (!out)
            REPORT_ERROR(
                ERR_IO_NOTOPEN,
                static_cast< std::string >("Matrix1D::write: File " + fn + " cannot be opened for output"));

        out << *this;
        out.close();
    }

    /** Edit with xmipp_editor.
     *
     * This function generates a random filename starting with PPP and
     * edits it with xmipp_editor. After closing the editor the file is
     * removed.
     */
    void edit()
    {
        FileName nam;
        nam.initRandom(15);

        nam = static_cast<std::string>("PPP" + nam + ".txt");
        write
        (nam);

        system(
            (static_cast<std::string>("xmipp_edit -i " + nam + " -remove &").c_str()));
    }
    //@}
};


 typedef Matrix1D<double> DVector;
 typedef Matrix1D<int> IVector;

/**@name Vector Related functions
 * These functions are not methods of Matrix1D
 */
//@{
/** Creates vector in R2.
 * After this function the vector is (x,y) in R2.
 *
 * @code
 * Matrix1D< double > v = vectorR2(1, 2);
 * @endcode
 */
Matrix1D<double> vectorR2(double x, double y);

/** Creates vector in R3.
 * After this function the vector is (x,y,z) in R3.
 *
 * @code
 * Matrix1D< double > v = vectorR2(1, 2, 1);
 * @endcode
 */
Matrix1D<double> vectorR3(double x, double y, double z);

/** Creates an integer vector in Z3.
 */
Matrix1D<int> vectorR3(int x, int y, int z);

/** Dot product.
 * Given any two vectors in Rn (n-dimensional vector), this function returns the
 * dot product of both. If the vectors are not of the same size or shape then an
 * exception is thrown. The dot product is defined as the sum of the component
 * by component multiplication.
 *
 * For the R3 vectors (V1x,V1y,V1z), (V2x, V2y, V2z) the result is V1x*V2x +
 * V1y*V2y + V1z*V2z.
 *
 * @code
 * Matrix1D< double > v1(1000);
 * v1.init_random(0, 10, "gaussian");
 * std::cout << "The power of this vector should be 100 and is " <<
 *     dotProduct(v1, v1) << std::endl;
 * @endcode
 */
template<typename T>
T dotProduct(const Matrix1D<T>& v1, const Matrix1D<T>& v2)
{
    if (!v1.sameShape(v2))
        REPORT_ERROR(ERR_MATRIX_SIZE,
                     "Dot product: vectors of different size or shape");

    T accumulate = 0;
    for (size_t j = 0; j < v1.vdim; j++)
        accumulate += v1.vdata[j] * v2.vdata[j];
    return accumulate;
}

/** Vector product in R3.
 * This function takes two R3 vectors and compute their vectorial product. For
 * two vectors (V1x,V1y,V1z), (V2x, V2y, V2z) the result is (V1y*V2z-V1z*v2y,
 * V1z*V2x-V1x*V2z, V1x*V2y-V1y*V2x). Pay attention that this operator is not
 * conmutative. An exception is thrown if the vectors are not of the same shape
 * or they don't belong to R3.
 *
 * @code
 * Matrix1D< T > X = vectorR3(1, 0, 0), Y = vector_R3(0, 1, 0);
 * std::cout << "X*Y=Z=" << vectorProduct(X,Y).transpose() << std::endl;
 * @endcode
 */
template<typename T>
Matrix1D<T> vectorProduct(const Matrix1D<T>& v1, const Matrix1D<T>& v2)
{
    if (v1.vdim != 3 || v2.vdim != 3)
        REPORT_ERROR(ERR_MATRIX_SIZE, "Vector_product: vectors are not in R3");

    if (v1.isRow() != v2.isRow())
        REPORT_ERROR(ERR_MATRIX_SIZE,
                     "Vector_product: vectors are of different shape");

    Matrix1D<T> result(3);
    vectorProduct(v1, v2, result);
    return result;
}

/** Vector product in R3.
 * This function computes the vector product of two R3 vectors.
 * No check is performed, it is assumed that the output vector
 * is already resized
 *
 */
template<typename T>
void vectorProduct(const Matrix1D<T>& v1, const Matrix1D<T>& v2,
                   Matrix1D<T> &result)
{
    XX(result) = YY(v1) * ZZ(v2) - ZZ(v1) * YY(v2);
    YY(result) = ZZ(v1) * XX(v2) - XX(v1) * ZZ(v2);
    ZZ(result) = XX(v1) * YY(v2) - YY(v1) * XX(v2);
}

/** Sort two vectors.
 * v1 and v2 must be of the same shape, if not an exception is thrown. After
 * calling this function all components in v1 are the minimum between the
 * corresponding components in v1 and v2, and all components in v2 are the
 * maximum.
 *
 * For instance, XX(v1)=MIN(XX(v1), XX(v2)), XX(v2)=MAX(XX(v1), XX(v2)). Notice
 * that both vectors are modified. This function is very useful for sorting two
 * corners. After calling it you can certainly perform a non-empty for (from
 * corner1 to corner2) loop.
 */
template<typename T>
void sortTwoVectors(Matrix1D<T>& v1, Matrix1D<T>& v2)
{
    T temp;
    if (!v1.sameShape(v2))
        REPORT_ERROR(ERR_MATRIX_SIZE,
                     "sortTwoVectors: vectors are not of the same shape");

    for (size_t j = 0; j < v1.vdim; j++)
    {
        temp = XMIPP_MIN(v1.vdata[j], v2.vdata[j]);
        v2.vdata[j] = XMIPP_MAX(v1.vdata[j], v2.vdata[j]);
        v1.vdata[j] = temp;
    }
}

/** Conversion from one type to another.
 * If we have an integer array and we need a double one, we can use this
 * function. The conversion is done through a type casting of each element
 * If n >= 0, only the nth volumes will be converted, otherwise all NSIZE volumes
 */
template<typename T1, typename T2>
void typeCast(const Matrix1D<T1>& v1, Matrix1D<T2>& v2)
{
    if (v1.vdim == 0)
    {
        v2.clear();
        return;
    }

    v2.resizeNoCopy(v1.vdim);
    for (size_t j = 0; j < v1.vdim; j++)
        v2.vdata[j] = static_cast<T2>(v1.vdata[j]);
}

/** Conversion from one type to another.
 * In some cases, the two types are the same. So a faster way is simply by assignment.
 */
template<typename T1>
void typeCast(const Matrix1D<T1>& v1, Matrix1D<T1>& v2)
{
    v2 = v1;
}

//@}
//@}
#endif /* MATRIX1D_H_ */
