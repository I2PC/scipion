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

#ifndef MATRIX_H
#define MATRIX_H

#include <fstream>
#include <complex>
#include <cmath>

#include "multidimensional_array.h"
#include "error.h"

template <typename T> class Matrix2D;

/// @defgroup Vectors Vectors
/// @ingroup MultidimensionalArrays

/** @defgroup VectorsSpeedUp Speed up macros
 * @ingroup Vectors
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

/** For all elements in the array
 * @ingroup VectorsSizeShape
 *
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
    for (int i=STARTINGX(v); i<=FINISHINGX(v); i++)

/** For all elements in the array between corners
 * @ingroup VectorsSizeShape
 *
 * This macro is used to generate loops for a vector in an easy manner. It needs
 * an externally defined Matrix1D< double > r(1). Then XX(r) ranges from
 * (int) XX(corner1) to (int) XX(corner2) (included limits) (notice that corner1
 * and corner2 need only to be Matrix1D).
 *
 * @code
 * Matrix1D< double > corner1(1), corner2(1), r(1);
 * XX(corner1) = -1;
 * XX(corner2) = 1;
 * FOR_ALL_ELEMENTS_IN_MATRIX1D_BETWEEN(corner1, corner2)
 * {
 *     std::cout << v(XX(r)) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_ELEMENTS_IN_MATRIX1D_BETWEEN(corner1, corner2) \
    for (XX(r)=(int) XX((corner1)); XX(r)<=(int) XX((corner2)); XX(r)++)

/** For all elements in common
 * @ingroup VectorsSizeShape
 *
 * This macro is used to generate loops for all the elements logically in common
 * between two vectors in an easy manner. Then i (locally defined) ranges from
 * MAX(STARTINGX(V1), STARTINGX(V2)) to MIN(FINISHINGX(V1), FINISHINGX(V2))
 * (included limits) respectively. You need to define SPEED_UP_temps.
 *
 * @code
 * Matrix2D< double > v1(10), v2(20);
 * v1.setXmippOrigin();
 * v2.setXmippOrigin();
 *
 * FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX1D(v1, v2)
 * {
 *     ...
 * }
 * @endcode
 */
#define FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX1D(v1, v2) \
    ispduptmp4 = XMIPP_MAX(STARTINGX(v1), STARTINGX(v2)); \
    ispduptmp5 = XMIPP_MIN(FINISHINGX(v1), FINISHINGX(v2)); \
    for (int i=ispduptmp4; i<=ispduptmp5; i++)

/** For all elements in the array, accessed physically
 * @ingroup VectorsSizeShape
 *
 * This macro is used to generate loops for the vector in an easy way using
 * physical indexes. It defines internal the index 'i' which ranges the vector
 * using its physical definition.
 *
 * @code
 * FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX1D(v)
 * {
 *     std::cout << DIRECT_MAT_ELEM(v, i) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX1D(v) \
    for (int i=0; i<v.xdim; i++)

/// @defgroup VectorsMemory Memory access
/// @ingroup VectorsSpeedUp

/** A short alias to previous function
 * @ingroup VectorsMemory
 */
#define dVi(v, i) DIRECT_VEC_ELEM(v, i)

/** Access to X component
 * @ingroup VectorsMemory
 *
 * @code
 * XX(v) = 1;
 * val = XX(v);
 * @endcode
 */
#define XX(v) DIRECT_VEC_ELEM(v, 0)

/** Access to Y component
 * @ingroup VectorsMemory
 *
 * @code
 * YY(v) = 1;
 * val = YY(v);
 * @endcode
 */
#define YY(v) DIRECT_VEC_ELEM(v, 1)

/** Access to Z component
 * @ingroup VectorsMemory
 *
 * @code
 * ZZ(v) = 1;
 * val = ZZ(v);
 * @endcode
 */
#define ZZ(v) DIRECT_VEC_ELEM(v, 2)

/** @defgroup VectorsGeometry Geometry vectors
 * @ingroup VectorsSpeedUp
 *
 * The vectors involved in these macros should be created with the correct size
 * before entering in them. These macros allow a fast creation of R2 and R3
 * vectors.
 */

/** Creates vector in R2
 * @ingroup VectorsGeometry
 *
 * The vector must be created beforehand to the correct size. After this macro
 * the vector is (x, y) in R2.
 *
 * @code
 * Matrix1D< double > v(2);
 * VECTOR_R2(v, 1, 2);
 * @endcode
 */
#define VECTOR_R2(v, x, y) { \
        XX(v) = x; YY(v) = y; }

/** Creates vector in R3
 * @ingroup VectorsGeometry
 *
 * The vector must be created beforehand to the correct size. After this macro
 * the vector is (x, y, z) in R3.
 *
 * @code
 * Matrix1D< double > v(3);
 * VECTOR_R2(v, 1, 2, 1);
 * @endcode
 */
#define VECTOR_R3(v, x, y, z) { \
        XX(v) = x; YY(v) = y; ZZ(v) = z;}

/** @defgroup VectorsArithmetic Arithmethic operations
 * @ingroup VectorsSpeedUp
 *
 * The vectors involved in these macros should be created with the correct size
 * before entering in them. These macros allow a fast operation on R2 and R3
 * vectors.
 */

/// @defgroup VectorsOperations2D R2 operations
/// @ingroup VectorsArithmetic

/** Adding two R2 vectors (a=b+c)
 * @ingroup VectorOperations2D
 *
 * @code
 * Matrix1D< double > a(2), b(2), c(2);
 * ...;
 * V2_PLUS_V2(a, b, c);
 * @endcode
 */
#define V2_PLUS_V2(a, b, c) { \
        XX(a) = XX(b) + XX(c); \
        YY(a) = YY(b) + YY(c); }

/** Substracting two R2 vectors (a=b-c)
 * @ingroup VectorsOperations2D
 *
 * @code
 * Matrix1D< double > a(2), b(2), c(2);
 * ...;
 * V2_MINUS_V2(a, b, c);
 * @endcode
 */
#define V2_MINUS_V2(a, b, c) { \
        XX(a) = XX(b) - XX(c); \
        YY(a) = YY(b) - YY(c); }

/** Adding/substracting a constant to a R2 vector (a=b-k).
 * @ingroup VectorsOperations2D
 *
 * @code
 * Matrix1D< double > a(2), b(2);
 * double k;
 * ...;
 * V2_PLUS_CT(a, b, k);
 *
 * Matrix1D< double > a(2), b(2);
 * double k;
 * ...;
 * V2_PLUS_CT(a, b, -k);
 * @endcode
 */
#define V2_PLUS_CT(a, b, k) { \
        XX(a) = XX(b) + (k); \
        YY(a) = YY(b) + (k); }

/** Multiplying/dividing by a constant a R2 vector (a=b*k)
 * @ingroup VectorsOperations2D
 *
 * @code
 * Matrix1D< double > a(2), b(2);
 * double k;
 * ...;
 * V2_BY_CT(a, b, k);
 *
 * Matrix1D< double > a(2), b(2);
 * double k;
 * ...;
 * V2_BY_CT(a, b, 1/k);
 * @endcode
 */
#define V2_BY_CT(a, b, k) { \
        XX(a) = XX(b) * (k); \
        YY(a) = YY(b) * (k); }

/// @defgroup VectorsOperations3D R3 operations
/// @ingroup VectorsArithmetic

/** Adding two R3 vectors (a=b+c)
 * @ingroup VectorOperations3D
 *
 * @code
 * Matrix1D< double > a(3), b(3), c(3);
 * ...;
 * V3_PLUS_V3(a, b, c);
 * @endcode
 */
#define V3_PLUS_V3(a, b, c) { \
        XX(a) = XX(b) + XX(c); \
        YY(a) = YY(b) + YY(c); \
        ZZ(a) = ZZ(b) + ZZ(c); }

/** Substracting two R3 vectors (a=b-c)
 * @ingroup VectorsOperations3D
 *
 * @code
 * Matrix1D< double > a(3), b(3), c(3);
 * ...;
 * V3_MINUS_V3(a, b, c);
 * @endcode
 */
#define V3_MINUS_V3(a, b, c) { \
        XX(a) = XX(b) - XX(c); \
        YY(a) = YY(b) - YY(c); \
        ZZ(a) = ZZ(b) - ZZ(c); }

/** Adding/substracting a constant to a R3 vector (a=b-k)
 * @ingroup VectorsOperations3D
 *
 * @code
 * Matrix1D< double > a(3), b(3);
 * double k;
 * ...;
 * V3_PLUS_CT(a, b, k);
 *
 * Matrix1D< double > a(3), b(3);
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
 * @ingroup VectorsOperations3D
 *
 * @code
 * Matrix1D< double > a(3), b(3);
 * double k;
 * ...;
 * V3_BY_CT(a, b, k);
 *
 * Matrix1D< double > a(3), b(3);
 * double k;
 * ...;
 * V3_BY_CT(a, b, 1/k);
 * @endcode
 */
#define V3_BY_CT(a, b, c) { \
        XX(a) = XX(b) * (c); \
        YY(a) = YY(b) * (c); \
        ZZ(a) = ZZ(b) * (c); }

/// Template class for Xmipp vectors
/// @ingroup Vectors
template<typename T>
class Matrix1D: public MultidimArray<T>
{
public:
#ifndef SWIG
    bool row; ///< 0=column vector (default), 1=row vector
#endif

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
    Matrix1D(bool column = true): MultidimArray<T>() 
    {
        row = ! column;
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
    Matrix1D(int dim, bool column = true): MultidimArray<T>()
    {
        row = ! column;
        MultidimArray<T>::resize(1,1,dim);
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
    Matrix1D(const Matrix1D<T>& v): MultidimArray<T>()
    {
        *this = v;
    }

    /** Clear.
     * @ingroup VectorsConstructors
     */
     void clear()
     {
        MultidimArray<T>::clear();
     }

    /// @defgroup VectorsInitialization Initialization
    /// @ingroup Vectors

    /** Linear initialization
     * @ingroup VectorsInitialization
     *
     * The vector is filled with values increasing/decreasing linearly within a
     * range or at given steps.
     *
     * Increment functionality: The default increment is 1, the initial point is
     * incremented by this value until the upper limit is reached. This is the
     * default working mode for the function.
     *
     * @code
     * v1.initLinear(1, 3); // v1=[1 2 3]
     * v1.initLinear(1.5, 3.1); // v1=[1.5 2.5]
     * v1.initLinear(0, 10, 3); // v1=[0 3 6 9]
     * v1.initLinear(0, 10, 3, "incr"); // v1=[0 3 6 9]
     * @endcode
     *
     * Step functionality: The given range is divided in as many points as
     * indicated (in the example 6 points).
     *
     * @code
     * v1.initLinear(0, 10, 6, "steps"); // v1=[0 2 4 6 8 10]
     * @endcode
     */
    void initLinear(T minF, T maxF, int n = 1, const std::string& mode = "incr")
    {
        double slope;
        int steps;

        if (mode == "incr")
        {
            steps = 1 + (int) FLOOR((double) ABS((maxF - minF)) / ((double) n));
            slope = n * SGN(maxF - minF);
        }
        else if (mode == "steps")
        {
            steps = n;
            slope = (maxF - minF) / (steps - 1);
        }
        else
            REPORT_ERROR(1005, "Init_linear: Mode not supported (" + mode +
                         ")");

        if (steps == 0)
            clear();
        else
        {
            resize(steps);
            for (int i = 0; i < steps; i++)
                VEC_ELEM(*this, i) = (T)((double) minF + slope * i);
        }
    }

    /** Zero initialisation with a new dimension
     * @ingroup VectorsInitialization
     *
     * @code
     * v1.initZeros(6);
     * @endcode
     */
    void initZeros(int dim)
    {
        resize(dim);
        initConstant((T) 0);
    }

    /** Zero initialisation with current dimension
     * @ingroup VectorsInitialization
     */
    void initZeros()
    {
        MultidimArray<T>::initZeros();
    }

    /** Zero initialisation with current dimension
     * @ingroup VectorsInitialization
     */
    template <typename T1>
    void initZeros(const Matrix1D<T1> &m)
    {
        MultidimArray<T>::initZeros(m);
    }

    /** @defgroup VectorsSize Size and shape
     * @ingroup Vectors
     *
     * The shape of a vector is defined by its origin, its size and if it is a
     * column or a row vector. The size is clear, and the origin is the logical
     * position of the first real position of the array. For instance, if we
     * have a vector of dimension 5 and origin -2, this means that the array is
     * representing the logical positions [-2 -1 0 1 2], we could access to any
     * of these positions (Ex: v(-2) = 3;) and actually any try to access to a
     * position related to 5 (Ex: v(4) = 3;), although it physically exists, is
     * not logically correct and hence it will throw an exception. The starting
     * and finishing positions for this sample vector are -2 and 2
     * respectively, and the "for" iterations through the vector should include
     * these 2 values if you want to cover the whole vector
     *
     * @code
     * for (int i=STARTINGX(v); i<=FINISHINGX(v); i++)
     *     VEC_ELEM(v, i) += 1;
     * @endcode
     */

    /** Print vector shape.
      * @ingroup VectorsSize */
    void printShape(std::ostream& out=std::cout) const
    {
        out << "Size: " << XSIZE(*this)
            << "i=[" << STARTINGX(*this) << ".." << FINISHINGX(*this) << "]";
    }

    /** Resize to a given size
     * @ingroup VectorsSize
     *
     * This function resize the actual array to the given size. The origin is
     * not modified. If the actual array is larger than the pattern then the
     * trailing values are lost, if it is smaller then 0's are added at the end.
     * An exception is thrown if there is no memory.
     *
     * @code
     * v1.resize(3);
     * @endcode
     */
    void resize(int Xdim)
    {
        MultidimArray<T>::resize(1,1,Xdim);
    }

    /** Resize taking the shape from another vector
     * @ingroup VectorsSize
     */
    template <typename T1>
    void resize(const Matrix1D<T1> &M)
    {
    	MultidimArray<T>::resize(M);
    }

#ifndef SWIG
    /** Produce a vector suitable for working with Numerical Recipes
     * @ingroup VectorsSize
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
        return MULTIDIM_ARRAY(*this) - 1;
    }

    /** Kill an array produced for Numerical Recipes.
     * @ingroup VectorsSize
     *
     * Nothing needs to be done in fact.
     *
     * This function is not ported to Python.
     */
    void killAdaptationForNumericalRecipes(T* m) const
        {}
#endif

    /** Outside
     * @ingroup VectorsSize
     *
     * TRUE if the logical index given is outside the definition region of this
     * array.
     */
    bool outside(int i) const
    {
        return (i < STARTINGX(*this) || i > FINISHINGX(*this));
    }

    /** Outside
     * @ingroup VectorsSize
     *
     * TRUE if the logical index given is outside the definition region of this
     * array.
     */
    bool outside(const Matrix1D<double> &r) const
    {
        if (XSIZE(r) < 1)
            REPORT_ERROR(1, "Outside: index vector has not got enough components");
    
        return (XX(r) < STARTINGX(*this) || XX(r) > FINISHINGX(*this));
    }

    /** True if vector is a row.
     * @ingroup VectorsSize
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
     * @ingroup VectorsSize
     *
     * @code
     * if (v.isCol())
     *     std::cout << "v is a column vector\n";
     * @endcode
     */
    int  isCol()  const
    {
        return !row;
    }

    /** Forces the vector to be a row vector
     * @ingroup VectorsSize
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
     * @ingroup VectorsSize
     *
     * @code
     * v.setCol();
     * @endcode
     */
    void setCol()
    {
        row = false;
    }

    /// @defgroup VectorsMemory Memory access
    /// @ingroup Vectors

    /** Vector element access
     * @ingroup VectorsMemory
     *
     * Returns the value of a vector logical position. In our example we could
     * access from v(-2) to v(2). The elements can be used either by value or by
     * reference. An exception is thrown if the index is outside the logical
     * range.
     *
     * @code
     * v(-2) = 1;
     * val = v(-2);
     * @endcode
     */
    T& operator()(int i) const
    {
        if (i > FINISHINGX(*this) || i < STARTINGX(*this))
            REPORT_ERROR(1003, static_cast< std::string >
	       ("Vector subscript not defined for this vector i=")+
	       integerToString(i));

        return VEC_ELEM(*this, i);
    }

    /** Interpolates the value of the 1D vector M at the point (x) knowing
     * that this vector is a set of B-spline coefficients
     * @ingroup VectorsMemory
     *
     * (x) is in logical coordinates
     *
     * To interpolate using splines you must first produce the Bspline
     * coefficients. An example to interpolate a vector at (0.5) using
     * splines would be:
     *
     * @code
     * Matrix1D< double > Bspline_coeffs;
     * myVector.produceSplineCoefficients(Bspline_coeffs, 3);
     * interpolated_value = Bspline_coeffs.interpolatedElementBSpline(0.5,3);
     * @endcode
     */
    T interpolatedElementBSpline(double x, int SplineDegree = 3) const
    {
        int SplineDegree_1 = SplineDegree - 1;

        // Logical to physical
        x -= STARTINGX(*this);

        int lmax = XSIZE(*this);
        int l1 = CEIL(x - SplineDegree_1);
        int l2 = l1 + SplineDegree;

        double sum = 0.0;
        for (int l = l1; l <= l2; l++)
        {
            double xminusl = x - (double) l;
	    int equivalent_l=l;
	    if      (l<0)             equivalent_l=-l-1;
	    else if (l>=XSIZE(*this)) equivalent_l=2*XSIZE(*this)-l-1;
            double Coeff = (double) DIRECT_VEC_ELEM(*this,equivalent_l);
            switch (SplineDegree)
            {
            case 2:
                sum += Coeff * Bspline02(xminusl);
                break;

            case 3:
                sum += Coeff * Bspline03(xminusl);
                break;

            case 4:
                sum += Coeff * Bspline04(xminusl);
                break;

            case 5:
                sum += Coeff * Bspline05(xminusl);
                break;

            case 6:
                sum += Coeff * Bspline06(xminusl);
                break;

            case 7:
                sum += Coeff * Bspline07(xminusl);
                break;

            case 8:
                sum += Coeff * Bspline08(xminusl);
                break;

            case 9:
                sum += Coeff * Bspline09(xminusl);
                break;
            }
        }
        return (T) sum;
    }

    /** Get element at i (logical access)
     * @ingroup VectorsMemory
     */
    T getElement(int i) const
    {
        return (*this)(i);
    }

    /** Set element at i (logical access)
     * @ingroup VectorsMemory
     */
    void setElement(int i, T val)
    {
        (*this)(i) = val;
    }

    /** Matrix element access via double vector
     * @ingroup VectorsMemory
     *
     * Returns the value of a matrix logical position, but this time the element
     * position is determined by a R1 vector. The elements can be used either by
     * value or by reference. No check is done on the validity of the vector.
     */
    T& operator()(const Matrix1D< double >& v) const
    {
        return VEC_ELEM(*this, ROUND(XX(v)));
    }

    /** Matrix element access via integer vector
     * @ingroup VectorsMemory
     */
    T& operator()(const Matrix1D< int >& v) const
    {
        return VEC_ELEM(*this, XX(v));
    }

    /** Access to X component
     * @ingroup VectorsMemory
     *
     * A special case of vectors are those of 2 or 3 components, they are
     * treated as vectors in R2 or R3 for which the X component is defined. This
     * function allows you to access in a more elegant way to this component
     * than the previous function, but you could always access to the first
     * position of the array. An exception is thrown if the vector doesn't
     * belong to R2 or R3.
     *
     * @code
     * v.X() = 1;
     * val = v.X();
     * @endcode
     */
    T& X() const
    {
        if (XSIZE(*this) > 3 || XSIZE(*this) <= 0)
            REPORT_ERROR(1003, "X: Subscript not defined for this dimension");
        return XX(*this);
    }

    /** Access to Y component
     * @ingroup VectorsMemory
     *
     * A special case of vectors are those of 2 or 3 components, they are
     * treated as vectors in R2 or R3 for which the Y component is defined. This
     * function allows you to access in a more elegant way to this component
     * than the previous function, but you could always access to the second
     * position of the array. An exception is thrown if the vector doesn't
     * belong to R2 or R3.
     *
     * @code
     * v.Y() = 1;
     * val = v.Y();
     * @endcode
     */
    T& Y() const
    {
        if (XSIZE(*this) > 3 || XSIZE(*this) <= 1)
            REPORT_ERROR(1003, "Y: Subscript not defined for this dimension");
        return YY(*this);
    }

    /** Access to Z component
     * @ingroup VectorsMemory
     *
     * A special case of vectors are those of 3 components, they are treated as
     * vectors in R3 for which the Z component is defined. This function allows
     * you to access in a more elegant way to this component than the previous
     * function, but you could always access to the second position of the
     * array. An exception is thrown if the vector doesn't belong to R3.
     *
     * @code
     * v.Z() = 1;
     * val = v.Z();
     * @endcode
     */
    T& Z() const
    {
        if (XSIZE(*this) != 3)
            REPORT_ERROR(1003, "Z: Subscript not defined for this dimension");
        return ZZ(*this);
    }

    /** Logical to physical index translation
     * @ingroup VectorsMemory
     *
     * This function returns the physical position of a logical one.
     *
     * @code
     * v.toPhysical(i_log, i_phys);
     * @endcode
     */
    void toPhysical(int i_log, int& i_phys) const
    {
        i_phys = i_log - STARTINGX(*this);
    }

    /** Physical to logical index translation.
     * @ingroup VectorsMemory
     *
     * This function returns the logical position of a physical one.
     *
     * @code
     * v.toLogical(i_phys, i_log);
     * @endcode
     */
    void toLogical(int i_phys, int& i_log) const
    {
        i_log = i_phys + STARTINGX(*this);
    }

    /// @defgroup VectorOperators Operators
    /// @ingroup Matrices

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
    Matrix1D<T>& operator=(const Matrix1D<T>& op1)
    {
	if (&op1 != this)
	{
            resize(op1);
            T* ptr=NULL;
	    unsigned long int n;
            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
	    	*ptr=DIRECT_MULTIDIM_ELEM(op1,n);
            row=op1.row;
	}

	return *this;
    }

    /** Unary minus.
     * @ingroup VectorOperators
     *
     * It is used to build arithmetic expressions. You can make a minus
     * of anything as long as it is correct semantically.
     */
    Matrix1D<T> operator-() const
    {
        Matrix1D<T> tmp(*this);
	T* ptr;
	unsigned long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(tmp,n,ptr)
            *ptr = -(*ptr);
        return tmp;
    }

    /** v3 = v1 + v2.
     * @ingroup VectorOperators
     */
    Matrix1D<T> operator+(const Matrix1D<T>& op1) const
    {
        Matrix1D<T> tmp;
        arrayByArray(*this, op1, tmp, '+');
        return tmp;
    }

    /** v3 = v1 - v2.
     * @ingroup VectorOperators
     */
    Matrix1D<T> operator-(const Matrix1D<T>& op1) const
    {
        Matrix1D<T> tmp;
        arrayByArray(*this, op1, tmp, '-');
        return tmp;
    }

    /** v3 = v1 * v2.
     * @ingroup VectorOperators
     */
    Matrix1D<T> operator*(const Matrix1D<T>& op1) const
    {
        Matrix1D<T> tmp;
        arrayByArray(*this, op1, tmp, '*');
        return tmp;
    }

    /** v3 = v1 / v2.
     * @ingroup VectorOperators
     */
    Matrix1D<T> operator/(const Matrix1D<T>& op1) const
    {
        Matrix1D<T> tmp;
        arrayByArray(*this, op1, tmp, '/');
        return tmp;
    }

    /** v3 += v2.
     * @ingroup VectorOperators
     */
    void operator+=(const Matrix1D<T>& op1)
    {
        arrayByArray(*this, op1, *this, '+');
    }

    /** v3 -= v2.
     * @ingroup VectorOperators
     */
    void operator-=(const Matrix1D<T>& op1)
    {
        arrayByArray(*this, op1, *this, '-');
    }

    /** v3 *= v2.
     * @ingroup VectorOperators
     */
    void operator*=(const Matrix1D<T>& op1)
    {
    	Matrix1D<T> tmp;
        arrayByArray(*this, op1, tmp, '*');
	*this=tmp;
    }

    /** v3 /= v2.
     * @ingroup VectorOperators
     */
    void operator/=(const Matrix1D<T>& op1)
    {
        arrayByArray(*this, op1, *this, '/');
    }

    /** v3 = v1 + k.
     * @ingroup VectorOperators
     */
    Matrix1D<T> operator+(T op1) const
    {
        Matrix1D<T> tmp;
        arrayByScalar(*this, op1, tmp, '+');
        return tmp;
    }

    /** v3 = v1 - k.
     * @ingroup VectorOperators
     */
    Matrix1D<T> operator-(T op1) const
    {
        Matrix1D<T> tmp;
        arrayByScalar(*this, op1, tmp, '-');
        return tmp;
    }

    /** v3 = v1 * k.
     * @ingroup VectorOperators
     */
    Matrix1D<T> operator*(T op1) const
    {
        Matrix1D<T> tmp;
        arrayByScalar(*this, op1, tmp, '*');
        return tmp;
    }

    /** v3 = v1 / k.
     * @ingroup VectorOperators
     */
    Matrix1D<T> operator/(T op1) const
    {
        Matrix1D<T> tmp;
        arrayByScalar(*this, op1, tmp, '/');
        return tmp;
    }

    /** v3 += k.
     * @ingroup VectorOperators
     *
     * This function is not ported to Python.
     */

    void operator+=(const T& op1)
    {
        arrayByScalar(*this, op1, *this, '+');
    }

    /** v3 -= k.
     * @ingroup VectorOperators
     *
     * This function is not ported to Python.
     */
    void operator-=(const T& op1)
    {
        arrayByScalar(*this, op1, *this, '-');
    }

    /** v3 *= k.
     * @ingroup VectorOperators
     *
     * This function is not ported to Python.
     */
    void operator*=(const T& op1)
    {
        arrayByScalar(*this, op1, *this, '*');
    }

    /** v3 /= k.
     * @ingroup VectorOperators
     *
     * This function is not ported to Python.
     */
    void operator/=(const T& op1)
    {
        arrayByScalar(*this, op1, *this, '/');
    }

    /** v3 = k + v2.
     * @ingroup VectorOperators
     */
    friend Matrix1D<T> operator+(T op1, const Matrix1D<T>& op2)
    {
        Matrix1D<T> tmp;
        scalarByArray(op1, op2, tmp, '+');
        return tmp;
    }

    /** v3 = k - v2.
     * @ingroup VectorOperators
     */
    friend Matrix1D<T> operator-(T op1, const Matrix1D<T>& op2)
    {
        Matrix1D<T> tmp;
        scalarByArray(op1, op2, tmp, '-');
        return tmp;
    }

    /** v3 = k * v2.
     * @ingroup VectorOperators
     */
    friend Matrix1D<T> operator*(T op1, const Matrix1D<T>& op2)
    {
        Matrix1D<T> tmp;
        scalarByArray(op1, op2, tmp, '*');
        return tmp;
    }

    /** v3 = k / v2
     * @ingroup VectorOperators
     */
    friend Matrix1D<T> operator/(T op1, const Matrix1D<T>& op2)
    {
        Matrix1D<T> tmp;
        scalarByArray(op1, op2, tmp, '/');
        return tmp;
    }

    /// @defgroup VectorsUtilities Utilities
    /// @ingroup Vectors

    /** Vector by matrix
     * @ingroup VectorsUtilities
     *
     * Algebraic vector by matrix multiplication. This function is actually
     * implemented in xmippMatrices2D
     */
    Matrix1D<T> operator*(const Matrix2D<T>& M);

    /** Algebraic transpose of vector
     * @ingroup VectorsUtilities
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
     * @ingroup VectorsUtilities
     *
     * The same as before but the result is stored in this same object.
     */
    void selfTranspose()
    {
        row = !row;
    }

    /** Reverse vector values, keep in this object
     * @ingroup VectorsUtilities
     */
    void selfReverse()
    {
        MultidimArray<T>::selfReverseX();
    }

    /** Module of the vector
     * @ingroup VectorsUtilities
     *
     * This module is defined as the square root of the sum of the squared
     * components. Euclidean norm of the vector.
     *
     * @code
     * double mod = v.module();
     * @endcode
     */
    double module() const
    {
        return sqrt(MultidimArray<T>::sum2());
    }

    /** Angle of the vector
     * @ingroup VectorsUtilities
     *
     * Supposing this vector is in R2 this function returns the angle of this
     * vector with X axis, ie, atan2(YY(v), XX(v))
     */
    double angle()
    {
        return atan2((double) YY(*this), (double) XX(*this));
    }

    /** Normalize this vector, store the result here
     * @ingroup VectorsUtilities
     */
    void selfNormalize()
    {
        double m = module();
        if (ABS(m) > XMIPP_EQUAL_ACCURACY)
            *this /= (T) m;
        else
            MultidimArray<T>::initZeros();
    }

    /** Compute the center of mass within a mask.
      * If no mask is to be used, supply NULL.
      * @ingroup VectorsUtilities */
    void centerOfMass(Matrix1D< double >& center, void* mask=NULL)
    {
	center.initZeros(1);
	double mass = 0;
	Matrix1D< int >* imask = (Matrix1D< int >*) mask;

	FOR_ALL_ELEMENTS_IN_MATRIX1D(*this)
	{
            if ((imask == NULL || VEC_ELEM(*imask, i)) &&
		VEC_ELEM(*this, i))
            {
        	XX(center) += i * VEC_ELEM(*this, i);
        	mass += VEC_ELEM(*this, i);
            }
	}
	if (mass != 0)
            center /= mass;
    }

    /** Sort vector elements
     * @ingroup VectorsUtilities
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
        Matrix1D< double > aux;

        if (XSIZE(*this) == 0)
            return temp;

        // Initialise data
        typeCast(*this, aux);

        // Sort
        double * aux_array = aux.adaptForNumericalRecipes();
        qcksrt(XSIZE(*this), aux_array);

        typeCast(aux, temp);
        return temp;
    }

    /** Gives a vector with the indexes for a sorted vector
     * @ingroup VectorsUtilities
     *
     * This function returns the indexes of a sorted vector. The input vector is
     * not modified at all. For instance, if the input vector is [3 2 -1 0] the
     * result of this function would be [2 3 1 0] meaning that the lowest value
     * is at index 2, then comes the element at index 3, ...
     *
     * @code
     * v2 = v1.indexSort();
     * @endcode
     */
    Matrix1D< int > indexSort() const
    {
        Matrix1D< int >   indx;
        Matrix1D< double > temp;

        if (XSIZE(*this) == 0)
            return indx;

        if (XSIZE(*this) == 1)
        {
            indx.resize(1);
            indx(0) = 1;
            return indx;
        }

        // Initialise data
        indx.resize(XSIZE(*this));
        typeCast(*this, temp);

        // Sort indexes
        double* temp_array = temp.adaptForNumericalRecipes();
        int* indx_array = indx.adaptForNumericalRecipes();
        indexx(XSIZE(*this), temp_array, indx_array);

        return indx;
    }

    /** Put a window to vector
     * @ingroup VectorsUtilities
     *
     * The vector is windowed within the two indexes given to this function.
     * Indexes always refer to logical indexes. If an index is outside the
     * actual vector range then the vector is padded winit_value. In the
     * following examples suppose that v1=[-2 -1 0 1 2] and that the origin is
     * -2.
     *
     * @code
     * v1.window(-1, 2); // v1=[-1 0 1 2]; v1.startingX() == -1
     *
     * v1.window(-3, 1); // v1=[0 -2 -1 0 1]; v1.startingX() == -3
     * @endcode
     */
    void window(int x0, int xF, T init_value = 0)
    {
        Matrix1D<T> result(xF - x0 + 1);
        STARTINGX(result) = x0;

        for (int j = x0; j <= xF; j++)
            if (j >= STARTINGX(*this) && j <= FINISHINGX(*this))
                VEC_ELEM(result, j) = VEC_ELEM(*this, j);
            else
                VEC_ELEM(result, j) = init_value;

        *this = result;
    }

    /** Maximum element
     * @ingroup VectorsUtilities
     *
     * This function returns the index of the maximum element of an array.
     * Returns -1 if the array is empty.
     */
    void maxIndex(int& imax) const
    {
        if (XSIZE(*this) == 0)
        {
            imax = -1;
            return;
        }

        imax = 0;
        T max = VEC_ELEM(*this, imax);

        FOR_ALL_ELEMENTS_IN_MATRIX1D(*this)
            if (VEC_ELEM(*this, i) > max)
            {
                max = VEC_ELEM(*this, i);
                imax = i;
            }
    }

    /** Minimum element
     * @ingroup VectorsUtilities
     *
     * This function returns the index of the minimum element of an array.
     * Returns -1 if the array is empty.
     */
    void minIndex(int& imin) const
    {
        if (XSIZE(*this) == 0)
        {
            imin = -1;
            return;
        }

        imin = 0;
        T min = VEC_ELEM(*this, imin);

        FOR_ALL_ELEMENTS_IN_MATRIX1D(*this)
            if (VEC_ELEM(*this, i) < min)
            {
                min = VEC_ELEM(*this, i);
                imin = i;
            }
    }

    /** Show using gnuplot
     * @ingroup VectorsUtilities
     *
     * This function uses gnuplot to plot this vector. You must supply the
     * xlabel, ylabel, and title.
     */
    void showWithGnuPlot(const std::string& xlabel, const std::string& title)
    {
        FileName fn_tmp;
        fn_tmp.init_random(10);
        MultidimArray<T>::write(static_cast<std::string>("PPP") +
            fn_tmp + ".txt");

        std::ofstream fh_gplot;
        fh_gplot.open((static_cast<std::string>("PPP") + fn_tmp +
            ".gpl").c_str());
        if (!fh_gplot)
            REPORT_ERROR(1,
            static_cast<std::string>("vector::showWithGnuPlot: Cannot open PPP")
                + fn_tmp + ".gpl for output");
        fh_gplot << "set xlabel \"" + xlabel + "\"\n";
        fh_gplot << "plot \"PPP" + fn_tmp + ".txt\" title \"" + title +
        "\" w l\n";
        fh_gplot << "pause 300 \"\"\n";
        fh_gplot.close();
        system((static_cast<std::string>("(gnuplot PPP") + fn_tmp +
            ".gpl; rm PPP" + fn_tmp + ".txt PPP" + fn_tmp + ".gpl) &").c_str());
    }
};

/**@defgroup VectorsRelated Related functions
 * @ingroup Vectors
 *
 * These functions are not methods of Matrix1D
 */

/// @defgroup VectorsGeom Geometry with vectors
/// @ingroup VectorsRelated

/** Creates vector in R2
 * @ingroup VectorsGeom
 *
 * After this function the vector is (x,y) in R2.
 *
 * @code
 * Matrix1D< double > v = vectorR2(1, 2);
 * @endcode
 */
Matrix1D< double > vectorR2(double x, double y);

/** Creates vector in R3
 * @ingroup VectorsGeom
 *
 * After this function the vector is (x,y,z) in R3.
 *
 * @code
 * Matrix1D< double > v = vectorR2(1, 2, 1);
 * @endcode
 */
Matrix1D< double > vectorR3(double x, double y, double z);

/** Creates an integer vector in Z3
 * @ingroup VectorsGeom
 */
Matrix1D< int > vectorR3(int x, int y, int z);

/** Dot product.
 * @ingroup VectorsGeom
 *
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
T dotProduct(const Matrix1D< T >& v1, const Matrix1D< T >& v2)
{
    if (!v1.sameShape(v2))
        REPORT_ERROR(1002, "Dot product: vectors of different size or shape");

    T accumulate = 0;
    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX1D(v1)
        accumulate += DIRECT_VEC_ELEM(v1,i) * DIRECT_VEC_ELEM(v2,i);

    return accumulate;
}

/** Vector product in R3
 * @ingroup VectorsUtilities
 *
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
Matrix1D< T > vectorProduct(const Matrix1D< T >& v1, const Matrix1D< T >& v2)
{
    if (v1.xdim != 3 || v2.xdim != 3)
        REPORT_ERROR(1002, "Vector_product: vectors are not in R3");

    if (v1.isRow() != v2.isRow())
        REPORT_ERROR(1007, "Vector_product: vectors are of different shape");

    Matrix1D< T > result(3);
    result.X() = v1.Y() * v2.Z() - v1.Z() * v2.Y();
    result.Y() = v1.Z() * v2.X() - v1.X() * v2.Z();
    result.Z() = v1.X() * v2.Y() - v1.Y() * v2.X();

    return result;
}

/// @defgroup VectorsMiscellaneous Miscellaneous
/// @ingroup VectorsRelated

/** Vector equality.
 * @ingroup VectorsMiscellaneous */
template<typename T>
bool operator==(const Matrix1D<T>& op1, const Matrix1D<T>& op2)
{
    return op1.equal(op2);
}

/** Vector inequality.
 * @ingroup VectorsMiscellaneous */
template<typename T>
bool operator!=(const Matrix1D<T>& op1, const Matrix1D<T>& op2)
{
    return !(op1==op2);
}

/** Reduce both vectors to a common size
 * @ingroup VectorsMiscellaneous
 *
 * Search the range of logical indexes for which both vectors have got valid
 * values, and cut both vectors to that size, the corresponding origin is
 * automatically computed.
 *
 * @code
 * Matrix1D< double > v1(5);
 * v1.startingX() = -2;
 *
 * Matrix1D< double > v2(4);
 * v2.startingX() = 0;
 *
 * cutToCommonSize(v1, v2);
 * // v1 and v2 range from 0 to 2
 * @endcode
 */
template<typename T>
void cutToCommonSize(Matrix1D<T>& V1, Matrix1D<T>& V2)
{
    int x0 = XMIPP_MAX(STARTINGX(V1), STARTINGX(V2));
    int xF = XMIPP_MIN(FINISHINGX(V1), FINISHINGX(V2));
    V1.window(x0, xF);
    V2.window(x0, xF);
}

/** Sort two vectors
 * @ingroup VectorsMiscellaneous
 *
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
        REPORT_ERROR(1007,
                     "sortTwoVectors: vectors are not of the same shape");

    FOR_ALL_ELEMENTS_IN_MATRIX1D(v1)
    {
        temp = XMIPP_MIN(VEC_ELEM(v1, i), VEC_ELEM(v2, i));
        VEC_ELEM(v2, i) = XMIPP_MAX(VEC_ELEM(v1, i), VEC_ELEM(v2, i));
        VEC_ELEM(v1, i) = temp;
    }
}

#ifndef SWIG
/** Optimize using Powell's method.
  * @ingroup VectorsMiscellaneous
  *
  * See Numerical Recipes Chapter 10.
  *
  * Problem: Minimize f(x) starting at point p. n is the dimension of x.
  * If changes are smaller than ftol then stop. The number of iterations is
  * returned in iter, fret contains the function value at the minimum and p
  * contains the minimum.
  *
  * Watch out that the evaluating function f must consider that x starts at
  * index 1, at least, and goes until n. i0 is used to allow optimization in
  * large spaces where only one part is going to be optimized. Thus, if in a
  * vector of dimension 20 you want to optimize the first 3 components then
  * i0=1, n=3; if you want to optimize it all, i0=1, n=20; and if you want to
  * optimize the last five components i0=15, n=5.
  *
  * The steps define the allowed steps on each variable. When a variable has
  * large ranges this step should be larger to allow bigger movements. The steps
  * should have only the dimension of the optimized part (3,20 and 5 in the
  * examples above).
  *
  * The option show forces the routine to show the convergence path
  *
  * Example of use:
  *
  * @code
  * Matrix1D<double> x(8), steps(8);
  * double fitness;
  * int iter;
  * steps.initConstant(1);
  * x,initZeros();
  * powellOptimizer(x,1,8,&wrapperFitness,0.01,fitness,iter,steps,true);
  * @endcode
  *
  */
void powellOptimizer(Matrix1D< double >& p,
                     int i0, int n,
                     double(*f)(double* x),
                     double ftol,
                     double& fret,
                     int& iter,
                     const Matrix1D< double >& steps,
                     bool show = false);
#endif

/** Show a vector.
  * @ingroup VectorsUtilities */
template<typename T>
std::ostream& operator<<(std::ostream& out, const Matrix1D<T>& v)
{
    if (XSIZE(v) == 0)
        out << "NULL vector\n";
    else
    {
        if (typeid(T) == typeid(std::complex<double>))
        {
            FOR_ALL_ELEMENTS_IN_MATRIX1D(v)
            {
                if (v.row)
                    out << VEC_ELEM(v, i) << ' ';
                else
                    out << VEC_ELEM(v, i) << '\n';
            }
        }
        else
        {
            // Look for the exponent
            Matrix1D<T> aux(v);
            aux.selfABSnD();
            int prec = bestPrecision(aux.computeMax(), 10);

            FOR_ALL_ELEMENTS_IN_MATRIX1D(v)
            {
                if (v.row)
                    out << floatToString((double) VEC_ELEM(v, i), 10, prec) << ' ';
                else
                    out << floatToString((double) VEC_ELEM(v, i), 10, prec) << '\n';
            }
        }
    }

    return out;
}
#endif
