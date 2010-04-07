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

#ifndef MATRIX_H
#define MATRIX_H

#include <fstream>
#include <complex>
#include <cmath>

#include "core_multidim_array.h"
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
class Matrix1D: public coreMultidimArray<T>
{
public:
    bool row; ///< 0=column vector (default), 1=row vector

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
    Matrix1D(bool column = true): coreMultidimArray<T>() 
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
    Matrix1D(int Xdim, bool column = true): coreMultidimArray<T>()
    {
        row = ! column;
        coreMultidimArray<T>::resize(1,1,1,Xdim);
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
    Matrix1D(const Matrix1D<T>& v): coreMultidimArray<T>()
    {
        *this = v;
    }

    /// @defgroup VectorsInitialization Initialization
    /// @ingroup Vectors

    /** Alias a multidimarray.
     * @ingroup VectorsInitialization
     *
     * Treat the multidimarray as if it were a vector. The data is not copied
     * into new memory, but a pointer to the multidimarray is copied.
     * You should not make any operation on this vector such that the
     * memory locations are changed
     */
     void alias(const coreMultidimArray<T> &m)
     {
         this->row=false;
         copyShape(m);
         this->data=m.data;
         this->destroyData=false;
     }

    /** Zero initialisation with current dimension
     * @ingroup VectorsInitialization
     */
    template <typename T1>
    void initZeros(const Matrix1D<T1> &m)
    {
        coreMultidimArray<T>::initZeros(m);
    }

    /** Zero initialisation with given dimension
     * @ingroup VectorsInitialization
     */
    void initZeros(int Xdim)
    {
        coreMultidimArray<T>::initZeros(Xdim);
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


    /** Resize with size
     * @ingroup VectorsSize
     * 
     * This function is needed to be able to resize inside this class
     */
    void resize(int Xdim)
    {
    	coreMultidimArray<T>::resize(Xdim);
    }

    /** Resize taking the shape from another vector
     * @ingroup VectorsSize
     */
    template <typename T1>
    void resize(const Matrix1D<T1> &M)
    {
    	coreMultidimArray<T>::resize(M);
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
        coreMultidimArray<T>::selfReverseX();
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
        return sqrt(coreMultidimArray<T>::sum2());
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
        {
            T im=(T) (1.0/m);
            *this *= im;
        }
        else
            coreMultidimArray<T>::initZeros();
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
        double * aux_array = aux.adaptForNumericalRecipes1D();
        qcksrt(XSIZE(*this), aux_array);

        typeCast(aux, temp);
        return temp;
    }

    /** Gives a vector with the indexes for a sorted vector
     * @ingroup VectorsUtilities
     *
     * This function returns the indexes of a sorted vector. The input vector is
     * not modified at all. For instance, if the input vector is [3 2 -1 0] the
     * result of this function would be [3 4 2 1] meaning that the lowest value
     * is at index 3, then comes the element at index 4, ... Note that
     * indexes start at 1.
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
        double* temp_array = temp.adaptForNumericalRecipes1D();
        int* indx_array = indx.adaptForNumericalRecipes1D();
        indexx(XSIZE(*this), temp_array, indx_array);

        return indx;
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
        coreMultidimArray<T>::write(static_cast<std::string>("PPP") +
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

    /** Compute numerical derivative
     * @ingroup VectorsUtilities
     *
     * The numerical derivative is of the same size as the input vector.
     * However, the first two and the last two samples are set to 0,
     * because the numerical method is not able to correctly estimate the
     * derivative there.
     */
    void numericalDerivative(Matrix1D<double> &result)
    {
        const double i12=1.0/12.0;
	result.initZeros(*this);
	for (int i=STARTINGX(*this)+2; i<=FINISHINGX(*this)-2; i++)
            result(i)=i12*(-(*this)(i+2)+8*(*this)(i+1)
                           -8*(*this)(i-1)+(*this)(i+2));
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
    XX(result) = YY(v1) * ZZ(v2) - ZZ(v1) * YY(v2);
    YY(result) = ZZ(v1) * XX(v2) - XX(v1) * ZZ(v2);
    ZZ(result) = XX(v1) * YY(v2) - YY(v1) * XX(v2);

    return result;
}

/** Vector product in R3
 * @ingroup VectorsUtilities
 *
 * This function computes the vector product of two R3 vectors.
 * No check is performed, it is assumed that the output vector
 * is already resized
 *
 */
template<typename T>
void vectorProduct(const Matrix1D< T >& v1, const Matrix1D< T >& v2,
   Matrix1D<T> &result)
{
    XX(result) = YY(v1) * ZZ(v2) - ZZ(v1) * YY(v2);
    YY(result) = ZZ(v1) * XX(v2) - XX(v1) * ZZ(v2);
    ZZ(result) = XX(v1) * YY(v2) - YY(v1) * XX(v2);
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

/** Generate random permutation
 * @ingroup VectorsMiscellaneous
 *
 * Generate a random permutation of the numbers between 0 and N-1
 */
void randomPermutation(int N, Matrix1D<int>& result);

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
  * If your function needs extra parameters you can pass them through 
  * the void pointer prm. If you don't need this feature set it to NULL.
  *
  * Example of use:
  *
  * @code
  * Matrix1D<double> x(8), steps(8);
  * double fitness;
  * int iter;
  * steps.initConstant(1);
  * x.initZeros();
  * powellOptimizer(x,1,8,&wrapperFitness,NULL,0.01,fitness,iter,steps,true);
  * @endcode
  *
  */
void powellOptimizer(Matrix1D< double >& p,
                     int i0, int n,
                     double(*f)(double* , void *),
                     void *prm, 
                     double ftol,
                     double& fret,
                     int& iter,
                     const Matrix1D< double >& steps,
                     bool show = false);

/** Gaussian interpolator
 * @ingroup VectorsMiscellaneous 
 *
 * This class helps to perform a quick evaluation of the N(0,1) Gaussian.
 * 1/sqrt(2*PI)*exp(-x^2)
 *
 * @code
 *  GaussianInterpolator GI;
 *  GI.initialize(6,60000);
 *  double g=GI.getValue(1.96);
 * @endcode
 */
class GaussianInterpolator {
    Matrix1D<double> v;
    double xstep;
    double xmax;
    double ixstep;
public:
    /** Constructor.
        xmax is the maximum value considered by the interpolator.
        N is the number of samples between 0 and xmax. 
        If normalize is set to true, then the factor 1/sqrt(2*PI)
        is introduced. */
    void initialize(double xmax, int N, bool normalize=true);
    
    /** Value */
    inline double getValue(double x) const {
        x=ABS(x);
        if (x>xmax) return 0;
        else
        {
            double aux=x*ixstep;
            int iaux=ROUND(aux);
            return DIRECT_VEC_ELEM(v,iaux);
        }
    }
};

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
