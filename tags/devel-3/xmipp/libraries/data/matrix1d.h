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

#include "funcs.h"

extern int bestPrecision(float F, int _width);
extern std::string floatToString(float F, int _width, int _prec);

/** @defgroup Vectors Vectors speed up macros
 * @ingroup MultidimArraysSpeedUp
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

/** Vector element: Physical access
 * @ingroup MultidimArraySizeShape
 *
 * Be careful because this is physical access, usually vectors follow the C
 * convention of starting index==0. This function should not be used as it goes
 * against the vector library philosophy unless you explicitly want to access
 * directly to any value in the vector without taking into account its logical
 * position.
 *
 * @code
 * DIRECT_VEC_ELEM(v, 0) = 1;
 * val = DIRECT_VEC_ELEM(v, 0);
 * @endcode
 */
#define DIRECT_VEC_ELEM(v, i) ((v).data[(i)])

/** A short alias to previous function
 * @ingroup MultidimArraySizeShape
 */
#define dVi(v, i) DIRECT_VEC_ELEM(v, i)


/** Access to X component
 * @ingroup Vectors
 *
 * @code
 * XX(v) = 1;
 * val = XX(v);
 * @endcode
 */
#define XX(v) (v).data[0]

/** Access to Y component
 * @ingroup Vectors
 *
 * @code
 * YY(v) = 1;
 * val = YY(v);
 * @endcode
 */
#define YY(v) (v).data[1]

/** Access to Z component
 * @ingroup Vectors
 *
 * @code
 * ZZ(v) = 1;
 * val = ZZ(v);
 * @endcode
 */
#define ZZ(v) (v).data[1]

/** Creates vector in R2
 * @ingroup Vectors
 *
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
 * @ingroup Vectors
 *
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
 * @ingroup Vector
 *
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
 * @ingroup Vectors
 *
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
 * @ingroup Vectors
 *
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
 * @ingroup Vectors
 *
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
 * @ingroup Vector
 *
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
 * @ingroup Vectors
 *
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
 * @ingroup Vectors
 *
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
 * @ingroup Vectors
 *
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



template<typename T>
class Matrix1D
{
public:
    /* The array itself.
    */
    T* data;

    // Destroy data
    bool destroyData;

    // Number of elements
    int xdim;

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
    Matrix1D(bool column = true)
    {
    	coreInit();
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
    Matrix1D(int dim, bool column = true)
    {
        row = ! column;
        resize(dim);
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
    Matrix1D(const Matrix1D<T>& v, bool column=true)
    {
        coreInit();
        *this = v;
        row!=column;
    }

    /** Destructor.
     * @ingroup MultidimArrayConstructors
     */
     ~Matrix1D()
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
        xdim=0;
        row=false;
        data=NULL;
        destroyData=true;
    }

    /** Core allocate.
     * @ingroup MultidimArrayCore
     */
    void coreAllocate(int _xdim)
    {
        if (_xdim<=0)
        {
            clear();
            return;
        }

        xdim=_xdim;
        data = new T [xdim];
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
    void resize(int Xdim)
    {
        if (Xdim == xdim)
            return;

        if (Xdim <= 0)
        {
            clear();
            return;
        }


        T * new_data;

        try
        {
        	new_data = new T [Xdim];
        }
        catch (std::bad_alloc &)
		{
			REPORT_ERROR(1001, "Allocate: No space left");
		}

		// Copy needed elements, fill with 0 if necessary
		for (int j = 0; j < Xdim; j++)
		{
			T val;
			if (j >= xdim)
				val = 0;
			else
				val = data[j];
			new_data[j] = val;
		}

		// deallocate old vector
		coreDeallocate();

		// assign *this vector to the newly created
		data = new_data;
		xdim = Xdim;

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
    	for (int j = 0; j < xdim; j++)
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
    void initZeros(int Xdim)
    {
        resize(Xdim);
        initZeros();
    }

    /** Same shape.
     * @ingroup MultidimSize
     *
     * Returns true if this object has got the same shape (origin and size)
     * than the argument
     */
    template <typename T1>
    bool sameShape(const Matrix1D<T1>& op) const
    {
        return (xdim == op.xdim);
    }

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
        return (*this).data - 1;
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
     * reference.
     *
     * @code
     * v(-2) = 1;
     * val = v(-2);
     * @endcode
     */
    T& operator()(int i) const
    {
        return data[i];
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
    Matrix1D<T>& operator=(const Matrix1D<T>& op1)
    {
        if (&op1 != this)
        {
            resize(op1);
            for (int i = 0; i < xdim; i++)
            	data[i] = op1.data[i];
            row=op1.row;
        }

        return *this;
    }

    /// @defgroup VectorsUtilities Utilities
    /// @ingroup Vectors

    /** Vector by matrix
     * @ingroup VectorsUtilities
     *
     * Algebraic vector by matrix multiplication. This function is actually
     * implemented in xmippMatrices2D
     */
    //Matrix1D<T> operator*(const Matrix2D<T>& M);

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

    /** Sum of squared vector values.
     * @ingroup VectorsUtilities
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
		for (int j = 0; j < xdim; j++)
		{
			sum += data[j] * data[j];
		}
		return sum;
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
        return sqrt(sum2());
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
            initZeros();
    }

    /** Reverse matrix values over X axis, keep in this object.
     * @ingroup Utilites
     *
     * Maybe better with an example:
     *
     * @code
     * slice 0
     * [01 02 03          [07 08 09
     *  04 05 06           04 05 06
     *  07 08 09]          01 02 03]
     *
     * ----->
     *
     * slice 1
     * [11 12 13          [17 18 19
     *  14 15 16           14 15 16
     *  17 18 19]          11 12 13]
     * @endcode
     *
     */
    void selfReverse()
    {
    	for (int j = 0; j <= (int)(xdim - 1) / 2; j++)
    	{
    		T aux;
    		SWAP(data[j], data[xdim-1-j], aux);
    	}
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

        if (xdim == 0)
            return temp;

        // Initialise data
        typeCast(*this, aux);

        // Sort
        double * aux_array = aux.adaptForNumericalRecipes();
        qcksrt(xdim, aux_array);

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

        if (xdim == 0)
            return indx;

        if (xdim == 1)
        {
            indx.resize(1);
            indx(0) = 1;
            return indx;
        }

        // Initialise data
        indx.resize(xdim);
        typeCast(*this, temp);

        // Sort indexes
        double* temp_array = temp.adaptForNumericalRecipes();
        int* indx_array = indx.adaptForNumericalRecipes();
        indexx(xdim, temp_array, indx_array);

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
        Matrix1D<T>::write(static_cast<std::string>("PPP") +
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

    /** Read from an ASCII file.
     * @ingroup Operators
     *
     * The array must be previously resized to the correct size.
     */
    void read(const FileName& fn)
    {
        std::ifstream in;
        in.open(fn.c_str(), std::ios::in);

        if (!in)
            REPORT_ERROR(1,
                         static_cast< std::string >("MultidimArray::read: File " +
                                                    fn + " not found"));

        in >> *this;
        in.close();
    }


    /** Write to an ASCII file.
     * @ingroup Operators
     */
    void write(const FileName& fn) const
    {
        std::ofstream out;
        out.open(fn.c_str(), std::ios::out);
        if (!out)
            REPORT_ERROR(1,
                         static_cast< std::string >("MultidimArray::write: File " + fn
                                                    + " cannot be opened for output"));

        out << *this;
        out.close();
    }

    /** Edit with xmipp_editor.
     * @ingroup Operators
     *
     *
     * This function generates a random filename starting with PPP and
     * edits it with xmipp_editor. After closing the editor the file is
     * removed.
     */
    // FIXME This is not good practice.. (system)
    void edit()
    {
        FileName nam;
        nam.init_random(15);

        nam = static_cast< std::string >("PPP" + nam + ".txt");
        write(nam);

        system((static_cast< std::string >("xmipp_edit -i " + nam +
                                           " -remove &").c_str()));
    }



};

 template<typename T>
 std::ostream& operator<<(std::ostream& ostrm, const Matrix1D<T>& v)
 {
     if (v.xdim == 0)
         ostrm << "NULL Array\n";
     else
         ostrm << std::endl;

     double max_val = ABS(v.data[0]);

     for (int j = 0; j < v.xdim; j++)
     {
    	 max_val = XMIPP_MAX(max_val, v.data[j]);
     }

     int prec = bestPrecision(max_val, 10);

     for (int j = 0; j < v.xdim; j++)
     {
    	 ostrm << floatToString((double) v.data[j], 10, prec)
    	 << std::endl;
     }
     return ostrm;
 }

 /**@defgroup VectorsRelated Related functions
  * @ingroup MultidimensionalArrays
  *
  * These functions are not methods of MultidimArray
  */

 /** Creates vector in R2
  * @ingroup VectorsRelated
  *
  * After this function the vector is (x,y) in R2.
  *
  * @code
  * Matrix1D< double > v = vectorR2(1, 2);
  * @endcode
  */
 Matrix1D< double > vectorR2(double x, double y);

 /** Creates vector in R3
  * @ingroup VectorsRelated
  *
  * After this function the vector is (x,y,z) in R3.
  *
  * @code
  * Matrix1D< double > v = vectorR2(1, 2, 1);
  * @endcode
  */
 Matrix1D< double > vectorR3(double x, double y, double z);

 /** Creates an integer vector in Z3
  * @ingroup VectorsRelated
  */
 Matrix1D< int > vectorR3(int x, int y, int z);

 /** Dot product.
  * @ingroup VectorsRelated
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
     for (int j = 0; j < v1.xdim; j++)
     {
    	 accumulate += v1.data[j] * v2.data[j];
     }
     return accumulate;
 }

 /** Vector product in R3
  * @ingroup VectorsRelated
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
  * @ingroup VectorsRelated
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

     for (int j = 0; j < v1.xdim; j++)
     {
    	 temp       = XMIPP_MIN(v1.data[j], v2.data[j]);
    	 v2.data[j] = XMIPP_MAX(v1.data[j], v2.data[j]);
    	 v1.data[j] = temp;
     }
 }

 /** Conversion from one type to another.
  * @ingroup MultidimFunctions
  *
  * If we have an integer array and we need a double one, we can use this
  * function. The conversion is done through a type casting of each element
  * If n >= 0, only the nth volumes will be converted, otherwise all NSIZE volumes
  */
 template<typename T1, typename T2>
 void typeCast(const Matrix1D<T1>& v1,  Matrix1D<T2>& v2)
 {
     if (v1.xdim == 0)
     {
         v2.clear();
         return;
     }

     v2.resize(v1.xdim);
     for (int j = 0; j < v1.xdim; j++)
     {
    	 v2.data[j] = static_cast< T2 > (v1.data[j]);
     }

 }

#endif /* MATRIX1D_H_ */
