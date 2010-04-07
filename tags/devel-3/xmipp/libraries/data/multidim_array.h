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

#ifndef MULTIDIM_ARRAY_H
#define MULTIDIM_ARRAY_H

#include <typeinfo>
#include <external/bilib/types/tsplinebasis.h>
#include <external/bilib/types/tboundaryconvention.h>
#include <external/bilib/headers/linearalgebra.h>
#include <external/bilib/headers/changebasis.h>
#include <external/bilib/headers/kernel.h>
#include <external/bilib/headers/pyramidtools.h>
#include "funcs.h"
#include "error.h"
#include "args.h"
#include "core_multidim_array.h"
#include "matrix1d.h"


// Forward declarations ====================================================
template<typename T> class MultidimArray;

template<typename T>
void coreArrayByArray(const MultidimArray<T>& op1, const MultidimArray<T>& op2,
    MultidimArray<T>& result, char operation);


/** Template class for Xmipp arrays.
  * @ingroup coreMultidimensionalArrays
  * This class provides physical and logical access.
*/
template<typename T>
class MultidimArray: public coreMultidimArray<T>
{
public:
    /// @defgroup MultidimArrayConstructors Constructors
    /// @ingroup MultidimensionalArrays
    /** Empty constructor.
     * @ingroup MultidimArrayConstructors
     *
     * The empty constructor creates an array with no memory associated, size=0.
     * 
     */
    MultidimArray(): coreMultidimArray<T>()
    {
    }
    
    /** Copy constructor
     * @ingroup MultidimArrayConstructors
     *
     * The created volume is a perfect copy of the input array but with a
     * different memory assignment.
     *
     * @code
     * Matrix3D< double > V2(V1);
     * @endcode
     */
    MultidimArray(const MultidimArray<T>& V)
    {
        // TODO Check this...
        //coreMultidimArray<T>::coreInit();
        *this = V;
    }

    /** Dimension constructor (4D)
     * @ingroup VectorsConstructors
     *
     * The dimension constructor creates a vector of the given size with memory associated
     *
     * @code
     * MultidimArray< double > v1(6);
     * @endcode
     */
    MultidimArray(unsigned long int Ndim, int Zdim, int Ydim, int Xdim): coreMultidimArray<T>()
    {
        coreMultidimArray<T>::resize(Ndim, Zdim, Ydim, Xdim);
    }

    /** Dimension constructor (3D)
     * @ingroup VectorsConstructors
     *
     * The dimension constructor creates a vector of the given size with memory associated
     *
     * @code
     * MultidimArray< double > v1(6);
     * @endcode
     */
    MultidimArray(int Zdim, int Ydim, int Xdim): coreMultidimArray<T>()
    {
        coreMultidimArray<T>::resize(1, Zdim, Ydim, Xdim);
    }

    /** Dimension constructor (2D)
     * @ingroup VectorsConstructors
     *
     * The dimension constructor creates a vector of the given size with memory associated
     *
     * @code
     * MultidimArray< double > v1(6);
     * @endcode
     */
    MultidimArray(int Ydim, int Xdim): coreMultidimArray<T>()
    {
        coreMultidimArray<T>::resize(1, 1, Ydim, Xdim);
    }

    /** Dimension constructor (1D)
     * @ingroup VectorsConstructors
     *
     * The dimension constructor creates a vector of the given size with memory associated
     *
     * @code
     * MultidimArray< double > v1(6);
     * @endcode
     */
    MultidimArray(int Xdim): coreMultidimArray<T>()
    {
        coreMultidimArray<T>::resize(1, 1, 1, Xdim);
    }

    /** Destructor.
     * @ingroup MultidimArrayConstructors
     */
     ~MultidimArray()
     {
        coreMultidimArray<T>::coreDeallocate();
     }

    /** Clear.
     * @ingroup MultidimArrayConstructors
     */
     void clear()
     {
         coreMultidimArray<T>::clear();
     }

    /** Alias a multidimarray.
     * @ingroup MultidimSize
     *
     * Treat the multidimarray as if it were a volume. The data is not copied
     * into new memory, but a pointer to the multidimarray is copied.
     * You should not make any operation on this volume such that the
     * memory locations are changed
     */
     void alias(const MultidimArray<T> &m)
     {
         coreMultidimArray<T>::copyShape(m);
         this->data=m.data;
         this->destroyData=false;
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
    void resize(const MultidimArray<T1> &v)
    {
        coreMultidimArray<T>::resize(v);
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
    void resize(unsigned long int Ndim, int Zdim, int Ydim, int Xdim)
    {
        coreMultidimArray<T>::resize(Ndim, Zdim, Ydim, Xdim);
    }
  
    /** Same shape.
     * @ingroup MultidimSize
     *
     * Returns true if this object has got the same shape (origin and size)
     * than the argument
     */
    template <typename T1>
    bool sameShape(const MultidimArray<T1>& op) const
    {
        return coreMultidimArray<T>::sameShape(op);
    }

    /** Initialize to zeros following a pattern.
     * @ingroup Initialization
     *
     * All values are set to 0, and the origin and size of the pattern are
     * adopted.
     *
     * @code
     * v2.initZeros(v1);
     * @endcode
     */
    template <typename T1>
    void initZeros(const MultidimArray<T1>& op)
    {
        coreMultidimArray<T>::initZeros(op);
    }

    /** Equality.
     * @ingroup Operators
     *
     * Returns true if this object has got the same shape (origin and size)
     * than the argument and the same values (within accuracy).
     */
    bool equal(const MultidimArray<T>& op,
    	double accuracy = XMIPP_EQUAL_ACCURACY) const
    {
        return coreMultidimArray<T>::equal(op, accuracy);
    }


//XXXXXXXXXXXXXXXXXXX UNTIL HERE CORE

    /** @defgroup MultidimSize MultidimArray Size
     * @ingroup MultidimensionalArrays
     */

    /** Put a 3D window to the nth volume
     * @ingroup MultidimSize
     *
     * The volume is windowed within the two positions given to this function.
     * Indexes always refer to logical indexes. If a position is outside the
     * actual matrix range then the matrix is padded init_value until the
     * new position is reached. In the following example suppose that m1
     * is the following and that the origin is (-1,-1,-1).
     *
     * @code
     * slice 0
     * [01 02 03          [
     *  04 05 06           04 05 06 0
     *  07 08 09]          07 08 09 0]
     *
     * ----->
     *
     * slice 1
     * [11 12 13          [
     *  14 15 16           14 15 16 0
     *  17 18 19]          17 18 19 0]
     * @endcode
     *
     * @code
     * V1.window(0, 0, -1, 1, 1, 2);
     * @endcode
     */
    void window(int z0, int y0, int x0, int zF, int yF, int xF, T init_value = 0, unsigned long n = 0)
    {
        MultidimArray<T> result(1, zF - z0 + 1, yF - y0 + 1, xF - x0 + 1);
        result.zinit = z0;
        result.yinit = y0;
        result.xinit = x0;

        for (int k = z0; k <= zF; k++)
            for (int i = y0; i <= yF; i++)
                for (int j = x0; j <= xF; j++)
                    if ((k >= STARTINGZ(*this) && k <= FINISHINGZ(*this)) &&
                        (i >= STARTINGY(*this) && i <= FINISHINGY(*this)) &&
                        (j >= STARTINGX(*this) && j <= FINISHINGX(*this)))
                        VOL_ELEM(result, k, i, j) = NZYX_ELEM(*this, n, k, i, j);
                    else
                        VOL_ELEM(result, k, i, j) = init_value;

        *this = result;
    }

    /** Put a 2D window to the nth matrix
     * @ingroup MultidimSize
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
    void window(int y0, int x0, int yF, int xF, T init_value = 0, unsigned long n = 0)
    {
        MultidimArray<T> result(1, 1, yF - y0 + 1, xF - x0 + 1);
        STARTINGY(result) = y0;
        STARTINGX(result) = x0;

        FOR_ALL_ELEMENTS_IN_MATRIX2D(result)
        if (j >= STARTINGX(*this) && j <= FINISHINGX(*this) &&
            i >= STARTINGY(*this) && i <= FINISHINGY(*this))
            MAT_ELEM(result, i, j) = NZYX_ELEM(*this, n, 0, i, j);
        else
            MAT_ELEM(result, i, j) = init_value;

        *this = result;
    }

    /** Put a 1D window to the nth vector
     * @ingroup MultidimSize
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
    void window(int x0, int xF, T init_value = 0, unsigned long n = 0)
    {
        MultidimArray<T> result(xF - x0 + 1);
        STARTINGX(result) = x0;

        for (int j = x0; j <= xF; j++)
            if (j >= STARTINGX(*this) && j <= FINISHINGX(*this))
                VEC_ELEM(result, j) = NZYX_ELEM(*this, n, 0, 0, j);
            else
                VEC_ELEM(result, j) = init_value;

        *this = result;
    }

    /** Print shape of multidimensional array.
     * @ingroup MultidimSize
     *
     * This function shows the size, starting and finishing indexes of the
     * given array. No end of line is printed neither at the beginning nor
     * the end.
     *
     * @code
     * v.printShape();
     *
     * std::ofstream fh;
     * ...;
     * v.printShape(fh);
     * @endcode
     */
    void printShape(std::ostream& out = std::cout) const
    {
        out << " Number of images = "<<NSIZE(*this)
            << " Size(Z,Y,X): " << ZSIZE(*this) << "x" << YSIZE(*this) << "x"
            << XSIZE(*this)
            << " k=[" << STARTINGZ(*this) << ".." << FINISHINGZ(*this) << "]"
            << " i=[" << STARTINGY(*this) << ".." << FINISHINGY(*this) << "]"
            << " j=[" << STARTINGX(*this) << ".." << FINISHINGX(*this) << "]";
    }

    /** Outside for 3D matrices
     * @ingroup MultidimSize
     *
     * TRUE if the logical index given is outside the definition region of this
     * array.
     */
    bool outside(int k, int i, int j) const
    {
        return (j < STARTINGX(*this) || j > FINISHINGX(*this) ||
                i < STARTINGY(*this) || i > FINISHINGY(*this) ||
                k < STARTINGZ(*this) || k > FINISHINGZ(*this));
    }

    /** Outside for 2D matrices
     * @ingroup MultidimSize
     *
     * TRUE if the logical index given is outside the definition region of this
     * array.
     */
    bool outside(int i, int j) const
    {
        return (j < STARTINGX(*this) || j > FINISHINGX(*this) ||
                i < STARTINGY(*this) || i > FINISHINGY(*this));
    }

    /** Outside for 1D matrices
     * @ingroup MultidimSize
     *
     * TRUE if the logical index given is outside the definition region of this
     * array.
     */
    bool outside(int i) const
    {
        return (i < STARTINGX(*this) || i > FINISHINGX(*this));
    }

    /** Outside
     * @ingroup MultidimSize
     *
     * TRUE if the logical index given is outside the definition region of this
     * array.
     */
    bool outside(const Matrix1D<double> &r) const
    {
        if (XSIZE(r) < 1)
        {
            REPORT_ERROR(1, "Outside: index vector has not got enough components");
        }
        else if (XSIZE(r)==1)
        {    
            return (XX(r) < STARTINGX(*this) || XX(r) > FINISHINGX(*this));
        }
        else if (XSIZE(r)==2)
        {
            return (XX(r) < STARTINGX(*this) || XX(r) > FINISHINGX(*this) ||
                    YY(r) < STARTINGY(*this) || YY(r) > FINISHINGY(*this));
        }
        else if (XSIZE(r)==3)
        {
            return (XX(r) < STARTINGX(*this) || XX(r) > FINISHINGX(*this) ||
                    YY(r) < STARTINGY(*this) || YY(r) > FINISHINGY(*this) ||
                    ZZ(r) < STARTINGZ(*this) || ZZ(r) > FINISHINGZ(*this));
        }
        else
            REPORT_ERROR(2,"Outside: index vector has too many components");
    }

    /** IsCorner (in 2D or 3D matrix)
     * @ingroup MultidimSize
     *
     * TRUE if the logical index given is a corner of the definition region of this
     * array.
     */
    bool isCorner(const Matrix1D< double >& v) const
    {
        if (XSIZE(v) < 2)
            REPORT_ERROR(1, "isCorner: index vector has got not enough components");

        else if (XSIZE(v)==2)
            return ((XX(v) == STARTINGX(*this)  && YY(v) == STARTINGY(*this))  ||
                    (XX(v) == STARTINGX(*this)  && YY(v) == FINISHINGY(*this)) ||
                    (XX(v) == FINISHINGX(*this) && YY(v) == STARTINGY(*this))  ||
                    (XX(v) == FINISHINGX(*this) && YY(v) == FINISHINGY(*this)));
        else if (XSIZE(v)==3)
            return ((XX(v) == STARTINGX(*this)  && YY(v) == STARTINGY(*this) && ZZ(v) == STARTINGZ(*this)) ||
                    (XX(v) == STARTINGX(*this)  && YY(v) == FINISHINGY(*this) && ZZ(v) == STARTINGZ(*this)) ||
                    (XX(v) == FINISHINGX(*this) && YY(v) == STARTINGY(*this) && ZZ(v) == STARTINGZ(*this))  ||
                    (XX(v) == FINISHINGX(*this) && YY(v) == FINISHINGY(*this) && ZZ(v) == STARTINGZ(*this)) ||
                    (XX(v) == STARTINGX(*this)  && YY(v) == STARTINGY(*this) && ZZ(v) == FINISHINGZ(*this)) ||
                    (XX(v) == STARTINGX(*this)  && YY(v) == FINISHINGY(*this) && ZZ(v) == FINISHINGZ(*this)) ||
                    (XX(v) == FINISHINGX(*this) && YY(v) == STARTINGY(*this) && ZZ(v) == FINISHINGZ(*this))  ||
                    (XX(v) == FINISHINGX(*this) && YY(v) == FINISHINGY(*this) && ZZ(v) == FINISHINGZ(*this)));
        else
            REPORT_ERROR(1, "isCorner: index vector has too many components");

    }

    /** @defgroup MultidimMemory Access to the pixel values
     * @ingroup MultidimensionalArrays
     */

    /** Volume element access via double vector.
     * @ingroup MultidimMemory
     *
     * Returns the value of a matrix logical position, but this time the
     * element position is determined by a R3 vector. The elements can be used
     * either by value or by reference. An exception is thrown if the index is
     * outside the logical range. Pay attention in the following example that
     * we are accessing the same element as in the previous function but, now
     * we have to give first the X position because we are building first a
     * vector of the form (x,y,z).
     *
     * @code
     * V(vectorR3(1, -2, 0)) = 1;
     * val = V(vectorR3(1, -2, 0));
     * @endcode
     */
    T& operator()(const Matrix1D< double >& v, unsigned long n = 0) const
    {
        v.resize(3);
        return NZYX_ELEM((*this), n, ROUND(ZZ(v)), ROUND(YY(v)), ROUND(XX(v)));
    }

    /** Volume element access via integer vector.
     * @ingroup MultidimMemory
     */
    T& operator()(const Matrix1D< int >& v, unsigned long n = 0) const
    {
        v.resize(3);
        return NZYX_ELEM((*this), n, ZZ(v), YY(v), XX(v));
    }

    /// @defgroup MultidimInit Initialization
    /// @ingroup MultidimensionalArrays

    /** Initialize with random values.
     * @ingroup MultidimInit
     *
     * This function allows you to initialize the array with a set of random
     * values picked from a uniform random distribution or a gaussian one. You
     * must choose two parameters for each, for the uniform distribution they
     * mean the range where to generate the random numbers, while in the
     * gaussian case they are the mean and the standard deviation. By default
     * the uniform distribution is selected. The size and origin of the array
     * are not modified.
     *
     * @code
     * v.initRandom(0, 1);
     * // uniform distribution between 0 and 1
     *
     * v.initRandom(0, 1, "uniform");
     * // the same
     *
     * v.initRandom(0, 1, "gaussian");
     * // gaussian distribution with 0 mean and stddev=1
     * @endcode
     */
    void initRandom(double op1, double op2, const std::string& mode = "uniform")
    {
        T* ptr=NULL;
	unsigned long int n;
        if (mode == "uniform")
            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
                *ptr = static_cast< T >(rnd_unif(op1, op2));
        else if (mode == "gaussian")
            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
                *ptr = static_cast< T >(rnd_gaus(op1, op2));
        else
            REPORT_ERROR(1005,
                         static_cast< std::string >("InitRandom: Mode not supported (" +
                                                    mode + ")"));
    }


    /** Linear initialization
     * @ingroup MultidimInit
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

        if ( YSIZE(*this) != 1 || YSIZE(*this) != 1)
            REPORT_ERROR(1,"MultidimArray ERROR: initLinear is only meant for 1D multidimArrays");

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

    /** Interpolates the value of the nth 3D matrix M at the point (x,y,z).
     * @ingroup MultidimMemory
     *
     * (x,y,z) are in logical coordinates.
     */
    T interpolatedElement3D(double x, double y, double z, T outside_value = (T) 0, int n = 0)
    {
        int x0 = FLOOR(x);
        double fx = x - x0;
        int x1 = x0 + 1;

        int y0 = FLOOR(y);
        double fy = y - y0;
        int y1 = y0 + 1;

        int z0 = FLOOR(z);
        double fz = z - z0;
        int z1 = z0 + 1;

        T d000 = (outside(z0, y0, x0)) ? outside_value : NZYX_ELEM(*this, n, z0, y0, x0);
        T d001 = (outside(z0, y0, x1)) ? outside_value : NZYX_ELEM(*this, n, z0, y0, x1);
        T d010 = (outside(z0, y1, x0)) ? outside_value : NZYX_ELEM(*this, n, z0, y1, x0);
        T d011 = (outside(z0, y1, x1)) ? outside_value : NZYX_ELEM(*this, n, z0, y1, x1);
        T d100 = (outside(z1, y0, x0)) ? outside_value : NZYX_ELEM(*this, n, z1, y0, x0);
        T d101 = (outside(z1, y0, x1)) ? outside_value : NZYX_ELEM(*this, n, z1, y0, x1);
        T d110 = (outside(z1, y1, x0)) ? outside_value : NZYX_ELEM(*this, n, z1, y1, x0);
        T d111 = (outside(z1, y1, x1)) ? outside_value : NZYX_ELEM(*this, n, z1, y1, x1);

        double dx00 = LIN_INTERP(fx, (double) d000, (double) d001);
        double dx01 = LIN_INTERP(fx, (double) d100, (double) d101);
        double dx10 = LIN_INTERP(fx, (double) d010, (double) d011);
        double dx11 = LIN_INTERP(fx, (double) d110, (double) d111);
        double dxy0 = LIN_INTERP(fy, (double) dx00, (double) dx10);
        double dxy1 = LIN_INTERP(fy, (double) dx01, (double) dx11);

        return (T) LIN_INTERP(fz, dxy0, dxy1);
    }

    /** Interpolates the value of the nth 2D matrix M at the point (x,y)
     * @ingroup MultidimMemory
     *
     * Bilinear interpolation. (x,y) are in logical coordinates.
     */
    inline T interpolatedElement2D(double x, double y, T outside_value = (T) 0, int n = 0) const
    {
        int x0 = FLOOR(x);
        double fx = x - x0;
        int x1 = x0 + 1;
        int y0 = FLOOR(y);
        double fy = y - y0;
        int y1 = y0 + 1;

        T d00 = outside(y0, x0) ? outside_value : NZYX_ELEM(*this, n, 0, y0, x0);
        T d10 = outside(y1, x0) ? outside_value : NZYX_ELEM(*this, n, 0, y1, x0);
        T d11 = outside(y1, x1) ? outside_value : NZYX_ELEM(*this, n, 0, y1, x1);
        T d01 = outside(y0, x1) ? outside_value : NZYX_ELEM(*this, n, 0, y0, x1);

        double d0 = (T) LIN_INTERP(fx, (double) d00, (double) d01);
        double d1 = (T) LIN_INTERP(fx, (double) d10, (double) d11);
        return (T) LIN_INTERP(fy, d0, d1);
    }

    /** Interpolates the value of the nth 3D matrix M at the point (x,y,z) knowing
     * that this image is a set of B-spline coefficients.
     * @ingroup MultidimMemory
     *
     * (x,y,z) are in logical coordinates.
     */
    T interpolatedElementBSpline3D(double x, double y, double z, 
                                   int SplineDegree = 3, unsigned long n = 0)
    {
        int SplineDegree_1 = SplineDegree - 1;

        // Logical to physical
        z -= STARTINGZ(*this);
        y -= STARTINGY(*this);
        x -= STARTINGX(*this);

        int lmax = XSIZE(*this);
        int mmax = YSIZE(*this);
        int nmax = ZSIZE(*this);

        int l1 = CEIL(x - SplineDegree_1);
        int l2 = l1 + SplineDegree;

        int m1 = CEIL(y - SplineDegree_1);
        int m2 = m1 + SplineDegree;

        int n1 = CEIL(z - SplineDegree_1);
        int n2 = n1 + SplineDegree;

        double zyxsum = 0.0;
        for (int nn = n1; nn <= n2; nn++) {
	    int equivalent_n=n;
	    if      (nn<0)             equivalent_n=-nn-1;
	    else if (nn>=ZSIZE(*this)) equivalent_n=2*ZSIZE(*this)-nn-1;
            double yxsum = 0.0;
            for (int m = m1; m <= m2; m++) {
		int equivalent_m=m;
		if      (m<0)             equivalent_m=-m-1;
		else if (m>=YSIZE(*this)) equivalent_m=2*YSIZE(*this)-m-1;
                double xsum = 0.0;
                for (int l = l1; l <= l2; l++)
                {
                    double xminusl = x - (double) l;
		    int equivalent_l=l;
		    if      (l<0)             equivalent_l=-l-1;
		    else if (l>=XSIZE(*this)) equivalent_l=2*XSIZE(*this)-l-1;
                    double Coeff = (double) DIRECT_NZYX_ELEM(*this, n, 
                        equivalent_n,equivalent_m,equivalent_l);
                    switch (SplineDegree)
                    {
                    case 2:
                        xsum += Coeff * Bspline02(xminusl);
                        break;
                    case 3:
                        xsum += Coeff * Bspline03(xminusl);
                        break;
                    case 4:
                        xsum += Coeff * Bspline04(xminusl);
                        break;
                    case 5:
                        xsum += Coeff * Bspline05(xminusl);
                        break;
                    case 6:
                        xsum += Coeff * Bspline06(xminusl);
                        break;
                    case 7:
                        xsum += Coeff * Bspline07(xminusl);
                        break;
                    case 8:
                        xsum += Coeff * Bspline08(xminusl);
                        break;
                    case 9:
                        xsum += Coeff * Bspline09(xminusl);
                        break;
                    }
                }

                double yminusm = y - (double) m;
                switch (SplineDegree)
                {
                case 2:
                    yxsum += xsum * Bspline02(yminusm);
                    break;
                case 3:
                    yxsum += xsum * Bspline03(yminusm);
                    break;
                case 4:
                    yxsum += xsum * Bspline04(yminusm);
                    break;
                case 5:
                    yxsum += xsum * Bspline05(yminusm);
                    break;
                case 6:
                    yxsum += xsum * Bspline06(yminusm);
                    break;
                case 7:
                    yxsum += xsum * Bspline07(yminusm);
                    break;
                case 8:
                    yxsum += xsum * Bspline08(yminusm);
                    break;
                case 9:
                    yxsum += xsum * Bspline09(yminusm);
                    break;
                }
	    }

            double zminusn = z - (double) n;
            switch (SplineDegree)
            {
            case 2:
                zyxsum += yxsum * Bspline02(zminusn);
                break;
            case 3:
                zyxsum += yxsum * Bspline03(zminusn);
                break;
            case 4:
                zyxsum += yxsum * Bspline04(zminusn);
                break;
            case 5:
                zyxsum += yxsum * Bspline05(zminusn);
                break;
            case 6:
                zyxsum += yxsum * Bspline06(zminusn);
                break;
            case 7:
                zyxsum += yxsum * Bspline07(zminusn);
                break;
            case 8:
                zyxsum += yxsum * Bspline08(zminusn);
                break;
            case 9:
                zyxsum += yxsum * Bspline09(zminusn);
                break;
            }
	}

        return (T) zyxsum;
    }

    /** Interpolates the value of the nth 2D matrix M at the point (x,y) knowing
     * that this image is a set of B-spline coefficients
     * @ingroup MultidimMemory
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
    inline T interpolatedElementBSpline2D(double x, double y, int SplineDegree = 3, 
                                          unsigned long n = 0) const
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
                double Coeff = (double) DIRECT_NZYX_ELEM(*this, n, 0, equivalent_m,equivalent_l);
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

    /** Interpolates the value of the nth 1D vector M at the point (x) knowing
     * that this vector is a set of B-spline coefficients
     * @ingroup MultidimMemory
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
    T interpolatedElementBSpline1D(double x, int SplineDegree = 3, unsigned long n = 0) const
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
            double Coeff = (double) DIRECT_NZYX_ELEM(*this, n, 0, 0, equivalent_l);
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


    /// @defgroup Statistics Statistics functions
    /// @ingroup MultidimensionalArrays
    /** Compute statistics within 2D region of 2D image.
     * @ingroup Statistics
     *
     * The 2D region is specified by two corners.
     * Note that this function only works for the 0th image in a multi-image array...
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

    /** Adjust the range of the array to a given one.
     * @ingroup Statistics
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
        if (NZYXSIZE(*this) <= 0)
            return;

        double min0, max0;
        coreMultidimArray<T>::computeDoubleMinMax(min0, max0);

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

    /** Adjust the range of the array to a given one within a mask.
     * @ingroup Statistics
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
    void rangeAdjust(T minF, T maxF, MultidimArray<int> &mask)
    {
        if (MULTIDIM_SIZE(*this) <= 0)
            return;

        double min0, max0;
        bool first=true;
        T* ptr=NULL;
        unsigned long int n;
        int * ptrMask=MULTIDIM_ARRAY(mask);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
        {
            if (*ptrMask)
            {
                T val= *ptr;
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
            ptrMask++;
       }

        // If max0==min0, it means that the vector is a constant one, so the
        // only possible transformation is to a fixed minF
        double slope;
        if (max0 != min0)
            slope = static_cast< double >(maxF - minF) /
                    static_cast< double >(max0 - min0);
        else
            slope = 0;

        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
            *ptr = minF + static_cast< T >(slope *
                static_cast< double >(*ptr - min0));
    }

    /** Adjust the range of the array to the range of another array in
        a least squares sense.
     * @ingroup Statistics
     *
     * A linear operation is performed on the values of the array such that
     * after it, the values of the self array are as similar as possible
     * (L2 sense) to the values of the array shown as sample
     */
    void rangeAdjust(const MultidimArray<T> &example,
        const MultidimArray<int> *mask=NULL)
    {
        if (NZYXSIZE(*this) <= 0)
            return;
        
        // y=a+bx
        double sumx=0, sumy=0, sumxy=0, sumx2=0;
        double* ptrExample=MULTIDIM_ARRAY(example);
        int* ptrMask=NULL;
        if (mask!=NULL) ptrMask=MULTIDIM_ARRAY(*mask);
        double* ptr=NULL;
        unsigned long int n;
        double N=0;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
        {
            bool process=true;
            if (mask!=NULL)
                if (*ptrMask==0) process=false;
            if (process)
            {
                double x=*ptr;
                double y=*ptrExample;
                sumy+=y;
                sumxy+=x*y;
                sumx+=x;
                sumx2+=x*x;
                N++;
            }
            ptrExample++;
            if (mask!=NULL) ptrMask++;
        }
        double a=(sumy*sumx2-sumx*sumxy)/(N*sumx2-sumx*sumx);
        double b=(N*sumxy-sumx*sumy)/(N*sumx2-sumx*sumx);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
            *ptr = static_cast< double >(a+b * static_cast< double > (*ptr));
    }

    /** Adjust the average and stddev of the array to given values.
     * @ingroup Statistics
     *
     * A linear operation is performed on the values of the array such
     * that after it, the average and standard deviation of the array
     * are the two values set. The actual array is modified itself
     *
     * @code
     * v.statisticsAdjust(0,1);
     * // The array has got now 0 mean and stddev=1
     * @endcode
     */
    // This function must be explictly implemented outside.
    void statisticsAdjust(double avgF, double stddevF)
    {
        double avg0, stddev0;
        double a, b;

        if (NZYXSIZE(*this) == 0)
            return;

        T minval, maxval;
        coreMultidimArray<T>::computeStats(avg0, stddev0, minval, maxval);

        if (stddev0 != 0)
            a = static_cast< double >(stddevF) / static_cast< double >(stddev0);
        else
            a = 0;

        b = static_cast< double >(avgF) - a * static_cast< double >(avg0);

        T* ptr=NULL;
	unsigned long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
            *ptr = static_cast< T >(a * static_cast< double > (*ptr) + b);
    }


    /** @defgroup MultidimUtilities Utilities
     *  @ingroup MultidimensionalArrays
     *
     * Here you have several easy functions to manage the values of
     * the array.
     */

    /** Add noise to actual values.
     * @ingroup MultidimUtilities
     *
     * This function add some noise to the actual values of the array according
     * to a certain random distribution. You must choose two parameters for
     * each, for the uniform distribution they mean the range where to generate
     * the random numbers, while in the gaussian case they are the mean and the
     * standard deviation. By default the uniform distribution is selected. The
     * size and origin of the array are not modified. The array itself is
     * modified.
     *
     * @code
     * v1.addNoise(0, 1);
     * // uniform distribution between 0 and 1
     *
     * v1.addNoise(0, 1, "uniform");
     * // the same
     *
     * v1.addNoise(0, 1, "gaussian");
     * // gaussian distribution with 0 mean and stddev=1
     *
     * v1.addNoise(0, 1, "student", 3);
     * // t-student distribution with 0 mean and stddev=1, and 3 degrees of freedom
     *

     * @endcode
     */
    void addNoise(double op1,
                  double op2,
                  const std::string& mode = "uniform",
                  double df = 3) const
    {
        T* ptr=NULL;
	unsigned long int n;
        if (mode == "uniform")
            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
                *ptr += static_cast< T >(rnd_unif(op1, op2));
        else if (mode == "gaussian")
            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
                *ptr += static_cast< T >(rnd_gaus(op1, op2));
        else if (mode == "student")
            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
                *ptr += static_cast< T >(rnd_student_t(df, op1, op2));
        else
            REPORT_ERROR(1005,
                         static_cast< std::string >("AddNoise: Mode not supported (" +
                                                    mode + ")"));
    }

    /** Computes the center of mass of the nth array
     * @ingroup MultidimUtilities
     */
    void centerOfMass(Matrix1D< double >& center, void * mask=NULL, unsigned long n = 0)
    {
	center.initZeros(3);
	double mass = 0;
	MultidimArray< int >* imask = (MultidimArray< int >*) mask;

	FOR_ALL_ELEMENTS_IN_MATRIX3D(*this)
	{
            if ((imask == NULL || NZYX_ELEM(*imask, n, k, i, j)) &&
		VOL_ELEM(*this, k, i, j) > 0)
            {
        	XX(center) += j * NZYX_ELEM(*this, n, k, i, j);
        	YY(center) += i * NZYX_ELEM(*this, n, k, i, j);
        	ZZ(center) += k * NZYX_ELEM(*this, n, k, i, j);

        	mass += NZYX_ELEM(*this, n, k, i, j);
            }
	}

	if (mass != 0)
            center /= mass;
    }

    /** Several thresholding.
     * @ingroup MultidimUtilities
     *
     * Apply a threshold to the array, the object is modified itself. There
     * are several kinds of thresholding and you must specify it, the values
     * given in the fuction have different meanings according to the threshold
     * applied.
     *
     * abs_above: if |x|>a => x=b
     * abs_below: if |x|<a => x=b
     * above:     if  x >a => x=b
     * below:     if  x <a => x=b
     * range:     if  x <a => x=a   and    if x>b => x=b
     *
     * @code
     * v.threshold("abs_above", 10, 10);
     * // any value whose absolute value is above 10 will be substituted by
     * // -10 (if it is negative) or 10 (if it is positive)
     *
     * v.threshold("abs_below", 0.1, 0);
     * // any value whose absolute value is below 0.1 will be substituted by
     * // -0 (if it is negative) or 0 (if it is positive)
     *
     * v.threshold("above", 10, 10);
     * // any value above 10 will be substituted by 10
     *
     * v.threshold("below", -10, -10);
     * // any value below -10 will be substituted by -10
     *
     * v.threshold("range", 0, 1);
     * // v is "saturated" by values 0 and 1, any value outside this range
     * // will be substituted by its nearest border
     * @endcode
     */
    void threshold(const std::string& type, 
	           T a, 
	           T b,
	           MultidimArray<int> * mask = NULL )
    {
        int mode;

        if (type == "abs_above")
            mode = 1;
        else if (type == "abs_below")
            mode = 2;
        else if (type == "above")
            mode = 3;
        else if (type == "below")
            mode = 4;
        else if (type == "range")
            mode = 5;
        else
            REPORT_ERROR(1005,
                         static_cast< std::string >("Threshold: mode not supported (" +
                                                    type + ")"));

        T* ptr=NULL;
	unsigned long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
        {
	    if (mask == NULL || DIRECT_MULTIDIM_ELEM(*mask,n) > 0 )
	    {
	        switch (mode)
	        {
	        case 1:
		    if (ABS(*ptr) > a)
		        *ptr = SGN(*ptr) * b;
		    break;
	        case 2:
		    if (ABS(*ptr) < a)
		        *ptr = SGN(*ptr) * b;
		    break;
	        case 3:
		    if (*ptr > a)
		        *ptr = b;
		    break;
	        case 4:
		    if (*ptr < a)
		        *ptr = b;
		    break;
	        case 5:
		    if (*ptr < a)
		        *ptr = a;
		    else if (*ptr > b)
		        *ptr = b;
		    break;
	        }
	    }
        }
    }

    /** Count with threshold.
     * @ingroup MultidimUtilities
     *
     * This function returns the number of elements meeting the threshold
     * condition.
     */
    unsigned long countThreshold(const std::string& type, 
		         T a, 
		         T b,
		         MultidimArray<int> * mask = NULL )
    {
        int mode;

        if (type == "abs_above")
            mode = 1;
        else if (type == "abs_below")
            mode = 2;
        else if (type == "above")
            mode = 3;
        else if (type == "below")
            mode = 4;
        else if (type == "range")
            mode = 5;
        else
            REPORT_ERROR(1005,
                         static_cast< std::string >("CountThreshold: mode not supported (" +
                                                    type + ")"));

        unsigned long ret = 0;

        T* ptr=NULL;
	unsigned long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
	    if (mask == NULL || DIRECT_MULTIDIM_ELEM(*mask,n) > 0 )
	    {
	        switch (mode)
	        {
	        case 1:
		    if (ABS(*ptr) > a)
		        ret++;
		    break;
	        case 2:
		    if (ABS(*ptr) < a)
		        ret++;
		    break;
	        case 3:
		    if (*ptr > a)
		        ret++;
		    break;
	        case 4:
		    if (*ptr < a)
		        ret++;
		    break;
	        case 5:
		    if (*ptr >= a && *ptr <= b)
		        ret++;
		    break;
	        }
	    }
        return ret;
    }

    /** Substitute a value by another.
     * @ingroup MultidimUtilities
     *
     * Substitute an old value by a new one. The accuracy is used to say if
     * the value in the array is equal to the old value. Set it to 0 for
     * perfect accuracy.
     */
    void substitute(T oldv,
                    T newv,
                    double accuracy = XMIPP_EQUAL_ACCURACY,
		    MultidimArray<int> * mask = NULL )
    {
        T* ptr=NULL;
	unsigned long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
	    if (mask == NULL || DIRECT_MULTIDIM_ELEM(*mask,n) > 0 )
	        if (ABS(*ptr - oldv) <= accuracy)
		    *ptr = newv;

    }

    /** Substitute a given value by a sample from a Gaussian distribution.
     * @ingroup MultidimUtilities
     *
     * Substitute  a given value by a sample from a Gaussian distribution.
     * The accuracy is used to say if the value in the array is equal 
     * to the old value.  Set it to 0 for perfect accuracy.
     */
    void randomSubstitute(T oldv,
		           T avgv,
		           T sigv,
		           double accuracy = XMIPP_EQUAL_ACCURACY,
		           MultidimArray<int> * mask = NULL )
    {
        T* ptr=NULL;
	unsigned long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
	    if (mask == NULL || DIRECT_MULTIDIM_ELEM(*mask,n) > 0 )
	        if (ABS(*ptr - oldv) <= accuracy)
		    *ptr = rnd_gaus(avgv, sigv);

    }

    /** Binarize.
     * @ingroup MultidimUtilities
     *
     * This functions substitutes all values in a volume which are greater
     * than val+accuracy by 1 and the rest are set to 0. Use threshold to get a
     * very powerful binarization.
     */
    void binarize(double val = 0, 
	          double accuracy = XMIPP_EQUAL_ACCURACY,
	          MultidimArray<int> * mask = NULL )
    {
        T* ptr=NULL;
	unsigned long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
	    if (mask == NULL || DIRECT_MULTIDIM_ELEM(*mask,n) > 0 )
	        if (*ptr <= val + accuracy)
		    *ptr = 0;
	        else
		    *ptr = 1;

    }

    /** Extracts the profile between two points
     * @ingroup MultidimUtilities
     *
     * Given two logical indexes, this function returns samples of the line that
     * joins them. This is done by bilinear interpolation. The number of samples
     * in the line is N.
     */
    void profile2D(int x0, int y0, int xF, int yF, int N,
                   Matrix1D< T >& profile) const
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


    /** Input from input stream.
     * @ingroup MultidimUtilities
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
    friend std::istream& operator>>(std::istream& in, MultidimArray<T>& v)
    {
        T* ptr;
	unsigned long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(v,n,ptr)
            in >> *ptr;
        return in;
    }


};

/// @defgroup MultidimFunctions Functions for all multidimensional arrays
/// @ingroup MultidimensionalArrays

/** Conversion from one type to another.
 * @ingroup MultidimFunctions
 *
 * If we have an integer array and we need a double one, we can use this
 * function. The conversion is done through a type casting of each element
 * If n >= 0, only the nth volumes will be converted, otherwise all NSIZE volumes
 */
template<typename T1, typename T2>
    void typeCast(const MultidimArray<T1>& v1, MultidimArray<T2>& v2, long n = -1)
{
    typecast(v1, v2, n);
    if (NZYXSIZE(v1) == 0)
    {
        v2.clear();
        return;
    }

    if (n < 0)
    {
        v2.resize(v1);
        T1* ptr1=NULL;
        unsigned long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(v1,n,ptr1)
            DIRECT_MULTIDIM_ELEM(v2,n) = static_cast< T2 >(*ptr1);
    }
    else
    {
        v2.resize(ZSIZE(v1),YSIZE(v1),XSIZE(v1));
        FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(v2)
            DIRECT_VOL_ELEM(v2,k,i,j) = static_cast< T2 >DIRECT_NZYX_ELEM(v1,n,k,i,j);
    }

}


/** MultidimArray equality.
 * @ingroup MultidimMisc */
template<typename T>
bool operator==(const MultidimArray<T>& op1, const MultidimArray<T>& op2)
{
    return op1.equal(op2);
}

/** MultidimArray inequality.
 * @ingroup MultidimMisc */
template<typename T>
bool operator!=(const MultidimArray<T>& op1, const MultidimArray<T>& op2)
{
    return !(op1==op2);
}


//XXXXXXXXXXXXXXX FROM HERE ON TO MULTIDIM_ARRAY

/** Reduce both volumes to a common size.
 * @ingroup MultidimMisc
 *
 * Search the range of logical indexes for which both volumes have got valid
 * values, and cut both to that size, the corresponding origin is automatically
 * computed.
 *
 * @code
 * Matrix3D< double > V1(4, 5, 3);
 * V1.startingX() = -2;
 * V1.startingY() = -2;
 * V1.startingZ() = -2;
 *
 * Matrix3D< double > V2(4, 2, 3);
 * V2.startingX() = 0;
 * V2.startingY() = 0;
 * V2.startingZ() = 0;
 *
 * // V1 and V2 range from (0,0,0)=(z,y,x) to (1,1,0)
 * cutToCommonSize(V1, V2);
 * @endcode
 */
template<typename T>
void cutToCommonSize(MultidimArray<T>& V1, MultidimArray<T>& V2)
{
    int z0 = XMIPP_MAX(STARTINGZ(V1), STARTINGZ(V2));
    int zF = XMIPP_MIN(FINISHINGZ(V1), FINISHINGZ(V2));
    int y0 = XMIPP_MAX(STARTINGY(V1), STARTINGY(V2));
    int yF = XMIPP_MIN(FINISHINGY(V1), FINISHINGY(V2));
    int x0 = XMIPP_MAX(STARTINGX(V1), STARTINGX(V2));
    int xF = XMIPP_MIN(FINISHINGX(V1), FINISHINGX(V2));

    V1.window(z0, y0, x0, zF, yF, xF);
    V2.window(z0, y0, x0, zF, yF, xF);
}

#endif
