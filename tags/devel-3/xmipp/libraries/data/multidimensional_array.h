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

#ifndef MULTIDIMENSIONALARRAY_H
#define MULTIDIMENSIONALARRAY_H

#include <typeinfo>

#include <external/bilib/types/tsplinebasis.h>
#include <external/bilib/types/tboundaryconvention.h>
#include <external/bilib/headers/linearalgebra.h>
#include <external/bilib/headers/changebasis.h>
#include <external/bilib/headers/kernel.h>
#include "funcs.h"
#include "error.h"
#include "args.h"

extern int bestPrecision(float F, int _width);
extern std::string floatToString(float F, int _width, int _prec);
extern std::string integerToString(int I, int _width, char fill_with);

/// @defgroup MultidimensionalArrays Multidimensional Arrays
/// @ingroup Arrays

/** @defgroup MultidimArraysSpeedUp Speed up macros
 * @ingroup MultidimensionalArrays
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

/// @defgroup MultidimArraysSizeShape Size and shape
/// @ingroup MultidimArraysSpeedUp

/** Returns the first X valid logical index
 * @ingroup MultidimArraysSizeShape
 */
#define STARTINGX(v) ((v).xinit)

/** Returns the last X valid logical index
 * @ingroup MultidimArraysSizeShape
 */
#define FINISHINGX(v) ((v).xinit + (v).xdim - 1)

/** Returns the first Y valid logical index
 * @ingroup MultidimArraysSizeShape
 */
#define STARTINGY(v) ((v).yinit)

/** Returns the last Y valid logical index
 * @ingroup MultidimArraysSizeShape
 */
#define FINISHINGY(v) ((v).yinit + (v).ydim - 1)

/** Returns the first Z valid logical index
 * @ingroup MultidimArraysSizeShape
 */
#define STARTINGZ(v) ((v).zinit)

/** Returns the last Z valid logical index
 * @ingroup MultidimArraysSizeShape
 */
#define FINISHINGZ(v) ((v).zinit + (v).zdim - 1)

/** Access to X dimension (size)
 * @ingroup MultidimArraysSizeShape
 */
#define XSIZE(v) ((v).xdim)

/** Access to Y dimension (size)
 * @ingroup MultidimArraysSizeShape
 */
#define YSIZE(v) ((v).ydim)

/** Access to Z dimension (size)
 * @ingroup MultidimArraysSizeShape
 */
#define ZSIZE(v) ((v).zdim)

/** Access to Z dimension (size)
 * @ingroup MultidimArraysSizeShape
 */
#define NSIZE(v) ((v).ndim)

/** Access to XY dimension (Ysize*Xsize)
 * @ingroup MultidimArraysSizeShape
 */
#define YXSIZE(v) ((v).yxdim)

/** Access to XYZ dimension (Zsize*Ysize*Xsize)
 * @ingroup MultidimArraysSizeShape
 */
#define ZYXSIZE(v) ((v).zyxdim)

/// FIXME MULTIDIM_SIZE will disappear
/** Access to XYZ dimension (Zsize*Ysize*Xsize)
 * @ingroup MultidimArraysSizeShape
 */
#define MULTIDIM_SIZE(v) ((v).nzyxdim)

/** Access to XYZN dimension (Nsize*Zsize*Ysize*Xsize)
 * @ingroup MultidimArraysSizeShape
 */
#define NZYXSIZE(v) ((v).nzyxdim)

/** Array access.
 * @ingroup MultidimArraysSizeShape
 *
 * This macro gives you access to the array (T **)
 */
#ifndef MULTIDIM_ARRAY
#define MULTIDIM_ARRAY(v) ((v).data)
#endif

/** Access to a direct element.
 * @ingroup MultidimArraysSizeShape.
 * v is the array, k is the slice (Z), i is the Y index and j is the X index.
 */
#define DIRECT_VOL_ELEM(v,k,i,j) ((v).data[(k)*YXSIZE(v)+((i)*XSIZE(v))+(j)])

/** Volume element: Logical access.
 * @ingroup VolumesMemory
 *
 * @code
 * VOL_ELEM(V, -1, -2, 1) = 1;
 * val = VOL_ELEM(V, -1, -2, 1);
 * @endcode
 */
#define VOL_ELEM(V, k, i, j) \
    DIRECT_VOL_ELEM((V),(k) - STARTINGZ(V), (i) - STARTINGY(V), (j) - STARTINGX(V))

/** Access to a direct element.
 * @ingroup MultidimArraysSizeShape.
 * v is the array, l is the image, k is the slice, i is the Y index and j is the X index.
 * i and j) within the slice.
 */
#define DIRECT_NZYX_ELEM(v, l, k, i, j) ((v).data[(l)*ZYXSIZE(v)+(k)*YXSIZE(v)+((i)*XSIZE(v))+(j)])

/** Multidim element: Logical access.
 * @ingroup MultidimArraysSizeShape.

 */
#define NZYX_ELEM(v, l, k, i, j)  \
    DIRECT_NZYX_ELEM((V), (l), (k) - STARTINGZ(V), (i) - STARTINGY(V), (j) - STARTINGX(V))

/** Access to a direct element.
 * @ingroup MultidimArraysSizeShape.
 * v is the array, k is the slice and n is the number of the pixel (combined
 * i and j) within the slice.
 */
#define DIRECT_MULTIDIM_ELEM(v,n) ((v).data[(n)])

/** Access to a direct element of a matrix.
 * @ingroup MultidimArraysSizeShape.
 * v is the array, i and j define the element v_ij.
 *
 * Be careful because this is physical access, usually matrices follow the C
 * convention of starting index==0 (X and Y). This function should not be used
 * as it goes against the vector library philosophy unless you explicitly want
 * to access directly to any value in the matrix without taking into account its
 * logical position
 *
 * @code
 * DIRECT_MAT_ELEM(m, 0, 0) = 1;
 * val = DIRECT_MAT_ELEM(m, 0, 0);
 * @endcode

 */
#define DIRECT_MAT_ELEM(v,i,j) ((v).data[(i)*(v).xdim+(j)])

/** Short alias for DIRECT_MAT_ELEM
 * @ingroup MultidimArraysSizeShape
 */
#define dMij(M, i, j) DIRECT_MAT_ELEM(M, i, j)

/** Matrix element: Logical access
 * @ingroup MatricesMemory
 *
 * @code
 * MAT_ELEM(m, -2, 1) = 1;
 * val = MAT_ELEM(m, -2, 1);
 * @endcode
 */
#define MAT_ELEM(v, i, j) \
    DIRECT_MAT_ELEM(v, (i) - STARTINGY(v), (j) - STARTINGX(v))

/** Vector element: Physical access
 * @ingroup MultidimArraysSizeShape
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

/** Vector element: Logical access
 * @ingroup MultidimArraysSizeShape
 *
 * @code
 * VEC_ELEM(v, -2) = 1;
 * val = VEC_ELEM(v, -2);
 * @endcode
 */
#define VEC_ELEM(v, i) DIRECT_VEC_ELEM(v, (i) - ((v).xinit))

/** For all direct elements in the array
 * @ingroup MultidimArraysSizeShape
 *
 * This macro is used to generate loops for the array in an easy manner. It
 * defines an internal index 'n' which goes over the slices and 'n' that
 * goes over the pixels in each slice.
 *
 * @code
 * FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(v)
 * {
 *     std::cout << DIRECT_MULTIDIM_ELEM(v,n) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(v) \
    for (unsigned long int n=0; n<NZYXSIZE(v); ++n)

/** For all direct elements in the array, pointer version
 * @ingroup MultidimArraysSizeShape
 *
 * This macro is used to generate loops for the array in an easy manner. It
 * defines an internal index 'k' which goes over the slices and 'n' that
 * goes over the pixels in each slice. Each element can be accessed through
 * an external pointer called ptr.
 *
 * @code
 * T* ptr=NULL;
 * unsigned long int n;
 * FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(v,n,ptr)
 * {
 *     std::cout << *ptr << " ";
 * }
 * @endcode
 */
#define FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(v,n,ptr) \
    for ((n)=0, (ptr)=(v).data; (n)<NZYXSIZE(v); ++(n), ++(ptr))

// Forward declarations ====================================================
template<typename T> class MultidimArray;
template<typename T> class Matrix1D;

template<typename T>
void coreArrayByScalar(const MultidimArray<T>& op1, const T& op2,
    MultidimArray<T>& result, char operation);

template<typename T>
void coreScalarByArray(const T& op1, const MultidimArray<T>& op2,
    MultidimArray<T>& result, char operation);

template<typename T>
void coreArrayByArray(const MultidimArray<T>& op1, const MultidimArray<T>& op2,
    MultidimArray<T>& result, char operation);

/** Template class for Xmipp arrays.
  * @ingroup MultidimensionalArrays
  * This class provides physical and logical access.
*/
template<typename T>
class MultidimArray
{
public:
    /* The array itself.
       The array is always a 3D array (Z,Y,X). For vectors the size of the array
       is (1,1,X) and for matrices (1,Y,X). The pixel (i,j) (y,x) is at the
       position data[i*Xdim+j] or data[y*Xdim+x]
    */
    T* data;

    // Destroy data
    bool destroyData;

    // Number of images
    int ndim;

    // Number of elements in Z
    int zdim;

    // Number of elements in Y
    int ydim;

    // Number of elements in X
    int xdim;

    // Number of elements in YX
    unsigned long int yxdim;

    // Number of elements in ZYX
    unsigned long int zyxdim;

    // Number of elements in NZYX
    unsigned long int nzyxdim;

    // Z init
    int zinit;

    // Y init
    int yinit;

    // X init
    int xinit;
public:
    /// @defgroup MultidimArrayConstructors Constructors
    /// @ingroup MultidimensionalArrays
    /** Empty constructor.
     * @ingroup MultidimArrayConstructors
     * The empty constructor creates an array with no memory associated,
     * size=0.
     */
    MultidimArray()
    {
        coreInit();
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
        coreInit();
        *this = V;
    }

    /** Destructor.
     * @ingroup MultidimArrayConstructors
     */
     ~MultidimArray()
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
        zdim=ydim=xdim=yxdim=zyxdim=nzyxdim=0;
        zinit=yinit=xinit=0;
        ndim=1;
        data=NULL;
        destroyData=true;
    }

    /** Core allocate.
     * @ingroup MultidimArrayCore
     */
    void coreAllocate(int _ndim, int _zdim, int _ydim, int _xdim)
    {
        if (_ndim <= 0 || _zdim <= 0 || _ydim<=0 || _xdim<=0)
        {
            clear();
            return;
        }

        ndim=_ndim;
        zdim=_zdim;
        ydim=_ydim;
        xdim=_xdim;
        yxdim=ydim*xdim;
        zyxdim=zdim*yxdim;
        nzyxdim=ndim*zyxdim;

        data = new T [nzyxdim];
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

    /// @defgroup MultidimSize Size
    /// @ingroup MultidimensionalArrays
    /** Copy the shape parameters
     * @ingroup MultidimSize
     *
     */
    void copyShape(const MultidimArray<T> &m)
    {
        ndim=m.ndim;
        zdim=m.zdim;
        ydim=m.ydim;
        xdim=m.xdim;
        yxdim=m.yxdim;
        zyxdim=m.zyxdim;
        nzyxdim=m.nzyxdim;
        zinit=m.zinit;
        yinit=m.yinit;
        xinit=m.xinit;
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
        if (NSIZE(*this) != NSIZE(v) || XSIZE(*this) != XSIZE(v) || 
            YSIZE(*this) != YSIZE(v) || ZSIZE(*this) != ZSIZE(v))
            resize(NSIZE(v), ZSIZE(v), YSIZE(v), XSIZE(v));

        STARTINGX(*this) = STARTINGX(v);
        STARTINGY(*this) = STARTINGY(v);
        STARTINGZ(*this) = STARTINGZ(v);
    }

    /** Resize a single 2D image
     * @ingroup MultidimSize
     *
     * This function assumes n and z are 1
     * @code
     * V1.resize(3, 2);
     * @endcode
     */
    void resize(int Ydim, int Xdim)
    {
        resize(1, 1, Ydim, Xdim);
    }

    /** Resize a single 3D image
     * @ingroup MultidimSize
     *
     * This function assumes n is 1
     * @code
     * V1.resize(3, 3, 2);
     * @endcode
     */
    void resize(int Zdim, int Ydim, int Xdim)
    {
        resize(1, Zdim, Ydim, Xdim);
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
#define FORCE_RESIZE true
    void resize(int Ndim, int Zdim, int Ydim, int Xdim)
    {
        if (Xdim == XSIZE(*this) && Ydim == YSIZE(*this) &&
            Zdim == ZSIZE(*this) && Ndim == NSIZE(*this) )
            return;

        if (Xdim <= 0 || Ydim <= 0 || Zdim <= 0 || Ndim <= 0)
        {
            clear();
            return;
        }

        // Ask for memory
        size_t YXdim=Ydim*Xdim;
        size_t ZYXdim=Zdim*YXdim;
        size_t NZYXdim=Ndim*ZYXdim;
	
	T * new_data;
	
	try
	{	
		new_data = new T [NZYXdim];
	}
	catch (std::bad_alloc &)
	{
		REPORT_ERROR(1001, "Allocate: No space left");
	}

        // Copy needed elements, fill with 0 if necessary
        for (int l = 0; l < Ndim; l++)
            for (int k = 0; k < Zdim; k++)
                for (int i = 0; i < Ydim; i++)
                    for (int j = 0; j < Xdim; j++)
                    {
                        T val;
                        if (k >= ZSIZE(*this))
                            val = 0;
                        else if (i >= YSIZE(*this))
                            val = 0;
                        else if (j >= XSIZE(*this))
                            val = 0;
                        else
                            val = DIRECT_VOL_ELEM(*this, k, i, j);
                        new_data[l*ZYXdim + k*YXdim+i*Xdim+j] = val;
                    }

        // deallocate old vector
        coreDeallocate();

        // assign *this vector to the newly created
        data = new_data;
        ndim = Ndim;
        xdim = Xdim;
        ydim = Ydim;
        zdim = Zdim;
        yxdim = Ydim * Xdim;
        zyxdim = Zdim * yxdim;
        nzyxdim = Ndim * zyxdim;

    }

    /** Get size.
     * @ingroup MultidimSize
     *
     * Returns the size of the object in a 4D vector. If the object is a matrix
     * or a vector, then the higher order dimensions will be set to 1, ie,
     * (Xdim, 1, 1) or (Xdim, Ydim, 1).
     *
     * This function is not ported to Python.
     */
    void getSize(int* size) const
    {
        size[0] = xdim;
        size[1] = ydim;
        size[2] = zdim;
        size[3] = ndim;
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

    /** Same shape.
     * @ingroup MultidimSize
     *
     * Returns true if this object has got the same shape (origin and size)
     * than the argument
     */
    template <typename T1>
    bool sameShape(const MultidimArray<T1>& op) const
    {
        return (NSIZE(*this) == NSIZE(op) &&
                XSIZE(*this) == XSIZE(op) &&
                YSIZE(*this) == YSIZE(op) &&
                ZSIZE(*this) == ZSIZE(op) &&
                STARTINGX(*this) == STARTINGX(op) &&
                STARTINGY(*this) == STARTINGY(op) &&
                STARTINGZ(*this) == STARTINGZ(op));
    }

    /** Set logical origin in Xmipp fashion.
     * @ingroup MultidimSize
     *
     * This function adjust the starting points in the array such that the
     * center of the array is defined in the Xmipp fashion.
     *
     * @code
     * V.setXmippOrigin();
     * @endcode
     */
    void setXmippOrigin()
    {
        zinit = FIRST_XMIPP_INDEX(zdim);
        yinit = FIRST_XMIPP_INDEX(ydim);
        xinit = FIRST_XMIPP_INDEX(xdim);
    }

    /** Move origin to.
     * @ingroup MultidimSize
     *
     * This function adjust logical indexes such that the Xmipp origin of the
     * array moves to the specified position. For instance, an array whose x
     * indexes go from -1 to 1, if we move the origin to 4, then the x indexes
     * go from 3 to 5. This is very useful for convolution operations where you
     * only need to move the logical starting of the array.
     *
     */
    void moveOriginTo(int k, int i, int j)
    {
        zinit = k + FIRST_XMIPP_INDEX(zdim);
        yinit = i + FIRST_XMIPP_INDEX(ydim);
        xinit = j + FIRST_XMIPP_INDEX(xdim);
    }

    /** Returns the first valid logical Z index.
     * @ingroup MultidimSize
     */
    int startingZ() const
    {
        return zinit;
    }

    /** Returns the last valid logical Z index.
     * @ingroup MultidimSize
     */
    int finishingZ() const
    {
        return zinit + zdim - 1;
    }

    /** Returns the first valid logical Y index.
     * @ingroup MultidimSize
     */
    int startingY() const
    {
        return yinit;
    }

    /** Returns the last valid logical Y index.
     * @ingroup MultidimSize
     */
    int finishingY() const
    {
        return yinit + ydim - 1;
    }

    /** Returns the first valid logical X index.
     * @ingroup MultidimSize
     */
    int startingX() const
    {
        return xinit;
    }

    /** Returns the last valid logical X index.
     * @ingroup MultidimSize
     */
    int finishingX() const
    {
        return xinit + xdim - 1;
    }

    /** Returns Z dimension.
     * @ingroup VolumesSizeShape
     */
    int sliceNumber() const
    {
        return zdim;
    }

    /** Returns Y dimension.
     * @ingroup VolumesSizeShape
     */
    int rowNumber() const
    {
        return ydim;
    }

    /** Returns X dimension.
     * @ingroup VolumesSizeShape
     */
    int colNumber() const
    {
        return xdim;
    }

    /// @defgroup Statistics Statistics functions
    /// @ingroup MultidimensionalArrays
    /** Print statistics in current line.
     * @ingroup Statistics
     *
     * No end of line character is written after this print out.
     *
     * @code
     * a.computeStats();
     * std::cout << "Statistics of variable a ";
     * a.printStats();
     * std::cout << std::endl;
     * @endcode
     */
    void printStats(std::ostream& out = std::cout) const
    {
        T minval, maxval;
        double avgval, devval;

        computeStats(avgval, devval, minval, maxval);

        out.setf(std::ios::showpoint);
        int old_prec = out.precision(7);

        out << " min= ";
        out.width(9);
        out << minval;
        out << " max= ";
        out.width(9);
        out << maxval;
        out << " avg= ";
        out.width(9);
        out << avgval;
        out << " dev= ";
        out.width(9);
        out << devval;

        out.precision(old_prec);
    }

    /** Maximum of the values in the array.
     * @ingroup Statistics
     *
     * The returned value is of the same type as the type of the array.
     */
    T computeMax() const
    {
        if (NZYXSIZE(*this) <= 0)
            return static_cast< T >(0);

        T maxval = data[0];

        T* ptr=NULL;
	unsigned long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
            if (*ptr > maxval)
                maxval = *ptr;

        return maxval;
    }

    /** Minimum of the values in the array.
     * @ingroup Statistics
     *
     * The returned value is of the same type as the type of the array.
     */
    T computeMin() const
    {
        if (NZYXSIZE(*this) <= 0)
            return static_cast< T >(0);

        T minval = data[0];

        T* ptr=NULL;
	unsigned long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
            if (*ptr < minval)
                minval = *ptr;

        return minval;
    }

    /** Minimum and maximum of the values in the array.
     * @ingroup Statistics
     *
     * As doubles.
     */
    void computeDoubleMinMax(double& minval, double& maxval) const
    {
        if (NZYXSIZE(*this) <= 0)
            return;

        minval = maxval = static_cast< double >(data[0]);

        T* ptr=NULL;
	unsigned long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
        {
            if (*ptr < minval)
                minval = static_cast< double >(*ptr);

            if (*ptr > maxval)
                maxval = static_cast< double >(*ptr);
        }
    }

    /** Average of the values in the array.
     * @ingroup Statistics
     *
     * The returned value is always double, independently of the type of the
     * array.
     */
    double computeAvg() const
    {
        if (NZYXSIZE(*this) <= 0)
            return 0;

        double sum = 0;

        T* ptr=NULL;
	unsigned long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
            sum += static_cast< double >(*ptr);

        return sum / NZYXSIZE(*this);
    }

    /** Standard deviation of the values in the array.
     * @ingroup Statistics
     *
     * Be careful that the standard deviation and NOT the variance is returned.
     * The returned value is always double, independently of the type of the
     * array.
     */
    double computeStddev() const
    {
        if (NZYXSIZE(*this) <= 1)
            return 0;

        double avg = 0, stddev = 0;

        T* ptr=NULL;
	unsigned long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
        {
            avg += static_cast< double >(*ptr);
            stddev += static_cast< double >(*ptr) *
                      static_cast< double >(*ptr);
        }

        avg /= NZYXSIZE(*this);
        stddev = stddev / NZYXSIZE(*this) - avg * avg;
        stddev *= NZYXSIZE(*this) / (NZYXSIZE(*this) - 1);

        // Foreseeing numerical instabilities
        stddev = sqrt(static_cast<double>((ABS(stddev))));

        return stddev;
    }

    /** Compute statistics.
     * @ingroup Statistics
     *
     * The average, standard deviation, minimum and maximum value are
     * returned.
     */
    void computeStats(double& avg, double& stddev, T& minval, T& maxval) const
    {
        if (NZYXSIZE(*this) <= 0)
            return;

        avg = 0;
        stddev = 0;

        minval = maxval = data[0];

        T* ptr=NULL;
	unsigned long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
        {
            avg += static_cast< double >(*ptr);
            stddev += static_cast< double >(*ptr) *
                      static_cast< double >(*ptr);

            if (*ptr > maxval)
                maxval = *ptr;

            if (*ptr < minval)
                minval = *ptr;
        }

        avg /= NZYXSIZE(*this);

        if (NZYXSIZE(*this) > 1)
        {
            stddev = stddev / NZYXSIZE(*this) - avg * avg;
            stddev *= NZYXSIZE(*this) / (NZYXSIZE(*this) - 1);

            // Foreseeing numerical instabilities
            stddev = sqrt(static_cast< double >(ABS(stddev)));
        }
        else
            stddev = 0;
    }

    /** Median
     * @ingroup Statistics
     *
     * Calculate the median element.
     *
     * @code
     * med = v1.computeMedian();
     * @endcode
     */
    double computeMedian() const
    {
        if (XSIZE(*this) == 0)
            return 0;

        if (XSIZE(*this) == 1)
            return DIRECT_MULTIDIM_ELEM(*this,0);

        // Initialise data
        MultidimArray< double > temp;
        temp.resize(*this);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(*this)
            DIRECT_MULTIDIM_ELEM(temp,n)=DIRECT_MULTIDIM_ELEM(*this,n);

        // Sort indexes
        double* temp_array = MULTIDIM_ARRAY(temp)-1;
        qcksrt(NZYXSIZE(*this), temp_array);

        // Get median
        if (NZYXSIZE(*this)%2==0)
            return 0.5*(DIRECT_MULTIDIM_ELEM(temp,NZYXSIZE(*this)/2-1)+
                        DIRECT_MULTIDIM_ELEM(temp,NZYXSIZE(*this)/2  ));
        else
            return DIRECT_MULTIDIM_ELEM(temp,NZYXSIZE(*this)/2);
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
        computeDoubleMinMax(min0, max0);

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
        computeStats(avg0, stddev0, minval, maxval);

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

    /// @defgroup Arithmethic Arithmethic operations
    /// @ingroup MultidimensionalArrays

    /** @defgroup ArrayByArray Array "by" array operations.
     * @ingroup Arithmethic
     *
     * These are operations that are performed between 2 arrays of the
     * SAME type (two integer vectors, two double matrices, ...). If they
     * are not of the same type you can convert one of the arrays to the
     * desired type using the function typeCast. The result must have been
     * defined to be of the same type as the operands.
     *
     * In this kind of operations each element of array 1 is operated with its
     * homologous in array 2, it is very important that both have got the
     * same size and starting origins. The result has also got the same
     * shape as the two operated arrays and its former content is lost.
     */

    /** Core array by array operation.
     * @ingroup ArrayByArray
     *
     * It assumes that the result is already resized.
     */
    friend void coreArrayByArray<>(const MultidimArray<T>& op1,
        const MultidimArray<T>& op2, MultidimArray<T>& result,
        char operation);
    
    /** Array by array
     * @ingroup ArrayByArray
     *
     * This function must take two vectors of the same size, and operate element
     * by element according to the operation required. This is the function
     * which really implements the operations. Simple calls to it perform much
     * faster than calls to the corresponding operators. Although it is supposed
     * to be a hidden function not useable by normal programmers.
     *
     * In the case of Matrix2D, the multiplication is performed as an algebraic
     * multiplication.
     */
    friend void arrayByArray(const MultidimArray<T>& op1,
        const MultidimArray<T>& op2, MultidimArray<T>& result,
        char operation)
    {
        if (!op1.sameShape(op2))
            REPORT_ERROR(1007,
                         (std::string) "Array_by_array: different shapes (" +
                         operation + ")");

        result.resize(op1);
        coreArrayByArray(op1, op2, result, operation);
    }

    /// FIXME: THIS SHOULD BE MOVED TO MATRIX2D!!!!
    /** Algebraic multiplication of two matrices
     * @ingroup ArrayByArray
     */
    friend void multiplyMatrix(const MultidimArray<T>& op1,
        const MultidimArray<T>& op2, MultidimArray<T>& result)
    {
        if (XSIZE(op1) != YSIZE(op2))
            REPORT_ERROR(1102, "Not compatible sizes in matrix multiplication");
        if (ZSIZE(op1) != 1 || ZSIZE(op2))
            REPORT_ERROR(1102, "Operation intended only for matrices");

        result.initZeros(1, YSIZE(op1), XSIZE(op2));
        for (int i = 0; i < YSIZE(op1); i++)
            for (int j = 0; j < XSIZE(op2); j++)
                for (int k = 0; k < XSIZE(op1); k++)
                    DIRECT_MAT_ELEM(result, i, j) += DIRECT_MAT_ELEM(op1, i, k) *
                                                     DIRECT_MAT_ELEM(op2, k, j);

        STARTINGZ(result) = STARTINGZ(op1);
        STARTINGY(result) = STARTINGY(op1);
        STARTINGX(result) = STARTINGX(op1);
    }

    /** v3 = v1 + v2.
     * @ingroup ArrayByArray
     */
    MultidimArray<T> operator+(const MultidimArray<T>& op1) const
    {
        MultidimArray<T> tmp;
        arrayByArray(*this, op1, tmp, '+');
        return tmp;
    }

    /** v3 = v1 - v2.
     * @ingroup ArrayByArray
     */
    MultidimArray<T> operator-(const MultidimArray<T>& op1) const
    {
        MultidimArray<T> tmp;
        arrayByArray(*this, op1, tmp, '-');
        return tmp;
    }

    /** v3 = v1 * v2.
     * @ingroup ArrayByArray
     */
    MultidimArray<T> operator*(const MultidimArray<T>& op1) const
    {
        MultidimArray<T> tmp;
        arrayByArray(*this, op1, tmp, '*');
        return tmp;
    }

    /** v3 = v1 / v2.
     * @ingroup ArrayByArray
     */
    MultidimArray<T> operator/(const MultidimArray<T>& op1) const
    {
        MultidimArray<T> tmp;
        arrayByArray(*this, op1, tmp, '/');
        return tmp;
    }

    /** v3 += v2.
     * @ingroup ArrayByArray
     */
    void operator+=(const MultidimArray<T>& op1)
    {
        arrayByArray(*this, op1, *this, '+');
    }

    /** v3 -= v2.
     * @ingroup ArrayByArray
     */
    void operator-=(const MultidimArray<T>& op1)
    {
        arrayByArray(*this, op1, *this, '-');
    }

    /** v3 *= v2.
     * @ingroup ArrayByArray
     */
    void operator*=(const MultidimArray<T>& op1)
    {
        arrayByArray(*this, op1, *this, '*');
    }

    /** v3 /= v2.
     * @ingroup ArrayByArray
     */
    void operator/=(const MultidimArray<T>& op1)
    {
        arrayByArray(*this, op1, *this, '/');
    }

    /** @defgroup ArrayByScalar Array "by" scalar operations
     * @ingroup Arithmethic
     *
     * These operations are between an array and a scalar (of the same type as
     * the array). The result must have been defined to be of the same type as
     * the operands.
     *
     * In this kind of operations each element of array 1 is operated with the
     * given constant. The result has also got the same shape as the input
     * array and its former content is lost
     */

    /** Array by scalar.
     * @ingroup ArrayByScalar
     *
     * This function must take one vector and a constant, and operate element
     * by element according to the operation required. This is the function
     * which really implements the operations. Simple calls to it perform much
     * faster than calls to the corresponding operators. Although it is
     * supposed to be a hidden function not useable by normal programmers.
     *
     * This function is not ported to Python.
     */
    friend void arrayByScalar(const MultidimArray<T>& op1,
                                T op2,
                                MultidimArray<T>& result,
                                char operation)
    {
        result.resize(op1);
        coreArrayByScalar(op1, op2, result, operation);
    }

    /** Core array by scalar operation.
     * @ingroup ArrayByScalar
     *
     * It assumes that the result is already resized.
     *
     * This function is not ported to Python.
     */
    friend void coreArrayByScalar<>(const MultidimArray<T>& op1,
                                    const T& op2,
                                    MultidimArray<T>& result,
                                    char operation);

    /** v3 = v1 + k.
     * @ingroup ArrayByScalar
     */
    MultidimArray<T> operator+(T op1) const
    {
        MultidimArray<T> tmp;
        arrayByScalar(*this, op1, tmp, '+');
        return tmp;
    }

    /** v3 = v1 - k.
     * @ingroup ArrayByScalar
     */
    MultidimArray<T> operator-(T op1) const
    {
        MultidimArray<T> tmp;
        arrayByScalar(*this, op1, tmp, '-');
        return tmp;
    }

    /** v3 = v1 * k.
     * @ingroup ArrayByScalar
     */
    MultidimArray<T> operator*(T op1) const
    {
        MultidimArray<T> tmp;
        arrayByScalar(*this, op1, tmp, '*');
        return tmp;
    }

    /** v3 = v1 / k.
     * @ingroup ArrayByScalar
     */
    MultidimArray<T> operator/(T op1) const
    {
        MultidimArray<T> tmp;
        arrayByScalar(*this, op1, tmp, '/');
        return tmp;
    }

    /** v3 += k.
     * @ingroup ArrayByScalar
     *
     * This function is not ported to Python.
     */
    void operator+=(const T& op1)
    {
        arrayByScalar(*this, op1, *this, '+');
    }

    /** v3 -= k.
     * @ingroup ArrayByScalar
     *
     * This function is not ported to Python.
     */
    void operator-=(const T& op1)
    {
        arrayByScalar(*this, op1, *this, '-');
    }

    /** v3 *= k.
     * @ingroup ArrayByScalar
     *
     * This function is not ported to Python.
     */
    void operator*=(const T& op1)
    {
        arrayByScalar(*this, op1, *this, '*');
    }

    /** v3 /= k.
     * @ingroup ArrayByScalar
     *
     * This function is not ported to Python.
     */
    void operator/=(const T& op1)
    {
        arrayByScalar(*this, op1, *this, '/');
    }

    /** @defgroup ScalarByArray Scalar "by" array operations
     * @ingroup Arithmethic
     *
     * These operations are between a scalar (of the same type as the array)
     * and an array. The result must have been defined to be of the same type
     * as the operand. The former content of the result array is lost after
     * the operation.
     *
     * In this kind of operations the constant is operated with each element
     * of array 2. The result has also got the same shape as the input array
     * and its former content is lost
     */

    /** Scalar by array.
     * @ingroup ScalarByArray
     *
     * This function must take one scalar and a vector, and operate element by
     * element according to the operation required. This is the function which
     * really implements the operations. Simple calls to it perform much faster
     * than calls to the corresponding operators. Although it is supposed to
     * be a hidden function not useable by normal programmers.
     *
     * This function is not ported to Python.
     */
    friend void scalarByArray(T op1,
                                const MultidimArray<T>& op2,
                                MultidimArray<T>& result,
                                char operation)
    {
        result.resize(op2);
        coreScalarByArray(op1, op2, result, operation);
    }

    /** Core array by scalar operation.
     * @ingroup ScalarByArray
     *
     * It assumes that the result is already resized.
     *
     * This function is not ported to Python.
     */
    friend void coreScalarByArray<>(const T& op1,
                                    const MultidimArray<T>& op2,
                                    MultidimArray<T>& result,
                                    char operation);

    /** v3 = k + v2.
     * @ingroup ScalarByArray
     */
    friend MultidimArray<T> operator+(T op1, const MultidimArray<T>& op2)
    {
        MultidimArray<T> tmp;
        scalarByArray(op1, op2, tmp, '+');
        return tmp;
    }

    /** v3 = k - v2.
     * @ingroup ScalarByArray
     */
    friend MultidimArray<T> operator-(T op1, const MultidimArray<T>& op2)
    {
        MultidimArray<T> tmp;
        scalarByArray(op1, op2, tmp, '-');
        return tmp;
    }

    /** v3 = k * v2.
     * @ingroup ScalarByArray
     */
    friend MultidimArray<T> operator*(T op1, const MultidimArray<T>& op2)
    {
        MultidimArray<T> tmp;
        scalarByArray(op1, op2, tmp, '*');
        return tmp;
    }

    /** v3 = k / v2
     * @ingroup ScalarByArray
     */
    friend MultidimArray<T> operator/(T op1, const MultidimArray<T>& op2)
    {
        MultidimArray<T> tmp;
        scalarByArray(op1, op2, tmp, '/');
        return tmp;
    }

    /// @defgroup Initialization Initialization
    /// @ingroup MultidimensionalArrays

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
        T* ptr=NULL;
	unsigned long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
            *ptr = val;
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
        resize(op);
        initConstant(static_cast< T >(0));
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
    void initZeros(int Zdim, int Ydim, int Xdim)
    {
        resize(Zdim,Ydim,Xdim);
        initConstant(static_cast< T >(0));
    }

    /** Initialize to zeros with a given size.
     * @ingroup Initialization
     */
    void initZeros(int Ndim, int Zdim, int Ydim, int Xdim)
    {
        resize(Ndim, Zdim,Ydim,Xdim);
        initConstant(static_cast< T >(0));
    }

    /** Initialize with random values.
     * @ingroup Initialization
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

    /** Add noise to actual values.
     * @ingroup Initialization
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

    /** @defgroup MultidimUtilities Utilities
     *  @ingroup MultidimensionalArrays
     *
     * Here you have several easy functions to manage the values of
     * the array.
     */

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

    /** ROUND n-dimensional.
     * @ingroup MultidimUtilities
     *
     * Applies a ROUND (look for the nearest integer) to each array element.
     */
    void selfROUNDnD()
    {
        T* ptr=NULL;
	unsigned long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
            *ptr = ROUND(*ptr);
    }

    /** CEILING n-dimensional.
     * @ingroup MultidimUtilities
     *
     * Applies a CEILING (look for the nearest larger integer) to each
     * array element.
     */
    void selfCEILnD()
    {
        T* ptr=NULL;
	unsigned long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
            *ptr = CEIL(*ptr);
    }

    /** FLOOR n-dimensional.
     * @ingroup MultidimUtilities
     *
     * Applies a FLOOR (look for the nearest larger integer) to each
     * array element.
     */
    void selfFLOORnD()
    {
        T* ptr=NULL;
	unsigned long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
            *ptr = FLOOR(*ptr);
    }

    /** ABS n-dimensional.
     * @ingroup MultidimUtilities
     *
     * Applies an ABS (absolute value) to each array element.
     */
    void selfABSnD()
    {
        T* ptr=NULL;
	unsigned long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
            *ptr = ABS(*ptr);
    }

    /** MAX n-dimensional.
     * @ingroup MultidimUtilities
     *
     * Each component of the result is the maximum of the correspoing
     * components of the two input arrays. They must have the same shape, if
     * not an exception is thrown
     */
    friend void MAXnD(const MultidimArray<T>& v1, const MultidimArray<T>& v2,
        MultidimArray<T>& result)
    {
        if (!v1.sameShape(v2))
            REPORT_ERROR(1007, "MAX: arrays of different shape");

        result.resize(v1);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(result)
            DIRECT_MULTIDIM_ELEM(result,n) = XMIPP_MAX(
                DIRECT_MULTIDIM_ELEM(v1,n),
                DIRECT_MULTIDIM_ELEM(v2,n));
    }

    /** MIN n-dimensional.
     * @ingroup MultidimUtilities
     *
     * Each component of the result is the minimum of the correspoing
     * components of the two input arrays. They must have the same shape, if
     * not an exception is thrown
     */
    friend void MINnD(const MultidimArray<T>& v1, const MultidimArray<T>& v2,
        MultidimArray<T>& result)
    {
        if (!v1.sameShape(v2))
            REPORT_ERROR(1007, "MIN: arrays of different shape");

        result.resize(v1);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(result)
            DIRECT_MULTIDIM_ELEM(result,n) = XMIPP_MIN(
                DIRECT_MULTIDIM_ELEM(v1,n),
                DIRECT_MULTIDIM_ELEM(v2,n));
    }

    /** Sqrt.
     * @ingroup MultidimUtilities
     *
     * Each component of the result is the square root of the original
     * component.
     */
    void selfSQRTnD()
    {
        T* ptr=NULL;
	unsigned long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
            *ptr = static_cast< T >(sqrt(static_cast< double >(*ptr)));
    }

    /** Sum of matrix values.
     * @ingroup MultidimUtilities
     *
     * This function returns the sum of all internal values.
     *
     * @code
     * double sum = m.sum();
     * @endcode
     */
    double sum() const
    {
        double sum = 0;
        T* ptr=NULL;
	unsigned long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
            sum += *ptr;
        return sum;
    }

    /** Sum of squared vector values.
     * @ingroup MultidimUtilities
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
        T* ptr=NULL;
	unsigned long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
            sum += *ptr * *ptr;
        return sum;
    }

    /** Log10.
     * @ingroup MultidimUtilities
     *
     * Each component of the result is the log10 of the original components.
     */
    void selfLog10()
    {
        T* ptr=NULL;
	unsigned long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
            *ptr = static_cast< T >(log10(static_cast< double >(*ptr)));
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
    void selfReverseX()
    {
        for (int l = 0; l < NSIZE(*this); l++)
            for (int k = 0; k < ZSIZE(*this); k++)
                for (int i = 0; i < YSIZE(*this); i++)
                    for (int j = 0; j <= (int)(XSIZE(*this) - 1) / 2; j++)
                    {
                        T aux;
                        SWAP(DIRECT_NZYX_ELEM(*this, l, k, i, j),
                             DIRECT_NZYX_ELEM(*this, l, k, i, XSIZE(*this) - 1 - j),
                             aux);
                    }

        STARTINGX(*this) = -FINISHINGX(*this);
    }

    /** Reverse matrix values over Y axis, keep in this object.
     * @ingroup Utilites
     *
     * Maybe better with an example:
     *
     * @code
     * slice 0
     * [01 02 03          [03 02 01
     *  04 05 06           06 05 04
     *  07 08 09]          09 08 07]
     *
     * ----->
     *
     * slice 1
     * [11 12 13          [13 12 11
     *  14 15 16           16 15 14
     *  17 18 19]          19 18 17]
     * @endcode
     *
     */
    void selfReverseY()
    {
        for (int l = 0; l < NSIZE(*this); l++)
            for (int k = 0; k < ZSIZE(*this); k++)
                for (int i = 0; i <= (int)(YSIZE(*this) - 1) / 2; i++)
                    for (int j = 0; j < XSIZE(*this); j++)
                    {
                        T aux;
                        SWAP(DIRECT_NZYX_ELEM(*this, l, k, i, j),
                             DIRECT_NZYX_ELEM(*this, l, k, YSIZE(*this) - 1 - i, j),
                             aux);
                    }

        STARTINGY(*this) = -FINISHINGY(*this);
    }

    /** Reverse matrix values over Z axis, keep in this object.
     * @ingroup Utilites
     *
     * Maybe better with an example:
     *
     * @code
     * slice 0
     * [01 02 03          [11 12 13
     *  04 05 06           14 15 16
     *  07 08 09]          17 18 19]
     *
     *  ----->
     *
     * slice 1
     * [11 12 13          [01 02 03
     *  14 15 16           04 05 06
     *  17 18 19]          07 08 09]
     * @endcode
     *
     */
    void selfReverseZ()
    {

        for (int l = 0; l < NSIZE(*this); l++)
            for (int k = 0; k <= (int)(ZSIZE(*this) - 1) / 2; k++)
                for (int i = 0; i < YSIZE(*this); i++)
                    for (int j = 0; j < XSIZE(*this); j++)
                    {
                        T aux;
                        SWAP(DIRECT_NZYX_ELEM(*this, l, k, i, j),
                             DIRECT_NZYX_ELEM(*this, l, ZSIZE(*this) - 1 - k, i, j),
                             aux);
                    }

        STARTINGZ(*this) = -FINISHINGZ(*this);
    }


    /// FIXME RETHINK FOR 4D!??

    /** Produce spline coefficients.
     * @ingroup Utilites
     */
#ifndef DBL_EPSILON
#define DBL_EPSILON 1e-50
#endif
    void produceSplineCoefficients(MultidimArray< double >& coeffs, int SplineDegree = 3)
    const
    {
        coeffs.initZeros(NSIZE(*this), ZSIZE(*this), YSIZE(*this), XSIZE(*this));
        STARTINGX(coeffs) = STARTINGX(*this);
        STARTINGY(coeffs) = STARTINGY(*this);
        STARTINGZ(coeffs) = STARTINGZ(*this);

        int Status;
        MultidimArray< double > aux;
        typeCast(*this, aux);

        ChangeBasisVolume(MULTIDIM_ARRAY(aux), MULTIDIM_ARRAY(coeffs),
                          XSIZE(*this), YSIZE(*this), ZSIZE(*this),
                          CardinalSpline, BasicSpline, SplineDegree,
                          MirrorOffBounds, DBL_EPSILON, &Status);
        if (Status)
            REPORT_ERROR(1, "Matrix3D::produceSplineCoefficients: Error");
    }

    /// FIXME RETHINK FOR 4D!??

    /** Produce image from B-spline coefficients.
     * @ingroup Utilites
     */
    void produceImageFromSplineCoefficients(
        MultidimArray< double >& img, int SplineDegree = 3) const
    {
        img.initZeros(ZSIZE(*this), YSIZE(*this), XSIZE(*this));
        STARTINGX(img) = STARTINGX(*this);
        STARTINGY(img) = STARTINGY(*this);
        STARTINGZ(img) = STARTINGZ(*this);

        int Status;
        MultidimArray< double > aux;
        typeCast(*this, aux);
        
        ChangeBasisVolume(MULTIDIM_ARRAY(aux), MULTIDIM_ARRAY(img),
                          XSIZE(*this), YSIZE(*this), ZSIZE(*this),
                          BasicSpline, CardinalSpline, SplineDegree,
                          MirrorOnBounds, DBL_EPSILON, &Status);
        if (Status)
            REPORT_ERROR(1, "Matrix3D::produce_spline_img: Error");
    }
#undef DBL_EPSILON

    /// @defgroup Operators Operators
    /// @ingroup MultidimensionalArrays

    /** Assignment.
     * @ingroup Operators
     *
     * You can build as complex assignment expressions as you like. Multiple
     * assignment is allowed.
     *
     * @code
     * v1 = v2 + v3;
     * v1 = v2 = v3;
     * @endcode
     *
     * This function is ported to Python as assign.
     */
    MultidimArray<T>& operator=(const MultidimArray<T>& op1)
    {
        if (&op1 != this)
        {
            resize(op1);
            T *ptr;
	    unsigned long int n;
            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(op1,n,ptr)
                DIRECT_MULTIDIM_ELEM(*this,n) = *ptr;
        }

        return *this;
    }

    /** Unary minus.
     * @ingroup Operators
     *
     * It is used to build arithmetic expressions. You can make a minus
     * of anything as long as it is correct semantically.
     *
     * @code
     * v1 = -v2;
     * v1 = -v2.transpose();
     * @endcode
     */
    MultidimArray<T> operator-() const
    {
        MultidimArray<T> tmp(*this);
        T* ptr;
	unsigned long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(tmp,n,ptr)
            *ptr = -(*ptr);
        return tmp;
    }

    /** Input from input stream.
     * @ingroup Operators
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

    /** Equality.
     * @ingroup Operators
     *
     * Returns true if this object has got the same shape (origin and size)
     * than the argument and the same values (within accuracy).
     */
    bool equal(const MultidimArray<T>& op,
    	double accuracy = XMIPP_EQUAL_ACCURACY) const
    {
        if (!sameShape(op))
            return false;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(*this)
            if (ABS(DIRECT_MULTIDIM_ELEM(*this,n) -
                    DIRECT_MULTIDIM_ELEM(op,n)) > accuracy)
                return false;
        return true;
    }


    //FIXME: THIS HAS TO BE RETHOUGHT!
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

    /** Read from a binary file.
     * @ingroup Operators
     *
     * The array must be previously resized to the correct size.
     */
    void readBinary(const FileName& fn)
    {
        std::ifstream in;
        in.open(fn.c_str(), std::ios::in | std::ios::binary);
        if (!in)
            REPORT_ERROR(1,
                         static_cast< std::string>("MultidimArray::read: File " + fn
                                                   + " not found"));

        T* ptr;
	unsigned long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
            in.read(reinterpret_cast< char* >(ptr), sizeof(T));

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

    /** Write to a binary file.
     * @ingroup Operators
     */
    void writeBinary(const FileName& fn) const
    {
        std::ofstream out;

        out.open(fn.c_str(), std::ios::out | std::ios::binary);
        if (!out)
            REPORT_ERROR(1,
                         static_cast< std::string >("MultidimArray::write: File " + fn
                                                    + " not found"));

        T* ptr;
	unsigned long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
            out.write(reinterpret_cast< char* >(ptr), sizeof(T));

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

/// @defgroup MultidimFunctions Functions for all multidimensional arrays
/// @ingroup MultidimensionalArrays

/** Conversion from one type to another.
 * @ingroup MultidimFunctions
 *
 * If we have an integer array and we need a double one, we can use this
 * function. The conversion is done through a type casting of each element
 */
template<typename T1, typename T2>
void typeCast(const MultidimArray<T1>& v1, MultidimArray<T2>& v2)
{
    if (NZYXSIZE(v1) == 0)
    {
        v2.clear();
        return;
    }

    v2.resize(v1);
    T1* ptr1=NULL;
    unsigned long int n;
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(v1,n,ptr1)
        DIRECT_MULTIDIM_ELEM(v2,n) = static_cast< T2 >(*ptr1);
}

template <typename T>
void coreArrayByArray(const MultidimArray<T>& op1,
    const MultidimArray<T>& op2, MultidimArray<T>& result,
    char operation)
{
    T* ptrResult=NULL;
    T* ptrOp1=NULL;
    T* ptrOp2=NULL;
    unsigned long int n;
    for (n=0, ptrResult=result.data, ptrOp1=op1.data,ptrOp2=op2.data;
	n<op1.zyxdim; ++n, ++ptrResult, ++ptrOp1, ++ptrOp2)
        switch (operation)
        {
        case '+':
            *ptrResult = *ptrOp1 + *ptrOp2;
            break;
        case '-':
            *ptrResult = *ptrOp1 - *ptrOp2;
            break;
        case '*':
            *ptrResult = *ptrOp1 * *ptrOp2;
            break;
        case '/':
            *ptrResult = *ptrOp1 / *ptrOp2;
            break;
        }
}

template <typename T>
void coreArrayByScalar(const MultidimArray<T>& op1,
                            const T& op2,
                            MultidimArray<T>& result,
                            char operation)
{
    T* ptrResult=NULL;
    T* ptrOp1=NULL;
    unsigned long int n;
    for (n=0, ptrResult=result.data, ptrOp1=op1.data;
         n<op1.zyxdim; ++n, ++ptrResult, ++ptrOp1)
        switch (operation)
        {
        case '+':
            *ptrResult = *ptrOp1 + op2;
            break;
        case '-':
            *ptrResult = *ptrOp1 - op2;
            break;
        case '*':
            *ptrResult = *ptrOp1 * op2;
            break;
        case '/':
            *ptrResult = *ptrOp1 / op2;
            break;
        }
}

template <typename T>
void coreScalarByArray(const T& op1,
                         const MultidimArray<T>& op2,
                         MultidimArray<T>& result,
                         char operation)
{
    T* ptrResult=NULL;
    T* ptrOp2=NULL;
    unsigned long int n;
    for (n=0, ptrResult=result.data, ptrOp2=op2.data;
         n<op2.zyxdim; ++n, ++ptrResult, ++ptrOp2)
        switch (operation)
        {
        case '+':
            *ptrResult = op1 + *ptrOp2;
            break;
        case '-':
            *ptrResult = op1 - *ptrOp2;
            break;
        case '*':
            *ptrResult = op1 * *ptrOp2;
            break;
        case '/':
            *ptrResult = op1 / *ptrOp2;
            break;
        }
}

template<>
inline void MultidimArray< std::complex< double > >::produceSplineCoefficients(
    MultidimArray< double >& coeffs, int SplineDegree) const
{
    // TODO Implement
    std::cerr << "Spline coefficients of a complex matrix is not implemented\n";
}

template<typename T>
std::ostream& operator<<(std::ostream& ostrm, const MultidimArray<T>& v)
{
    if (v.xdim == 0)
        ostrm << "NULL Array\n";
    else
        ostrm << std::endl;

    double max_val = ABS(DIRECT_VOL_ELEM(v , 0, 0, 0));

    T* ptr;
    unsigned long int n;
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(v,n,ptr)
    	max_val = XMIPP_MAX(max_val, ABS(*ptr));

    int prec = bestPrecision(max_val, 10);

    if (YSIZE(v)==1 && ZSIZE(v)==1)
    {
        for (int j = STARTINGX(v); j <= FINISHINGX(v); j++)
            ostrm << floatToString((double) VOL_ELEM(v, 0, 0, j), 10, prec)
                  << std::endl;
    }
    else
    {
        for (int l = 0; l < NSIZE(v); l++)
        {
            if (NSIZE(v)>1) ostrm << "Image No. " << l << std::endl;
            for (int k = STARTINGZ(v); k <= FINISHINGZ(v); k++)
            {
                if (ZSIZE(v)>1) ostrm << "Slice No. " << k << std::endl;
                for (int i = STARTINGY(v); i <= FINISHINGY(v); i++)
                {
                    for (int j = STARTINGX(v); j <= FINISHINGX(v); j++)
                    {
                        ostrm << floatToString((double) VOL_ELEM(v, k, i, j), 10, prec) << ' ';
                    }
                    ostrm << std::endl;
                }
            }
        }
    }

    return ostrm;
}

template<>
std::ostream& operator<<(std::ostream& ostrm,
    const MultidimArray< std::complex<double> >& v);
#endif
