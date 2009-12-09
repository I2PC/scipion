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

#ifndef MATRIX3D_H
#define MATRIX3D_H

/* FIXME When adding functions to this module, don't forget to add them in the
 * MultidimInstantiation.cc, too
 */

#include <iostream>
#include <string>
#include <complex>

#include <external/bilib/headers/pyramidtools.h>
#include "multidimensional_array.h"
#include "matrix1d.h"
#include "matrix2d.h"
#include "error.h"

#ifndef SWIG
template<typename T> class Matrix3D;

template<typename T>
void applyGeometry(Matrix3D<T>& V2, const Matrix2D< double > &A,
    const Matrix3D<T>& V1, bool inv, bool wrap,
    T outside = 0);

template<typename T>
void applyGeometryBSpline(Matrix3D<T>& V2, const Matrix2D< double > &A,
    const Matrix3D<T>& V1, int Splinedegree, bool inv, bool wrap,
    T outside = 0);
#endif

/// @defgroup Volumes Volumes
/// @ingroup MultidimensionalArrays

/** @defgroup Matrix3dSpeedUp Speed up macros
 *  @ingroup Volumes
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

/** @defgroup VolumesSizeShape Size and shape
 * @ingroup Matrix3dSpeedUp
 */

/** For all elements in the array.
 * @ingroup VolumesSizeShape
 *
 * This macro is used to generate loops for the volume in an easy way. It
 * defines internal indexes 'k','i' and 'j' which ranges the volume using its
 * mathematical definition (ie, logical access).
 *
 * @code
 * FOR_ALL_ELEMENTS_IN_MATRIX3D(V)
 * {
 *     std::cout << V(k, i, j) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_ELEMENTS_IN_MATRIX3D(V) \
    for (size_t k=STARTINGZ(V); k<=FINISHINGZ(V); k++) \
        for (size_t i=STARTINGY(V); i<=FINISHINGY(V); i++) \
            for (size_t j=STARTINGX(V); j<=FINISHINGX(V); j++)

/** For all elements in the array between corners.
 * @ingroup VolumesSizeShape
 *
 * This macro is used to generate loops for a volume in an easy manner. Then
 *  ZZ(r), YY(r) and XX(r) range from
 *
 * (int) ZZ(corner1) to (int)ZZ(corner2),
 * (int) YY(corner1) to (int)YY(corner2),
 * (int) XX(corner1) to (int) XX(corner2) (included limits) respectively.
 *
 * Notice that corner1 and corner2 need only be Matrix1D.
 *
 * @code
 * Matrix1D< double > corner1(3), corner2(3), r(3);
 * XX(corner1) = -1; XX(corner2) = 1;
 * YY(corner1) = -2; YY(corner2) = 2;
 * ZZ(corner1) = -3; ZZ(corner2) = 3;
 *
 * FOR_ALL_ELEMENTS_IN_MATRIX3D_BETWEEN(corner1, corner2)
 * {
 *     std::cout << v(r) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_ELEMENTS_IN_MATRIX3D_BETWEEN(corner1, corner2) \
    for (ZZ(r)=ZZ((corner1)); ZZ(r)<=ZZ((corner2)); ZZ(r)++) \
        for (YY(r)=YY((corner1)); YY(r)<=YY((corner2)); YY(r)++) \
            for (XX(r)=XX((corner1)); XX(r)<=XX((corner2)); XX(r)++)

/** For all elements in common.
 * @ingroup VolumesSizeShape
 *
 * This macro is used to generate loops for all the elements logically in common
 * between two volumes in an easy manner. Then k, i and j (locally defined)
 * range from
 *
 * MAX(STARTINGZ(V1),STARTINGZ(V2)) to MIN(FINISHINGZ(V1),FINISHINGZ(V2)),
 * MAX(STARTINGY(V1),STARTINGY(V2)) to MIN(FINISHINGY(V1),FINISHINGY(V2)),
 * MAX(STARTINGX(V1),STARTINGX(V2)) to MIN(FINISHINGX(V1),FINISHINGX(V2))
 *
 * (included limits) respectively. You need to define SPEED_UP_temps.
 *
 * @code
 * SPEED_UP_temps; 
 * Matrix3D< double > V1(10, 10, 10), V2(20, 20, 20);
 * V1.setXmippOrigin();
 * V2.setXmippOrigin();
 *
 * FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX3D(V1, V2)
 * {
 *    // ...
 * }
 * @endcode
 */
#define FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX3D(V1, V2) \
    ispduptmp0 = XMIPP_MAX(STARTINGZ(V1), STARTINGZ(V2)); \
    ispduptmp1 = XMIPP_MIN(FINISHINGZ(V1),FINISHINGZ(V2)); \
    ispduptmp2 = XMIPP_MAX(STARTINGY(V1), STARTINGY(V2)); \
    ispduptmp3 = XMIPP_MIN(FINISHINGY(V1),FINISHINGY(V2)); \
    ispduptmp4 = XMIPP_MAX(STARTINGX(V1), STARTINGX(V2)); \
    ispduptmp5 = XMIPP_MIN(FINISHINGX(V1),FINISHINGX(V2)); \
    for (int k=ispduptmp0; k<=ispduptmp1; k++) \
        for (int i=ispduptmp2; i<=ispduptmp3; i++) \
            for (int j=ispduptmp4; j<=ispduptmp5; j++)

/** For all direct elements in the array.
 * @ingroup VolumeSizeShape
 *
 * This macro is used to generate loops for the volume in an easy way. It
 * defines internal indexes 'k','i' and 'j' which ranges the volume using its
 * physical definition.
 *
 * @code
 * FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(V)
 * {
 *     std::cout << DIRECT_VOL_ELEM(m, k, i, j) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(V) \
    for (int k=0; k<ZSIZE(V); k++) \
        for (int i=0; i<YSIZE(V); i++) \
            for (int j=0; j<XSIZE(V); j++)

/** @defgroup VolumesMemory Memory access
 * @ingroup Matrix3dSpeedUp
 */

/** A short alias for the previous function.
 * @ingroup VolumesMemory
 *
 * To avoid writing so much
 */
#define dVkij(V, k, i, j) DIRECT_VOL_ELEM(V, k, i, j)

/** Template class for Xmipp volumes
 */
template<typename T>
class Matrix3D: public MultidimArray<T>
{
public:
    /// @defgroup VolumeConstructors Constructors
    /// @ingroup Matrices

    /** Empty constructor
     * @ingroup VolumeConstructors
     */
    Matrix3D(): MultidimArray<T>()
    {
    }

    /** Dimension constructor
     * @ingroup VolumeConstructors
     *
     * The dimension constructor creates a matrix with memory associated (but
     * not assigned to anything, could be full of garbage) origin=0, size=the
     * given one. e careful that first number is the Y dimension (number of
     * rows), and the second the X dimension (number of columns).
     *
     * @code
     * MAtrix3D< double > v1(6, 3);
     * @endcode
     */
    Matrix3D(int Zdim, int Ydim, int Xdim): MultidimArray<T>()
    {
    	resize(Zdim, Ydim,Xdim);
    }

    /** Copy constructor
     * @ingroup VolumeConstructors
     *
     * The created matrix is a perfect copy of the input matrix but with a
     * different memory assignment.
     *
     * @code
     * MAtrix3D< double > m2(m1);
     * @endcode
     */
    Matrix3D(const Matrix3D<T>& v)
    {
        *this = v;
    }

    /** Clear.
     * @ingroup VolumeConstructors
     */
     void clear()
     {
        MultidimArray<T>::clear();
     }

    /** @defgroup VolumesInitialization Initialisation
     * @ingroup Volumes
     */

    /** Alias a multidimarray.
     * @ingroup VolumesInitialization
     *
     * Treat the multidimarray as if it were a volume. The data is not copied
     * into new memory, but a pointer to the multidimarray is copied.
     * You should not make any operation on this volume such that the
     * memory locations are changed
     */
     void alias(const MultidimArray<T> &m)
     {
         copyShape(m);
         this->data=m.data;
         this->destroyData=false;
     }

    /** Zero initialisation with a new dimension.
     * @ingroup VolumesInitialization
     *
     * Be careful to the size order (Zdim, Ydim, Xdim).
     *
     * @code
     * v1.initZeros(6, 3);
     * @endcode
     */
    void initZeros(int Zdim, int Ydim, int Xdim)
    {
        resize(Zdim, Ydim, Xdim);
        initConstant(static_cast<T>(0));
    }

    /** Zero initialisation with current dimension
     * @ingroup VolumesInitialization
     */
    void initZeros()
    {
        MultidimArray<T>::initZeros();
    }

    /** Zero initialisation with current dimension
     * @ingroup VolumesInitialization
     */
    template <typename T1>
    void initZeros(const Matrix3D<T1> &m)
    {
        MultidimArray<T>::initZeros(m);
    }

    /** @defgroup VolumesSizeShape Size and shape
     * @ingroup Volumes
     *
     * The shape of a volume is defined by its origin and its size. The size is
     * clear, and the origin is the logical position of the first real position
     * of the array. For instance, if we have a matrix of dimension (2,5,3) =
     * (Zdim,Ydim, Xdim) and origin (0,-2,-1), this means that the array is
     * representing the logical positions
     *
     * @code
     * Slice 0
     * [(0,-2,-1) (0,-2,0) (0,-2,1)
     *  (0,-1,-1) (0,-1,0) (0,-1,1)
     *  (0, 0,-1) (0, 0,0) (0, 0,1)
     *  (0, 1,-1) (0, 1,0) (0, 1,1)
     *  (0, 2,-1) (0, 2,0) (0, 2,1)]
     *
     * Slice 1
     * [(1,-2,-1) (1,-2,0) (1,-2,1)
     *  (1,-1,-1) (1,-1,0) (1,-1,1)
     *  (1, 0,-1) (1, 0,0) (1, 0,1)
     *  (1, 1,-1) (1, 1,0) (1, 1,1)
     *  (1, 2,-1) (1, 2,0) (1, 2,1)]
     * @endcode
     *
     * we could access to any of these positions (Ex: v(0, -2, 1) = 3;) and
     * actually any try to access to a position related to 5 (Ex: v(0, 4, 1) =
     * 3;), although it physically exists, is not logically correct and hence it
     * will throw an exception. The startingX and finishingX positions for this
     * sample vector are -1 and 1 respectively, for Y are -2 and 2 and for Z are
     * 0 and 1. The "for" iterations through the matrix should include these two
     * values if you want to cover the whole matrix.
     *
     * @code
     * for (int k=STARTINGZ(V); k<=FINISHINGZ(V); k++)
     *     for (int i=STARTINGY(V); i<=FINISHINGY(V); i++)
     *         for (int j=STARTINGX(V); j<=FINISHINGX(V); j++)
     *             VOL_ELEM(V, k, i, j) += 1;
     * @endcode
     */

    /** Produce an array suitable for working with Numerical Recipes.
     * @ingroup VolumesSizeShape
     *
     * This function must be used only as a preparation for routines which need
     * that the first physical index is 1 and not 0 as it usually is in C. New
     * memory is needed to hold the new double pointer array.
     */
    T*** adaptForNumericalRecipes() const
    {
        T*** m = NULL;
        ask_Tvolume(m, 1, ZSIZE(*this), 1, YSIZE(*this), 1, XSIZE(*this));

        FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(*this)
        m[k+1][i+1][j+1] = DIRECT_VOL_ELEM(*this, k, i, j);

        return m;
    }

    /** Kill an array produced for numerical recipes.
     * @ingroup VolumesSizeShape
     */
    void killAdaptationForNumericalRecipes(T*** m) const
    {
        free_Tvolume(m, 1, ZSIZE(*this), 1, YSIZE(*this), 1, XSIZE(*this));
    }

    /** Resize to a given size
     * @ingroup VolumesSizeShape
     */
    void resize(int Zdim, int Ydim, int Xdim)
    {
    	MultidimArray<T>::resize(Zdim, Ydim, Xdim);
    }

    /** Resize taking the shape from another volume
     * @ingroup VolumesSizeShape
     */
    template <typename T1>
    void resize(const Matrix3D<T1> &M)
    {
    	MultidimArray<T>::resize(M);
    }

    /** Resize taking the shape from another volume which
        is given as a MultidimArray
     * @ingroup VolumesSizeShape
     */
    template <typename T1>
    void resize(const MultidimArray<T1> &M)
    {
    	MultidimArray<T>::resize(M);
    }

    /** Outside.
     * @ingroup VolumesSizeShape
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

    /** Outside
     * @ingroup VolumesSizeShape
     *
     * TRUE if the logical index given is outside the definition region of this
     * array.
     */
    bool outside(const Matrix1D<double> &r) const
    {
        if (XSIZE(r) < 3)
            REPORT_ERROR(1, "Outside: index vector has not got enough components");
    
        return (XX(r) < STARTINGX(*this) || XX(r) > FINISHINGX(*this) ||
                YY(r) < STARTINGY(*this) || YY(r) > FINISHINGY(*this) ||
                ZZ(r) < STARTINGZ(*this) || ZZ(r) > FINISHINGZ(*this));
    }

    /** Returns the volume dimension.
     * @ingroup VolumesSizeShape
     *
     * Pay attention to the dimension order (Z,Y,X).
     *
     * @code
     * V.getDimension(Zdim, Ydim, Xdim);
     * @endcode
     */
    void getDimension(int& Ydim, int& Xdim, int& Zdim) const
    {
        Xdim = XSIZE(*this);
        Ydim = YSIZE(*this);
        Zdim = ZSIZE(*this);
    }

    /** @defgroup VolumesMemory Memory access
     * @ingroup Volumes
     *
     * This functions allows you to access to the matrix elements.
     */

    /** Volume element access via index.
     * @ingroup VolumesMemory
     *
     * Returns the value of a matrix logical position. In our example we could
     * access from v(0,-2,-1) to v(1,2,1). The elements can be used either by
     * value or by reference. An exception is thrown if the index is outside
     * the logical range. Be careful that the argument order is (Z,Y,X).
     *
     * @code
     * V(0, -2, 1) = 1;
     * val = V(0, -2, 1);
     * @endcode
     */
    T& operator()(int k, int i, int j) const
    {
        if (k < STARTINGZ(*this) || k > FINISHINGZ(*this))
            REPORT_ERROR(1203,static_cast< std::string > (
                         "Matrix3D::operator (): Matrix3D subscript (k) out of range k=")+
                         integerToString(k));

        if (i < STARTINGY(*this) || i > FINISHINGY(*this))
            REPORT_ERROR(1203,static_cast< std::string > (
                         "Matrix3D::operator (): Matrix3D subscript (i) out of range i=")+
                         integerToString(i));
        if (j < STARTINGX(*this) || j > FINISHINGX(*this))
            REPORT_ERROR(1203,static_cast< std::string > (
                         "Matrix3D::operator (): Matrix3D subscript (j) out of range j=")+
                         integerToString(j));

        return VOL_ELEM(*this, k, i, j);
    }

    /** Get the voxel at (k,i,j) (logical access).
     * @ingroup VolumesMemory
     */
    T getVoxel(int k, int i , int j) const
    {
        return (*this)(k, i, j);
    }

    /** Set the voxel at (k,i,j) (logical access).
     * @ingroup VolumesMemory
     */
    void setVoxel(int k, int i, int j, T val)
    {
        (*this)(k, i, j) = val;
    }

    /** Volume element access via double vector.
     * @ingroup VolumesMemory
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
    T& operator()(const Matrix1D< double >& v) const
    {
        return VOL_ELEM((*this), ROUND(ZZ(v)),
                        ROUND(YY(v)), ROUND(XX(v)));
    }

    /** Volume element access via integer vector.
     * @ingroup VolumesMemory
     */
    T& operator()(const Matrix1D< int >& v) const
    {
        return VOL_ELEM((*this), ZZ(v), YY(v), XX(v));
    }

    /** Interpolates the value of the 3D matrix M at the point (x,y,z).
     * @ingroup VolumesMemory
     *
     * (x,y,z) are in logical coordinates.
     */
    T interpolatedElement(double x, double y, double z, T outside_value = (T) 0)
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

        T d000 = (outside(z0, y0, x0)) ? outside_value : VOL_ELEM(*this, z0, y0, x0);
        T d001 = (outside(z0, y0, x1)) ? outside_value : VOL_ELEM(*this, z0, y0, x1);
        T d010 = (outside(z0, y1, x0)) ? outside_value : VOL_ELEM(*this, z0, y1, x0);
        T d011 = (outside(z0, y1, x1)) ? outside_value : VOL_ELEM(*this, z0, y1, x1);
        T d100 = (outside(z1, y0, x0)) ? outside_value : VOL_ELEM(*this, z1, y0, x0);
        T d101 = (outside(z1, y0, x1)) ? outside_value : VOL_ELEM(*this, z1, y0, x1);
        T d110 = (outside(z1, y1, x0)) ? outside_value : VOL_ELEM(*this, z1, y1, x0);
        T d111 = (outside(z1, y1, x1)) ? outside_value : VOL_ELEM(*this, z1, y1, x1);

        double dx00 = LIN_INTERP(fx, (double) d000, (double) d001);
        double dx01 = LIN_INTERP(fx, (double) d100, (double) d101);
        double dx10 = LIN_INTERP(fx, (double) d010, (double) d011);
        double dx11 = LIN_INTERP(fx, (double) d110, (double) d111);
        double dxy0 = LIN_INTERP(fy, (double) dx00, (double) dx10);
        double dxy1 = LIN_INTERP(fy, (double) dx01, (double) dx11);

        return (T) LIN_INTERP(fz, dxy0, dxy1);
    }

    /** Interpolates the value of the 3D matrix M at the point (x,y,z) knowing
     * that this image is a set of B-spline coefficients.
     * @ingroup VolumesMemory
     *
     * (x,y,z) are in logical coordinates.
     */
    T interpolatedElementBSpline(double x, double y, double z,
                                   int SplineDegree = 3)
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
        for (int n = n1; n <= n2; n++) {
	    int equivalent_n=n;
	    if      (n<0)             equivalent_n=-n-1;
	    else if (n>=ZSIZE(*this)) equivalent_n=2*ZSIZE(*this)-n-1;
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
                    double Coeff = (double) DIRECT_VOL_ELEM(*this,
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

    /** Slice access for reading.
     * @ingroup VolumesMemory
     *
     * This function returns a slice (a matrix) corresponding to the choosen
     * slice inside matrix, the numbering of the slices is also logical not
     * physical. This function differs from the previous one in that this one
     * cuts and assign in a single step instead of in two steps, as in
     * the previous example.
     *
     * @code
     * V.slice(0, m);
     * @endcode
     */
    void getSlice(int k, Matrix2D<T>& M, char axis = 'Z') const
    {
        if (XSIZE(*this) == 0)
        {
            M.clear();
            return;
        }

        switch (axis)
        {
        case 'Z':
            if (k < STARTINGZ(*this) || k > FINISHINGZ(*this))
                REPORT_ERROR(1203,
                             "Slice: Matrix3D subscript (k) out of range");

            k = k - STARTINGZ(*this);
            M.resize(YSIZE(*this), XSIZE(*this));
            FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(M)
                DIRECT_MAT_ELEM(M, i, j) = DIRECT_VOL_ELEM(*this, k, i, j);
            STARTINGX(M) = STARTINGX(*this);
            STARTINGY(M) = STARTINGY(*this);
            break;
        case 'Y':
            if (k < STARTINGY(*this) || k > FINISHINGY(*this))
                REPORT_ERROR(1203,
                             "Slice: Matrix3D subscript (i) out of range");

            k = k - STARTINGY(*this);
            M.resize(ZSIZE(*this), XSIZE(*this));
            FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(M)
                DIRECT_MAT_ELEM(M, i, j) = DIRECT_VOL_ELEM(*this, i, k, j);
            STARTINGX(M) = STARTINGX(*this);
            STARTINGY(M) = STARTINGZ(*this);
            break;
        case 'X':
            if (k < STARTINGX(*this) || k > FINISHINGX(*this))
                REPORT_ERROR(1203,
                             "Slice: Matrix3D subscript (j) out of range");

            k = k - STARTINGX(*this);
            M.resize(ZSIZE(*this), YSIZE(*this));
            FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(M)
                DIRECT_MAT_ELEM(M, i, j) = DIRECT_VOL_ELEM(*this, i, j, k);
            STARTINGX(M) = STARTINGY(*this);
            STARTINGY(M) = STARTINGZ(*this);
            break;
        default:
            REPORT_ERROR(1205,
                         (std::string) "Slice: not supported axis " + axis);
        }
    }

    /** Slice access for writing.
     * @ingroup VolumesMemory
     *
     * This function sets a matrix corresponding to the choosen slice inside
     * volume, the numbering of the slices is also logical not physical.
     *
     * @code
     * // Copies slice 0 in slice 1
     * V.setSlice(1, (V.slice(0)));
     * @endcode
     */
    void setSlice(int k, const Matrix2D<T>& v)
    {
        if (XSIZE(*this) == 0)
            return;

        if (k < STARTINGZ(*this) || k > FINISHINGZ(*this))
            REPORT_ERROR(1203,
                         "setSlice: Matrix3D subscript (k) out of range");

        if (v.rowNumber() != YSIZE(*this) || v.colNumber() != XSIZE(*this))
            REPORT_ERROR(1202,
                         "setSlice: Matrix3D dimensions different from the matrix ones");

        k = k - STARTINGZ(*this);

        FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(v)
            DIRECT_VOL_ELEM(*this, k, i, j) = DIRECT_MAT_ELEM(v, i, j);
    }

    /** Logical to physical index translation.
     * @ingroup VolumesMemory
     *
     * This function returns the physical position of a logical one.
     *
     * @code
     * m.toPhysical(k_log, i_log, j_log, k_phys, i_phys, j_phys);
     * @endcode
     */
    void toPhysical(int k_log, int i_log, int j_log,
                          int& k_phys, int& i_phys, int& j_phys) const
    {
        k_phys = k_log - STARTINGZ(*this);
        i_phys = i_log - STARTINGY(*this);
        j_phys = j_log - STARTINGX(*this);
    }

    /** Physical to logical index translation.
     * @ingroup VolumesMemory
     *
     * This function returns the logical position of a physical one.
     *
     * @code
     * m.toLogical(i_phys, j_phys, i_log, j_log);
     * @endcode
     */
    void toLogical(int k_phys, int i_phys, int j_phys,
                          int& k_log, int& i_log, int& j_log) const
    {
        k_log = k_phys + STARTINGZ(*this);
        i_log = i_phys + STARTINGY(*this);
        j_log = j_phys + STARTINGX(*this);
    }

    /// @defgroup VolumeOperators Operators
    /// @ingroup Matrices

    /** Assignment.
     * @ingroup VolumeOperators
     *
     * You can build as complex assignment expressions as you like. Multiple
     * assignment is allowed.
     *
     * @code
     * v1 = v2 + v3;
     * v1 = v2 = v3;
     * @endcode
     */
    Matrix3D<T>& operator=(const Matrix3D<T>& op1)
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
     * @ingroup VolumeOperators
     *
     * It is used to build arithmetic expressions. You can make a minus
     * of anything as long as it is correct semantically.
     */
    Matrix3D<T> operator-() const
    {
        Matrix3D<T> tmp(*this);
	T* ptr;
	unsigned long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(tmp,n,ptr)
            *ptr = -(*ptr);
        return tmp;
    }

    /** v3 = v1 + v2.
     * @ingroup VolumeOperators
     */
    Matrix3D<T> operator+(const Matrix3D<T>& op1) const
    {
        Matrix3D<T> tmp;
        arrayByArray(*this, op1, tmp, '+');
        return tmp;
    }

    /** v3 = v1 - v2.
     * @ingroup VolumeOperators
     */
    Matrix3D<T> operator-(const Matrix3D<T>& op1) const
    {
        Matrix3D<T> tmp;
        arrayByArray(*this, op1, tmp, '-');
        return tmp;
    }

    /** v3 = v1 * v2.
     * @ingroup VolumeOperators
     */
    Matrix3D<T> operator*(const Matrix3D<T>& op1) const
    {
        Matrix3D<T> tmp;
        arrayByArray(*this, op1, tmp, '*');
        return tmp;
    }

    /** v3 = v1 / v2.
     * @ingroup VolumeOperators
     */
    Matrix3D<T> operator/(const Matrix3D<T>& op1) const
    {
        Matrix3D<T> tmp;
        arrayByArray(*this, op1, tmp, '/');
        return tmp;
    }

    /** v3 += v2.
     * @ingroup VolumeOperators
     */
    void operator+=(const Matrix3D<T>& op1)
    {
        arrayByArray(*this, op1, *this, '+');
    }

    /** v3 -= v2.
     * @ingroup VolumeOperators
     */
    void operator-=(const Matrix3D<T>& op1)
    {
        arrayByArray(*this, op1, *this, '-');
    }

    /** v3 *= v2.
     * @ingroup VolumeOperators
     */
    void operator*=(const Matrix3D<T>& op1)
    {
        arrayByArray(*this, op1, *this, '*');
    }

    /** v3 /= v2.
     * @ingroup VolumeOperators
     */
    void operator/=(const Matrix3D<T>& op1)
    {
        arrayByArray(*this, op1, *this, '/');
    }

    /** v3 = v1 + k.
     * @ingroup VolumeOperators
     */
    Matrix3D<T> operator+(T op1) const
    {
        Matrix3D<T> tmp;
        arrayByScalar(*this, op1, tmp, '+');
        return tmp;
    }

    /** v3 = v1 - k.
     * @ingroup VolumeOperators
     */
    Matrix3D<T> operator-(T op1) const
    {
        Matrix3D<T> tmp;
        arrayByScalar(*this, op1, tmp, '-');
        return tmp;
    }

    /** v3 = v1 * k.
     * @ingroup VolumeOperators
     */
    Matrix3D<T> operator*(T op1) const
    {
        Matrix3D<T> tmp;
        arrayByScalar(*this, op1, tmp, '*');
        return tmp;
    }

    /** v3 = v1 / k.
     * @ingroup VolumeOperators
     */
    Matrix3D<T> operator/(T op1) const
    {
        Matrix3D<T> tmp;
        arrayByScalar(*this, op1, tmp, '/');
        return tmp;
    }

    /** v3 += k.
     * @ingroup VolumeOperators
     *
     * This function is not ported to Python.
     */

    void operator+=(const T& op1)
    {
        arrayByScalar(*this, op1, *this, '+');
    }

    /** v3 -= k.
     * @ingroup VolumeOperators
     *
     * This function is not ported to Python.
     */
    void operator-=(const T& op1)
    {
        arrayByScalar(*this, op1, *this, '-');
    }

    /** v3 *= k.
     * @ingroup VolumeOperators
     *
     * This function is not ported to Python.
     */
    void operator*=(const T& op1)
    {
        arrayByScalar(*this, op1, *this, '*');
    }

    /** v3 /= k.
     * @ingroup VolumeOperators
     *
     * This function is not ported to Python.
     */
    void operator/=(const T& op1)
    {
        arrayByScalar(*this, op1, *this, '/');
    }

    /** v3 = k + v2.
     * @ingroup VolumeOperators
     */
    friend Matrix3D<T> operator+(T op1, const Matrix3D<T>& op2)
    {
        Matrix3D<T> tmp;
        scalarByArray(op1, op2, tmp, '+');
        return tmp;
    }

    /** v3 = k - v2.
     * @ingroup VolumeOperators
     */
    friend Matrix3D<T> operator-(T op1, const Matrix3D<T>& op2)
    {
        Matrix3D<T> tmp;
        scalarByArray(op1, op2, tmp, '-');
        return tmp;
    }

    /** v3 = k * v2.
     * @ingroup VolumeOperators
     */
    friend Matrix3D<T> operator*(T op1, const Matrix3D<T>& op2)
    {
        Matrix3D<T> tmp;
        scalarByArray(op1, op2, tmp, '*');
        return tmp;
    }

    /** v3 = k / v2
     * @ingroup VolumeOperators
     */
    friend Matrix3D<T> operator/(T op1, const Matrix3D<T>& op2)
    {
        Matrix3D<T> tmp;
        scalarByArray(op1, op2, tmp, '/');
        return tmp;
    }

    /** @defgroup VolumesUtilites Utilities
     * @ingroup Volumes
     */

    /** Put a window to volume.
     * @ingroup VolumesUtilites
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
    void window(int z0, int y0, int x0, int zF, int yF, int xF, T init_value = 0)
    {
        Matrix3D<T> result(zF - z0 + 1, yF - y0 + 1, xF - x0 + 1);
        result.zinit = z0;
        result.yinit = y0;
        result.xinit = x0;

        for (int k = z0; k <= zF; k++)
            for (int i = y0; i <= yF; i++)
                for (int j = x0; j <= xF; j++)
                    if ((k >= STARTINGZ(*this) && k <= FINISHINGZ(*this)) &&
                        (i >= STARTINGY(*this) && i <= FINISHINGY(*this)) &&
                        (j >= STARTINGX(*this) && j <= FINISHINGX(*this)))
                        VOL_ELEM(result, k, i, j) = VOL_ELEM(*this, k, i, j);
                    else
                        VOL_ELEM(result, k, i, j) = init_value;

        *this = result;
    }

    /** Computes the center of mass of a volume.
     * @ingroup VolumesUtilites
     */
    void centerOfMass(Matrix1D< double >& center, void * mask=NULL)
    {
	center.initZeros(3);
	double mass = 0;
	Matrix3D< int >* imask = (Matrix3D< int >*) mask;

	FOR_ALL_ELEMENTS_IN_MATRIX3D(*this)
	{
            if ((imask == NULL || VOL_ELEM(*imask, k, i, j)) &&
		VOL_ELEM(*this, k, i, j) > 0)
            {
        	XX(center) += j * VOL_ELEM(*this, k, i, j);
        	YY(center) += i * VOL_ELEM(*this, k, i, j);
        	ZZ(center) += k * VOL_ELEM(*this, k, i, j);

        	mass += VOL_ELEM(*this, k, i, j);
            }
	}

	if (mass != 0)
            center /= mass;
    }

    /** Adjust the range of the array to a given one within a mask.
     * @ingroup VolumesUtilites
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
    void rangeAdjust(T minF, T maxF, Matrix3D<int> &mask)
    {
        if (MULTIDIM_SIZE(*this) <= 0)
            return;

        double min0, max0;
		bool first=true;
		FOR_ALL_ELEMENTS_IN_MATRIX3D(*this)
		{
			if (mask(k,i,j))
			{
				T val=(*this)(k,i,j);
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
     * @ingroup VolumesUtilites
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

    /** Adjust the range of the array to a given one.
     * @ingroup VolumesUtilites
     *
     * A linear operation is performed on the values of the array such that
     * after it, the values of the self array are as similar as possible
     * (L2 sense) to the values of the array shown as sample
     */
    void rangeAdjust(const Matrix3D<T> &example,
        const Matrix3D<int> *mask=NULL)
    {
        MultidimArray<T>::rangeAdjust(example,mask);
    }

    /** @defgroup VolumesGeometrical Geometrical Transformations
     * @ingroup Volumes
     *
     * In all geometrical transformations a periodic extension of the volume
     * is supposed, ie, if a voxel goes out on the left, it is entering on
     * the right, ...
     */

    /** Applies a geometrical transformation.
     * @ingroup VolumesGeometrical
     *
     * Any geometrical transformation defined by the matrix A (double (4x4)!!
     * ie, in homogeneous R3 coordinates) is applied to the volume V1.
     * The result is stored in V2 (it cannot be the same as the input volume).
     * An exception is thrown if the transformation matrix is not 4x4.
     *
     * Structure of the transformation matrix: It should have the following
     * components
     *
     * r11 r12 r13 x
     * r21 r22 r23 y
     * r31 r32 r33 z
     * 0   0   0   1
     *
     * where (x,y,z) is the translation desired, and Rij are the components of
     * the rotation matrix R. If you want to apply a scaling factor to the
     * transformation, then multiply r11, r22 and r33 by it.
     *
     * The result volume is resized to the same dimensions as V1 if V2 is empty
     * (0x0) at the beginning, if it is not, ie, if V2 has got some size
     * then only those values in the volume are filled, this is very
     * useful for resizing the volume, then you manually resize the output
     * volume to the desired size and then call this routine.
     *
     * The relationship between the output coordinates and the input ones are
     *
     * @code
     * out = A * in
     * (x, y, z) = A * (x', y', z')
     * @endcode
     *
     * This function works independently from the logical indexing of each
     * matrix, it sets the logical center and the physical center of the image
     * and work with these 2 coordinate spaces. At the end the original logical
     * indexing of each matrix is kept.
     *
     * The procedure followed goes from coordinates in the output volume
     * to the ones in the input one, so the inverse of the A matrix is
     * needed. There is a flag telling if the given matrix is already
     * the inverse one or the normal one. If it is the normal one internally
     * the matrix is inversed. If you are to do many "rotations" then
     * some time is spent in inverting the matrix. Normally the matrix is the
     * normal one.
     *
     * There is something else to tell about the geometrical tranformation.
     * The value of the voxel in the output volume is computed via
     * bilinear interpolation in the input volume. If any of the voxels
     * participating in the interpolation falls outside the input volume,
     * then automatically the corresponding output voxel is set to 0, unless
     * that the wrap flag has been set to 1. In this case if the voxel
     * falls out by the right hand then it is "wrapped" and the corresponding
     * voxel in the left hand is used. The same is appliable to top-bottom.
     * Usually wrap mode is off. Wrap mode is interesting for translations
     * but not for rotations, for example.
     *
     * The inverse mode and wrapping mode should be taken by default by the
     * routine, g++ seems to have problems with template functions outside
     * a class with default parameters. So, I'm sorry, you will have to
     * put them always. The usual combination is
     *
     * applyGeometry(..., IS_NOT_INV, DONT_WRAP).
     *
     * Although you can also use the constants IS_INV, or WRAP.
     *
     * @code
     * Matrix2D< double > A(4,4);
     * A.initIdentity;
     * applyGeometry(V2, A, V1);
     * @endcode
     */
    friend void applyGeometry<>(Matrix3D<T>& V2, const Matrix2D< double > &A,
                             const Matrix3D<T>& V1, bool inv, bool wrap, T outside);

    /** Apply geom with B-spline interpolation.
     * @ingroup VolumesGeometrical
     */
    friend void applyGeometryBSpline<>(Matrix3D<T>& V2, const Matrix2D< double > &A,
                                     const Matrix3D<T>& V1, int Splinedegree,
                                     bool inv, bool wrap, T outside);

    /** Self apply geom.
     * @ingroup VolumesGeometrical
     *
     * As apply geometry, but the result is kept in this object
     */
    void selfApplyGeometry(const Matrix2D< double > &A, bool inv, bool wrap, T outside = 0)
    {
        Matrix3D<T> aux;
        applyGeometry(aux, A, *this, inv, wrap, outside);
        *this = aux;
    }

    /** Self apply geom Bspline.
     * @ingroup VolumesGeometrical
     */
    void selfApplyGeometryBSpline(const Matrix2D< double > &A, int SplineDegree,
                                 bool inv, bool wrap, T outside = 0)
    {
        Matrix3D<T> aux;
        applyGeometryBSpline(aux, A, *this, SplineDegree, inv, wrap, outside);
        *this = aux;
    }

    /** Rotate a volume around system axis.
     * @ingroup VolumesGeometrical
     *
     * The rotation angle is in degrees, and the rotational axis is either
     * 'X', 'Y' or 'Z'. An exception is thrown if the axis given is not one of
     * these.
     *
     * @code
     * V2 = V1.rotate(60);
     * @endcode
     */
    void rotate(double ang, char axis, Matrix3D<T>& result,
        bool wrap = DONT_WRAP, T outside = 0) const
    {
        Matrix2D< double > tmp = rotation3DMatrix(ang, axis);
        applyGeometry(result, tmp, *this, IS_NOT_INV, wrap, outside);
    }

    /** Rotate a volume arounf system axis (BSpline).
     * @ingroup VolumesGeometrical
     */
    void rotateBSpline(int Splinedegree, double ang, char axis,
        Matrix3D<T>& result, bool wrap = DONT_WRAP, T outside = 0) const
    {
        Matrix2D< double > temp = rotation3DMatrix(ang, axis);
        applyGeometryBSpline(result, temp, *this, IS_NOT_INV, wrap, outside);
    }

    /** Rotate a volume around system axis, keep in this object.
     * @ingroup VolumesGeometrical
     */
    void selfRotate(double ang, char axis, bool wrap = DONT_WRAP)
    {
        Matrix3D<T> aux;
        rotate(ang, axis, aux, wrap);
        *this = aux;
    }

    /** Rotate a volume around system axis, keep in this object (Bspline).
     * @ingroup VolumesGeometrical
     */
    void selfRotateBSpline(int Splinedegree, double ang, char axis,
                             bool wrap = DONT_WRAP)
    {
        Matrix3D<T> aux;
        rotateBSpline(Splinedegree, ang, axis, aux, wrap);
        *this = aux;
    }

    /** Rotate a volume around any axis.
     * @ingroup VolumesGeometrical
     *
     * The rotation angle is in degrees, and the rotational axis is given as a
     * R3 vector. An exception is thrown if the axis is not a R3 vector. The
     * axis needs not to be unitary.
     *
     * @code
     * V2 = V1.rotate(60, vectorR3(1, 1, 1));
     * @endcode
     */
    void rotate(double ang, const Matrix1D< double >& axis, Matrix3D<T>& result,
                bool wrap = DONT_WRAP, T outside = 0) const
    {
        Matrix2D< double > tmp = rotation3DMatrix(ang, axis);
        applyGeometry(result, tmp, *this, IS_NOT_INV, wrap, outside);
    }

    /** Rotate a volume around any axis (Bspline).
     * @ingroup VolumesGeometrical
     */
    void rotateBSpline(int Splinedegree, double ang,
                        const Matrix1D< double >& axis, Matrix3D<T>& result,
                        bool wrap = DONT_WRAP, T outside = 0) const
    {
        Matrix2D< double > tmp = rotation3DMatrix(ang, axis);
        applyGeometryBSpline(result, tmp, *this, Splinedegree, IS_NOT_INV,
                           wrap, outside);
    }

    /** Rotate a volume around any axis, keep in this object.
     * @ingroup VolumesGeometrical
     */
    void selfRotate(double ang, const Matrix1D< double >& v,
                     bool wrap = DONT_WRAP)
    {
        Matrix3D<T> aux;
        rotate(ang, v, aux, wrap);
        *this = aux;
    }

    /** Rotate a volume around any axis, keep in this object (Bspline).
     * @ingroup VolumesGeometrical
     */
    void selfRotateBSpline(int Splinedegree, double ang,
                             const Matrix1D< double >& v, bool wrap = DONT_WRAP)
    {
        Matrix3D<T> aux;
        rotateBSpline(Splinedegree, ang, v, aux, wrap);
        *this = aux;
    }

    /** Translate a volume.
     * @ingroup VolumesGeometrical
     *
     * The shift is given as a R3 vector (shift_X, shift_Y, shift_Z).
     * An exception is thrown if the displacement is not a R3 vector.
     *
     * @code
     * // Displacement of 2 pixels down
     * V2 = V1.translate(vectorR3(0, 0, 2));
     * @endcode
     */
    void translate(const Matrix1D< double >& v, Matrix3D<T>& result,
        bool wrap = WRAP, T outside = 0)
    const
    {
        Matrix2D< double > tmp = translation3DMatrix(v);
        applyGeometry(result, tmp, *this, IS_NOT_INV, wrap, outside);
    }

    /** Translate a volume (Bspline).
     * @ingroup VolumesGeometrical
     */
    void translateBSpline(int Splinedegree, const Matrix1D< double >& v,
                           Matrix3D<T>& result, bool wrap = WRAP, T outside = 0) const
    {
        Matrix2D< double > tmp = translation3DMatrix(v);
        applyGeometryBSpline(result, tmp, *this, Splinedegree, IS_NOT_INV, wrap, outside);
    }

    /** Translate a volume, keep in this object.
     * @ingroup VolumesGeometrical
     */
    void selfTranslate(const Matrix1D< double >& v, bool wrap = WRAP, T outside = 0)
    {
        Matrix3D<T> aux;
        translate(v, aux, wrap, outside);
        *this = aux;
    }

    /** Translate a volume, keep in this object (Bspline).
     * @ingroup VolumesGeometrical
     */
    void selfTranslateBSpline(int Splinedegree, const Matrix1D< double >& v,
                                bool wrap = WRAP, T outside = 0)
    {
        Matrix3D<T> aux;
        translateBSpline(Splinedegree, v, aux, wrap, outside);
        *this = aux;
    }

    /** Translate center of mass to center.
     * @ingroup VolumesGeometrical
     *
     * If the input has very high values, sometimes it is better to rescale it
     * to be between 0 and 1.
     */
    void selfTranslateCenterOfMassToCenter(bool wrap = WRAP)
    {
        MultidimArray<T>::setXmippOrigin();
        Matrix1D< double > center;
        centerOfMass(center);
        center *= -1;
        selfTranslate(center, wrap);
    }

    /** Translate center of mass to center (Bspline).
     * @ingroup VolumesGeometrical
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

    /** Scales to a new size.
     * @ingroup VolumesGeometrical
     *
     * The volume is scaled (resampled) to fill a new size. It is not the
     * same as "window" in this same class. The size can be larger or smaller
     * than the actual one.
     *
     * @code
     * V2 = V1.scaleToSize(128, 128, 128);
     * @endcode
     */
    void scaleToSize(int Zdim, int Ydim, int Xdim, Matrix3D<T>& result) const
    {
        Matrix2D< double > tmp(4, 4);
        tmp.initIdentity();

        DIRECT_MAT_ELEM(tmp, 0, 0) = (double) Xdim / (double) XSIZE(*this);
        DIRECT_MAT_ELEM(tmp, 1, 1) = (double) Ydim / (double) YSIZE(*this);
        DIRECT_MAT_ELEM(tmp, 2, 2) = (double) Zdim / (double) ZSIZE(*this);

        result.resize(Zdim, Ydim, Xdim);

        applyGeometry(result, tmp, *this, IS_NOT_INV, WRAP);
    }

    /** Scales to a new size (Bspline).
     * @ingroup VolumesGeometrical
     */
    void scaleToSizeBSpline(int Splinedegree, int Zdim, int Ydim, int Xdim,
                               Matrix3D<T>& result) const
    {
        Matrix2D< double > tmp(4, 4);
        tmp.initIdentity();

        DIRECT_MAT_ELEM(tmp, 0, 0) = (double) Xdim / (double) XSIZE(*this);
        DIRECT_MAT_ELEM(tmp, 1, 1) = (double) Ydim / (double) YSIZE(*this);
        DIRECT_MAT_ELEM(tmp, 2, 2) = (double) Zdim / (double) ZSIZE(*this);

        result.resize(Zdim, Ydim, Xdim);

        applyGeometryBSpline(result, tmp, *this, Splinedegree, IS_NOT_INV, WRAP);
    }

    /** Scales to a new size, keep in this object.
     * @ingroup VolumesGeometrical
     */
    void selfScaleToSize(int Zdim, int Ydim, int Xdim)
    {
        Matrix3D<T> aux;
        scaleToSize(Zdim, Ydim, Xdim, aux);
        *this = aux;
    }

    /** Scales to a new size, keep in this object (Bspline).
     * @ingroup VolumesGeometrical
     */
    void selfScaleToSizeBSpline(int Splinedegree,
                                    int Zdim, int Ydim, int Xdim)
    {
        Matrix3D<T> aux;
        scaleToSizeBSpline(Splinedegree, Zdim, Ydim, Xdim, aux);
        *this = aux;
    }

    /** Reduce the image by 2 using a BSpline pyramid.
     * @ingroup VolumesGeometrical
     */
    void pyramidReduce(Matrix3D< double >& result, int levels = 1) const
    {
        Matrix3D< double > aux, aux2;
        MultidimArray<T>::produceSplineCoefficients(aux, 3);

        for (int i = 0; i < levels; i++)
        {
            aux.reduceBSpline(aux2, 3);
            aux = aux2;
        }

        aux2.produceImageFromSplineCoefficients(result, 3);
    }

    /** Expand the image by 2 using a BSpline pyramid.
     * @ingroup VolumesGeometrical
     */
    void pyramidExpand(Matrix3D< double >& result, int levels = 1) const
    {
        Matrix3D< double > aux, aux2;
        MultidimArray<T>::produceSplineCoefficients(aux, 3);

        for (int i = 0; i < levels; i++)
        {
            aux.expandBSpline(aux2, 3);
            aux = aux2;
        }

        aux2.produceImageFromSplineCoefficients(result, 3);
    }

    /** Expand a set of B-spline coefficients.
     * @ingroup VolumesGeometrical
     *
     * Knowing that this matrix is a set of B-spline coefficients, produce the
     * expanded set of B-spline coefficients using the two-scale relationship.
     */
    void expandBSpline(Matrix3D< double >& expanded, int SplineDegree = 3) const
    {
        double g[200]; // Coefficients of the reduce filter
        long ng; // Number of coefficients of the reduce filter
        double h[200]; // Coefficients of the expansion filter
        long nh; // Number of coefficients of the expansion filter
        short IsCentered; // Equal TRUE if the filter is a centered spline, FALSE otherwise */

        // Get the filter
        if (GetPyramidFilter("Centered Spline", SplineDegree, g, &ng, h, &nh,
                             &IsCentered))
            REPORT_ERROR(1, "Unable to load the filter coefficients");

        Matrix3D< double > aux;
        typeCast(*this, aux);
        expanded.resize(2 * ZSIZE(aux), 2 * YSIZE(aux), 2 * XSIZE(aux));
        Expand_3D(MULTIDIM_ARRAY(aux), XSIZE(aux), YSIZE(aux), ZSIZE(aux),
                  MULTIDIM_ARRAY(expanded), h, nh, IsCentered);
    }

    /** Reduce a set of B-spline coefficients.
     * @ingroup VolumesGeometrical
     *
     * Knowing that this matrix is a set of B-spline coefficients, produce the
     * reduced set of B-spline coefficients using the two-scale relationship.
     */
    void reduceBSpline(Matrix3D< double >& reduced, int SplineDegree = 3) const
    {
        double g[200]; // Coefficients of the reduce filter
        long ng; // Number of coefficients of the reduce filter
        double h[200]; // Coefficients of the expansion filter
        long nh; // Number of coefficients of the expansion filter
        short IsCentered; // Equal TRUE if the filter is a centered spline, FALSE otherwise

        // Get the filter
        if (GetPyramidFilter("Centered Spline", SplineDegree, g, &ng, h, &nh,
                             &IsCentered))
            REPORT_ERROR(1, "Unable to load the filter coefficients");

        Matrix3D< double > aux;
        typeCast(*this, aux);

        if (XSIZE(aux) % 2 != 0 && YSIZE(aux) % 2 != 0 && ZSIZE(aux) % 2 != 0)
            aux.resize(ZSIZE(aux - 1), YSIZE(aux) - 1, XSIZE(aux) - 1);
        else if (XSIZE(aux) % 2 != 0 && YSIZE(aux) % 2 != 0 && ZSIZE(aux) % 2
                 == 0)
            aux.resize(ZSIZE(aux), YSIZE(aux) - 1, XSIZE(aux) - 1);
        else if (XSIZE(aux) % 2 != 0 && YSIZE(aux) % 2 == 0 && ZSIZE(aux) % 2
                 != 0)
            aux.resize(ZSIZE(aux) - 1, YSIZE(aux), XSIZE(aux) - 1);
        else if (XSIZE(aux) % 2 != 0 && YSIZE(aux) % 2 == 0 && ZSIZE(aux) % 2
                 == 0)
            aux.resize(ZSIZE(aux), YSIZE(aux), XSIZE(aux) - 1);
        else if (XSIZE(aux) % 2 == 0 && YSIZE(aux) % 2 != 0 && ZSIZE(aux) % 2
                 != 0)
            aux.resize(ZSIZE(aux) - 1, YSIZE(aux) - 1, XSIZE(aux));
        else if (XSIZE(aux) % 2 == 0 && YSIZE(aux) % 2 != 0 && ZSIZE(aux) % 2
                 == 0)
            aux.resize(ZSIZE(aux), YSIZE(aux) - 1, XSIZE(aux));
        else if (XSIZE(aux) % 2 == 0 && YSIZE(aux) % 2 == 0 && ZSIZE(aux) % 2
                 != 0)
            aux.resize(ZSIZE(aux) - 1, YSIZE(aux), XSIZE(aux));

        reduced.resize(ZSIZE(aux) / 2, YSIZE(aux) / 2, XSIZE(aux) / 2);

        Reduce_3D(MULTIDIM_ARRAY(aux), XSIZE(aux), YSIZE(aux), ZSIZE(aux),
                  MULTIDIM_ARRAY(reduced), g, ng, IsCentered);
    }

    /** Maximum element.
     * @ingroup VolumesGeometrical
     *
     * This function returns the index of the maximum element of an array.
     * array(k,i,j). Returns -1 if the array is empty
     */
    void maxIndex(int& kmax, int& imax, int& jmax) const
    {
        if (XSIZE(*this) == 0)
        {
            kmax = imax = jmax = -1;
            return;
        }

        kmax = STARTINGZ(*this);
        imax = STARTINGY(*this);
        jmax = STARTINGX(*this);
        T maxval = VOL_ELEM(*this, kmax, imax, jmax);

        FOR_ALL_ELEMENTS_IN_MATRIX3D(*this)
            if (VOL_ELEM(*this, k, i, j) > maxval)
            {
                maxval = VOL_ELEM(*this, k, i, j);
                kmax = k;
                imax = i;
                jmax = j;
            }
    }

    /** Minimum element.
     * @ingroup VolumesGeometrical
     *
     * This function returns the index of the minimum element of an array.
     * array(k,i,j). Returns -1 if the array is empty
     */
    void minIndex(int& kmin, int& imin, int& jmin) const
    {
        if (XSIZE(*this) == 0)
        {
            kmin = imin = jmin = -1;
            return;
        }

        kmin = STARTINGZ(*this);
        imin = STARTINGY(*this);
        jmin = STARTINGX(*this);
        T minval = VOL_ELEM(*this, kmin, imin, jmin);

        FOR_ALL_ELEMENTS_IN_MATRIX3D(*this)
            if (VOL_ELEM(*this, k, i, j) > minval)
            {
                minval = VOL_ELEM(*this, k, i, j);
                kmin = k;
                imin = i;
                jmin = j;
            }
    }
};

/** @defgroup VolumesRelated Related functions
 * @ingroup Volumes
 *
 * These functions are not methods of Matrix1D
 */

/** @defgroup VolumesMisc Miscellaneous
 * @ingroup VolumesRelated
 */

/** Volume equality.
 * @ingroup VolumesMisc */
template<typename T>
bool operator==(const Matrix3D<T>& op1, const Matrix3D<T>& op2)
{
    return op1.equal(op2);
}

/** Volume inequality.
 * @ingroup VolumesMisc */
template<typename T>
bool operator!=(const Matrix3D<T>& op1, const Matrix3D<T>& op2)
{
    return !(op1==op2);
}

/** Reduce both volumes to a common size.
 * @ingroup VolumesMisc
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
void cutToCommonSize(Matrix3D<T>& V1, Matrix3D<T>& V2)
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


/** Does a radial average of a volume, around the voxel where is the origin.
 * @ingroup VolumesMisc
 *
 * A vector radial_mean is returned where:
 * - the first element is the mean of the voxels whose
 *   distance to the origin is (0-1),
 * - the second element is the mean of the voxels
 *   whose distance to the origin is (1-2)
 * - and so on.
 *
 * A second vector radial_count is returned containing the number of voxels
 * over which each radial average was calculated.
 *
 * Sjors nov2003: if rounding=true, element=round(distance);
 * - so the first element is the mean of the voxels whose distance to the
 *   origin is (0.5-1.5),
 * - the second element is the mean of the voxels whose distance to the origin
 *   is (1.5-2.5)
 * - and so on.
 */
template<typename T>
void radialAverage(const Matrix3D< T >& m,
                    const Matrix1D< int >& center_of_rot,
                    Matrix1D< T >& radial_mean,
                    Matrix1D< int >& radial_count,
                    const bool& rounding = false)
{
    Matrix1D< double > idx(3);

    // First determine the maximum distance that one should expect, to set the
    // dimension of the radial average vector
    Matrix1D< int > distances(8);

    double z = STARTINGZ(m) - ZZ(center_of_rot);
    double y = STARTINGY(m) - YY(center_of_rot);
    double x = STARTINGX(m) - XX(center_of_rot);

    distances(0) = (int) floor(sqrt(x * x + y * y + z * z));
    x = FINISHINGX(m) - XX(center_of_rot);

    distances(1) = (int) floor(sqrt(x * x + y * y + z * z));
    y = FINISHINGY(m) - YY(center_of_rot);

    distances(2) = (int) floor(sqrt(x * x + y * y + z * z));
    x = STARTINGX(m) - XX(center_of_rot);

    distances(3) = (int) floor(sqrt(x * x + y * y + z * z));
    z = FINISHINGZ(m) - ZZ(center_of_rot);

    distances(4) = (int) floor(sqrt(x * x + y * y + z * z));
    x = FINISHINGX(m) - XX(center_of_rot);

    distances(5) = (int) floor(sqrt(x * x + y * y + z * z));
    y = STARTINGY(m) - YY(center_of_rot);

    distances(6) = (int) floor(sqrt(x * x + y * y + z * z));
    x = STARTINGX(m) - XX(center_of_rot);

    distances(7) = (int) floor(sqrt(x * x + y * y + z * z));

    int dim = (int) CEIL(distances.computeMax()) + 1;
    if (rounding)
        dim++;

    // Define the vectors
    radial_mean.resize(dim);
    radial_mean.initZeros();
    radial_count.resize(dim);
    radial_count.initZeros();

    // Perform the radial sum and count pixels that contribute to every
    // distance
    FOR_ALL_ELEMENTS_IN_MATRIX3D(m)
    {
        ZZ(idx) = k - ZZ(center_of_rot);
        YY(idx) = i - YY(center_of_rot);
        XX(idx) = j - XX(center_of_rot);

        // Determine distance to the center
        int distance;
        if (rounding)
            distance = (int) ROUND(idx.module());
        else
            distance = (int) floor(idx.module());

        // Sum te value to the pixels with the same distance
        radial_mean(distance) += m(k, i, j);

        // Count the pixel
        radial_count(distance)++;
    }

    // Perform the mean
    FOR_ALL_ELEMENTS_IN_MATRIX1D(radial_mean)
        radial_mean(i) /= (T) radial_count(i);
}

template<typename T>
std::ostream& operator<<(std::ostream& ostrm, const Matrix3D<T>& v)
{
    if (v.xdim == 0)
        ostrm << "NULL Matrix3D\n";
    else
        ostrm << std::endl;

    double max_val;
    int prec;
    if (typeid(T)!=typeid(std::complex<double>))
    {
        max_val = std::abs(DIRECT_VOL_ELEM(v, 0, 0, 0));
    	T* ptr;
	unsigned long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(v,n,ptr)
            max_val = XMIPP_MAX(max_val, std::abs(*ptr));
        prec = bestPrecision(max_val, 10);
    }

    for (int k = STARTINGZ(v); k <= FINISHINGZ(v); k++)
    {
        ostrm << "Slice No. " << k << std::endl;
        for (int i = STARTINGY(v); i <= FINISHINGY(v); i++)
        {
            for (int j = STARTINGX(v); j <= FINISHINGX(v); j++)
            {
                if (typeid(T) == typeid(std::complex<double>))
                    ostrm << VOL_ELEM(v, k, i, j) << ' ';
                else
                    ostrm << floatToString((double) VOL_ELEM(v, k, i, j), 10, prec) << ' ';
            }
            ostrm << std::endl;
        }
    }

    return ostrm;
}

// Specialization for complex matrices
std::ostream& operator<<(std::ostream& ostrm,
    const Matrix3D< std::complex<double> >& v);

//#define DEBUG
template<typename T>
void applyGeometry(Matrix3D<T>& V2, const Matrix2D< double > &A,
    const Matrix3D<T>& V1, bool inv, bool wrap, T outside)
{
    int m1, n1, o1, m2, n2, o2;
    double x, y, z, xp, yp, zp;
    double minxp, minyp, maxxp, maxyp, minzp, maxzp;
    int cen_x, cen_y, cen_z, cen_xp, cen_yp, cen_zp;

    // Weights in X,Y,Z directions for bilinear interpolation
    double wx, wy, wz;

    if ((XSIZE(A) != 4) || (YSIZE(A) != 4))
        REPORT_ERROR(1102,
                     "Apply_geom3D: geometrical transformation is not 4x4");

    if (A.isIdentity())
    {
        V2 = V1;
        return;
    }

    if (XSIZE(V1) == 0)
    {
        V2.clear();
        return;
    }

    Matrix2D<double> Ainv;
    const Matrix2D<double> * Aptr=&A;
    if (!inv)
    {
        Ainv = A.inv();
        Aptr=&Ainv;
    }
    const Matrix2D<double> &Aref=*Aptr;

    if (XSIZE(V2) == 0)
        V2.resize(V1);

    // For scalings the output matrix is resized outside to the final
    // size instead of being resized inside the routine with the
    // same size as the input matrix
    if (XSIZE(V2) == 0)
        V2.resize(V1);

    // Find center of Matrix3D
    cen_z = (int)(V2.zdim / 2);
    cen_y = (int)(V2.ydim / 2);
    cen_x = (int)(V2.xdim / 2);
    cen_zp = (int)(V1.zdim / 2);
    cen_yp = (int)(V1.ydim / 2);
    cen_xp = (int)(V1.xdim / 2);
    minxp = -cen_xp;
    minyp = -cen_yp;
    minzp = -cen_zp;
    maxxp = V1.xdim - cen_xp - 1;
    maxyp = V1.ydim - cen_yp - 1;
    maxzp = V1.zdim - cen_zp - 1;

#ifdef DEBUG
    std::cout << "Geometry 2 center=("
              << cen_z  << "," << cen_y  << "," << cen_x  << ")\n"
              << "Geometry 1 center=("
              << cen_zp << "," << cen_yp << "," << cen_xp << ")\n"
              << "           min=("
              << minzp  << "," << minyp  << "," << minxp  << ")\n"
              << "           max=("
              << maxzp  << "," << maxyp  << "," << maxxp  << ")\n"
    ;
#endif

    // Now we go from the output Matrix3D to the input Matrix3D, ie, for any
    // voxel in the output Matrix3D we calculate which are the corresponding
    // ones in the original Matrix3D, make an interpolation with them and put
    // this value at the output voxel

    // V2 is not initialised to 0 because all its pixels are rewritten
    for (int k = 0; k < V2.zdim; k++)
        for (int i = 0; i < V2.ydim; i++)
        {
            // Calculate position of the beginning of the row in the output
            // Matrix3D
            x = -cen_x;
            y = i - cen_y;
            z = k - cen_z;

            // Calculate this position in the input image according to the
            // geometrical transformation they are related by
            // coords_output(=x,y) = A * coords_input (=xp,yp)
            xp = x * dMij(Aref, 0, 0) + y * dMij(Aref, 0, 1) + z * dMij(Aref, 0, 2)
                 + dMij(Aref, 0, 3);
            yp = x * dMij(Aref, 1, 0) + y * dMij(Aref, 1, 1) + z * dMij(Aref, 1, 2)
                 + dMij(Aref, 1, 3);
            zp = x * dMij(Aref, 2, 0) + y * dMij(Aref, 2, 1) + z * dMij(Aref, 2, 2)
                 + dMij(Aref, 2, 3);

            for (int j = 0; j < V2.xdim; j++)
            {
                bool interp;
                T tmp;

#ifdef DEBUG
                bool show_debug = false;
                if ((i == 0 && j == 0 && k == 0) ||
                    (i == V2.ydim - 1 && j == V2.xdim - 1 && k == V2.zdim - 1))
                    show_debug = true;

                if (show_debug)
                    std::cout << "(x,y,z)-->(xp,yp,zp)= "
                    << "(" << x  << "," << y  << "," << z  << ") "
                    << "(" << xp << "," << yp << "," << zp << ")\n";
#endif

                // If the point is outside the volume, apply a periodic
                // extension of the volume, what exits by one side enters by
                // the other
                interp  = true;
                if (wrap)
                {
                    if (xp < minxp - XMIPP_EQUAL_ACCURACY ||
                        xp > maxxp + XMIPP_EQUAL_ACCURACY)
                        xp = realWRAP(xp, minxp - 0.5, maxxp + 0.5);

                    if (yp < minyp - XMIPP_EQUAL_ACCURACY ||
                        yp > maxyp + XMIPP_EQUAL_ACCURACY)
                        yp = realWRAP(yp, minyp - 0.5, maxyp + 0.5);

                    if (zp < minzp - XMIPP_EQUAL_ACCURACY ||
                        zp > maxzp + XMIPP_EQUAL_ACCURACY)
                        zp = realWRAP(zp, minzp - 0.5, maxzp + 0.5);
                }
                else
                {
                    if (xp < minxp - XMIPP_EQUAL_ACCURACY ||
                        xp > maxxp + XMIPP_EQUAL_ACCURACY)
                        interp = false;

                    if (yp < minyp - XMIPP_EQUAL_ACCURACY ||
                        yp > maxyp + XMIPP_EQUAL_ACCURACY)
                        interp = false;

                    if (zp < minzp - XMIPP_EQUAL_ACCURACY ||
                        zp > maxzp + XMIPP_EQUAL_ACCURACY)
                        interp = false;
                }

                if (interp)
                {
                    // Calculate the integer position in input volume, be
                    // careful that it is not the nearest but the one at the
                    // top left corner of the interpolation square. Ie,
                    // (0.7,0.7) would give (0,0)
                    // Calculate also weights for point m1+1,n1+1
                    wx = xp + cen_xp;
                    m1 = (int) wx;
                    wx = wx - m1;
                    m2 = m1 + 1;
                    wy = yp + cen_yp;
                    n1 = (int) wy;
                    wy = wy - n1;
                    n2 = n1 + 1;
                    wz = zp + cen_zp;
                    o1 = (int) wz;
                    wz = wz - o1;
                    o2 = o1 + 1;

#ifdef DEBUG
                    if (show_debug)
                    {
                        std::cout << "After wrapping(xp,yp,zp)= "
                                  << "(" << xp << "," << yp << "," << zp << ")\n";
                        std::cout << "(m1,n1,o1)-->(m2,n2,o2)="
                                  << "(" << m1 << "," << n1 << "," << o1 << ") "
                                  << "(" << m2 << "," << n2 << "," << o2 << ")\n";
                        std::cout << "(wx,wy,wz)="
                                  << "(" << wx << "," << wy << "," << wz << ")\n";
                    }
#endif

                    // Perform interpolation
                    // if wx == 0 means that the rightest point is useless for
                    // this interpolation, and even it might not be defined if
                    // m1=xdim-1
                    // The same can be said for wy.
                    tmp  = (T)((1 - wz) * (1 - wy) * (1 - wx) * dVkij(V1, o1, n1,
                               m1));

                    if (wx != 0 && m2 < V1.xdim)
                        tmp += (T)((1 - wz) * (1 - wy) * wx * dVkij(V1, o1, n1,
                                   m2));

                    if (wy != 0 && n2 < V1.ydim)
                    {
                        tmp += (T)((1 - wz) * wy * (1 - wx) * dVkij(V1, o1, n2,
                                   m1));
                        if (wx != 0 && m2 < V1.xdim)
                            tmp += (T)((1 - wz) * wy * wx * dVkij(V1, o1, n2,
                                                                  m2));
                    }

                    if (wz != 0 && o2 < V1.zdim)
                    {
                        tmp += (T)(wz * (1 - wy) * (1 - wx) * dVkij(V1, o2, n1,
                                   m1));
                        if (wx != 0 && m2 < V1.xdim)
                            tmp += (T)(wz * (1 - wy) * wx * dVkij(V1, o2, n1,
                                                                  m2));
                        if (wy != 0 && n2 < V1.ydim)
                        {
                            tmp += (T)(wz * wy * (1 - wx) * dVkij(V1, o2, n2,
                                                                  m1));
                            if (wx != 0 && m2 < V1.xdim)
                                tmp += (T)(wz * wy * wx * dVkij(V1, o2, n2,
                                                                m2));
                        }
                    }

                    dVkij(V2 , k, i, j) = tmp;
#ifdef DEBUG
                    if (show_debug)
                        std::cout <<
                        "tmp1=" << dVkij(V1, o1, n1, m1) << " " << (T)((1 - wz)
                                *(1 - wy) *(1 - wx) * dVkij(V1, o1, n1, m1)) <<
                        std::endl <<
                        "tmp2=" << dVkij(V1, o1, n1, m2) << " " << (T)((1 - wz)
                                *(1 - wy) * wx * dVkij(V1, o1, n1, m2)) <<
                        std::endl <<
                        "tmp3=" << dVkij(V1, o1, n2, m1) << " " << (T)((1 - wz)
                                * wy *(1 - wx) * dVkij(V1, o1, n2, m1)) <<
                        std::endl <<
                        "tmp4=" << dVkij(V1, o1, n2, m2) << " " << (T)((1 - wz)
                                * wy * wx * dVkij(V1, o1, n2, m2)) <<
                        std::endl <<
                        "tmp5=" << dVkij(V1, o2, n1, m1) << " " << (T)(wz *
                                (1 - wy) *(1 - wx) * dVkij(V1, o2, n1, m1))
                        << std::endl <<
                        "tmp6=" << dVkij(V1, o2, n1, m2) << " " << (T)(wz *
                                (1 - wy) * wx * dVkij(V1, o2, n1, m2)) <<
                        std::endl <<
                        "tmp7=" << dVkij(V1, o2, n2, m1) << " " << (T)(wz *
                                wy *(1 - wx) * dVkij(V1, o2, n2, m1)) <<
                        std::endl <<
                        "tmp8=" << dVkij(V1, o2, n2, m2) << " " << (T)(wz *
                                wy * wx * dVkij(V1, o2, n2, m2)) <<
                        std::endl <<
                        "tmp= " << tmp << std::endl;
#endif
                }
                else
                    dVkij(V2, k, i, j) = outside;


                // Compute new point inside input image
                xp += dMij(Aref, 0, 0);
                yp += dMij(Aref, 1, 0);
                zp += dMij(Aref, 2, 0);
            }
        }
}
#undef DEBUG

//#define DEBUG
template<typename T>
void applyGeometryBSpline(Matrix3D<T>& V2, const Matrix2D< double > &A,
    const Matrix3D<T>& V1, int Splinedegree, bool inv, bool wrap, T outside)
{
    int m1, n1, o1, m2, n2, o2;
    double x, y, z, xp, yp, zp;
    double minxp, minyp, maxxp, maxyp, minzp, maxzp;
    int   cen_x, cen_y, cen_z, cen_xp, cen_yp, cen_zp;

    if ((XSIZE(A) != 4) || (YSIZE(A) != 4))
        REPORT_ERROR(1102, "Apply_geom3D: geometrical transformation is not 4x4");

    if (A.isIdentity())
    {
        V2 = V1;
        return;
    }

    if (XSIZE(V1) == 0)
    {
        V2.clear();
        return;
    }

    Matrix2D<double> Ainv;
    const Matrix2D<double> * Aptr=&A;
    if (!inv)
    {
        Ainv = A.inv();
        Aptr=&Ainv;
    }
    const Matrix2D<double> &Aref=*Aptr;

    if (XSIZE(V2) == 0)
        V2.resize(V1);

    // For scalings the output matrix is resized outside to the final
    // size instead of being resized inside the routine with the
    // same size as the input matrix
    if (XSIZE(V2) == 0)
        V2.resize(V1);

    // Find center of Matrix3D
    cen_z = (int)(V2.zdim / 2);
    cen_y = (int)(V2.ydim / 2);
    cen_x = (int)(V2.xdim / 2);
    cen_zp = (int)(V1.zdim / 2);
    cen_yp = (int)(V1.ydim / 2);
    cen_xp = (int)(V1.xdim / 2);
    minxp = -cen_xp;
    minyp = -cen_yp;
    minzp = -cen_zp;
    maxxp = V1.xdim - cen_xp - 1;
    maxyp = V1.ydim - cen_yp - 1;
    maxzp = V1.zdim - cen_zp - 1;
#ifdef DEBUG
    std::cout << "Geometry 2 center=("
              << cen_z  << "," << cen_y  << "," << cen_x  << ")\n"
              << "Geometry 1 center=("
              << cen_zp << "," << cen_yp << "," << cen_xp << ")\n"
              << "           min=("
              << minzp  << "," << minyp  << "," << minxp  << ")\n"
              << "           max=("
              << maxzp  << "," << maxyp  << "," << maxxp  << ")\n"
    ;
#endif

    // Build the B-spline coefficients
    Matrix3D< double > Bcoeffs;
    V1.produceSplineCoefficients(Bcoeffs, Splinedegree);
    STARTINGX(Bcoeffs) = (int) minxp;
    STARTINGY(Bcoeffs) = (int) minyp;
    STARTINGZ(Bcoeffs) = (int) minzp;

    // Now we go from the output Matrix3D to the input Matrix3D, ie, for any
    // voxel in the output Matrix3D we calculate which are the corresponding
    // ones in the original Matrix3D, make an interpolation with them and put
    // this value at the output voxel

    // V2 is not initialised to 0 because all its pixels are rewritten
    for (int k = 0; k < V2.zdim; k++)
        for (int i = 0; i < V2.ydim; i++)
        {
            // Calculate position of the beginning of the row in the output
            // Matrix3D
            x = -cen_x;
            y = i - cen_y;
            z = k - cen_z;

            // Calculate this position in the input image according to the
            // geometrical transformation they are related by
            // coords_output(=x,y) = A * coords_input (=xp,yp)
            xp = x * dMij(Aref, 0, 0) + y * dMij(Aref, 0, 1) + z * dMij(Aref, 0, 2)
                 + dMij(Aref, 0, 3);
            yp = x * dMij(Aref, 1, 0) + y * dMij(Aref, 1, 1) + z * dMij(Aref, 1, 2)
                 + dMij(Aref, 1, 3);
            zp = x * dMij(Aref, 2, 0) + y * dMij(Aref, 2, 1) + z * dMij(Aref, 2, 2)
                 + dMij(Aref, 2, 3);

            for (int j = 0; j < V2.xdim; j++)
            {
                bool interp;
                T tmp;

#ifdef DEBUG
                bool show_debug = false;
                if ((i == 0 && j == 0 && k == 0) ||
                    (i == V2.ydim - 1 && j == V2.xdim - 1 && k == V2.zdim - 1))
                    show_debug = true;

                if (show_debug)
                    std::cout << "(x,y,z)-->(xp,yp,zp)= "
                    << "(" << x  << "," << y  << "," << z  << ") "
                    << "(" << xp << "," << yp << "," << zp << ")\n";
#endif

                // If the point is outside the volume, apply a periodic
                // extension of the volume, what exits by one side enters by
                // the other
                interp = true;
                if (wrap)
                {
                    if (xp < minxp - XMIPP_EQUAL_ACCURACY ||
                        xp > maxxp + XMIPP_EQUAL_ACCURACY)
                        xp = realWRAP(xp, minxp - 0.5, maxxp + 0.5);

                    if (yp < minyp - XMIPP_EQUAL_ACCURACY ||
                        yp > maxyp + XMIPP_EQUAL_ACCURACY)
                        yp = realWRAP(yp, minyp - 0.5, maxyp + 0.5);

                    if (zp < minzp - XMIPP_EQUAL_ACCURACY ||
                        zp > maxzp + XMIPP_EQUAL_ACCURACY)
                        zp = realWRAP(zp, minzp - 0.5, maxzp + 0.5);
                }
                else
                {
                    if (xp < minxp - XMIPP_EQUAL_ACCURACY ||
                        xp > maxxp + XMIPP_EQUAL_ACCURACY)
                        interp = false;

                    if (yp < minyp - XMIPP_EQUAL_ACCURACY ||
                        yp > maxyp + XMIPP_EQUAL_ACCURACY)
                        interp = false;

                    if (zp < minzp - XMIPP_EQUAL_ACCURACY ||
                        zp > maxzp + XMIPP_EQUAL_ACCURACY)
                        interp = false;
                }

                if (interp)
                {
                    dVkij(V2, k, i, j) =
                        (T) Bcoeffs.interpolatedElementBSpline(xp, yp, zp,
                                Splinedegree);
                }
                else
                    dVkij(V2, k, i, j) = outside;

                // Compute new point inside input image
                xp += dMij(Aref, 0, 0);
                yp += dMij(Aref, 1, 0);
                zp += dMij(Aref, 2, 0);
            }
        }
}

// Specific instantiations for complexes
template<>
std::complex<double> Matrix3D< std::complex< double> >::interpolatedElement(double x,
        double y, double z, std::complex< double > outside_value);

#endif
