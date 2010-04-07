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

#ifndef CORE_MULTIDIM_ARRAY_H
#define CORE_MULTIDIM_ARRAY_H

#include <typeinfo>
#include "funcs.h"
#include "error.h"
#include "args.h"

extern int bestPrecision(float F, int _width);
extern std::string floatToString(float F, int _width, int _prec);
extern std::string integerToString(int I, int _width, char fill_with);

/// @defgroup coreMultidimensionalArrays core Multidimensional Arrays
/// @ingroup Arrays

/** @defgroup coreMultidimArraysSpeedUp Speed up macros
 * @ingroup coreMultidimensionalArrays
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

/// @defgroup coreMultidimArraysSizeShape Size and shape
/// @ingroup coreMultidimArraysSpeedUp

/** Returns the first X valid logical index
 * @ingroup coreMultidimArraysSizeShape
 */
#define STARTINGX(v) ((v).xinit)

/** Returns the last X valid logical index
 * @ingroup coreMultidimArraysSizeShape
 */
#define FINISHINGX(v) ((v).xinit + (v).xdim - 1)

/** Returns the first Y valid logical index
 * @ingroup coreMultidimArraysSizeShape
 */
#define STARTINGY(v) ((v).yinit)

/** Returns the last Y valid logical index
 * @ingroup coreMultidimArraysSizeShape
 */
#define FINISHINGY(v) ((v).yinit + (v).ydim - 1)

/** Returns the first Z valid logical index
 * @ingroup coreMultidimArraysSizeShape
 */
#define STARTINGZ(v) ((v).zinit)

/** Returns the last Z valid logical index
 * @ingroup coreMultidimArraysSizeShape
 */
#define FINISHINGZ(v) ((v).zinit + (v).zdim - 1)

/** Access to X dimension (size)
 * @ingroup coreMultidimArraysSizeShape
 */
#define XSIZE(v) ((v).xdim)

/** Access to Y dimension (size)
 * @ingroup coreMultidimArraysSizeShape
 */
#define YSIZE(v) ((v).ydim)

/** Access to Z dimension (size)
 * @ingroup coreMultidimArraysSizeShape
 */
#define ZSIZE(v) ((v).zdim)

/** Access to N dimension (size)
 * @ingroup coreMultidimArraysSizeShape
 */
#define NSIZE(v) ((v).ndim)

/** Access to XY dimension (Ysize*Xsize)
 * @ingroup coreMultidimArraysSizeShape
 */
#define YXSIZE(v) ((v).yxdim)

/** Access to XYZ dimension (Zsize*Ysize*Xsize)
 * @ingroup coreMultidimArraysSizeShape
 */
#define ZYXSIZE(v) ((v).zyxdim)

/// FIXME MULTIDIM_SIZE will disappear
/** Access to XYZ dimension (Zsize*Ysize*Xsize)
 * @ingroup coreMultidimArraysSizeShape
 */
#define MULTIDIM_SIZE(v) ((v).nzyxdim)

/** Access to XYZN dimension (Nsize*Zsize*Ysize*Xsize)
 * @ingroup coreMultidimArraysSizeShape
 */
#define NZYXSIZE(v) ((v).nzyxdim)

/** Array access.
 * @ingroup coreMultidimArraysSizeShape
 *
 * This macro gives you access to the array (T **)
 */
#ifndef MULTIDIM_ARRAY
#define MULTIDIM_ARRAY(v) ((v).data)
#endif

/** Access to a direct element.
 * @ingroup coreMultidimArraysSizeShape.
 * v is the array, l is the image, k is the slice, i is the Y index and j is the X index.
 * i and j) within the slice.
 */
#define DIRECT_NZYX_ELEM(v, l, k, i, j) ((v).data[(l)*ZYXSIZE(v)+(k)*YXSIZE(v)+((i)*XSIZE(v))+(j)])

/** Multidim element: Logical access.
 * @ingroup coreMultidimArraysSizeShape.

 */
#define NZYX_ELEM(v, l, k, i, j)  \
    DIRECT_NZYX_ELEM((v), (l), (k) - STARTINGZ(v), (i) - STARTINGY(v), (j) - STARTINGX(v))

/** Access to a direct element.
 * @ingroup coreMultidimArraysSizeShape.
 * v is the array, k is the slice and n is the number of the pixel (combined
 * i and j) within the slice.
 */
#define DIRECT_MULTIDIM_ELEM(v,n) ((v).data[(n)])

/** Access to X component of a Matrix1D
 * @ingroup coreMultidimArraysSizeShape.
 *
 * @code
 * XX(v) = 1;
 * val = XX(v);
 * @endcode
 */
#define XX(v) DIRECT_VEC_ELEM(v, 0)

/** Access to Y component of a Matrix1D
 * @ingroup coreMultidimArraysSizeShape.
 *
 * @code
 * YY(v) = 1;
 * val = YY(v);
 * @endcode
 */
#define YY(v) DIRECT_VEC_ELEM(v, 1)

/** Access to Z component of a Matrix1D
 * @ingroup coreMultidimArraysSizeShape.
 *
 * @code
 * ZZ(v) = 1;
 * val = ZZ(v);
 * @endcode
 */
#define ZZ(v) DIRECT_VEC_ELEM(v, 2)

/** For all direct elements in the array
 * @ingroup coreMultidimArraysSizeShape
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

/** For all direct elements in the array
 * @ingroup coreMultidimArraysSizeShape
 *
 * This macro is used to generate loops for the array in an easy
 * manner. It defines internal indexes 'l', 'k','i' and 'j' which
 * ranges over the n volume using its physical definition.
 *
 * @code
 * FOR_ALL_DIRECT_NZYX_ELEMENTS_IN_MULTIDIMARRAY(v)
 * {
 *     std::cout << DIRECT_NZYX_ELEM(v,l, k, i, j) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_DIRECT_NZYX_ELEMENTS_IN_MULTIDIMARRAY(V) \
    for (int l=0; l<NSIZE(V); l++) \
        for (int k=0; k<ZSIZE(V); k++) \
            for (int i=0; i<YSIZE(V); i++)      \
                for (int j=0; j<XSIZE(V); j++)

/** For all direct elements in the array
 * @ingroup coreMultidimArraysSizeShape
 *
 * This macro is used to generate loops for the array in an easy
 * manner. It defines internal indexes 'l', 'k','i' and 'j' which
 * ranges over the n volume using its logical definition.
 *
 * @code
 * FOR_ALL_NZYX_ELEMENTS_IN_MULTIDIMARRAY(v)
 * {
 *     std::cout << NZYX_ELEM(v,l, k, i, j) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_NZYX_ELEMENTS_IN_MULTIDIMARRAY(V) \
    for (int l=0; l<NSIZE(V); l++) \
        for (int k=STARTINGZ(V); k<=FINISHINGZ(V); k++) \
            for (int i=STARTINGY(V); i<=FINISHINGY(V); i++)     \
                for (int j=STARTINGX(V); j<=FINISHINGX(V); j++)

/** For all direct elements in the array, pointer version
 * @ingroup coreMultidimArraysSizeShape
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



/// @defgroup VolumesSizeShape 3D matrices size and shape (assume ndim=1)
/// @ingroup coreMultidimArraysSizeShape.

/** Access to a direct element.
 * @ingroup coreMultidimArraysSizeShape.
 * v is the array, k is the slice (Z), i is the Y index and j is the X index.
 */
#define DIRECT_VOL_ELEM(v,k,i,j) ((v).data[(k)*YXSIZE(v)+((i)*XSIZE(v))+(j)])

/** A short alias for the previous function.
 * @ingroup coreMultidimArraysSizeShape.
 *
 */
#define dVkij(V, k, i, j) DIRECT_VOL_ELEM(V, k, i, j)

/** Volume element: Logical access.
 * @ingroup coreMultidimArraysSizeShape.
 *
 * @code
 * VOL_ELEM(V, -1, -2, 1) = 1;
 * val = VOL_ELEM(V, -1, -2, 1);
 * @endcode
 */
#define VOL_ELEM(V, k, i, j) \
    DIRECT_VOL_ELEM((V),(k) - STARTINGZ(V), (i) - STARTINGY(V), (j) - STARTINGX(V))

/** For all elements in the array.
 * @ingroup coreMultidimArraysSizeShape.
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
    for (int k=STARTINGZ(V); k<=FINISHINGZ(V); k++) \
        for (int i=STARTINGY(V); i<=FINISHINGY(V); i++) \
            for (int j=STARTINGX(V); j<=FINISHINGX(V); j++)

/** For all elements in the array between corners.
 * @ingroup coreMultidimArraysSizeShape.
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
 * @ingroup coreMultidimArraysSizeShape
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


/// @defgroup coreMatricesSizeShape 2D Matrices size and shape
/// @ingroup coreMultidimArraysSpeedUp

/** Access to a direct element of a matrix.
 * @ingroup coreMatricesSizeShape
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
 * @ingroup coreMatricesSizeShape
 */
#define dMij(M, i, j) DIRECT_MAT_ELEM(M, i, j)

/** Matrix element: Logical access
 * @ingroup coreMatricesSizeShape
 *
 * @code
 * MAT_ELEM(m, -2, 1) = 1;
 * val = MAT_ELEM(m, -2, 1);
 * @endcode
 */
#define MAT_ELEM(v, i, j) \
    DIRECT_MAT_ELEM(v, (i) - STARTINGY(v), (j) - STARTINGX(v))

/** TRUE if both arrays have the same shape
 * @ingroup coreMatricesSizeShape
 *
 * Two arrays have the same shape if they have the same size and the same
 * starting point. Be aware that this is a macro which simplifies to a boolean.
 */
#define SAME_SHAPE2D(v1, v2) \
    (XSIZE(v1) == XSIZE(v2) && \
     YSIZE(v1) == YSIZE(v2) && \
     STARTINGX(v1) == STARTINGX(v2) && \
     STARTINGY(v1) == STARTINGY(v2))

/** For all elements in the array
 * @ingroup coreMatricesSizeShape
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
    for (int i=STARTINGY(m); i<=FINISHINGY(m); i++) \
        for (int j=STARTINGX(m); j<=FINISHINGX(m); j++)

/** For all elements in the array between corners
 * @ingroup coreMatricesSizeShape
 *
 * This macro is used to generate loops for a matrix in an easy manner. It needs
 * an externally defined Matrix1D< double > r(2). Then YY(r) and XX(r) range
 * from (int) YY(corner1) to (int)YY(corner2), (int) XX(corner1) to (int)
 * XX(corner2) (included limits) respectively. Notice that corner1 and corner2
 * need only be Matrix1D.
 *
 * @code
 * Matrix1D< double > corner1(2), corner2(2);
 * Matrix1D< int > r(2);
 * XX(corner1) = -1;
 * XX(corner2) = 1;
 * YY(corner1) = -2;
 * YY(corner2) = 2;
 *
 * FOR_ALL_ELEMENTS_IN_MATRIX2D_BETWEEN(corner1, corner2)
 * {
 *     std::cout << v(r) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_ELEMENTS_IN_MATRIX2D_BETWEEN(corner1, corner2) \
    for (YY(r)=YY((corner1)); YY(r)<=YY((corner2)); YY(r)++) \
        for (XX(r)=XX((corner1)); XX(r)<=XX((corner2)); XX(r)++)

/** For all elements in common
 * @ingroup coreMatricesSizeShape
 *
 * This macro is used to generate loops for all the elements logically in common
 * between two images in an easy manner. Then i and j (locally defined) range
 * from MAX(STARTINGY(V1), STARTINGY(V2)) to MIN(FINISHINGY(V1),
 * FINISHINGY(V2)), MAX(STARTINGX(V1), STARTINGX(V2)) to MIN(FINISHINGX(V1),
 * FINISHINGX(V2)) (included limits) respectively. You need to define
 * SPEED_UP_temps.
 *
 * @code
 * Matrix2D< double > m1(10, 10), m2(20, 20);
 * m1.setXmippOrigin();
 * m2.setXmippOrigin();
 *
 * FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX2D(m1, m2)
 * {
 *     ...
 * }
 * @endcode
 */
#define FOR_ALL_ELEMENTS_IN_COMMON_IN_MATRIX2D(m1, m2) \
    ispduptmp2 = XMIPP_MAX(STARTINGY(m1), STARTINGY(m2)); \
    ispduptmp3 = XMIPP_MIN(FINISHINGY(m1), FINISHINGY(m2)); \
    ispduptmp4 = XMIPP_MAX(STARTINGX(m1), STARTINGX(m2)); \
    ispduptmp5 = XMIPP_MIN(FINISHINGX(m1), FINISHINGX(m2)); \
    for (int i=ispduptmp2; i<=ispduptmp3; i++) \
        for (int j=ispduptmp4; j<=ispduptmp5; j++)

/** For all elements in the array, accessed physically
 * @ingroup coreMatricesSizeShape
 *
 * This macro is used to generate loops for the matrix in an easy way using
 * physical indexes. It defines internal indexes 'i' and 'j' which ranges the
 * matrix using its physical definition.
 *
 * @code
 * FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(m)
 * {
 *     std::cout << DIRECT_MAT_ELEM(m, i, j) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(m) \
    for (int i=0; i<YSIZE(m); i++) \
        for (int j=0; j<XSIZE(m); j++)

/// @defgroup VectorsSizeShape 1D Matrices size and shape
/// @ingroup coreMultidimArraysSpeedUp

/** Vector element: Physical access
 * @ingroup VectorsSizeShape
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
 * @ingroup VectorsSizeShape
 */
#define dVi(v, i) DIRECT_VEC_ELEM(v, i)

/** Vector element: Logical access
 * @ingroup VectorsSizeShape
 *
 * @code
 * VEC_ELEM(v, -2) = 1;
 * val = VEC_ELEM(v, -2);
 * @endcode
 */
#define VEC_ELEM(v, i) DIRECT_VEC_ELEM(v, (i) - ((v).xinit))

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


// Forward declarations ====================================================
template<typename T> class coreMultidimArray;

template<typename T>
void coreArrayByScalar(const coreMultidimArray<T>& op1, const T& op2,
    coreMultidimArray<T>& result, char operation);

template<typename T>
void coreScalarByArray(const T& op1, const coreMultidimArray<T>& op2,
    coreMultidimArray<T>& result, char operation);

template<typename T>
void coreArrayByArray(const coreMultidimArray<T>& op1, const coreMultidimArray<T>& op2,
    coreMultidimArray<T>& result, char operation);


/** Template class for Xmipp arrays.
  * @ingroup coreMultidimensionalArrays
  * This class provides physical and logical access.
*/
template<typename T>
class coreMultidimArray
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
    unsigned long ndim;

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
    /// @defgroup coreMultidimArrayConstructors Constructors
    /// @ingroup coreMultidimensionalArrays
    /** Empty constructor.
     * @ingroup coreMultidimArrayConstructors
     * The empty constructor creates an array with no memory associated,
     * size=0.
     */
    coreMultidimArray()
    {
        coreInit();
    }
    
    /** Copy constructor
     * @ingroup coreMultidimArrayConstructors
     *
     * The created volume is a perfect copy of the input array but with a
     * different memory assignment.
     *
     * @code
     * Matrix3D< double > V2(V1);
     * @endcode
     */
    coreMultidimArray(const coreMultidimArray<T>& V)
    {
        coreInit();
        *this = V;
    }

    /** Destructor.
     * @ingroup coreMultidimArrayConstructors
     */
     ~coreMultidimArray()
     {
        coreDeallocate();
     }

    /** Clear.
     * @ingroup coreMultidimArrayConstructors
     */
     void clear()
     {
        coreDeallocate();
        coreInit();
     }

    /// @defgroup coreMultidimArrayCore Core memory operations for MultidimArrays
    /// @ingroup coreMultidimensionalArrays
    /** Core init.
     * @ingroup coreMultidimArrayCore
     * Initialize everything to 0
     */
    void coreInit()
    {
        xdim=yxdim=zyxdim=nzyxdim=0;
        ydim=zdim=ndim=1;
        zinit=yinit=xinit=0;
        data=NULL;
        destroyData=true;
    }

    /** Core allocate.
     * @ingroup coreMultidimArrayCore
     */
    void coreAllocate(unsigned long int _ndim, int _zdim, int _ydim, int _xdim)
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
     * @ingroup coreMultidimArrayCore
     * Free all data.
     */
    void coreDeallocate()
    {
        if (data != NULL && destroyData)
            delete[] data;
        data=NULL;
    }

    /// @defgroup coreMultidimSize Size
    /// @ingroup coreMultidimensionalArrays

    /** Copy the shape parameters
     * @ingroup coreMultidimSize
     *
     */
    void copyShape(const coreMultidimArray<T> &m)
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

    /** Resize to a given size
     * @ingroup coreMultidimSize
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
        for (unsigned long int l = 0; l < Ndim; l++)
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

    /** Resize a single 3D image
     * @ingroup coreMultidimSize
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

    /** Resize a single 2D image
     * @ingroup coreMultidimSize
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

    /** Resize a single 1D image
     * @ingroup coreMultidimSize
     *
     * This function assumes n and z and y are 1
     * @code
     * V1.resize(2);
     * @endcode
     */
    void resize(int Xdim)
    {
        resize(1, 1, 1, Xdim);
    }

    /** Resize according to a pattern.
     * @ingroup coreMultidimSize
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
    void resize(const coreMultidimArray<T1> &v)
    {
        if (NSIZE(*this) != NSIZE(v) || XSIZE(*this) != XSIZE(v) || 
            YSIZE(*this) != YSIZE(v) || ZSIZE(*this) != ZSIZE(v))
            resize(NSIZE(v), ZSIZE(v), YSIZE(v), XSIZE(v));

        STARTINGX(*this) = STARTINGX(v);
        STARTINGY(*this) = STARTINGY(v);
        STARTINGZ(*this) = STARTINGZ(v);
    }

    /** Returns the multidimArray dimension.
     * @ingroup coreMultidimSize
     *
     * Pay attention to the dimension order (N,Z,Y,X).
     *
     * @code
     * V.getDimension(Ndim, Zdim, Ydim, Xdim);
     * @endcode
     */
    void getDimension(unsigned long int &Ndim, int& Zdim, int& Ydim, int& Xdim) const
    {
        Xdim = XSIZE(*this);
        Ydim = YSIZE(*this);
        Zdim = ZSIZE(*this);
        Ndim = NSIZE(*this);
    }

    /** Get size.
     * @ingroup coreMultidimSize
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

    /** Same shape.
     * @ingroup coreMultidimSize
     *
     * Returns true if this object has got the same shape (origin and size)
     * than the argument
     */
    template <typename T1>
    bool sameShape(const coreMultidimArray<T1>& op) const
    {
        return (NSIZE(*this) == NSIZE(op) &&
                XSIZE(*this) == XSIZE(op) &&
                YSIZE(*this) == YSIZE(op) &&
                ZSIZE(*this) == ZSIZE(op) &&
                STARTINGX(*this) == STARTINGX(op) &&
                STARTINGY(*this) == STARTINGY(op) &&
                STARTINGZ(*this) == STARTINGZ(op));
    }


    /** 3D Logical to physical index translation.
     * @ingroup coreMultidimMemory
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

    /** 3D Physical to logical index translation.
     * @ingroup coreMultidimMemory
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

    /** 2D Logical to physical index translation
     * @ingroup coreMultidimMemory
     *
     * This function returns the physical position of a logical one.
     *
     * @code
     * m.toPhysical(i_log, j_log, i_phys, j_phys);
     * @endcode
     */
    void toPhysical(int i_log, int j_log, int& i_phys, int& j_phys) const
    {
        i_phys = i_log - STARTINGY(*this);
        j_phys = j_log - STARTINGX(*this);
    }

    /** 2D Physical to logical index translation
     * @ingroup coreMultidimMemory
     *
     * This function returns the logical position of a physical one.
     *
     * @code
     * m.toLogical(i_phys, j_phys, i_log, j_log);
     * @endcode
     */
    void toLogical(int i_phys, int j_phys, int &i_log, int& j_log) const
    {
        i_log = i_phys + STARTINGY(*this);
        j_log = j_phys + STARTINGX(*this);
    }

    /** 1D Logical to physical index translation
     * @ingroup coreMultidimMemory
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

    /** 1D Physical to logical index translation.
     * @ingroup coreMultidimMemory
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


    /** Set logical origin in Xmipp fashion.
     * @ingroup coreMultidimSize
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

    /** Returns the first valid logical Z index.
     * @ingroup coreMultidimSize
     */
    int startingZ() const
    {
        return zinit;
    }

    /** Returns the last valid logical Z index.
     * @ingroup coreMultidimSize
     */
    int finishingZ() const
    {
        return zinit + zdim - 1;
    }

    /** Returns the first valid logical Y index.
     * @ingroup coreMultidimSize
     */
    int startingY() const
    {
        return yinit;
    }

    /** Returns the last valid logical Y index.
     * @ingroup coreMultidimSize
     */
    int finishingY() const
    {
        return yinit + ydim - 1;
    }

    /** Returns the first valid logical X index.
     * @ingroup coreMultidimSize
     */
    int startingX() const
    {
        return xinit;
    }

    /** Returns the last valid logical X index.
     * @ingroup coreMultidimSize
     */
    int finishingX() const
    {
        return xinit + xdim - 1;
    }

    /** Returns Z dimension.
     * @ingroup coreMultidimSize
     */
    int sliceNumber() const
    {
        return zdim;
    }

    /** Returns Y dimension.
     * @ingroup coreMultidimSize
     */
    int rowNumber() const
    {
        return ydim;
    }

    /** Returns X dimension.
     * @ingroup coreMultidimSize
     */
    int colNumber() const
    {
        return xdim;
    }

    /// @defgroup Initialization Initialization
    /// @ingroup coreMultidimensionalArrays

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
    void initZeros(const coreMultidimArray<T1>& op)
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
    void initZeros(unsigned long int Ndim, int Zdim, int Ydim, int Xdim)
    {
        resize(Ndim, Zdim,Ydim,Xdim);
        initZeros();
    }

    /** Initialize to zeros with a given size.
     * @ingroup Initialization
     */
    void initZeros(int Xdim)
    {
        initZeros(1, 1, 1, Xdim);
    }

    /** Initialize to zeros with a given size.
     * @ingroup Initialization
     */
    void initZeros(int Ydim, int Xdim)
    {
        initZeros(1, 1, Ydim, Xdim);
    }

    /** Initialize to zeros with a given size.
     * @ingroup Initialization
     */
    void initZeros(int Zdim, int Ydim, int Xdim)
    {
        initZeros(1, Zdim, Ydim, Xdim);
    }


    /** Produce a 3D array suitable for working with Numerical Recipes.
     * @ingroup coreMultidimSize
     *
     * This function must be used only as a preparation for routines which need
     * that the first physical index is 1 and not 0 as it usually is in C. New
     * memory is needed to hold the new double pointer array.
     */
    T*** adaptForNumericalRecipes3D(unsigned long n = 0) const
    {
        T*** m = NULL;
        ask_Tvolume(m, 1, ZSIZE(*this), 1, YSIZE(*this), 1, XSIZE(*this));

        FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(*this)
            m[k+1][i+1][j+1] = DIRECT_NZYX_ELEM(*this, n, k, i, j);

        return m;
    }

    /** Kill a 3D array produced for numerical recipes.
     * @ingroup coreMultidimSize
     */
    void killAdaptationForNumericalRecipes3D(T*** m) const
    {
        free_Tvolume(m, 1, ZSIZE(*this), 1, YSIZE(*this), 1, XSIZE(*this));
    }

    /** Produce a 2D array suitable for working with Numerical Recipes
     * @ingroup coreMultidimSize
     *
     * This function must be used only as a preparation for routines which need
     * that the first physical index is 1 and not 0 as it usually is in C. New
     * memory is needed to hold the new double pointer array.
     */
    T** adaptForNumericalRecipes2D(unsigned long n = 0) const
    {
        T** m = NULL;
        ask_Tmatrix(m, 1, YSIZE(*this), 1, XSIZE(*this));

        FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(*this)
            m[i+1][j+1] = DIRECT_NZYX_ELEM(*this, n, 0, i, j);

        return m;
    }

    /** Produce a 1D pointer suitable for working with Numerical Recipes (2)
     * @ingroup coreMultidimSize
     *
     * This function meets the same goal as the one before, however this one
     * work with 2D arrays as a single pointer. The first element of the array
     * is pointed by result[1*Xdim+1], and in general result[i*Xdim+j]
     */
    T* adaptForNumericalRecipes22D() const
    {
        return MULTIDIM_ARRAY(*this) - 1 - XSIZE(*this);
    }

    /** Load 2D array from numerical recipes result.
     * @ingroup coreMultidimSize
     */
    void loadFromNumericalRecipes2D(T** m, int Ydim, int Xdim)
    {
        resize(Ydim, Xdim);

        for (int i = 1; i <= Ydim; i++)
            for (int j = 1; j <= Xdim; j++)
                (*this)(i - 1, j - 1) = m[i][j];
    }

    /** Kill a 2D array produced for numerical recipes
     * @ingroup coreMultidimSize
     *
     * The allocated memory is freed.
     */
    void killAdaptationForNumericalRecipes2D(T** m) const
    {
        free_Tmatrix(m, 1, YSIZE(*this), 1, XSIZE(*this));
    }

    /** Kill a 2D array produced for numerical recipes, 2.
     * @ingroup coreMultidimSize
     *
     * Nothing needs to be done.
     */
    void killAdaptationForNumericalRecipes22D(T** m) const
        {}

    /** Produce a 1D array suitable for working with Numerical Recipes
     * @ingroup coreMultidimSize
     *
     * This function must be used only as a preparation for routines which need
     * that the first physical index is 1 and not 0 as it usually is in C. In
     * fact the vector provided for Numerical recipes is exactly this same one
     * but with the indexes changed.
     *
     * This function is not ported to Python.
     */
    T* adaptForNumericalRecipes1D() const
    {
        return MULTIDIM_ARRAY(*this) - 1;
    }

    /** Kill a 1D array produced for Numerical Recipes.
     * @ingroup coreMultidimSize
     *
     * Nothing needs to be done in fact.
     *
     * This function is not ported to Python.
     */
    void killAdaptationForNumericalRecipes1D(T* m) const
        {}


    /** @defgroup coreMultidimUtilities Utilities
     *  @ingroup coreMultidimensionalArrays
     *
     */

    /** Reverse matrix values over X axis, keep in this object.
     * @ingroup coreMultidimUtilites
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
     * @ingroup coreMultidimUtilites
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
     * @ingroup coreMultidimUtilites
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

    /** Sqrt.
     * @ingroup coreMultidimUtilities
     *
     * Each component of the result is the square root of the original
     * component.
     */
    void selfSQRT()
    {
        T* ptr=NULL;
	unsigned long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
            *ptr = static_cast< T >(sqrt(static_cast< double >(*ptr)));
    }

    /** Log10.
     * @ingroup coreMultidimUtilities
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

    /** ROUND
     * @ingroup coreMultidimUtilities
     *
     * Applies a ROUND (look for the nearest integer) to each array element.
     */
    void selfROUND()
    {
        T* ptr=NULL;
	unsigned long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
            *ptr = ROUND(*ptr);
    }

    /** CEILING
     * @ingroup coreMultidimUtilities
     *
     * Applies a CEILING (look for the nearest larger integer) to each
     * array element.
     */
    void selfCEIL()
    {
        T* ptr=NULL;
	unsigned long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
            *ptr = CEIL(*ptr);
    }

    /** FLOOR 
     * @ingroup coreMultidimUtilities
     *
     * Applies a FLOOR (look for the nearest larger integer) to each
     * array element.
     */
    void selfFLOOR()
    {
        T* ptr=NULL;
	unsigned long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
            *ptr = FLOOR(*ptr);
    }

    /** ABS
     * @ingroup coreMultidimUtilities
     *
     * Applies an ABS (absolute value) to each array element.
     */
    void selfABS()
    {
        T* ptr=NULL;
	unsigned long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
            *ptr = ABS(*ptr);
    }

    /** MAX
     * @ingroup coreMultidimUtilities
     *
     * Each component of the result is the maximum of the correspoing
     * components of the two input arrays. They must have the same shape, if
     * not an exception is thrown
     */
    friend void MAX(const coreMultidimArray<T>& v1, const coreMultidimArray<T>& v2,
        coreMultidimArray<T>& result)
    {
        if (!v1.sameShape(v2))
            REPORT_ERROR(1007, "MAX: arrays of different shape");

        result.resize(v1);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(result)
            DIRECT_MULTIDIM_ELEM(result,n) = XMIPP_MAX(
                DIRECT_MULTIDIM_ELEM(v1,n),
                DIRECT_MULTIDIM_ELEM(v2,n));
    }

    /** MIN
     * @ingroup MultidimUtilities
     *
     * Each component of the result is the minimum of the correspoing
     * components of the two input arrays. They must have the same shape, if
     * not an exception is thrown
     */
    friend void MIN(const coreMultidimArray<T>& v1, const coreMultidimArray<T>& v2,
        coreMultidimArray<T>& result)
    {
        if (!v1.sameShape(v2))
            REPORT_ERROR(1007, "MIN: arrays of different shape");

        result.resize(v1);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(result)
            DIRECT_MULTIDIM_ELEM(result,n) = XMIPP_MIN(
                DIRECT_MULTIDIM_ELEM(v1,n),
                DIRECT_MULTIDIM_ELEM(v2,n));
    }

    /// @defgroup coreStatistics Statistics functions
    /// @ingroup coreMultidimensionalArrays

    /** Sum of array values.
     * @ingroup coreStatistics
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

    /** Sum of squared array values.
     * @ingroup coreStatistics
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

    /** Maximum of the values in the array.
     * @ingroup coreStatistics
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
     * @ingroup coreStatistics
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

    /** Print statistics in current line.
     * @ingroup coreStatistics
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

    /** 4D Indices for the minimum element.
     * @ingroup coreStatistics
     *
     * This function returns the index of the minimum element of an array.
     * array(l,k,i,j). Returns -1 if the array is empty
     */
    void minIndex(int &lmin, int& kmin, int& imin, int& jmin) const
    {
        if (XSIZE(*this) == 0)
        {
            lmin = kmin = imin = jmin = -1;
            return;
        }

        kmin = STARTINGZ(*this);
        imin = STARTINGY(*this);
        jmin = STARTINGX(*this);
        lmin = 0;
        T minval = NZYX_ELEM(*this, lmin, kmin, imin, jmin);


        FOR_ALL_NZYX_ELEMENTS_IN_MULTIDIMARRAY(*this)
            if (NZYX_ELEM(*this, l, k, i, j) > minval)
            {
                minval = NZYX_ELEM(*this, l, k, i, j);
                lmin = l;
                kmin = k;
                imin = i;
                jmin = j;
            }

    }

    /** 3D Indices for the minimum element.
     * @ingroup coreStatistics
     *
     * This function just calls to the 4D function
     */
    void minIndex(int& kmin, int& imin, int& jmin) const
    {
        minIndex(0,kmin,imin,jmin);
    }

    /** 2D Indices for the minimum element.
     * @ingroup coreStatistics
     *
     * This function just calls to the 4D function
     */
    void minIndex(int& imin, int& jmin) const
    {
        minIndex(0,0,imin,jmin);
    }

    /** 1D Indices for the minimum element.
     * @ingroup coreStatistics
     *
     * This function just calls to the 4D function
     */
    void minIndex(int& jmin) const
    {
        minIndex(0,0,0,jmin);
    }

    /** 4D Indices for the maximum element.
     * @ingroup coreStatistics
     *
     * This function returns the index of the maximum element of an array.
     * array(l,k,i,j). Returns -1 if the array is empty
     */
    void maxIndex(int &lmax, int& kmax, int& imax, int& jmax) const
    {
        if (XSIZE(*this) == 0)
        {
            lmax = kmax = imax = jmax = -1;
            return;
        }

        kmax = STARTINGZ(*this);
        imax = STARTINGY(*this);
        jmax = STARTINGX(*this);
        lmax = 0;
        T maxval = NZYX_ELEM(*this, lmax, kmax, imax, jmax);

        FOR_ALL_NZYX_ELEMENTS_IN_MULTIDIMARRAY(*this)
            if (NZYX_ELEM(*this, l, k, i, j) > maxval)
            {
                maxval = NZYX_ELEM(*this, l, k, i, j);
                lmax = l;
                kmax = k;
                imax = i;
                jmax = j;
            }
    }

    /** 3D Indices for the maximum element.
     * @ingroup coreStatistics
     *
     * This function just calls to the 4D function
     */
    void maxIndex(int& kmax, int& imax, int& jmax) const
    {
        maxIndex(0,kmax,imax,jmax);
    }

    /** 2D Indices for the maximum element.
     * @ingroup coreStatistics
     *
     * This function just calls to the 4D function
     */
    void maxIndex(int& imax, int& jmax) const
    {
        maxIndex(0,0,imax,jmax);
    }

    /** 1D Indices for the maximum element.
     * @ingroup coreStatistics
     *
     * This function just calls to the 4D function
     */
    void maxIndex(int& jmax) const
    {
        maxIndex(0,0,0,jmax);
    }

    /** Minimum and maximum of the values in the array.
     * @ingroup coreStatistics
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
     * @ingroup coreStatistics
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
     * @ingroup coreStatistics
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
     * @ingroup coreStatistics
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
     * @ingroup coreStatistics
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
        coreMultidimArray< double > temp;
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

     /** Volume element access via index.
     * @ingroup coreMultidimMemory
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
    T& operator()(const int k, const int i, const int j) const
    {
        return VOL_ELEM(*this, k, i, j);
    }

    /** Matrix element access via index
     * @ingroup coreMultidimMemory
     *
     * Returns the value of a matrix logical position. In our example we could
     * access from v(-2,-1) to v(2,1). The elements can be used either by value
     * or by reference. An exception is thrown if the index is outside the
     * logical range. The first argument is the Y position and the second the X
     * position.
     *
     * @code
     * m(-2, 1) = 1;
     * val = m(-2, 1);
     * @endcode
     */
    T& operator()(const int i, const int j) const
    {
        return MAT_ELEM(*this, i, j);
    }

    /** Vector element access
     * @ingroup coreMultidimMemory
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
    T& operator()(const int i) const
    {
       return VEC_ELEM(*this, i);
    }

    /** Get a single 1,2 or 3D image from a multi-image array
     *  @ingroup coreMultidimMemory
     * 
     * This function extracts a single-image array from a multi-image one.
     * @code
     * V.getImage(0, m);
     * @endcode
     */
    void getImage(unsigned long n, coreMultidimArray<T>& M) const
    {
        if (XSIZE(*this) == 0)
        {
            M.clear();
            return;
        }

        if (n > NSIZE(*this))
            REPORT_ERROR(1," Multidimarray getImage: n larger than NSIZE");

        M.resize(1, ZSIZE(*this), YSIZE(*this), XSIZE(*this));
        FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(M)
            DIRECT_MAT_ELEM(M, i, j) = DIRECT_NZYX_ELEM(*this, n, k, i, j);
        
        STARTINGX(M) = STARTINGX(*this);
        STARTINGY(M) = STARTINGY(*this);
        STARTINGZ(M) = STARTINGZ(*this);

    }

    /** 2D Slice access for reading.
     * @ingroup coreMultidimMemory
     *
     * This function returns a slice (a 2D matrix) corresponding to the choosen
     * slice inside the nth 3D matrix, the numbering of the slices is also logical not
     * physical. This function differs from the previous one in that this one
     * cuts and assign in a single step instead of in two steps, as in
     * the previous example.
     *
     * @code
     * V.slice(0, m);
     * @endcode
     */
    void getSlice(int k, coreMultidimArray<T>& M, char axis = 'Z', unsigned long n = 0) const
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
                             "Slice: Multidim subscript (k) out of range");

            k = k - STARTINGZ(*this);
            M.resize(1, 1, YSIZE(*this), XSIZE(*this));
            FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(M)
                DIRECT_MAT_ELEM(M, i, j) = DIRECT_NZYX_ELEM(*this, n, k, i, j);
            STARTINGX(M) = STARTINGX(*this);
            STARTINGY(M) = STARTINGY(*this);
            break;
        case 'Y':
            if (k < STARTINGY(*this) || k > FINISHINGY(*this))
                REPORT_ERROR(1203,
                             "Slice: Multidim subscript (i) out of range");

            k = k - STARTINGY(*this);
            M.resize(ZSIZE(*this), XSIZE(*this));
            FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(M)
                DIRECT_MAT_ELEM(M, i, j) = DIRECT_NZYX_ELEM(*this, n, i, k, j);
            STARTINGX(M) = STARTINGX(*this);
            STARTINGY(M) = STARTINGZ(*this);
            break;
        case 'X':
            if (k < STARTINGX(*this) || k > FINISHINGX(*this))
                REPORT_ERROR(1203,
                             "Slice: Multidim subscript (j) out of range");

            k = k - STARTINGX(*this);
            M.resize(ZSIZE(*this), YSIZE(*this));
            FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(M)
                DIRECT_MAT_ELEM(M, i, j) = DIRECT_NZYX_ELEM(*this, n, i, j, k);
            STARTINGX(M) = STARTINGY(*this);
            STARTINGY(M) = STARTINGZ(*this);
            break;
        default:
            REPORT_ERROR(1205,
                         (std::string) "Slice: not supported axis " + axis);
        }
    }

    /** Slice access for writing.
     * @ingroup coreMultidimMemory
     *
     * This function sets a 2D matrix corresponding to the choosen slice inside the nth
     * volume, the numbering of the slices is also logical not physical.
     *
     * @code
     * // Copies slice 0 in slice 1
     * V.setSlice(1, (V.slice(0)));
     * @endcode
     */
    void setSlice(int k, const coreMultidimArray<T>& v, unsigned long n = 0)
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
            DIRECT_NZYX_ELEM(*this, n, k, i, j) = DIRECT_MAT_ELEM(v, i, j);
    }

    /** Get row
     * @ingroup coreMultidimMemory
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
    void getRow(int i, coreMultidimArray<T>& v, unsigned long n = 0) const
    {
        if (XSIZE(*this) == 0 || YSIZE(*this) == 0)
        {
            v.clear();
            return;
        }

        if (i < STARTINGY(*this) || i > FINISHINGY(*this))
            REPORT_ERROR(1103, "getRow: Matrix subscript (i) greater than matrix dimension");

        v.resize(XSIZE(*this));
        STARTINGX(v) = STARTINGX(*this);

        for (int j = STARTINGX(*this); j <= FINISHINGX(*this); j++)
            VEC_ELEM(v, j) = NZYX_ELEM(*this, n, 0, i, j);

        v.setRow();
    }

    /** Return row. The same as previous.
     * @ingroup coreMultidimMemory
      */
    coreMultidimArray<T> Row(int i, unsigned long n = 0) const
    {
        coreMultidimArray<T> aux;
        getRow(i, aux, n);
        return aux;
    }

    /** Get Column
     * @ingroup coreMultidimMemory
     *
     * This function returns a column vector corresponding to the
     * choosen column inside the nth 2D matrix, the numbering of the
     * column is also logical not physical.
     *
     * @code
     * std::vector< double > v;
     * m.getCol(-1, v);
     * @endcode
     */
    void getCol(int j, coreMultidimArray<T>& v, unsigned long n = 0) const
    {
        if (XSIZE(*this) == 0 || YSIZE(*this) == 0)
        {
            v.clear();
            return;
        }

        if (j < STARTINGX(*this) || j > FINISHINGX(*this))
            REPORT_ERROR(1103,"getCol: Matrix subscript (j) greater than matrix dimension");

        v.resize(YSIZE(*this));
        STARTINGX(v)  = STARTINGY(*this);

        for (int i = STARTINGY(*this); i <= FINISHINGY(*this); i++)
            VEC_ELEM(v, i) = NZYX_ELEM(*this, n, 0, i, j);

        // FIXME: THINK WHETHER THIS IS NECESSARY... FOR NOW IT CONFLICTS WITH setCol below.
        //v.setCol();
    }

    /** Return Column.
     * @ingroup coreMultidimMemory
     *
     *  The same as the previous function.
     */
    coreMultidimArray<T> Col(int i, unsigned long n = 0) const
    {
        coreMultidimArray<T> aux;
        getCol(i, aux, n);
        return aux;
    }

    /** Set Row
     * @ingroup coreMultidimMemory
     *
     * This function sets a row vector corresponding to the choosen
     * row inside the nth 2D matrix, the numbering of the rows is also
     * logical not physical.
     *
     * @code
     * m.setRow(-2, m.row(1)); // Copies row 1 in row -2
     * @endcode
     */
    void setRow(int i, const coreMultidimArray<T>& v, unsigned long n = 0)
    {
        if (XSIZE(*this) == 0 || YSIZE(*this) == 0)
            REPORT_ERROR(1, "setRow: Target matrix is empty");

        if (i < STARTINGY(*this) || i > FINISHINGY(*this))
            REPORT_ERROR(1103, "setRow: Matrix subscript (i) out of range");

        if (XSIZE(v) != XSIZE(*this))
            REPORT_ERROR(1102,
                         "setRow: Vector dimension different from matrix one");

        //FIXME....
        //if (!v.isRow())
        //    REPORT_ERROR(1107, "setRow: Not a row vector in assignment");

        i = i - STARTINGY(*this);
        for (int j = 0; j < XSIZE(*this); j++)
            DIRECT_NZYX_ELEM(*this, n, 0, i, j) = DIRECT_VEC_ELEM(v, j);
    }

    /** Set Column
     * @ingroup coreMultidimMemory
     *
     * This function sets a column vector corresponding to the choosen column
     * inside matrix, the numbering of the column is also logical not physical.
     *
     * @code
     * m.setCol(-1, (m.row(1)).transpose()); // Copies row 1 in column -1
     * @endcode
     */
    void setCol(int j, const coreMultidimArray<T>& v, unsigned long n = 0)
    {
        if (XSIZE(*this) == 0 || YSIZE(*this) == 0)
            REPORT_ERROR(1, "setCol: Target matrix is empty");

        if (j < STARTINGX(*this) || j > FINISHINGX(*this))
            REPORT_ERROR(1103, "setCol: Matrix subscript (j) out of range");

        if (XSIZE(v) != YSIZE(*this))
            REPORT_ERROR(1102,
                         "setCol: Vector dimension different from matrix one");
        
        // FIXME...
        //if (!v.isCol())
        //    REPORT_ERROR(1107, "setCol: Not a column vector in assignment");

        j = j - STARTINGX(*this);
        for (int i = 0; i < YSIZE(*this); i++)
            DIRECT_NZYX_ELEM(*this, n, 0, i, j) = DIRECT_VEC_ELEM(v, i);
    }



    /// @defgroup coreArithmethic Arithmethic operations
    /// @ingroup coreMultidimensionalArrays

    /// @defgroup Operators Operators
    /// @ingroup coreMultidimensionalArrays

    /** Equality.
     * @ingroup Operators
     *
     * Returns true if this object has got the same shape (origin and size)
     * than the argument and the same values (within accuracy).
     */
    bool equal(const coreMultidimArray<T>& op,
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

    /** Assignment.
     * @ingroup coreArithmethic
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
    coreMultidimArray<T>& operator=(const coreMultidimArray<T>& op1)
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
     * @ingroup coreArithmethic
     *
     * It is used to build arithmetic expressions. You can make a minus
     * of anything as long as it is correct semantically.
     *
     * @code
     * v1 = -v2;
     * v1 = -v2.transpose();
     * @endcode
     */
    coreMultidimArray<T> operator-() const
    {
        coreMultidimArray<T> tmp(*this);
        T* ptr;
	unsigned long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(tmp,n,ptr)
            *ptr = -(*ptr);
        return tmp;
    }

    /** @defgroup coreArrayByScalar Array "by" scalar operations
     * @ingroup coreArithmethic
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
     * @ingroup coreArrayByScalar
     *
     * This function must take one vector and a constant, and operate element
     * by element according to the operation required. This is the function
     * which really implements the operations. Simple calls to it perform much
     * faster than calls to the corresponding operators. Although it is
     * supposed to be a hidden function not useable by normal programmers.
     *
     * This function is not ported to Python.
     */
    friend void arrayByScalar(const coreMultidimArray<T>& op1,
                                T op2,
                                coreMultidimArray<T>& result,
                                char operation)
    {
        result.resize(op1);
        coreArrayByScalar(op1, op2, result, operation);
    }

    /** Core array by scalar operation.
     * @ingroup coreArrayByScalar
     *
     * It assumes that the result is already resized.
     *
     * This function is not ported to Python.
     */
    friend void coreArrayByScalar<>(const coreMultidimArray<T>& op1,
                                    const T& op2,
                                    coreMultidimArray<T>& result,
                                    char operation);

    /** v3 = v1 + k.
     * @ingroup coreArrayByScalar
     */
    coreMultidimArray<T> operator+(T op1) const
    {
        coreMultidimArray<T> tmp;
        arrayByScalar(*this, op1, tmp, '+');
        return tmp;
    }

    /** v3 = v1 - k.
     * @ingroup coreArrayByScalar
     */
    coreMultidimArray<T> operator-(T op1) const
    {
        coreMultidimArray<T> tmp;
        arrayByScalar(*this, op1, tmp, '-');
        return tmp;
    }

    /** v3 = v1 * k.
     * @ingroup coreArrayByScalar
     */
    coreMultidimArray<T> operator*(T op1) const
    {
        coreMultidimArray<T> tmp;
        arrayByScalar(*this, op1, tmp, '*');
        return tmp;
    }

    /** v3 = v1 / k.
     * @ingroup coreArrayByScalar
     */
    coreMultidimArray<T> operator/(T op1) const
    {
        coreMultidimArray<T> tmp;
        arrayByScalar(*this, op1, tmp, '/');
        return tmp;
    }

    /** v3 += k.
     * @ingroup coreArrayByScalar
     *
     * This function is not ported to Python.
     */
    void operator+=(const T& op1)
    {
        arrayByScalar(*this, op1, *this, '+');
    }

    /** v3 -= k.
     * @ingroup coreArrayByScalar
     *
     * This function is not ported to Python.
     */
    void operator-=(const T& op1)
    {
        arrayByScalar(*this, op1, *this, '-');
    }

    /** v3 *= k.
     * @ingroup coreArrayByScalar
     *
     * This function is not ported to Python.
     */
    void operator*=(const T& op1)
    {
        arrayByScalar(*this, op1, *this, '*');
    }

    /** v3 /= k.
     * @ingroup coreArrayByScalar
     *
     * This function is not ported to Python.
     */
    void operator/=(const T& op1)
    {
        arrayByScalar(*this, op1, *this, '/');
    }

    /** @defgroup coreScalarByArray Scalar "by" array operations
     * @ingroup coreArithmethic
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
     * @ingroup coreScalarByArray
     *
     * This function must take one scalar and a vector, and operate element by
     * element according to the operation required. This is the function which
     * really implements the operations. Simple calls to it perform much faster
     * than calls to the corresponding operators. Although it is supposed to
     * be a hidden function not useable by normal programmers.
     *
     * This function is not ported to Python.
     */

    /// FIXME: DO I REALLY NEED FRIENDS HERE AND BELOW??

    friend void scalarByArray(T op1,
                                const coreMultidimArray<T>& op2,
                                coreMultidimArray<T>& result,
                                char operation)
    {
        result.resize(op2);
        coreScalarByArray(op1, op2, result, operation);
    }

    /** Core array by scalar operation.
     * @ingroup coreScalarByArray
     *
     * It assumes that the result is already resized.
     *
     * This function is not ported to Python.
     */
    friend void coreScalarByArray<>(const T& op1,
                                    const coreMultidimArray<T>& op2,
                                    coreMultidimArray<T>& result,
                                    char operation);

    /** v3 = k + v2.
     * @ingroup coreScalarByArray
     */
    friend coreMultidimArray<T> operator+(T op1, const coreMultidimArray<T>& op2)
    {
        coreMultidimArray<T> tmp;
        scalarByArray(op1, op2, tmp, '+');
        return tmp;
    }

    /** v3 = k - v2.
     * @ingroup coreScalarByArray
     */
    friend coreMultidimArray<T> operator-(T op1, const coreMultidimArray<T>& op2)
    {
        coreMultidimArray<T> tmp;
        scalarByArray(op1, op2, tmp, '-');
        return tmp;
    }

    /** v3 = k * v2.
     * @ingroup coreScalarByArray
     */
    friend coreMultidimArray<T> operator*(T op1, const coreMultidimArray<T>& op2)
    {
        coreMultidimArray<T> tmp;
        scalarByArray(op1, op2, tmp, '*');
        return tmp;
    }

    /** v3 = k / v2
     * @ingroup coreScalarByArray
     */
    friend coreMultidimArray<T> operator/(T op1, const coreMultidimArray<T>& op2)
    {
        coreMultidimArray<T> tmp;
        scalarByArray(op1, op2, tmp, '/');
        return tmp;
    }

    /** @defgroup coreArrayByArray  Array "by" array arithmethic operations
     * @ingroup coreArithmethic
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
    friend void coreArrayByArray<>(const coreMultidimArray<T>& op1,
        const coreMultidimArray<T>& op2, coreMultidimArray<T>& result,
        char operation);
    
    /** Array by array
     * @ingroup coreArrayByArray
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
    friend void arrayByArray(const coreMultidimArray<T>& op1,
        const coreMultidimArray<T>& op2, coreMultidimArray<T>& result,
        char operation)
    {
        if (!op1.sameShape(op2))
            REPORT_ERROR(1007,
                         (std::string) "Array_by_array: different shapes (" +
                         operation + ")");

        result.resize(op1);
        coreArrayByArray(op1, op2, result, operation);
    }

    /** v3 = v1 + v2.
     * @ingroup coreArrayByArray
     */
    coreMultidimArray<T> operator+(const coreMultidimArray<T>& op1) const
    {
        coreMultidimArray<T> tmp;
        arrayByArray(*this, op1, tmp, '+');
        return tmp;
    }

    /** v3 = v1 - v2.
     * @ingroup coreArrayByArray
     */
    coreMultidimArray<T> operator-(const coreMultidimArray<T>& op1) const
    {
        coreMultidimArray<T> tmp;
        arrayByArray(*this, op1, tmp, '-');
        return tmp;
    }

    /** v3 = v1 * v2.
     * @ingroup coreArrayByArray
     */
    coreMultidimArray<T> operator*(const coreMultidimArray<T>& op1) const
    {
        coreMultidimArray<T> tmp;
        arrayByArray(*this, op1, tmp, '*');
        return tmp;
    }

    /** v3 = v1 / v2.
     * @ingroup coreArrayByArray
     */
    coreMultidimArray<T> operator/(const coreMultidimArray<T>& op1) const
    {
        coreMultidimArray<T> tmp;
        arrayByArray(*this, op1, tmp, '/');
        return tmp;
    }

    /** v3 += v2.
     * @ingroup coreArrayByArray
     */
    void operator+=(const coreMultidimArray<T>& op1)
    {
        arrayByArray(*this, op1, *this, '+');
    }

    /** v3 -= v2.
     * @ingroup coreArrayByArray
     */
    void operator-=(const coreMultidimArray<T>& op1)
    {
        arrayByArray(*this, op1, *this, '-');
    }

    /** v3 *= v2.
     * @ingroup coreArrayByArray
     */
    void operator*=(const coreMultidimArray<T>& op1)
    {
        arrayByArray(*this, op1, *this, '*');
    }

    /** v3 /= v2.
     * @ingroup coreArrayByArray
     */
    void operator/=(const coreMultidimArray<T>& op1)
    {
        arrayByArray(*this, op1, *this, '/');
    }

    /** @defgroup rwCoreMultidimArray "by" array operations.
     * @ingroup coreMultidimensionalArrays

    /** Read from an ASCII file.
     * @ingroup rwCoreMultidimArray
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
     * @ingroup rwCoreMultidimArray
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
     * @ingroup rwCoreMultidimArray
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
     * @ingroup rwCoreMultidimArray
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
     * @ingroup rwCoreMultidimArray
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

/// @defgroup coreMultidimFunctions Functions for all multidimensional arrays
/// @ingroup coreMultidimensionalArrays

/** Conversion from one type to another.
 * @ingroup coreMultidimFunctions
 *
 * If we have an integer array and we need a double one, we can use this
 * function. The conversion is done through a type casting of each element
 * If n >= 0, only the nth volumes will be converted, otherwise all NSIZE volumes
 */
template<typename T1, typename T2>
    void typeCast(const coreMultidimArray<T1>& v1, coreMultidimArray<T2>& v2, long n = -1)
{
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

template <typename T>
void coreArrayByScalar(const coreMultidimArray<T>& op1,
                            const T& op2,
                            coreMultidimArray<T>& result,
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
                         const coreMultidimArray<T>& op2,
                         coreMultidimArray<T>& result,
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

template <typename T>
void coreArrayByArray(const coreMultidimArray<T>& op1,
    const coreMultidimArray<T>& op2, coreMultidimArray<T>& result,
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

/** MultidimArray equality.
 * @ingroup coreMultidimMisc */
template<typename T>
bool operator==(const coreMultidimArray<T>& op1, const coreMultidimArray<T>& op2)
{
    return op1.equal(op2);
}

/** MultidimArray inequality.
 * @ingroup coreMultidimMisc */
template<typename T>
bool operator!=(const coreMultidimArray<T>& op1, const coreMultidimArray<T>& op2)
{
    return !(op1==op2);
}


template<typename T>
std::ostream& operator<<(std::ostream& ostrm, const coreMultidimArray<T>& v)
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

// Specializations case for complex numbers
template<>
std::ostream& operator<<(std::ostream& ostrm, const coreMultidimArray< std::complex<double> >& v);


#endif
