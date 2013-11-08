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

#ifdef XMIPP_MMAP
#include <sys/mman.h>
#endif
#include "../../external/bilib/types/tsplinebasis.h"
#include "../../external/bilib/headers/kernel.h"
#include "xmipp_strings.h"
#include "matrix1d.h"
#include "matrix2d.h"

extern int bestPrecision(float F, int _width);
extern String floatToString(float F, int _width, int _prec);

/// @defgroup MultidimensionalArrays Multidimensional Arrays
/// @ingroup DataLibrary
//@{
/** @name MultidimArraysSpeedUp Speed up macros
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
/** Returns the first X valid logical index
 */
#define STARTINGX(v) ((v).xinit)

/** Returns the last X valid logical index
 */
#define FINISHINGX(v) ((v).xinit + (int)(v).xdim - 1)

/** Returns the first Y valid logical index
 */
#define STARTINGY(v) ((v).yinit)

/** Returns the last Y valid logical index
 */
#define FINISHINGY(v) ((v).yinit + (int)(v).ydim - 1)

/** Returns the first Z valid logical index
 */
#define STARTINGZ(v) ((v).zinit)

/** Returns the last Z valid logical index
 */
#define FINISHINGZ(v) ((v).zinit + (int)(v).zdim - 1)

/** Check if x is inside logical bounds
 */
#define INSIDEX(v, x) ((x) >= STARTINGX(v) && (x) <= FINISHINGX(v))
/** Check if y is inside logical bounds
 */
#define INSIDEY(v, y) ((y) >= STARTINGY(v) && (y) <= FINISHINGY(v))
/** Check if z is inside logical bounds
 */
#define INSIDEZ(v, z) ((z) >= STARTINGZ(v) && (z) <= FINISHINGZ(v))
/** Check if a position x, y is inside the logical index bounds
 */
#define INSIDEXY(v, x, y) (INSIDEX(v, x) && INSIDEY(v, y))
/** Check if a position x, y is inside the logical index bounds
 */
#define INSIDEXYZ(v, x, y, z) (INSIDEX(v, x) && INSIDEY(v, y) && INSIDEZ(v,z))

/** Access to X dimension (size)
 */
#define XSIZE(v) ((v).xdim)

/** Access to Y dimension (size)
 */
#define YSIZE(v) ((v).ydim)

/** Access to Z dimension (size)
 */
#define ZSIZE(v) ((v).zdim)

/** Access to N dimension (size)
 */
#define NSIZE(v) ((v).ndim)

/** Access to XY dimension (Ysize*Xsize)
 */
#define YXSIZE(v) ((v).yxdim)

/** Access to XYZ dimension (Zsize*Ysize*Xsize)
 */
#define ZYXSIZE(v) ((v).zyxdim)

/** Access to XYZN dimension (Nsize*Zsize*Ysize*Xsize)
 */
#define MULTIDIM_SIZE(v) ((v).nzyxdim)

/** Access to XYZN dimension (Nsize*Zsize*Ysize*Xsize)
 */
#define NZYXSIZE(v) ((v).nzyxdim)

/** Array access.
 *
 * This macro gives you access to the array (T **)
 */
#ifndef MULTIDIM_ARRAY
#define MULTIDIM_ARRAY(v) ((v).data)
#endif

/** Access to a direct element.
 * v is the array, l is the image, k is the slice, i is the Y index and j is the X index.
 * i and j) within the slice.
 */
#define DIRECT_NZYX_ELEM(v, l, k, i, j) ((v).data[(l)*ZYXSIZE(v)+(k)*YXSIZE(v)+((i)*XSIZE(v))+(j)])

/** Access to a direct element.
 * v is the array, l is the image, k is the slice, i is the Y index and j is the X index.
 * i and j) within the slice.
 */
#define DIRECT_ZYX_ELEM(v, k, i, j) ((v).data[(k)*YXSIZE(v)+((i)*XSIZE(v))+(j)])

/** Access to a direct element.
 * v is the array, l is the image, k =0, i is the Y index and j is the X index.
 * i and j) within the slice.
 */
#define DIRECT_N_YX_ELEM(v, l, k, i, j) ((v).data[(l)*ZYXSIZE(v)             +((i)*XSIZE(v))+(j)])
/** Access to a direct element.
 * v is the array, l is the image, k =0, i = 0 and j is the X index.
 * i and j) within the slice.
 */
#define DIRECT_N__X_ELEM(v, l, k, i, j) ((v).data[(l)*ZYXSIZE(v)+(j)])

/** Multidim element: Logical access.
 */
#define NZYX_ELEM(v, l, k, i, j)  \
    DIRECT_NZYX_ELEM((v), (l), (k) - STARTINGZ(v), (i) - STARTINGY(v), (j) - STARTINGX(v))

/** Access to a direct element.
 * v is the array, k is the slice and n is the number of the pixel (combined i and j)
 * within the slice.
 */
#define DIRECT_MULTIDIM_ELEM(v,n) ((v).data[(n)])

/** For all direct elements in the array
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
    for (size_t n=0; n<NZYXSIZE(v); ++n)

/** For all direct elements in the array
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
    for (size_t l=0; l<NSIZE(V); ++l) \
        for (size_t k=0; k<ZSIZE(V); ++k) \
            for (size_t i=0; i<YSIZE(V); ++i)      \
                for (size_t j=0; j<XSIZE(V); ++j)

/** For all direct elements in the array
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
    for (size_t l=0; l<NSIZE(V); ++l) \
        for (int k=STARTINGZ(V); k<=FINISHINGZ(V); ++k) \
            for (int i=STARTINGY(V); i<=FINISHINGY(V); ++i)     \
                for (int j=STARTINGX(V); j<=FINISHINGX(V); ++j)

/** For all direct elements in the array, pointer version
 *
 * This macro is used to generate loops for the array in an easy manner. It
 * defines an internal index 'k' which goes over the slices and 'n' that
 * goes over the pixels in each slice. Each element can be accessed through
 * an external pointer called ptr.
 *
 * @code
 * T* ptr=NULL;
 * size_t n;
 * FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(v,n,ptr)
 * {
 *     std::cout << *ptr << " ";
 * }
 * @endcode
 */
#define FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(v,n,ptr) \
    for ((n)=0, (ptr)=(v).data; (n)<NZYXSIZE(v); ++(n), ++(ptr))

/** Access to a direct element.
 * v is the array, k is the slice (Z), i is the Y index and j is the X index.
 */
#define DIRECT_A3D_ELEM(v,k,i,j) ((v).data[YXSIZE(v)*(k)+((i)*XSIZE(v))+(j)])

/** A short alias for the previous function.
 *
 */
#define dAkij(V, k, i, j) DIRECT_A3D_ELEM(V, k, i, j)

/** Volume element: Logical access.
 *
 * @code
 * A3D_ELEM(V, -1, -2, 1) = 1;
 * val = A3D_ELEM(V, -1, -2, 1);
 * @endcode
 */
#define A3D_ELEM(V, k, i, j) \
    DIRECT_A3D_ELEM((V),(k) - STARTINGZ(V), (i) - STARTINGY(V), (j) - STARTINGX(V))

/** For all elements in the array.
 *
 * This macro is used to generate loops for the volume in an easy way. It
 * defines internal indexes 'k','i' and 'j' which ranges the volume using its
 * mathematical definition (ie, logical access).
 *
 * @code
 * FOR_ALL_ELEMENTS_IN_ARRAY3D(V)
 * {
 *     std::cout << V(k, i, j) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_ELEMENTS_IN_ARRAY3D(V) \
    for (int k=STARTINGZ(V); k<=FINISHINGZ(V); ++k) \
        for (int i=STARTINGY(V); i<=FINISHINGY(V); ++i) \
            for (int j=STARTINGX(V); j<=FINISHINGX(V); ++j)

/** For all elements in common.
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
 * MultidimArray< double > V1(10, 10, 10), V2(20, 20, 20);
 * V1.setXmippOrigin();
 * V2.setXmippOrigin();
 *
 * FOR_ALL_ELEMENTS_IN_COMMON_IN_ARRAY3D(V1, V2)
 * {
 *    // ...
 * }
 * @endcode
 */
#define FOR_ALL_ELEMENTS_IN_COMMON_IN_ARRAY3D(V1, V2) \
    ispduptmp0 = XMIPP_MAX(STARTINGZ(V1), STARTINGZ(V2)); \
    ispduptmp1 = XMIPP_MIN(FINISHINGZ(V1),FINISHINGZ(V2)); \
    ispduptmp2 = XMIPP_MAX(STARTINGY(V1), STARTINGY(V2)); \
    ispduptmp3 = XMIPP_MIN(FINISHINGY(V1),FINISHINGY(V2)); \
    ispduptmp4 = XMIPP_MAX(STARTINGX(V1), STARTINGX(V2)); \
    ispduptmp5 = XMIPP_MIN(FINISHINGX(V1),FINISHINGX(V2)); \
    for (int k=ispduptmp0; k<=ispduptmp1; ++k) \
        for (int i=ispduptmp2; i<=ispduptmp3; ++i) \
            for (int j=ispduptmp4; j<=ispduptmp5; ++j)

/** For all direct elements in the array.
 *
 * This macro is used to generate loops for the volume in an easy way. It
 * defines internal indexes 'k','i' and 'j' which ranges the volume using its
 * physical definition.
 *
 * @code
 * FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(V)
 * {
 *     std::cout << DIRECT_A3D_ELEM(m, k, i, j) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(V) \
    for (size_t k=0; k<ZSIZE(V); ++k) \
        for (size_t i=0; i<YSIZE(V); ++i) \
            for (size_t j=0; j<XSIZE(V); ++j)

/** Access to a direct element of a matrix.
 * v is the array, i and j define the element v_ij.
 *
 * Be careful because this is physical access, usually matrices follow the C
 * convention of starting index==0 (X and Y). This function should not be used
 * as it goes against the vector library philosophy unless you explicitly want
 * to access directly to any value in the matrix without taking into account its
 * logical position
 *
 * @code
 * DIRECT_A2D_ELEM(m, 0, 0) = 1;
 * val = DIRECT_A2D_ELEM(m, 0, 0);
 * @endcode
 */
#define DIRECT_A2D_ELEM(v,i,j) ((v).data[(i)*(v).xdim+(j)])

/** Short alias for DIRECT_A2D_ELEM
 */
#define dAij(M, i, j) DIRECT_A2D_ELEM(M, i, j)

/** Matrix element: Logical access
 *
 * @code
 * A2D_ELEM(m, -2, 1) = 1;
 * val = A2D_ELEM(m, -2, 1);
 * @endcode
 */
#define A2D_ELEM(v, i, j) \
    DIRECT_A2D_ELEM(v, (i) - STARTINGY(v), (j) - STARTINGX(v))

/** TRUE if both arrays have the same shape
 *
 * Two arrays have the same shape if they have the same size and the same
 * starting point. Be aware that this is a macro which simplifies to a boolean.
 */
#define SAME_SHAPE2D(v1, v2) \
    (XSIZE(v1) == XSIZE(v2) && \
     YSIZE(v1) == YSIZE(v2) && \
     STARTINGX(v1) == STARTINGX(v2) && \
     STARTINGY(v1) == STARTINGY(v2))

#define SAME_SHAPE3D(v1, v2) \
    (XSIZE(v1) == XSIZE(v2) && \
     YSIZE(v1) == YSIZE(v2) && \
     ZSIZE(v1) == ZSIZE(v2) && \
     STARTINGX(v1) == STARTINGX(v2) && \
     STARTINGY(v1) == STARTINGY(v2) && \
     STARTINGZ(v1) == STARTINGZ(v2))


/** For all elements in the array
 *
 * This macro is used to generate loops for the matrix in an easy way. It
 * defines internal indexes 'i' and 'j' which ranges the matrix using its
 * mathematical definition (ie, logical access).
 *
 * @code
 * FOR_ALL_ELEMENTS_IN_ARRAY2D(m)
 * {
 *     std::cout << m(i, j) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_ELEMENTS_IN_ARRAY2D(m) \
    for (int i=STARTINGY(m); i<=FINISHINGY(m); ++i) \
        for (int j=STARTINGX(m); j<=FINISHINGX(m); ++j)

/** For all elements in common
 *
 * This macro is used to generate loops for all the elements logically in common
 * between two images in an easy manner. Then i and j (locally defined) range
 * from MAX(STARTINGY(V1), STARTINGY(V2)) to MIN(FINISHINGY(V1),
 * FINISHINGY(V2)), MAX(STARTINGX(V1), STARTINGX(V2)) to MIN(FINISHINGX(V1),
 * FINISHINGX(V2)) (included limits) respectively. You need to define
 * SPEED_UP_temps.
 *
 * @code
 * MultidimArray< double > m1(10, 10), m2(20, 20);
 * m1.setXmippOrigin();
 * m2.setXmippOrigin();
 *
 * FOR_ALL_ELEMENTS_IN_COMMON_IN_ARRAY2D(m1, m2)
 * {
 *     ...
 * }
 * @endcode
 */
#define FOR_ALL_ELEMENTS_IN_COMMON_IN_ARRAY2D(m1, m2) \
    ispduptmp2 = XMIPP_MAX(STARTINGY(m1), STARTINGY(m2)); \
    ispduptmp3 = XMIPP_MIN(FINISHINGY(m1), FINISHINGY(m2)); \
    ispduptmp4 = XMIPP_MAX(STARTINGX(m1), STARTINGX(m2)); \
    ispduptmp5 = XMIPP_MIN(FINISHINGX(m1), FINISHINGX(m2)); \
    for (int i=ispduptmp2; i<=ispduptmp3; ++i) \
        for (int j=ispduptmp4; j<=ispduptmp5; ++j)

/** For all elements in the array between corners.
 *
 * This macro is used to generate loops for a volume in an easy manner. Then
 *  ZZ(r), YY(r) and XX(r) range from
 *
 * (int) ZZ(corner1) to (int)ZZ(corner2),
 * (int) YY(corner1) to (int)YY(corner2),
 * (int) XX(corner1) to (int) XX(corner2) (included limits) respectively.
 *
 * Notice that corner1 and corner2 need only be MultidimArray.
 *
 * @code
 * MultidimArray< double > corner1(3), corner2(3), r(3);
 * XX(corner1) = -1; XX(corner2) = 1;
 * YY(corner1) = -2; YY(corner2) = 2;
 * ZZ(corner1) = -3; ZZ(corner2) = 3;
 *
 * FOR_ALL_ELEMENTS_IN_ARRAY3D_BETWEEN(corner1, corner2)
 * {
 *     std::cout << v(r) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_ELEMENTS_IN_ARRAY3D_BETWEEN(corner1, corner2) \
    for (ZZ(r)=ZZ((corner1)); ZZ(r)<=ZZ((corner2)); ++ZZ(r)) \
        for (YY(r)=YY((corner1)); YY(r)<=YY((corner2)); ++YY(r)) \
            for (XX(r)=XX((corner1)); XX(r)<=XX((corner2)); ++XX(r))

/** For all elements in the array between corners
 *
 * This macro is used to generate loops for a matrix in an easy manner. It needs
 * an externally defined MultidimArray< double > r(2). Then YY(r) and XX(r) range
 * from (int) YY(corner1) to (int)YY(corner2), (int) XX(corner1) to (int)
 * XX(corner2) (included limits) respectively. Notice that corner1 and corner2
 * need only be MultidimArray.
 *
 * @code
 * MultidimArray< double > corner1(2), corner2(2);
 * MultidimArray< int > r(2);
 * XX(corner1) = -1;
 * XX(corner2) = 1;
 * YY(corner1) = -2;
 * YY(corner2) = 2;
 *
 * FOR_ALL_ELEMENTS_IN_ARRAY2D_BETWEEN(corner1, corner2)
 * {
 *     std::cout << v(r) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_ELEMENTS_IN_ARRAY2D_BETWEEN(corner1, corner2) \
    for (YY(r)=YY((corner1)); YY(r)<=YY((corner2)); ++YY(r)) \
        for (XX(r)=XX((corner1)); XX(r)<=XX((corner2)); ++XX(r))

/** For all elements in the array between corners
 *
 * This macro is used to generate loops for a vector in an easy manner. It needs
 * an externally defined MultidimArray< double > r(1). Then XX(r) ranges from
 * (int) XX(corner1) to (int) XX(corner2) (included limits) (notice that corner1
 * and corner2 need only to be MultidimArray).
 *
 * @code
 * MultidimArray< double > corner1(1), corner2(1), r(1);
 * XX(corner1) = -1;
 * XX(corner2) = 1;
 * FOR_ALL_ELEMENTS_IN_ARRAY1D_BETWEEN(corner1, corner2)
 * {
 *     std::cout << v(XX(r)) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_ELEMENTS_IN_ARRAY1D_BETWEEN(corner1, corner2) \
    for (XX(r)=(int) XX((corner1)); XX(r)<=(int) XX((corner2)); ++XX(r))

/** For all elements in the array, accessed physically
 *
 * This macro is used to generate loops for the matrix in an easy way using
 * physical indexes. It defines internal indexes 'i' and 'j' which ranges the
 * matrix using its physical definition.
 *
 * @code
 * FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(m)
 * {
 *     std::cout << DIRECT_A2D_ELEM(m, i, j) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(m) \
    for (size_t i=0; i<YSIZE(m); ++i) \
        for (size_t j=0; j<XSIZE(m); ++j)

/** Vector element: Physical access
 *
 * Be careful because this is physical access, usually vectors follow the C
 * convention of starting index==0. This function should not be used as it goes
 * against the vector library philosophy unless you explicitly want to access
 * directly to any value in the vector without taking into account its logical
 * position.
 *
 * @code
 * DIRECT_A1D_ELEM(v, 0) = 1;
 * val = DIRECT_A1D_ELEM(v, 0);
 * @endcode
 */
#define DIRECT_A1D_ELEM(v, i) ((v).data[(i)])

/** A short alias to previous function
 */
#define dAi(v, i) DIRECT_A1D_ELEM(v, i)

/** Vector element: Logical access
 *
 * @code
 * A1D_ELEM(v, -2) = 1;
 * val = A1D_ELEM(v, -2);
 * @endcode
 */
#define A1D_ELEM(v, i) DIRECT_A1D_ELEM(v, (i) - ((v).xinit))

/** For all elements in the array
 *
 * This macro is used to generate loops for the vector in an easy manner. It
 * defines an internal index 'i' which ranges the vector using its mathematical
 * definition (ie, logical access).
 *
 * @code
 * FOR_ALL_ELEMENTS_IN_ARRAY1D(v)
 * {
 *     std::cout << v(i) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_ELEMENTS_IN_ARRAY1D(v) \
    for (int i=STARTINGX(v); i<=FINISHINGX(v); ++i)

/** For all elements in common
 *
 * This macro is used to generate loops for all the elements logically in common
 * between two vectors in an easy manner. Then i (locally defined) ranges from
 * MAX(STARTINGX(V1), STARTINGX(V2)) to MIN(FINISHINGX(V1), FINISHINGX(V2))
 * (included limits) respectively. You need to define SPEED_UP_temps.
 *
 * @code
 * MultidimArray< double > v1(10), v2(20);
 * v1.setXmippOrigin();
 * v2.setXmippOrigin();
 *
 * FOR_ALL_ELEMENTS_IN_COMMON_IN_ARRAY1D(v1, v2)
 * {
 *     ...
 * }
 * @endcode
 */
#define FOR_ALL_ELEMENTS_IN_COMMON_IN_ARRAY1D(v1, v2) \
    ispduptmp4 = XMIPP_MAX(STARTINGX(v1), STARTINGX(v2)); \
    ispduptmp5 = XMIPP_MIN(FINISHINGX(v1), FINISHINGX(v2)); \
    for (int i=ispduptmp4; i<=ispduptmp5; ++i)

/** For all elements in the array, accessed physically
 *
 * This macro is used to generate loops for the vector in an easy way using
 * physical indexes. It defines internal the index 'i' which ranges the vector
 * using its physical definition.
 *
 * @code
 * FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(v)
 * {
 *     std::cout << DIRECT_A2D_ELEM(v, i) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(v) \
    for (size_t i=0; i<v.xdim; ++i)

/** Macro to check whether a point is inside or outside a given matrix. */
#define OUTSIDE(i,j) \
            ((j) < STARTINGX(*this) || (j) > FINISHINGX(*this) || \
             (i) < STARTINGY(*this) || (i) > FINISHINGY(*this))

/** Macro to check whether a point is inside or outside a given matrix. */
#define OUTSIDE3D(k, i,j) \
            ((j) < STARTINGX(*this) || (j) > FINISHINGX(*this) || \
             (i) < STARTINGY(*this) || (i) > FINISHINGY(*this) || \
             (k) < STARTINGZ(*this) || (k) > FINISHINGZ(*this))
//@}

// Forward declarations ====================================================
template<typename T>
class MultidimArray;

template<typename T>
void coreArrayByScalar(const MultidimArray<T>& op1, const T& op2,
                       MultidimArray<T>& result, char operation);

template<typename T>
void coreScalarByArray(const T& op1, const MultidimArray<T>& op2,
                       MultidimArray<T>& result, char operation);

template<typename T>
void coreArrayByArray(const MultidimArray<T>& op1, const MultidimArray<T>& op2,
                      MultidimArray<T>& result, char operation);

/**
 *  Structure with the dimensions information of an image
 */
struct ArrayDim
{
    // Number of images
    size_t ndim;
    // Number of elements in Z
    size_t zdim;
    // Number of elements in Y
    size_t ydim;
    // Number of elements in X
    size_t xdim;
    // Number of elements in YX
    size_t yxdim;
    // Number of elements in ZYX
    size_t zyxdim;
    // Number of elements in NZYX
    size_t nzyxdim;

    ArrayDim()
    {
    	ndim = 0;
    	zdim = 0;
    	ydim = 0;
    	xdim = 0;
    	yxdim = 0;
    	zyxdim = 0;
    	nzyxdim = 0;
    }

    bool operator==(ArrayDim &adim)
    {
        return (this->ndim == adim.ndim &&
                this->zdim == adim.zdim &&
                this->ydim == adim.ydim &&
                this->xdim == adim.xdim );
    }
}
;

/**
 *  Structure with a set of coordinates in an image
 */
struct ArrayCoord
{
    // Number of images
    size_t n;
    // Number of elements in Z
    int z;
    // Number of elements in Y
    int y;
    // Number of elements in X
    int x;
} ;

/**
 * Possible views for 3D MuldimArray
 */
typedef enum
{
    VIEW_Z_NEG,    // Front view (Z negative)
    VIEW_Z_POS,    //  Z positve
    VIEW_Y_NEG,    // Align -Y axis to Z axis, rotating 90 degrees around X axis");
    VIEW_Y_POS, // Align Y axis to Z axis, rotating -90 degrees around X axis");
    VIEW_X_NEG,   // Align -X axis to Z axis, rotating -90 degrees around Y axis");
    VIEW_X_POS   // Align X axis to Z axis, rotating 90 degrees around Y axis");
} AxisView;

/** Template class for Xmipp arrays.
  * This class provides physical and logical access.
*/
class MultidimArrayBase
{
public:
    // Destroy data
    bool destroyData;

    // Number of images
    size_t ndim;

    // Number of elements in Z
    size_t zdim;

    // Number of elements in Y
    size_t ydim;

    // Number of elements in X
    size_t xdim;

    // Number of elements in YX
    size_t yxdim;

    // Number of elements in ZYX
    size_t zyxdim;

    // Number of elements in NZYX
    size_t nzyxdim;

    // Z init
    int zinit;

    // Y init
    int yinit;

    // X init
    int xinit;

    //Alloc memory or map to a file
    bool     mmapOn;
    //Mapped file handler
    FILE*      mFd;
    // Number of elements in NZYX in allocated memory
    size_t nzyxdimAlloc;
public:
    virtual ~MultidimArrayBase()
    {}

    // Virtual declarations to be used from MultidimArrayGeneric
    virtual void clear() = 0;
    virtual void selfReverseX() = 0;
    virtual void selfReverseY() = 0;
    virtual void selfReverseZ() = 0;
    virtual double computeAvg() const = 0;
    virtual void computeDoubleMinMaxRange(double& minval, double& maxval,size_t offset, size_t size) const = 0;
    virtual void maxIndex(int &lmax, int& kmax, int& imax, int& jmax) const = 0;
    virtual void coreAllocateReuse() = 0;
    virtual void coreDeallocate()= 0;

    /* return the value of the data pointer
     */
    virtual void * getArrayPointer() const = 0;

    /// @name Size
    //@{

    /** Sets new N dimension.
        *
        *  Note that the dataArray is NOT resized. This should be done separately with coreAllocate()
        */
    void setNdim(int Ndim);

    /** Sets new Z dimension.
     *
     *  Note that the dataArray is NOT resized. This should be done separately with coreAllocate()
     *
     */
    void setZdim(int Zdim);

    /** Sets new Y dimension.
     *
     *  Note that the dataArray is NOT resized. This should be done separately with coreAllocate()
     *
     */
    void setYdim(int Ydim);

    /** Sets new X dimension.
      *
      *  Note that the dataArray is NOT resized. This should be done separately with coreAllocate()
      *
      */
    void setXdim(int Xdim);

    /** Sets new 4D dimensions.
     *  Note that the dataArray is NOT resized. This should be done separately with coreAllocate()
     */
    void setDimensions(int Xdim, int Ydim, int Zdim, int Ndim);

    /** Sets new 4D dimensions.
     *  Note that the dataArray is NOT resized. This should be done separately with coreAllocate()
     */
    void setDimensions(ArrayDim &newDim);

    /** Get the array dimensions.
     */
    void getDimensions(size_t& Xdim, size_t& Ydim, size_t& Zdim, size_t &Ndim) const;
    void getDimensions(ArrayDim &idim) const;
    ArrayDim getDimensions() const;


    /** Get dimensions.
     *
     * Returns the size of the object in a 4D vector. If the object is a matrix
     * or a vector, then the higher order dimensions will be set to 1, ie,
     * (Xdim, 1, 1) or (Xdim, Ydim, 1).
     *
     * This function is not ported to Python.
     */
    void getDimensions(int* size) const;

    /** Returns the total size of the multidimArray
     *
     * @code
     * if (V.getSize() > 1) ...
     * @endcode
     */
    size_t getSize() const;

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
    virtual void resize(size_t Ndim, size_t Zdim, size_t Ydim, size_t Xdim, bool copy=true) = 0;

    /** Resize a single 3D image
     *
     * This function assumes n is 1
     * @code
     * V1.resize(3, 3, 2);
     * @endcode
     */
    void resize(size_t Zdim, size_t Ydim, size_t Xdim)
    {
        resize(1, Zdim, Ydim, Xdim);
    }

    /** Resize a single 2D image
     *
     * This function assumes n and z are 1
     * @code
     * V1.resize(3, 2);
     * @endcode
     */
    void resize(size_t Ydim, size_t Xdim)
    {
        resize(1, 1, Ydim, Xdim);
    }

    /** Resize a single 1D image
     *
     * This function assumes n and z and y are 1
     * @code
     * V1.resize(2);
     * @endcode
     */
    void resize(size_t Xdim)
    {
        resize(1, 1, 1, Xdim);
    }

    /** Resize an image using the dimensions
     *  from an ArrayDim structure.
     */
    void resize(ArrayDim &adim, bool copy=true);

    /** Resize with no copy a single 3D image
        */
    void resizeNoCopy(size_t Ndim, size_t Zdim, size_t Ydim, size_t Xdim)
    {
        resize(Ndim, Zdim, Ydim, Xdim, false);
    }

    /** Resize with no copy a single 3D image
        */
    void resizeNoCopy(size_t Zdim, size_t Ydim, size_t Xdim)
    {
        resize(1, Zdim, Ydim, Xdim, false);
    }

    /** Resize a single 2D image with no copy
     */
    void resizeNoCopy(size_t Ydim, size_t Xdim)
    {
        resize(1, 1, Ydim, Xdim, false);
    }

    /** Resize a single 1D image with no copy
     */
    void resizeNoCopy(size_t Xdim)
    {
        resize(1, 1, 1, Xdim, false);
    }

    /** Returns Y dimension.
       */
    inline size_t rowNumber() const
    {
        return ydim;
    }

    /** Returns X dimension.
     */
    inline size_t colNumber() const
    {
        return xdim;
    }

    /** Copy the shape parameters
      *
      */
    void copyShape(const MultidimArrayBase &m);

    /** Same shape.
     *
     * Returns true if this object has got the same shape (origin and size)
     * than the argument
     */
    inline bool sameShape(const MultidimArrayBase &op) const
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
      *
      * This function adjust the starting points in the array such that the
      * center of the array is defined in the Xmipp fashion.
      *
      * @code
      * V.setXmippOrigin();
      * @endcode
      */
    void setXmippOrigin();

    /** Reset logical origin to zeros.
     *
     * This function adjust the starting points in the array such
     * that upper left corner begins in zero.
     *
     * @code
     * V.resetOrigin();
     * @endcode
     */
    void resetOrigin();

    /** Move origin to.
      *
      * This function adjust logical indexes such that the Xmipp origin of the
      * array moves to the specified position. For instance, an array whose x
      * indexes go from -1 to 1, if we move the origin to 4, then the x indexes
      * go from 3 to 5. This is very useful for convolution operations where you
      * only need to move the logical starting of the array.
      *
      */
    void moveOriginTo(int k, int i, int j);

    /** Move origin to.
      *
      * This function adjust logical indexes such that the Xmipp origin of the
      * array moves to the specified position. For instance, an array whose x
      * indexes go from -1 to 1, if we move the origin to 4, then the x indexes
      * go from 3 to 5. This is very useful for convolution operations where you
      * only need to move the logical starting of the array.
      *
      */
    void moveOriginTo(int i, int j);

    /** Returns the first valid logical Z index.
      */
    inline int startingZ() const
    {
        return zinit;
    }

    /** Returns the last valid logical Z index.
     */
    inline int finishingZ() const
    {
        return zinit + zdim - 1;
    }

    /** Returns the first valid logical Y index.
     */
    inline int startingY() const
    {
        return yinit;
    }

    /** Returns the last valid logical Y index.
     */
    inline int finishingY() const
    {
        return yinit + ydim - 1;
    }

    /** Returns the first valid logical X index.
     */
    inline int startingX() const
    {
        return xinit;
    }

    /** Returns the last valid logical X index.
     */
    inline int finishingX() const
    {
        return xinit + xdim - 1;
    }

    /** IsCorner (in 2D or 3D matrix)
         *
         * TRUE if the logical index given is a corner of the definition region of this
         * array.
         */
    bool isCorner(const Matrix1D< double >& v) const;

    /** Outside for 3D matrices
       *
       * TRUE if the logical index given is outside the definition region of this
       * array.
       */
    inline bool outside(int k, int i, int j) const
    {
        return (j < STARTINGX(*this) || j > FINISHINGX(*this) ||
                i < STARTINGY(*this) || i > FINISHINGY(*this) ||
                k < STARTINGZ(*this) || k > FINISHINGZ(*this));
    }

    /** Outside for 2D matrices
     *
     * TRUE if the logical index given is outside the definition region of this
     * array.
     */
    inline bool outside(int i, int j) const
    {
        return (j < STARTINGX(*this) || j > FINISHINGX(*this) ||
                i < STARTINGY(*this) || i > FINISHINGY(*this));
    }

    /** Outside for 1D matrices
     *
     * TRUE if the logical index given is outside the definition region of this
     * array.
     */
    inline bool outside(int i) const
    {
        return (i < STARTINGX(*this) || i > FINISHINGX(*this));
    }

    /** Outside
     *
     * TRUE if the logical index given is outside the definition region of this
     * array.
     */
    bool outside(const Matrix1D<double> &r) const;

    //@}


    /** Returns the multidimArray dimension.
     *
     * @code
     * int dim = V.getDim();
     * @endcode
     */
    inline int getDim() const
    {
        if (NZYXSIZE(*this) < 1)
            return 0;
        if (NSIZE(*this) > 1)
            return 4;
        if (ZSIZE(*this) > 1)
            return 3;
        if (YSIZE(*this) > 1)
            return 2;
        return 1;
    }

    /** Sets mmap.
     *
     * Sets on/off mmap flag to allocate memory in a file.
     *
     */
    void setMmap(bool mmap)
    {
        coreDeallocate();
        mmapOn = mmap;
    }

    void maxIndex(ArrayCoord &pos) const
    {
        maxIndex((int&)pos.n, pos.z, pos.y, pos.x);
    }

    /** 3D Indices for the maximum element.
      *
      * This function just calls to the 4D function
      */
    void maxIndex(int& kmax, int& imax, int& jmax) const
    {
        int dum;
        maxIndex(dum, kmax, imax, jmax);
    }

    /** 2D Indices for the maximum element.
     *
     * This function just calls to the 4D function
     */
    void maxIndex(int& imax, int& jmax) const
    {
        int dum;
        maxIndex(dum, dum, imax, jmax);
    }

    /** 1D Indices for the maximum element.
     *
     * This function just calls to the 4D function
     */
    void maxIndex(int& jmax) const
    {
        int dum;
        maxIndex(dum, dum, dum, jmax);
    }

    /** Print shape of multidimensional array.
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
    void printShape(std::ostream& out = std::cout) const;
};

template<typename T>
class MultidimArray: public MultidimArrayBase
{
public:
    /* The array itself.
       The array is always a 3D array (Z,Y,X). For vectors the size of the array
       is (1,1,X) and for matrices (1,Y,X). The pixel (i,j) (y,x) is at the
       position data[i*Xdim+j] or data[y*Xdim+x]
    */
    T* data;

public:
    /// @name Constructors
    //@{

    /** Empty constructor.
     * The empty constructor creates an array with no memory associated,
     * size=0.
     */
    MultidimArray()
    {
        coreInit();
    }

    /** Size constructor with 4D size.
     * The Size constructor creates an array with memory associated,
     * and fills it with zeros.
     */
    MultidimArray(size_t Ndim, int Zdim, int Ydim, int Xdim)
    {
        coreInit();
        coreAllocate(Ndim, Zdim, Ydim, Xdim);
    }

    /** Size constructor with 3D size.
     * The Size constructor creates an array with memory associated,
     * and fills it with zeros.
     */
    MultidimArray( int Zdim, int Ydim, int Xdim)
    {
        coreInit();
        coreAllocate(1UL, Zdim, Ydim, Xdim);
    }

    /** Size constructor with 2D size.
     * The Size constructor creates an array with memory associated,
     * and fills it with zeros.
     */
    MultidimArray(int Ydim, int Xdim)
    {
        coreInit();
        coreAllocate(1UL, 1, Ydim, Xdim);
    }

    /** Size constructor with 1D size.
     * The Size constructor creates an array with memory associated,
     * and fills it with zeros.
     */
    MultidimArray(int Xdim)
    {
        coreInit();
        coreAllocate(1UL, 1, 1, Xdim);
    }

    /** Copy constructor
     *
     * The created volume is a perfect copy of the input array but with a
     * different memory assignment.
     *
     * @code
     * MultidimArray< double > V2(V1);
     * @endcode
     */
    MultidimArray(const MultidimArray<T>& V)
    {
        coreInit();
        *this = V;
    }

    /** Copy constructor from a Matrix1D.
     * The Size constructor creates an array with memory associated,
     * and fills it with zeros.
     */
    MultidimArray(const Matrix1D<T>& V)
    {
        coreInit();
        coreAllocate(1, 1, 1, V.size());
        for (size_t i = 0; i < V.size(); i++)
            DIRECT_A1D_ELEM(*this,i) = VEC_ELEM(V,i);
    }

    /** Constructor from vector 1D
     * This will create a MultidimArray 1D
     * the size and elements will be copied from
     * the std::vector
     */
    MultidimArray(const std::vector<T> &vector)
    {
        coreInit();
        coreAllocate(1, 1, 1, vector.size());
        for (size_t i = 0; i < vector.size(); i++)
            DIRECT_A1D_ELEM(*this,i) = vector[i];
    }

    /** Destructor.
     */
    virtual ~MultidimArray()
    {
        coreDeallocate();
    }

    /** Clear.
     */
    void clear()
    {
        coreDeallocate();
        coreInit();
    }
    //@}

    /// @name Core memory operations
    //@{

    /** Core init.
     * Initialize everything to 0
     */
    void coreInit()
    {
        xdim = yxdim = zyxdim = nzyxdim = 0;
        ydim = zdim = ndim = 0;
        zinit = yinit = xinit = 0;
        data = NULL;
        nzyxdimAlloc = 0;
        destroyData = true;
        mmapOn = false;
        mFd = NULL;
    }

    /** Core allocate with dimensions.
     */
    void coreAllocate(size_t _ndim, int _zdim, int _ydim, int _xdim)
    {
        if (_ndim <= 0 || _zdim <= 0 || _ydim<=0 || _xdim<=0)
        {
            clear();
            return;
        }
        if(data!=NULL)
            REPORT_ERROR(ERR_MEM_NOTDEALLOC, "do not allocate space for an image if you have not deallocate it first");

        ndim=_ndim;
        zdim=_zdim;
        ydim=_ydim;
        xdim=_xdim;
        yxdim=(size_t)ydim*xdim;
        zyxdim=yxdim*zdim;
        nzyxdim=zyxdim*ndim;

        coreAllocate();
    }

    /** Core allocate without dimensions.
     *
     * It is supposed the dimensions are set previously with setXdim(x), setYdim(y)
     * setZdim(z), setNdim(n) or with setDimensions(Xdim, Ydim, Zdim, Ndim);
     *
     */
    void coreAllocate()
    {
        if(data!=NULL)
            REPORT_ERROR(ERR_MEM_NOTDEALLOC, "do not allocate space for an image if you have not deallocate it first");
        if (nzyxdim < 0)
            REPORT_ERROR(ERR_MEM_BADREQUEST,"coreAllocate:Cannot allocate a negative number of bytes");

        if (mmapOn)
            mFd = mmapFile(data, nzyxdim);
        else
        {
            try
            {
                data = new T [nzyxdim];
                if (data == NULL)
                {
                    setMmap(true);
                    mFd = mmapFile(data, nzyxdim);
                }
            }
            catch (std::bad_alloc &)
            {
                setMmap(true);
                mFd = mmapFile(data, nzyxdim);
            }
        }
        memset(data,0,nzyxdim*sizeof(T));
        nzyxdimAlloc = nzyxdim;
    }

    /** Core allocate without dimensions.
     *
     * It is supposed the dimensions are set previously with setXdim(x), setYdim(y)
     * setZdim(z), setNdim(n) or with setDimensions(Xdim, Ydim, Zdim, Ndim);
     *
     */
    void coreAllocateReuse()
    {
        if(data != NULL && nzyxdim <= nzyxdimAlloc)
            return;
        else if (nzyxdim > nzyxdimAlloc)
            coreDeallocate();

        if (nzyxdim < 0)
            REPORT_ERROR(ERR_MEM_BADREQUEST,"coreAllocateReuse:Cannot allocate a negative number of bytes");

        if (mmapOn)
            mFd = mmapFile(data, nzyxdim);
        else
        {
            data = new T [nzyxdim];
            if (data == NULL)
                REPORT_ERROR(ERR_MEM_NOTENOUGH, "Allocate: No space left");
        }
        memset(data,0,nzyxdim*sizeof(T));
        nzyxdimAlloc = nzyxdim;
    }

    /* Create a temporary file to mmap the data.
     */
    FILE* mmapFile(T* &_data, size_t nzyxDim) const
    {
#ifdef XMIPP_MMAP
        FILE* fMap = tmpfile();
        int Fd = fileno(fMap);

        if ((lseek(Fd, nzyxDim*sizeof(T)-1, SEEK_SET) == -1) || (::write(Fd,"",1) == -1))
        {
            fclose(fMap);
            REPORT_ERROR(ERR_IO_NOWRITE,"MultidimArray::resize: Error 'stretching' the map file.");
        }
        if ( (_data = (T*) mmap(0,nzyxDim*sizeof(T), PROT_READ | PROT_WRITE, MAP_SHARED, Fd, 0)) == (void*) MAP_FAILED )
            REPORT_ERROR(ERR_MMAP_NOTADDR,formatString("MultidimArray::resize: mmap failed. Error %s", strerror(errno)));

        return fMap;
#else

        REPORT_ERROR(ERR_MMAP,"Mapping not supported in Windows");
#endif

    }

    /** Core deallocate.
     * Free all data.
     */
    void coreDeallocate()
    {
        if (data != NULL && destroyData)
        {
            if (mmapOn)
            {
#ifdef XMIPP_MMAP
                munmap(data,nzyxdimAlloc*sizeof(T));
                fclose(mFd);
#else

                REPORT_ERROR(ERR_MMAP,"Mapping not supported in Windows");
#endif

            }
            else
                delete[] data;
        }
        data = NULL;
        destroyData = true;
        nzyxdimAlloc = 0;
    }

    /** Alias a multidimarray.
     *
     * Treat the multidimarray as if it were a volume. The data is not copied
     * into new memory, but a pointer to the multidimarray is copied.
     * You should not make any operation on this volume such that the
     * memory locations are changed
     */
    void alias(const MultidimArray<T> &m)
    {
        coreDeallocate();
        copyShape(m);
        this->data=m.data;
        this->nzyxdimAlloc = this->nzyxdim;
        this->destroyData = false;
    }

    /** Alias a row in an image.
         *
         * Treat the multidimarray as if it were a single slice. The data is not copied
         * into new memory, but a pointer to the selected slice in the multidimarray is copied.
         * You should not make any operation on this volume such that the
         * memory locations are changed.
         * Select_slice starts at 0 towards Zsize.
         */
    void aliasRow(const MultidimArray<T> &m, size_t select_row)
    {
        if (select_row >= YSIZE(m))
            REPORT_ERROR(ERR_MULTIDIM_SIZE, "aliasRow: Selected row cannot be higher than Y size.");

        coreDeallocate();
        setDimensions(XSIZE(m),1, 1, 1);
        this->data = m.data + XSIZE(m)*select_row;
        this->nzyxdimAlloc = this->nzyxdim;
        this->destroyData = false;
    }

    /** Alias a slice in a multidimarray.
         *
         * Treat the multidimarray as if it were a single slice. The data is not copied
         * into new memory, but a pointer to the selected slice in the multidimarray is copied.
         * You should not make any operation on this volume such that the
         * memory locations are changed.
         * Select_slice starts at 0 towards Zsize.
         */
    void aliasSlice(const MultidimArray<T> &m, size_t select_slice)
    {
        if (select_slice >= ZSIZE(m))
            REPORT_ERROR(ERR_MULTIDIM_SIZE, "aliasSlice: Selected slice cannot be higher than Z size.");

        coreDeallocate();
        setDimensions(XSIZE(m), YSIZE(m), 1, 1);
        this->data = m.data + XSIZE(m)*YSIZE(m)*(select_slice);
        this->nzyxdimAlloc = this->nzyxdim;
        this->destroyData = false;
    }

    /** Alias an image in a stack.
         *
         * Treat the multidimarray as if it were a single slice. The data is not copied
         * into new memory, but a pointer to the selected image in the multidimarray is copied.
         * You should not make any operation on this volume such that the
         * memory locations are changed.
         * Select_slice starts at 0 towards Nsize.
         */
    void aliasImageInStack(const MultidimArray<T> &m, size_t select_image)
    {
        if (select_image >= NSIZE(m))
            REPORT_ERROR(ERR_MULTIDIM_SIZE, "aliasImageInStack: Selected image cannot be higher than N size.");
        if (ZSIZE(m)!=1)
            REPORT_ERROR(ERR_MULTIDIM_SIZE, "aliasImageInStack: This function is not meant for volumes");

        coreDeallocate();
        setDimensions(XSIZE(m), YSIZE(m), 1, 1);
        this->data = m.data + XSIZE(m)*YSIZE(m)*(select_image);
        this->nzyxdimAlloc = this->nzyxdim;
        this->destroyData = false;
    }



    //@}

    /// @name Size
    //@{

    /* These "using" declarations must be done due to c++ cannot overload methods that have been
     * declared in base class.
     */
    using MultidimArrayBase::resizeNoCopy;
    using MultidimArrayBase::resize;
    using MultidimArrayBase::maxIndex;

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
    void resize(size_t Ndim, size_t Zdim, size_t Ydim, size_t Xdim, bool copy=true)
    {
        if (Ndim*Zdim*Ydim*Xdim == nzyxdimAlloc && data != NULL)
        {
            ndim = Ndim;
            xdim = Xdim;
            ydim = Ydim;
            zdim = Zdim;
            yxdim = Ydim * Xdim;
            zyxdim = Zdim * yxdim;
            nzyxdim = Ndim * zyxdim;
            return;
        }
        else if (!destroyData)
            REPORT_ERROR(ERR_MULTIDIM_SIZE, "Cannot resize array when accessing through alias.");

        if (Xdim <= 0 || Ydim <= 0 || Zdim <= 0 || Ndim <= 0)
        {
            clear();
            return;
        }

        // data can be NULL while xdim etc are set to non-zero values
        // (This can happen for reading of images...)
        // In that case, initialize data to zeros.
        if (NZYXSIZE(*this) > 0 && data == NULL)
        {
            ndim = Ndim;
            xdim = Xdim;
            ydim = Ydim;
            zdim = Zdim;
            yxdim = Ydim * Xdim;
            zyxdim = Zdim * yxdim;
            nzyxdim = Ndim * zyxdim;

            coreAllocate();
            return;
        }

        // Ask for memory
        size_t YXdim=(size_t)Ydim*Xdim;
        size_t ZYXdim=YXdim*Zdim;
        size_t NZYXdim=ZYXdim*Ndim;
        FILE*  new_mFd=NULL;

        T * new_data=NULL;

        try
        {
            if (mmapOn)
                new_mFd = mmapFile(new_data, NZYXdim);
            else
                new_data = new T [NZYXdim];

            memset(new_data,0,NZYXdim*sizeof(T));
        }
        catch (std::bad_alloc &)
        {
            if (!mmapOn)
            {
                setMmap(true);
                resize(Ndim, Zdim, Ydim, Xdim, copy);
                return;
            }
            else
            {
                std::ostringstream sstream;
                sstream << "Allocate: No space left to allocate ";
                sstream << (NZYXdim * sizeof(T)/1024/1024/1024) ;
                sstream << "Gb." ;
                REPORT_ERROR(ERR_MEM_NOTENOUGH, sstream.str());
            }
        }
        // Copy needed elements, fill with 0 if necessary
        if (copy)
        {
            T zero=0; // Very useful for complex matrices
            T *val=NULL;
            for (size_t l = 0; l < Ndim; l++)
                for (size_t k = 0; k < Zdim; k++)
                    for (size_t i = 0; i < Ydim; i++)
                        for (size_t j = 0; j < Xdim; j++)
                        {
                            if (l >= NSIZE(*this))
                                val = &zero;
                            else if (k >= ZSIZE(*this))
                                val = &zero;
                            else if (i >= YSIZE(*this))
                                val = &zero;
                            else if (j >= XSIZE(*this))
                                val = &zero;
                            else
                                val = &DIRECT_NZYX_ELEM(*this, l, k, i, j);
                            new_data[l*ZYXdim + k*YXdim+i*Xdim+j] = *val;
                        }
        }

        // deallocate old array
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
        mFd = new_mFd;
        nzyxdimAlloc = nzyxdim;
    }

    /** Resize according to a pattern.
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
    void resize(const MultidimArray<T1> &v, bool copy = true)
    {
        if (NSIZE(*this) != NSIZE(v) || XSIZE(*this) != XSIZE(v) ||
            YSIZE(*this) != YSIZE(v) || ZSIZE(*this) != ZSIZE(v) || data==NULL)
            resize(NSIZE(v), ZSIZE(v), YSIZE(v), XSIZE(v), copy);

        STARTINGX(*this) = STARTINGX(v);
        STARTINGY(*this) = STARTINGY(v);
        STARTINGZ(*this) = STARTINGZ(v);
    }

    /** Resize according to a pattern with no copy.
     */
    template<typename T1>
    void resizeNoCopy(const MultidimArray<T1> &v)
    {
        if (NSIZE(*this) != NSIZE(v) || XSIZE(*this) != XSIZE(v) ||
            YSIZE(*this) != YSIZE(v) || ZSIZE(*this) != ZSIZE(v) || data==NULL)
            resize(NSIZE(v), ZSIZE(v), YSIZE(v), XSIZE(v), false);

        STARTINGX(*this) = STARTINGX(v);
        STARTINGY(*this) = STARTINGY(v);
        STARTINGZ(*this) = STARTINGZ(v);
    }

    /** Check dimension.
     *
     * returns true if the dimension is equal to the argument and false otherwise
     * It also prints an error message in the latter case.
     */
#define checkDimension(dim) checkDimensionWithDebug(dim,__FILE__,__LINE__)
    void checkDimensionWithDebug(int dim, const char *file, int line) const
    {
        if (getDim() != dim)
        {
            std::cerr<<" Check for dimension: "  << dim <<std::endl;
            std::cerr << "MultidimArray shape: ";
            printShape(std::cerr);
            std::cerr << std::endl;
            std::cerr << "Check called from file "<<file<<" line "<<line<<std::endl;
            exit(1);
        }
    }

    /** Generic selfWindow routine (dim independent)
     *
     * This function will call to 3D,2D or 1D specific selfWindow routines
     */
    void selfWindow(int n0,int z0, int y0, int x0,
                    int nF,int zF, int yF, int xF,
                    T init_value = 0)
    {
        if (this->ndim >1)
            REPORT_ERROR(ERR_MULTIDIM_DIM,"stack windowing not implemented");
        if (this->zdim >1)
        {//call 3Dwindow
            selfWindow( z0,  y0,  x0,
                        zF,  yF,  xF,
                        init_value);
        }
        else if (this->ydim >1)
        {//call 2Dwindow
            selfWindow( y0,  x0,
                        yF,  xF,
                        init_value);

        }
        else if (this->xdim >1)
        {//call 1Dwindow
            selfWindow( x0, xF, init_value);
        }
    }

    template <class T1>
    void window(MultidimArray<T1> &result, int n0,int z0, int y0, int x0,
                int nF,int zF, int yF, int xF,
                T1 init_value = 0) const
    {
        if (this->ndim >1)
            REPORT_ERROR(ERR_MULTIDIM_DIM,"stack windowing not implemented");
        if (this->zdim >1)
        {//call 3Dwindow
            window(result, z0,  y0,  x0,
                   zF,  yF,  xF,
                   init_value);
        }
        else if (this->ydim >1)
        {//call 2Dwindow
            window(result, y0,  x0,
                   yF,  xF,
                   init_value);
        }
        else if (this->xdim >1)
        {//call 1Dwindow
            window(result, x0, xF, init_value);
        }
    }

    /** Put a 3D selfWindow to the nth volume
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
     * V1.selfWindow(0, 0, -1, 1, 1, 2);
     * @endcode
     */
    template <class T1>
    void window(MultidimArray<T1> &result, int z0, int y0, int x0, int zF, int yF, int xF,
                T1 init_value = 0) const
    {
        result.resizeNoCopy(zF - z0 + 1, yF - y0 + 1, xF - x0 + 1);
        STARTINGZ(result) = z0;
        STARTINGY(result) = y0;
        STARTINGX(result) = x0;

        FOR_ALL_ELEMENTS_IN_ARRAY3D(result)
        if ((k >= STARTINGZ(*this) && k <= FINISHINGZ(*this)) &&
            (i >= STARTINGY(*this) && i <= FINISHINGY(*this)) &&
            (j >= STARTINGX(*this) && j <= FINISHINGX(*this)))
            A3D_ELEM(result, k, i, j) = (T1) A3D_ELEM(*this, k, i, j);
        else
            A3D_ELEM(result, k, i, j) = init_value;
    }

    /** 3D Self window */
    void selfWindow(int z0, int y0, int x0, int zF, int yF, int xF,
                    T init_value = 0)
    {
        if (z0 == STARTINGZ(*this) && zF == FINISHINGZ(*this) &&
            y0 == STARTINGY(*this) && yF == FINISHINGY(*this) &&
            x0 == STARTINGX(*this) && xF == FINISHINGX(*this))
            return;

        MultidimArray<T> result;
        window(result,z0,y0,x0,zF,yF,xF,init_value);
        *this=result;
    }

    /** Put a 2D selfWindow to the nth matrix
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
     * m1.selfWindow(-1, -1, 1, 2);
     * @endcode
     */
    template <class T1>
    void window(MultidimArray<T1> &result, int y0, int x0, int yF, int xF,
                T1 init_value = 0) const
    {
        result.resizeNoCopy(yF - y0 + 1, xF - x0 + 1);
        STARTINGY(result) = y0;
        STARTINGX(result) = x0;

        FOR_ALL_ELEMENTS_IN_ARRAY2D(result)
        if (j >= STARTINGX(*this) && j <= FINISHINGX(*this) &&
            i >= STARTINGY(*this) && i <= FINISHINGY(*this))
            A2D_ELEM(result, i, j) = A2D_ELEM(*this, i, j);
        else
            A2D_ELEM(result, i, j) = init_value;
    }

    /** 2D Self window */
    void selfWindow(int y0, int x0, int yF, int xF,
                    T init_value = 0)
    {
        MultidimArray<T> result;
        window(result,y0,x0,yF,xF,init_value);
        *this=result;
    }

    /** Put a 1D selfWindow to the nth vector
     *
     * The vector is windowed within the two indexes given to this function.
     * Indexes always refer to logical indexes. If an index is outside the
     * actual vector range then the vector is padded winit_value. In the
     * following examples suppose that v1=[-2 -1 0 1 2] and that the origin is
     * -2.
     *
     * @code
     * v1.selfWindow(-1, 2); // v1=[-1 0 1 2]; v1.startingX() == -1
     *
     * v1.selfWindow(-3, 1); // v1=[0 -2 -1 0 1]; v1.startingX() == -3
     * @endcode
     */
    template <class T1>
    void window(MultidimArray<T1> &result, int x0, int xF, T1 init_value = 0) const
    {
        result.resizeNoCopy(xF - x0 + 1);
        STARTINGX(result) = x0;

        FOR_ALL_ELEMENTS_IN_ARRAY1D(result)
        if (i >= STARTINGX(*this) && i <= FINISHINGX(*this))
            A1D_ELEM(result, i) = A1D_ELEM(*this, i);
        else
            A1D_ELEM(result, i) = init_value;
    }

    /** 1D Self window */
    void selfWindow(int x0, int xF,
                    T init_value = 0)
    {
        MultidimArray<T> result;
        window(result,x0,xF,init_value);
        *this=result;
    }

    /** Make a patch with the input array in the given positions */
    void patch(MultidimArray<T> patchArray, int x, int y)
    {
        int n = XSIZE(patchArray)*sizeof(T);

        for (size_t i=0; i < YSIZE(patchArray); ++i)
            memcpy(&dAij(*this, y+i, x), &dAij(patchArray, i, 0), n);
    }

    //@}

    ///@name Access to the pixel values
    //@{

    /** Volume element access via double vector.
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
        switch (VEC_XSIZE(v))
        {
        case 1:
            return A1D_ELEM((*this), ROUND(XX(v)));
        case 2:
            return A2D_ELEM((*this), ROUND(YY(v)), ROUND(XX(v)));
        case 3:
            return A3D_ELEM((*this), ROUND(ZZ(v)), ROUND(YY(v)), ROUND(XX(v)));
        default:
            REPORT_ERROR(ERR_ARG_INCORRECT,"Cannot handle indexes with dimension larger than 3");
        }
    }

    /** Volume element access via integer vector.
     */
    T& operator()(const Matrix1D< int >& v) const
    {
        switch (VEC_XSIZE(v))
        {
        case 1:
            return A1D_ELEM((*this), XX(v));
        case 2:
            return A2D_ELEM((*this), YY(v), XX(v));
        case 3:
            return A3D_ELEM((*this), ZZ(v), YY(v), XX(v));
        default:
            REPORT_ERROR(ERR_ARG_INCORRECT,"Cannot handle indexes with dimension larger than 3");
        }
    }

    /** 4D element access via index.
    *
    * Returns the value of a matrix logical position. In our example we could
    * access from v(0, 0,-2,-1) to v(0, 1,2,1). The elements can be used either by
    * value or by reference. An exception is thrown if the index is outside
    * the logical range. Be careful that the argument order is (Z,Y,X).
    *
    * @code
    * V(0, 0, -2, 1) = 1;
    * val = V(0, 0, -2, 1);
    * @endcode
    */
    inline T& operator()(size_t n, int k, int i, int j) const
    {
        return NZYX_ELEM(*this, n, k, i, j);
    }

    /** 3D element access via index.
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
    inline T& operator()(int k, int i, int j) const
    {
        return A3D_ELEM(*this, k, i, j);
    }

    /** Matrix element access via index
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
    inline T& operator()(int i, int j) const
    {
        return A2D_ELEM(*this, i, j);
    }

    /** Vector element access
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
    inline T& operator()(int i) const
    {
        return A1D_ELEM(*this, i);
    }


    /** Return the void pointer to the internal data array
     */
    void* getArrayPointer() const
    {
        return (void*) data;

    }

    /** Copy an image from a stack to another
     *
     * Copy image image n from this MDA to image n2 in MDA M.
     * @code
     * V.getImage(0, m, 3);
     * @endcode
     */
    void getImage(size_t n, MultidimArray<T>& M, size_t n2 = 0) const
    {
        if (XSIZE(*this) == 0)
        {
            M.clear();
            return;
        }

        if (n > NSIZE(*this))
            REPORT_ERROR(ERR_INDEX_OUTOFBOUNDS," Multidimarray getImage: n larger than NSIZE");

        if (ZSIZE(*this) != ZSIZE(M) || YSIZE(*this) != YSIZE(M) || XSIZE(*this) != XSIZE(M))
        {
            if (n2 == 0)
                M.resizeNoCopy(*this);
            else
                REPORT_ERROR(ERR_MULTIDIM_SIZE, "MultidimArray::getImage: Target dimensions do not match source dimensions.");
        }

        if (n2 > NSIZE(M))
            REPORT_ERROR(ERR_INDEX_OUTOFBOUNDS," Multidimarray getImage: n larger than MultidimArray target NSIZE");


        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(M)
        DIRECT_NZYX_ELEM(M, n2, k, i, j) = DIRECT_NZYX_ELEM(*this, n, k, i, j);

        STARTINGX(M) = STARTINGX(*this);
        STARTINGY(M) = STARTINGY(*this);
        STARTINGZ(M) = STARTINGZ(*this);
    }

    /** 2D Slice access for reading.
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
    template <typename T1>
    void getSlice(int k, MultidimArray<T1>& M, char axis = 'Z', bool reverse = false, size_t n = 0) const
    {
        if (XSIZE(*this) == 0)
        {
            M.clear();
            return;
        }

        T* ptr=NULL;
        switch (axis)
        {
        case 'Z':
            if (k < STARTINGZ(*this) || k > FINISHINGZ(*this))
                REPORT_ERROR(ERR_INDEX_OUTOFBOUNDS,
                             "Slice: Multidim subscript (k) out of range");

            M.resize(1, 1, YSIZE(*this), XSIZE(*this),false);
            k = k - STARTINGZ(*this);

            if (reverse)
            {
                int zEnd = ZSIZE(*this) - 1;
                FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(M)
                DIRECT_A2D_ELEM(M, i, j) = (T1) DIRECT_NZYX_ELEM(*this, n, k,  zEnd-i, j);
            }
            else
            {
                ptr=&(DIRECT_NZYX_ELEM(*this, n, k, 0, 0));
                FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(M)
                DIRECT_MULTIDIM_ELEM(M, n) = (T1) *(ptr++);
            }

            STARTINGX(M) = STARTINGX(*this);
            STARTINGY(M) = STARTINGY(*this);
            break;
        case 'Y':
            if (k < STARTINGY(*this) || k > FINISHINGY(*this))
                REPORT_ERROR(ERR_INDEX_OUTOFBOUNDS,
                             "Slice: Multidim subscript (i) out of range");

            k = k - STARTINGY(*this);
            M.resizeNoCopy(ZSIZE(*this), XSIZE(*this));

            if (reverse)
            {
                int zEnd = ZSIZE(*this) - 1;
                FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(M)
                DIRECT_A2D_ELEM(M, i, j) = (T1) DIRECT_NZYX_ELEM(*this, n, zEnd-i, k, j);
            }
            else
            {
                FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(M)
                DIRECT_A2D_ELEM(M, i, j) = (T1) DIRECT_NZYX_ELEM(*this, n, i, k, j);
            }

            STARTINGX(M) = STARTINGX(*this);
            STARTINGY(M) = STARTINGZ(*this);
            break;
        case 'X':
            if (k < STARTINGX(*this) || k > FINISHINGX(*this))
                REPORT_ERROR(ERR_INDEX_OUTOFBOUNDS,
                             "Slice: Multidim subscript (j) out of range");

            k = k - STARTINGX(*this);
            M.resizeNoCopy(YSIZE(*this), ZSIZE(*this));

            if (reverse)
            {
                int zEnd = ZSIZE(*this) - 1;
                FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(M)
                DIRECT_A2D_ELEM(M, i, j) = (T1) DIRECT_NZYX_ELEM(*this, n, zEnd-j, i, k);
            }
            else
            {
                FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(M)
                DIRECT_A2D_ELEM(M, i, j) = (T1) DIRECT_NZYX_ELEM(*this, n, j, i, k);
            }

            STARTINGX(M) = STARTINGY(*this);
            STARTINGY(M) = STARTINGZ(*this);
            break;
        default:
            REPORT_ERROR(ERR_VALUE_INCORRECT,
                         formatString("Slice: not supported axis %c", axis));
        }
    }

    /** Get Z slice as matrix */
    void getSliceAsMatrix(size_t k, Matrix2D<T> &m) const
    {
        m.resizeNoCopy(YSIZE(*this),XSIZE(*this));
        memcpy(&MAT_ELEM(m,0,0),&A3D_ELEM(*this,k,0,0),YSIZE(*this),XSIZE(*this)*sizeof(double));
    }


    /** Return the data aliased as a row vector in a Matrix1D
     */
    void getAliasAsRowVector(Matrix1D<T> &m) const
    {
        m.vdim = NZYXSIZE(*this);
        m.destroyData = false;
        m.row = true;
        m.vdata = data;
    }

    /** Slice access for writing.
     *
     * This function sets a 2D matrix corresponding to the chosen slice inside the nth
     * volume, the numbering of the slices is also logical not physical.
     *
     * @code
     * // Copies slice 0 in slice 1
     * V.setSlice(1, (V.slice(0)));
     * @endcode
     */
    template <typename T1>
    void setSlice(int k, const MultidimArray<T1>& v, size_t n=0)
    {
        if (xdim == 0)
            return;

        if (k < STARTINGZ(*this) || k > FINISHINGZ(*this))
            REPORT_ERROR(ERR_INDEX_OUTOFBOUNDS,formatString(
                             "setSlice: MultidimArray subscript (k=%d) out of range [%d, %d]",
                             k,STARTINGZ(*this),FINISHINGZ(*this)));

        if (v.rowNumber() != YSIZE(*this) || v.colNumber() != XSIZE(*this))
            REPORT_ERROR(ERR_MULTIDIM_DIM,
                         "setSlice: MultidimArray dimensions different from the matrix ones");

        k-=STARTINGZ(*this);
        T *ptr=&(DIRECT_NZYX_ELEM(*this, n, k, 0, 0));
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(v)
        *(ptr++) = (T) DIRECT_MULTIDIM_ELEM(v, n);
    }

    /** Reslice the volume aliging any X or Y direction with Z axis
     *
     * @param face Select the face to become the new Z direction
     * @param out  The resliced volume is returned
     * @param flip Invert the positions of Z planes, keeping the X-Y orientation
     * @param n Select the number of image in case of stacks
     */

    template <typename T1>
    void reslice(MultidimArray<T1>& out, AxisView face, bool flip = false, size_t n = 0) const
    {
        ArrayDim aDim, aDimOut;
        getDimensions(aDim);

        char axis;
        bool reverse;

        aDimOut = aDim;

        if (face == VIEW_Y_NEG || face == VIEW_Y_POS)
        {
            axis = 'Y';
            aDimOut.ydim = aDim.zdim;
            aDimOut.zdim = aDim.ydim;
            reverse = (face == VIEW_Y_NEG);
        }
        else if (face == VIEW_X_NEG || face == VIEW_X_POS)
        {
            axis = 'X';
            aDimOut.xdim = aDim.zdim;
            aDimOut.zdim = aDim.xdim;
            reverse = (face == VIEW_X_NEG);
        }

        flip = flip^reverse;

        out.resize(aDimOut, false);

        MultidimArray<T1> imTemp;

        int index;

        for (size_t k = 0; k < aDimOut.zdim; k++)
        {
            imTemp.aliasSlice(out, k);
            index = k + (aDimOut.zdim - 1 - 2*k) * (int)flip;
            this->getSlice(index, imTemp, axis, !reverse);
        }
    }

    /** Reslice the current volume
     *
     * @param face Select the face to become the new Z direction
     * @param flip Invert the positions of Z planes, keeping the X-Y orientation
     * @param n Select the number of image in case of stacks
     */
    void reslice(AxisView face, bool flip = false, size_t n = 0)
    {
        MultidimArray<T> mTemp;
        reslice(mTemp, face, flip, n);
        *this = mTemp;
    }

    /** Get Column
     *
     * This function returns a column vector corresponding to the
     * choosen column.
     *
     * @code
     * std::vector< double > v;
     * m.getCol(-1, v);
     * @endcode
     */
    void getCol(size_t j, MultidimArray<T>& v) const
    {
        if (xdim == 0 || ydim == 0)
        {
            v.clear();
            return;
        }

        if (j < 0 || j >= xdim)
            REPORT_ERROR(ERR_INDEX_OUTOFBOUNDS,"getCol: Matrix subscript (j) greater than matrix dimension");

        v.resizeNoCopy(ydim);
        for (size_t i = 0; i < ydim; i++)
            DIRECT_A1D_ELEM(v,i) = DIRECT_A2D_ELEM(*this,i, j);
    }

    /** Set Column
     *
     * This function sets a column vector corresponding to the choosen column
     * inside matrix.
     *
     * @code
     * m.setCol(0, (m.row(1)).transpose()); // Copies row 1 in column 0
     * @endcode
     */
    void setCol(size_t j, const MultidimArray<T>& v)
    {
        if (xdim == 0 || ydim == 0)
            REPORT_ERROR(ERR_MULTIDIM_EMPTY, "setCol: Target matrix is empty");

        if (j < 0 || j>= xdim)
            REPORT_ERROR(ERR_INDEX_OUTOFBOUNDS, "setCol: Matrix subscript (j) out of range");

        if (v.xdim != ydim)
            REPORT_ERROR(ERR_MULTIDIM_SIZE,
                         "setCol: Vector dimension different from matrix one");

        for (size_t i = 0; i < ydim; i++)
            DIRECT_A2D_ELEM(*this,i, j) = DIRECT_A1D_ELEM(v,i);
    }

    /** Get row
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
    void getRow(size_t i, MultidimArray<T>& v) const
    {
        if (xdim == 0 || ydim == 0)
        {
            v.clear();
            return;
        }

        if (i < 0 || i >= ydim)
            REPORT_ERROR(ERR_INDEX_OUTOFBOUNDS, "getRow: Matrix subscript (i) greater than matrix dimension");

        v.resizeNoCopy(xdim);
        for (size_t j = 0; j < xdim; j++)
            DIRECT_A1D_ELEM(v,j) = DIRECT_A2D_ELEM(*this,i, j);
    }

    /** Set Row
     *
     * This function sets a row vector corresponding to the choosen row in the 2D Matrix
     *
     * @code
     * m.setRow(-2, m.row(1)); // Copies row 1 in row -2
     * @endcode
     */
    void setRow(int i, const MultidimArray<T>& v)
    {
        if (xdim == 0 || ydim == 0)
            REPORT_ERROR(ERR_MULTIDIM_EMPTY, "setRow: Target matrix is empty");

        if (i < 0 || i >= ydim)
            REPORT_ERROR(ERR_INDEX_OUTOFBOUNDS, "setRow: Matrix subscript (i) out of range");

        if (v.xdim != xdim)
            REPORT_ERROR(ERR_MULTIDIM_SIZE,
                         "setRow: Vector dimension different from matrix one");

        memcpy(&A2D_ELEM(*this,i,0),&A1D_ELEM(v,0),xdim*sizeof(double));
    }

    void getReal(MultidimArray< double > & realImg) const
    {
        REPORT_ERROR(ERR_TYPE_INCORRECT, "MultidimArray: Non complex datatype.");
    }

    void getImag(MultidimArray< double > & imagImg) const
    {
        REPORT_ERROR(ERR_TYPE_INCORRECT, "MultidimArray: Non complex datatype.");
    }


    /** 3D Logical to physical index translation.
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

    /** Interpolates the value of the nth 3D matrix M at the point (x,y,z).
    *
    * (x,y,z) are in logical coordinates.
    */
    T interpolatedElement3D(double x, double y, double z, T outside_value = (T) 0) const
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

        double doutside_value=outside_value;
        double d000 = (OUTSIDE3D(z0, y0, x0)) ? doutside_value : A3D_ELEM(*this, z0, y0, x0);
        double d001 = (OUTSIDE3D(z0, y0, x1)) ? doutside_value : A3D_ELEM(*this, z0, y0, x1);
        double d010 = (OUTSIDE3D(z0, y1, x0)) ? doutside_value : A3D_ELEM(*this, z0, y1, x0);
        double d011 = (OUTSIDE3D(z0, y1, x1)) ? doutside_value : A3D_ELEM(*this, z0, y1, x1);
        double d100 = (OUTSIDE3D(z1, y0, x0)) ? doutside_value : A3D_ELEM(*this, z1, y0, x0);
        double d101 = (OUTSIDE3D(z1, y0, x1)) ? doutside_value : A3D_ELEM(*this, z1, y0, x1);
        double d110 = (OUTSIDE3D(z1, y1, x0)) ? doutside_value : A3D_ELEM(*this, z1, y1, x0);
        double d111 = (OUTSIDE3D(z1, y1, x1)) ? doutside_value : A3D_ELEM(*this, z1, y1, x1);

        double dx00 = LIN_INTERP(fx, d000, d001);
        double dx01 = LIN_INTERP(fx, d100, d101);
        double dx10 = LIN_INTERP(fx, d010, d011);
        double dx11 = LIN_INTERP(fx, d110, d111);
        double dxy0 = LIN_INTERP(fy, dx00, dx10);
        double dxy1 = LIN_INTERP(fy, dx01, dx11);

        return (T) LIN_INTERP(fz, dxy0, dxy1);
    }

    /** Interpolates the value of the nth 2D matrix M at the point (x,y)
     *
     * Bilinear interpolation. (x,y) are in logical coordinates.
     */
    T interpolatedElement2D(double x, double y, T outside_value = (T) 0) const
    {
        int x0 = floor(x);
        double fx = x - x0;
        int x1 = x0 + 1;
        int y0 = floor(y);
        double fy = y - y0;
        int y1 = y0 + 1;

        int i0=STARTINGY(*this);
        int j0=STARTINGX(*this);
        int iF=FINISHINGY(*this);
        int jF=FINISHINGX(*this);

#define ASSIGNVAL(d,i,j) \
     if ((j) < j0 || (j) > jF || (i) < i0 || (i) > iF) \
      d=outside_value;\
        else \
         d=A2D_ELEM(*this, i, j);

        double d00, d10, d11, d01;
        ASSIGNVAL(d00,y0,x0);
        ASSIGNVAL(d01,y0,x1);
        ASSIGNVAL(d10,y1,x0);
        ASSIGNVAL(d11,y1,x1);

        double d0 = LIN_INTERP(fx, d00, d01);
        double d1 = LIN_INTERP(fx, d10, d11);
        return (T) LIN_INTERP(fy, d0, d1);
    }

    /** Interpolates the value of the nth 3D matrix M at the point (x,y,z) knowing
     * that this image is a set of B-spline coefficients.
     *
     * (x,y,z) are in logical coordinates.
     */
    T interpolatedElementBSpline3D(double x, double y, double z,
                                   int SplineDegree = 3) const
    {
        int SplineDegree_1 = SplineDegree - 1;

        // Logical to physical
        z -= STARTINGZ(*this);
        y -= STARTINGY(*this);
        x -= STARTINGX(*this);

        int l1 = (int)ceil(x - SplineDegree_1);
        int l2 = l1 + SplineDegree;

        int m1 = (int)ceil(y - SplineDegree_1);
        int m2 = m1 + SplineDegree;

        int n1 = (int)ceil(z - SplineDegree_1);
        int n2 = n1 + SplineDegree;

        double zyxsum = 0.0;
        double aux;
        int Xdim=(int)XSIZE(*this);
        int Ydim=(int)YSIZE(*this);
        int Zdim=(int)ZSIZE(*this);
        for (int nn = n1; nn <= n2; nn++)
        {
            int equivalent_nn=nn;
            if      (nn<0)
                equivalent_nn=-nn-1;
            else if (nn>=Zdim)
                equivalent_nn=2*Zdim-nn-1;
            double yxsum = 0.0;
            for (int m = m1; m <= m2; m++)
            {
                int equivalent_m=m;
                if      (m<0)
                    equivalent_m=-m-1;
                else if (m>=Ydim)
                    equivalent_m=2*Ydim-m-1;
                double xsum = 0.0;
                for (int l = l1; l <= l2; l++)
                {
                    double xminusl = x - (double) l;
                    int equivalent_l=l;
                    if      (l<0)
                        equivalent_l=-l-1;
                    else if (l>=Xdim)
                        equivalent_l=2*Xdim-l-1;
                    double Coeff = (double) DIRECT_A3D_ELEM(*this,
                                                            equivalent_nn,equivalent_m,equivalent_l);
                    switch (SplineDegree)
                    {
                    case 2:
                        xsum += Coeff * Bspline02(xminusl);
                        break;
                    case 3:
                        BSPLINE03(aux,xminusl);
                        xsum += Coeff * aux;
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
                    BSPLINE03(aux,yminusm);
                    yxsum += xsum * aux;
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

            double zminusn = z - (double) nn;
            switch (SplineDegree)
            {
            case 2:
                zyxsum += yxsum * Bspline02(zminusn);
                break;
            case 3:
                BSPLINE03(aux,zminusn);
                zyxsum += yxsum * aux;
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
     *
     * (x,y) are in logical coordinates
     *
     * To interpolate using splines you must first produce the Bspline
     * coefficients. An example to interpolate an image at (0.5,0.5) using
     * splines would be:
     *
     * @code
     * MultidimArray< double > Bspline_coeffs;
     * myImage.produceSplineCoefficients(Bspline_coeffs, 3);
     * interpolated_value = Bspline_coeffs.interpolatedElementBSpline(0.5,
     * 0.5,3);
     * @endcode
     */
    inline T interpolatedElementBSpline2D(double x, double y, int SplineDegree = 3) const
    {
        int SplineDegree_1 = SplineDegree - 1;

        // Logical to physical
        y -= STARTINGY(*this);
        x -= STARTINGX(*this);

        int l1 = (int)ceil(x - SplineDegree_1);
        int l2 = l1 + SplineDegree;
        int m1 = (int)ceil(y - SplineDegree_1);
        int m2 = m1 + SplineDegree;

        double columns = 0.0;
        double aux;
        int Ydim=(int)YSIZE(*this);
        int Xdim=(int)XSIZE(*this);
        for (int m = m1; m <= m2; m++)
        {
            int equivalent_m=m;
            if      (m<0)
                equivalent_m=-m-1;
            else if (m>=Ydim)
                equivalent_m=2*Ydim-m-1;
            double rows = 0.0;
            for (int l = l1; l <= l2; l++)
            {
                double xminusl = x - (double) l;
                int equivalent_l=l;
                if      (l<0)
                    equivalent_l=-l-1;
                else if (l>=Xdim)
                    equivalent_l=2*Xdim-l-1;
                double Coeff = DIRECT_A2D_ELEM(*this, equivalent_m,equivalent_l);
                switch (SplineDegree)
                {
                case 2:
                    rows += Coeff * Bspline02(xminusl);
                    break;

                case 3:
                    BSPLINE03(aux,xminusl);
                    rows += Coeff * aux;
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
                BSPLINE03(aux,yminusm);
                columns += rows * aux;
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
     *
     * (x) is in logical coordinates
     *
     * To interpolate using splines you must first produce the Bspline
     * coefficients. An example to interpolate a vector at (0.5) using
     * splines would be:
     *
     * @code
     * MultidimArray< double > Bspline_coeffs;
     * myVector.produceSplineCoefficients(Bspline_coeffs, 3);
     * interpolated_value = Bspline_coeffs.interpolatedElementBSpline(0.5,3);
     * @endcode
     */
    T interpolatedElementBSpline1D(double x, int SplineDegree = 3) const
    {
        int SplineDegree_1 = SplineDegree - 1;

        // Logical to physical
        x -= STARTINGX(*this);

        int l1 = (int)ceil(x - SplineDegree_1);
        int l2 = l1 + SplineDegree;
        int Xdim=(int)XSIZE(*this);
        double sum = 0.0;
        for (int l = l1; l <= l2; l++)
        {
            double xminusl = x - (double) l;
            int equivalent_l=l;
            if      (l<0)
                equivalent_l=-l-1;
            else if (l>=Xdim)
                equivalent_l=2*Xdim-l-1;
            double Coeff = (double) DIRECT_A1D_ELEM(*this, equivalent_l);
            double aux;
            switch (SplineDegree)
            {
            case 2:
                sum += Coeff * Bspline02(xminusl);
                break;

            case 3:
                BSPLINE03(aux,xminusl);
                sum += Coeff * aux;
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
    //@}

    /// @name Statistics functions
    //@{

    /** Print statistics in current line.
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
     *
     * The returned value is of the same type as the type of the array.
     */
    T computeMax() const
    {
        if (NZYXSIZE(*this) <= 0)
            return static_cast< T >(0);

        T maxval = data[0];

        T* ptr=NULL;
        size_t n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
        if (*ptr > maxval)
            maxval = *ptr;

        return maxval;
    }

    /** 1D Indices for the maximum element.
     *
     * This function just calls to the 4D function
     */
    void maxIndex(int& jmax) const
    {
        int zeroInt=0;
        maxIndex(zeroInt,zeroInt,zeroInt,jmax);
    }

    /** Minimum of the values in the array.
     *
     * The returned value is of the same type as the type of the array.
     */
    T computeMin() const
    {
        if (NZYXSIZE(*this) <= 0)
            return static_cast< T >(0);

        T minval = data[0];

        T* ptr=NULL;
        size_t n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
        if (*ptr < minval)
            minval = *ptr;

        return minval;
    }

    /** 4D Indices for the minimum element.
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
        size_t n=0;
        T minval = DIRECT_MULTIDIM_ELEM(*this, n);

        FOR_ALL_NZYX_ELEMENTS_IN_MULTIDIMARRAY(*this)
        {
            T val=DIRECT_MULTIDIM_ELEM(*this,n);
            if (val > minval)
            {
                minval = val;
                lmin = l;
                kmin = k;
                imin = i;
                jmin = j;
            }
            ++n;
        }
    }

    /** 3D Indices for the minimum element.
     *
     * This function just calls to the 4D function
     */
    void minIndex(int& kmin, int& imin, int& jmin) const
    {
        int zeroInt=0;
        minIndex(zeroInt,kmin,imin,jmin);
    }

    /** 2D Indices for the minimum element.
     *
     * This function just calls to the 4D function
     */
    void minIndex(int& imin, int& jmin) const
    {
        int zeroInt=0;
        minIndex(zeroInt,zeroInt,imin,jmin);
    }

    /** 1D Indices for the minimum element.
     *
     * This function just calls to the 4D function
     */
    void minIndex(int& jmin) const
    {
        int zeroInt=0;
        minIndex(zeroInt,zeroInt,zeroInt,jmin);
    }

    /** 4D Indices for the maximum element.
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
        size_t n=0;
        T maxval = DIRECT_MULTIDIM_ELEM(*this, n);

        FOR_ALL_NZYX_ELEMENTS_IN_MULTIDIMARRAY(*this)
        {
            T val=DIRECT_MULTIDIM_ELEM(*this, n);
            if (val > maxval)
            {
                maxval = val;
                lmax = l;
                kmax = k;
                imax = i;
                jmax = j;
            }
            ++n;
        }
    }

    /** Minimum and maximum of the values in the array.
     *
     * As doubles.
     */
    void computeDoubleMinMax(double& minval, double& maxval) const
    {
        if (NZYXSIZE(*this) <= 0)
            return;

        T* ptr=NULL;
        size_t n;
        T Tmin=DIRECT_MULTIDIM_ELEM(*this,0);
        T Tmax=Tmin;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
        {
            T val=*ptr;
            if (val < Tmin)
                Tmin = val;
            else if (val > Tmax)
                Tmax = val;
        }
        minval = static_cast< double >(Tmin);
        maxval = static_cast< double >(Tmax);
    }

    /** Minimum and maximum of the values in the array.
     *
     * As doubles.
     */
    void computeDoubleMinMaxRange(double& minval, double& maxval,size_t offset, size_t size) const
    {
        if (NZYXSIZE(*this) <= 0)
            return;

        minval = maxval = static_cast< double >(data[offset]);

        T* ptr=NULL;
        T val;
        size_t n;

        for (n=offset,ptr=data+offset; n<size; ++n, ++ptr)
        {
            val = *ptr;
            if (val < minval)
                minval = static_cast< double >(val);
            else if (val > maxval)
                maxval = static_cast< double >(val);
        }
    }

    /** Average of the values in the array.
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
        size_t n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
        sum += static_cast< double >(*ptr);

        return sum / NZYXSIZE(*this);
    }

    /** Standard deviation of the values in the array.
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
        size_t n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
        {
            double val=static_cast< double >(*ptr);
            avg += val;
            stddev += val * val;
        }

        avg /= NZYXSIZE(*this);
        stddev = stddev / NZYXSIZE(*this) - avg * avg;
        stddev *= NZYXSIZE(*this) / (NZYXSIZE(*this) - 1);

        // Foreseeing numerical instabilities
        stddev = sqrt(static_cast<double>((ABS(stddev))));

        return stddev;
    }

    /** Compute statistics.
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
        size_t n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
        {
            T Tval=*ptr;
            double val=Tval;
            avg += val;
            stddev += val * val;

            if (Tval > maxval)
                maxval = Tval;
            else if (Tval < minval)
                minval = Tval;
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

    /** Compute statistics.
     *
     * The average, standard deviation, minimum and maximum value are
     * returned.
     */
    void computeAvgStdev(double& avg, double& stddev) const
    {
        if (NZYXSIZE(*this) <= 0)
            return;

        avg = 0;
        stddev = 0;

        T* ptr=NULL;
        size_t n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
        {
            T Tval=*ptr;
            double val=Tval;
            avg += val;
            stddev += val * val;
        }

        avg /= NZYXSIZE(*this);

        if (NZYXSIZE(*this) > 1)
        {
            stddev = stddev / NZYXSIZE(*this) - avg * avg;
            stddev *= NZYXSIZE(*this) / (NZYXSIZE(*this) - 1);

            // Foreseeing numerical instabilities
            stddev = sqrt(fabs(stddev));
        }
        else
            stddev = 0;
    }

    /** Compute statistics in the active area
     *
     * Only the statistics for values in the overlapping between the mask and the
     * volume for those the mask is not 0 are computed.
     */
    void computeAvgStdev_within_binary_mask(const MultidimArray< int >& mask,
                                            double& avg, double& stddev) const
    {
        SPEED_UP_tempsInt;
        double sum1 = 0;
        double sum2 = 0;
        int N = 0;

        FOR_ALL_ELEMENTS_IN_COMMON_IN_ARRAY3D(mask, *this)
        {
            if (A3D_ELEM(mask, k, i, j) != 0)
            {
                ++N;
                double aux=A3D_ELEM(*this, k, i, j);
                sum1 += aux;
                sum2 += aux*aux;
            }
        }

        // average and standard deviation
        avg  = sum1 / (double) N;
        if (N > 1)
            stddev = sqrt(fabs(sum2 / N - avg * avg) * N / (N - 1));
        else
            stddev = 0;
    }

    /** Compute statistics within 2D region of 2D image.
     *
     * The 2D region is specified by two corners.
     * Note that this function only works for the 0th image in a multi-image array...
     */
    void computeStats(double& avg,
                      double& stddev,
                      T& min_val,
                      T& max_val,
                      Matrix1D< int >& corner1,
                      Matrix1D< int >& corner2,
                      size_t n = 0)
    {
        (*this).checkDimension(2);
        min_val = max_val = NZYX_ELEM((*this), n, 0, YY(corner1), XX(corner1));

        Matrix1D< double > r(3);
        double N = 0, sum = 0, sum2 = 0;

        FOR_ALL_ELEMENTS_IN_ARRAY2D_BETWEEN(corner1, corner2)
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

    /** Median
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
        MultidimArray< double > temp(*this);

        // Sort indexes
        double* temp_array = MULTIDIM_ARRAY(temp)-1;
        qcksrt(NZYXSIZE(*this), temp_array); //FIXME: Valgrind:: Invalid read of size 8: qcksrt(int, double*) (numerical_recipes.cpp:266)


        // Get median
        if (NZYXSIZE(*this)%2==0)
            return 0.5*(DIRECT_MULTIDIM_ELEM(temp,NZYXSIZE(*this)/2-1)+
                        DIRECT_MULTIDIM_ELEM(temp,NZYXSIZE(*this)/2  ));
        else
            return DIRECT_MULTIDIM_ELEM(temp,NZYXSIZE(*this)/2);
    }

    /** Adjust the range of the array to a given one.
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
        size_t n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
        *ptr = minF + static_cast< T >(slope *
                                       static_cast< double >(*ptr - min0));
    }

    /** Adjust the range of the array to a given one within a mask.
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
        size_t n;
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
     *
     * A linear operation is performed on the values of the array such that
     * after it, the values of the self array are as similar as possible
     * (L2 sense) to the values of the array shown as sample
     */

    //As written this will only work for T=double
    //nevertheless since this is used is better
    //to use T than double or will create problem for int multidim arrays
    void rangeAdjust(const MultidimArray<T> &example,
                     const MultidimArray<int> *mask=NULL)
    {
        if (NZYXSIZE(*this) <= 0)
            return;

        double avgExample, stddevExample, avgThis, stddevThis;
        if (mask!=NULL)
        {
            example.computeAvgStdev_within_binary_mask(*mask,avgExample,stddevExample);
            computeAvgStdev_within_binary_mask(*mask,avgThis,stddevThis);
        }
        else
        {
            computeAvgStdev(avgThis,stddevThis);
            example.computeAvgStdev(avgExample,stddevExample);
        }

        // y=a+bx
        double b=stddevThis>0? stddevExample/stddevThis:0;
        double a=avgExample-avgThis*b;

        size_t n;
        T *ptr=NULL;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
        *ptr = static_cast< T >(a+b * static_cast< double > (*ptr));
    }

    /** Adjust the average and stddev of the array to given values.
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
        double avg0=0.0, stddev0=0.0;
        double a, b;

        if (NZYXSIZE(*this) == 0)
            return;

        computeAvgStdev(avg0, stddev0);

        if (stddev0 != 0)
            a = stddevF / stddev0;
        else
            a = 0;

        b = avgF - a * avg0;

        T* ptr=&DIRECT_MULTIDIM_ELEM(*this,0);
        size_t nmax=(nzyxdim/4)*4;
        for (size_t n=0; n<nmax; n+=4, ptr+=4)
        {
            *(ptr  )= static_cast< T >(a * (*(ptr  )) + b);
            *(ptr+1)= static_cast< T >(a * (*(ptr+1)) + b);
            *(ptr+2)= static_cast< T >(a * (*(ptr+2)) + b);
            *(ptr+3)= static_cast< T >(a * (*(ptr+3)) + b);
        }
        for (size_t n=nmax; n<nzyxdim; ++n, ptr+=1)
            *(ptr  )= static_cast< T >(a * (*(ptr  )) + b);
    }
    //@}

    /** @name Array "by" array operations.
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
    //@{

    /** Core array by array operation.
     *
     * It assumes that the result is already resized.
     */
    inline friend void coreArrayByArray(const MultidimArray<T>& op1,
                                        const MultidimArray<T>& op2, MultidimArray<T>& result,
                                        char operation)
    {
        T* ptrResult=NULL;
        T* ptrOp1=NULL;
        T* ptrOp2=NULL;
        size_t n;
        // Loop unrolling
        const size_t unroll=4;
        size_t nmax=unroll*(op1.zyxdim/unroll);
        switch (operation)
        {
        case '+':
            for (n=0, ptrResult=result.data, ptrOp1=op1.data,ptrOp2=op2.data;
                 n<nmax; n+=unroll, ptrResult+=unroll, ptrOp1+=unroll, ptrOp2+=unroll)
            {
                *ptrResult = *ptrOp1 + *ptrOp2;
                *(ptrResult+1) = *(ptrOp1+1) + *(ptrOp2+1);
                *(ptrResult+2) = *(ptrOp1+2) + *(ptrOp2+2);
                *(ptrResult+3) = *(ptrOp1+3) + *(ptrOp2+3);
            }
            for (n=nmax, ptrResult=result.data+nmax, ptrOp1=op1.data+nmax, ptrOp2=op2.data+nmax;
                 n<op1.zyxdim; ++n, ++ptrResult, ++ptrOp1, ++ptrOp2)
                *ptrResult = *ptrOp1 + *ptrOp2;
            break;
        case '-':
                for (n=0, ptrResult=result.data, ptrOp1=op1.data,ptrOp2=op2.data;
                     n<nmax; n+=unroll, ptrResult+=unroll, ptrOp1+=unroll, ptrOp2+=unroll)
                {
                    *ptrResult = *ptrOp1 - *ptrOp2;
                    *(ptrResult+1) = *(ptrOp1+1) - *(ptrOp2+1);
                    *(ptrResult+2) = *(ptrOp1+2) - *(ptrOp2+2);
                    *(ptrResult+3) = *(ptrOp1+3) - *(ptrOp2+3);
                }
            for (n=nmax, ptrResult=result.data+nmax, ptrOp1=op1.data+nmax, ptrOp2=op2.data+nmax;
                 n<op1.zyxdim; ++n, ++ptrResult, ++ptrOp1, ++ptrOp2)
                *ptrResult = *ptrOp1 - *ptrOp2;
            break;
        case '*':
                for (n=0, ptrResult=result.data, ptrOp1=op1.data,ptrOp2=op2.data;
                     n<nmax; n+=unroll, ptrResult+=unroll, ptrOp1+=unroll, ptrOp2+=unroll)
                {
                    *ptrResult = *ptrOp1 * *ptrOp2;
                    *(ptrResult+1) = *(ptrOp1+1) * *(ptrOp2+1);
                    *(ptrResult+2) = *(ptrOp1+2) * *(ptrOp2+2);
                    *(ptrResult+3) = *(ptrOp1+3) * *(ptrOp2+3);
                }
            for (n=nmax, ptrResult=result.data+nmax, ptrOp1=op1.data+nmax, ptrOp2=op2.data+nmax;
                 n<op1.zyxdim; ++n, ++ptrResult, ++ptrOp1, ++ptrOp2)
                *ptrResult = *ptrOp1 * *ptrOp2;
            break;
        case '/':
                for (n=0, ptrResult=result.data, ptrOp1=op1.data,ptrOp2=op2.data;
                     n<nmax; n+=unroll, ptrResult+=unroll, ptrOp1+=unroll, ptrOp2+=unroll)
                {
                    *ptrResult = *ptrOp1 / *ptrOp2;
                    *(ptrResult+1) = *(ptrOp1+1) / *(ptrOp2+1);
                    *(ptrResult+2) = *(ptrOp1+2) / *(ptrOp2+2);
                    *(ptrResult+3) = *(ptrOp1+3) / *(ptrOp2+3);
                }
            for (n=nmax, ptrResult=result.data+nmax, ptrOp1=op1.data+nmax, ptrOp2=op2.data+nmax;
                 n<op1.zyxdim; ++n, ++ptrResult, ++ptrOp1, ++ptrOp2)
                *ptrResult = *ptrOp1 / *ptrOp2;
            break;
        }
    }

    /** Array by array
     *
     * This function must take two vectors of the same size, and operate element
     * by element according to the operation required. This is the function
     * which really implements the operations. Simple calls to it perform much
     * faster than calls to the corresponding operators. Although it is supposed
     * to be a hidden function not useable by normal programmers.
     *
     */
    inline friend void arrayByArray(const MultidimArray<T>& op1,
                                    const MultidimArray<T>& op2, MultidimArray<T>& result,
                                    char operation)
    {
        if (!op1.sameShape(op2))
            REPORT_ERROR(ERR_MULTIDIM_SIZE,
                         formatString("Array_by_array: different shapes (%c)", operation));
        if (result.data == NULL || !result.sameShape(op1))
            result.resizeNoCopy(op1);
        coreArrayByArray(op1, op2, result, operation);
    }

    /** v3 = v1 + v2.
     */
    MultidimArray<T> operator+(const MultidimArray<T>& op1) const
    {
        MultidimArray<T> tmp;
        arrayByArray(*this, op1, tmp, '+');
        return tmp;
    }

    /** v3 = v1 - v2.
     */
    MultidimArray<T> operator-(const MultidimArray<T>& op1) const
    {
        MultidimArray<T> tmp;
        arrayByArray(*this, op1, tmp, '-');
        return tmp;
    }

    /** v3 = v1 * v2.
     */
    MultidimArray<T> operator*(const MultidimArray<T>& op1) const
    {
        MultidimArray<T> tmp;
        arrayByArray(*this, op1, tmp, '*');
        return tmp;
    }

    /** v3 = v1 / v2.
     */
    MultidimArray<T> operator/(const MultidimArray<T>& op1) const
    {
        MultidimArray<T> tmp;
        arrayByArray(*this, op1, tmp, '/');
        return tmp;
    }

    /** v3 += v2.
     */
    void operator+=(const MultidimArray<T>& op1)
    {
        arrayByArray(*this, op1, *this, '+');
    }

    /** v3 -= v2.
     */
    void operator-=(const MultidimArray<T>& op1)
    {
        arrayByArray(*this, op1, *this, '-');
    }

    /** v3 *= v2.
     */
    void operator*=(const MultidimArray<T>& op1)
    {
        arrayByArray(*this, op1, *this, '*');
    }

    /** v3 /= v2.
     */
    void operator/=(const MultidimArray<T>& op1)
    {
        arrayByArray(*this, op1, *this, '/');
    }

    /** Dot product */
    double dotProduct(const MultidimArray<T>& op1)
    {
        if (!sameShape(op1))
            REPORT_ERROR(ERR_MULTIDIM_SIZE,"The two arrays for dot product are not of the same shape");
        double dot=0;
        size_t n;
        T* ptrOp1=NULL;
        T* ptrOp2=NULL;
        const size_t unroll=4;
        size_t nmax=unroll*(MULTIDIM_SIZE(*this)/unroll);
        for (n=0, ptrOp1=data,ptrOp2=op1.data;
             n<nmax; n+=unroll, ptrOp1+=unroll, ptrOp2+=unroll)
        {
            dot += *ptrOp1 * *ptrOp2;
            dot += *(ptrOp1+1) * *(ptrOp2+1);
            dot += *(ptrOp1+2) * *(ptrOp2+2);
            dot += *(ptrOp1+3) * *(ptrOp2+3);
        }
        for (n=nmax, ptrOp1=data+nmax, ptrOp2=op1.data+nmax;
             n<MULTIDIM_SIZE(*this); ++n, ++ptrOp1, ++ptrOp2)
            dot += *ptrOp1 * *ptrOp2;
        return dot;
    }
    //@}

    /** @name Array "by" scalar operations
     *
     * These operations are between an array and a scalar (of the same type as
     * the array). The result must have been defined to be of the same type as
     * the operands.
     *
     * In this kind of operations each element of array 1 is operated with the
     * given constant. The result has also got the same shape as the input
     * array and its former content is lost
     */
    //@{

    /** Core array by scalar operation.
     *
     * It assumes that the result is already resized.
     *
     * This function is not ported to Python.
     */
    inline friend void coreArrayByScalar(const MultidimArray<T>& op1,
                                         const T& op2,
                                         MultidimArray<T>& result,
                                         char operation)
    {
        T* ptrResult=NULL;
        T* ptrOp1=NULL;
        size_t n;
        switch (operation)
        {
        case '+':
            for (n=0, ptrResult=result.data, ptrOp1=op1.data;
                 n<op1.zyxdim; ++n, ++ptrResult, ++ptrOp1)
                *ptrResult = *ptrOp1 + op2;
            break;
        case '-':
                for (n=0, ptrResult=result.data, ptrOp1=op1.data;
                     n<op1.zyxdim; ++n, ++ptrResult, ++ptrOp1)
                    *ptrResult = *ptrOp1 - op2;
            break;
        case '*':
                for (n=0, ptrResult=result.data, ptrOp1=op1.data;
                     n<op1.zyxdim; ++n, ++ptrResult, ++ptrOp1)
                    *ptrResult = *ptrOp1 * op2;
            break;
        case '/':
                for (n=0, ptrResult=result.data, ptrOp1=op1.data;
                     n<op1.zyxdim; ++n, ++ptrResult, ++ptrOp1)
                    *ptrResult = *ptrOp1 / op2;
            break;
        case '=':
                for (n=0, ptrResult=result.data, ptrOp1=op1.data;
                     n<op1.zyxdim; ++n, ++ptrResult, ++ptrOp1)
                    *ptrResult = *ptrOp1 == op2;
            break;
        }
    }

    /** Array by scalar.
     *
     * This function must take one vector and a constant, and operate element
     * by element according to the operation required. This is the function
     * which really implements the operations. Simple calls to it perform much
     * faster than calls to the corresponding operators. Although it is
     * supposed to be a hidden function not useable by normal programmers.
     *
     * This function is not ported to Python.
     */
    inline friend void arrayByScalar(const MultidimArray<T>& op1,
                                     T op2,
                                     MultidimArray<T>& result,
                                     char operation)
    {
        if (result.data == NULL || !result.sameShape(op1))
            result.resizeNoCopy(op1);
        coreArrayByScalar(op1, op2, result, operation);
    }

    /** v3 = v1 + k.
     */
    MultidimArray<T> operator+(T op1) const
    {
        MultidimArray<T> tmp;
        arrayByScalar(*this, op1, tmp, '+');
        return tmp;
    }

    /** v3 = v1 - k.
     */
    MultidimArray<T> operator-(T op1) const
    {
        MultidimArray<T> tmp;
        arrayByScalar(*this, op1, tmp, '-');
        return tmp;
    }

    /** v3 = v1 * k.
     */
    MultidimArray<T> operator*(T op1) const
    {
        MultidimArray<T> tmp;
        arrayByScalar(*this, op1, tmp, '*');
        return tmp;
    }

    /** v3 = v1 / k.
     */
    MultidimArray<T> operator/(T op1) const
    {
        MultidimArray<T> tmp;
        arrayByScalar(*this, op1, tmp, '/');
        return tmp;
    }


    /** v3 = (v1 == k).
     */
    void equal(T op1, MultidimArray<char> &result) const
    {
        result.resizeNoCopy(*this);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(result)
        DIRECT_MULTIDIM_ELEM(result,n) = DIRECT_MULTIDIM_ELEM(*this,n) == op1;
    }

    /** v3 += k.
     *
     * This function is not ported to Python.
     */
    void operator+=(const T& op1)
    {
        arrayByScalar(*this, op1, *this, '+');
    }

    /** v3 -= k.
     *
     * This function is not ported to Python.
     */
    void operator-=(const T& op1)
    {
        arrayByScalar(*this, op1, *this, '-');
    }

    /** v3 *= k.
     *
     * This function is not ported to Python.
     */
    void operator*=(const T& op1)
    {
        arrayByScalar(*this, op1, *this, '*');
    }

    /** v3 /= k.
     *
     * This function is not ported to Python.
     */
    void operator/=(const T& op1)
    {
        arrayByScalar(*this, op1, *this, '/');
    }
    //@}

    /** @name Scalar "by" array operations
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
    //@{

    /** Core array by scalar operation.
     *
     * It assumes that the result is already resized.
     *
     * This function is not ported to Python.
     */
    inline friend void coreScalarByArray(const T& op1,
                                         const MultidimArray<T>& op2,
                                         MultidimArray<T>& result,
                                         char operation)
    {
        T* ptrResult=NULL;
        T* ptrOp2=NULL;
        size_t n;
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

    /** Scalar by array.
     *
     * This function must take one scalar and a vector, and operate element by
     * element according to the operation required. This is the function which
     * really implements the operations. Simple calls to it perform much faster
     * than calls to the corresponding operators. Although it is supposed to
     * be a hidden function not useable by normal programmers.
     *
     * This function is not ported to Python.
     */
    inline friend void scalarByArray(T op1,
                                     const MultidimArray<T>& op2,
                                     MultidimArray<T>& result,
                                     char operation)
    {
        if (result.data == NULL || !result.sameShape(op2))
            result.resizeNoCopy(op2);
        coreScalarByArray(op1, op2, result, operation);
    }

    /** v3 = k + v2.
     */
    friend MultidimArray<T> operator+(T op1, const MultidimArray<T>& op2)
    {
        MultidimArray<T> tmp;
        scalarByArray(op1, op2, tmp, '+');
        return tmp;
    }

    /** v3 = k - v2.
     */
    friend MultidimArray<T> operator-(T op1, const MultidimArray<T>& op2)
    {
        MultidimArray<T> tmp;
        scalarByArray(op1, op2, tmp, '-');
        return tmp;
    }

    /** v3 = k * v2.
     */
    friend MultidimArray<T> operator*(T op1, const MultidimArray<T>& op2)
    {
        MultidimArray<T> tmp;
        scalarByArray(op1, op2, tmp, '*');
        return tmp;
    }

    /** v3 = k / v2
     */
    friend MultidimArray<T> operator/(T op1, const MultidimArray<T>& op2)
    {
        MultidimArray<T> tmp;
        scalarByArray(op1, op2, tmp, '/');
        return tmp;
    }
    //@}

    /// @name Initialization
    /// @{

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
        T* ptr=NULL;
        size_t n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
        *ptr = val;
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
    template <typename T1>
    void initZeros(const MultidimArray<T1>& op)
    {
        if (data == NULL || !sameShape(op))
            resizeNoCopy(op);
        memset(data,0,nzyxdim*sizeof(T));
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
    inline void initZeros()
    {
        memset(data,0,nzyxdim*sizeof(T));
    }

    /** Initialize to zeros with a given size.
     */
    inline void initZeros(size_t Ndim, size_t Zdim, size_t Ydim, size_t Xdim)
    {
        if (xdim!=Xdim || ydim!=Ydim || zdim!=Zdim || ndim!=Ndim)
            resize(Ndim, Zdim,Ydim,Xdim,false);
        memset(data,0,nzyxdim*sizeof(T));
    }

    /** Initialize to zeros with a given size.
     */
    void initZeros(int Xdim)
    {
        initZeros(1, 1, 1, Xdim);
    }

    /** Initialize to zeros with a given size.
     */
    void initZeros(int Ydim, int Xdim)
    {
        initZeros(1, 1, Ydim, Xdim);
    }

    /** Initialize to zeros with a given size.
     */
    void initZeros(int Zdim, int Ydim, int Xdim)
    {
        initZeros(1, Zdim, Ydim, Xdim);
    }

    /** Linear initialization (only for 1D)
     *
     * The 1D vector is filled with values increasing/decreasing linearly within a
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
    void initLinear(T minF, T maxF, int n = 1, const String& mode = "incr")
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
            REPORT_ERROR(ERR_VALUE_INCORRECT, "Init_linear: Mode not supported (" + mode +
                         ")");

        if (steps == 0)
            clear();
        else
        {
            resizeNoCopy(steps);
            for (int i = 0; i < steps; i++)
                A1D_ELEM(*this, i) = (T)((double) minF + slope * i);
        }
    }

    /** Initialize with random values.
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
    void initRandom(double op1, double op2, RandomMode mode = RND_UNIFORM)
    {
        T* ptr=NULL;
        size_t n;
        if (mode == RND_UNIFORM)
            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
            *ptr = static_cast< T >(rnd_unif(op1, op2));
        else if (mode == RND_GAUSSIAN)
            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
            *ptr = static_cast< T >(rnd_gaus(op1, op2));
        else
            REPORT_ERROR(ERR_VALUE_INCORRECT,
                         formatString("InitRandom: Mode not supported"));
    }

    /** Add noise to actual values.
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
                  const String& mode = "uniform",
                  double df = 3.) const
    {
        T* ptr=NULL;
        size_t n;
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
            REPORT_ERROR(ERR_VALUE_INCORRECT,
                         formatString("AddNoise: Mode not supported (%s)", mode.c_str()));
    }
    //@}

    /** @name Utilities
     *
     * Here you have several easy functions to manage the values of
     * the array.
     */
    //@{

    /** Produce a 3D array suitable for working with Numerical Recipes.
     *
     * This function must be used only as a preparation for routines which need
     * that the first physical index is 1 and not 0 as it usually is in C. New
     * memory is needed to hold the new double pointer array.
     */
    T*** adaptForNumericalRecipes3D(size_t n = 0) const
    {
        T*** m = NULL;
        ask_Tvolume(m, 1, ZSIZE(*this), 1, YSIZE(*this), 1, XSIZE(*this));

        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(*this)
        m[k+1][i+1][j+1] = DIRECT_NZYX_ELEM(*this, n, k, i, j);

        return m;
    }

    /** Kill a 3D array produced for numerical recipes.
     */
    void killAdaptationForNumericalRecipes3D(T*** m) const
    {
        free_Tvolume(m, 1, ZSIZE(*this), 1, YSIZE(*this), 1, XSIZE(*this));
    }

    /** Produce a 2D array suitable for working with Numerical Recipes
     *
     * This function must be used only as a preparation for routines which need
     * that the first physical index is 1 and not 0 as it usually is in C. New
     * memory is needed to hold the new double pointer array.
     */
    T** adaptForNumericalRecipes2D(size_t n = 0) const
    {
        T** m = NULL;
        ask_Tmatrix(m, 1, YSIZE(*this), 1, XSIZE(*this));

        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(*this)
        m[i+1][j+1] = DIRECT_NZYX_ELEM(*this, n, 0, i, j);

        return m;
    }

    /** Produce a 1D pointer suitable for working with Numerical Recipes (2)
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
     */
    void loadFromNumericalRecipes2D(T** m, int Ydim, int Xdim)
    {
        resizeNoCopy(Ydim, Xdim);

        for (int i = 1; i <= Ydim; i++)
            for (int j = 1; j <= Xdim; j++)
                DIRECT_A2D_ELEM(*this,i - 1, j - 1) = m[i][j];
    }

    /** Kill a 2D array produced for numerical recipes
     *
     * The allocated memory is freed.
     */
    void killAdaptationForNumericalRecipes2D(T** m) const
    {
        free_Tmatrix(m, 1, YSIZE(*this), 1, XSIZE(*this));
    }

    /** Kill a 2D array produced for numerical recipes, 2.
     *
     * Nothing needs to be done.
     */
    void killAdaptationForNumericalRecipes22D(T** m) const
        {}

    /** Produce a 1D array suitable for working with Numerical Recipes
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
     *
     * Nothing needs to be done in fact.
     *
     * This function is not ported to Python.
     */
    void killAdaptationForNumericalRecipes1D(T* m) const
        {}

    /** Computes the center of mass of the nth array
     */
    void centerOfMass(Matrix1D< double >& center, void * mask=NULL, size_t n = 0)
    {
        center.initZeros(3);
        double mass = 0;
        MultidimArray< int >* imask = (MultidimArray< int >*) mask;

        FOR_ALL_ELEMENTS_IN_ARRAY3D(*this)
        {
            if ((imask == NULL || NZYX_ELEM(*imask, n, k, i, j)) &&
                A3D_ELEM(*this, k, i, j) > 0)
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

    /** Sort 1D vector elements
     *
     * Sort in ascending order the vector elements. You can use the "reverse"
     * function to sort in descending order.
     *
     * @code
     * v2 = v1.sort();
     * @endcode
     */
    void sort(MultidimArray<T> &result) const
    {
        checkDimension(1);

        MultidimArray<T> temp;
        MultidimArray< double > aux;

        if (xdim == 0)
        {
            result.clear();
            return;
        }

        // Initialise data
        typeCast(*this, aux);

        // Sort
        double * aux_array = aux.adaptForNumericalRecipes1D();
        qcksrt(xdim, aux_array);

        typeCast(aux, result);
    }

    /** Gives a vector with the indexes for a sorted vector
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
    void indexSort(MultidimArray< int > &indx) const
    {
        checkDimension(1);

        MultidimArray< double > temp;
        indx.clear();

        if (xdim == 0)
            return;

        if (xdim == 1)
        {
            indx.resizeNoCopy(1);
            DIRECT_A1D_ELEM(indx,0) = 1;
            return;
        }

        // Initialise data
        indx.resizeNoCopy(xdim);
        typeCast(*this, temp);

        // Sort indexes
        double* temp_array = temp.adaptForNumericalRecipes1D();
        int* indx_array = indx.adaptForNumericalRecipes1D();
        indexx(XSIZE(*this), temp_array, indx_array);
    }

    /** Cumulative Density Function.
     * For each entry in the array, give what is the probability of having a value smaller or equal than this entry.*/
    void cumlativeDensityFunction(MultidimArray<double> &cdf)
    {
        MultidimArray<int> indx;
        indexSort(indx);
        cdf.resizeNoCopy(*this);
        double iSize=1.0/XSIZE(indx);
        FOR_ALL_ELEMENTS_IN_ARRAY1D(indx)
        {
            int ii=A1D_ELEM(indx,i)-1;
            A1D_ELEM(cdf,ii)=(i+1)*iSize;
        }
    }

    /** Several thresholding.
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
    void threshold(const String& type,
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
            REPORT_ERROR(ERR_VALUE_INCORRECT,
                         formatString("Threshold: mode not supported (%s)", type.c_str() ));

        T* ptr=NULL;
        size_t n;
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
     *
     * This function returns the number of elements meeting the threshold
     * condition.
     */
    size_t countThreshold(const String& type,
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
            REPORT_ERROR(ERR_VALUE_INCORRECT,
                         formatString("CountThreshold: mode not supported (%s)", type.c_str()));

        size_t ret = 0;

        T* ptr=NULL;
        size_t n;
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
        size_t n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
        if (mask == NULL || DIRECT_MULTIDIM_ELEM(*mask,n) > 0 )
            if (ABS(*ptr - oldv) <= accuracy)
                *ptr = newv;
    }

    /** Substitute a given value by a sample from a Gaussian distribution.
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
        size_t n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
        if (mask == NULL || DIRECT_MULTIDIM_ELEM(*mask,n) > 0 )
            if (ABS(*ptr - oldv) <= accuracy)
                *ptr = rnd_gaus(avgv, sigv);
    }

    /** Binarize.
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
        size_t n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
        if (mask == NULL || DIRECT_MULTIDIM_ELEM(*mask,n) > 0 )
        {
            if (*ptr <= val + accuracy)
                *ptr = 0;
            else
                *ptr = 1;
        }
    }

    /** Binarize using a range
     *
     * This functions substitutes all values in a volume which are in the range between valMin
     * and valMax by 1 and the rest are set to 0.
     */
    void binarizeRange(double valMin = 0, double valMax = 255,
                       MultidimArray<int> * mask = NULL )
    {
        T* ptr=NULL;
        size_t n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
        if (mask == NULL || DIRECT_MULTIDIM_ELEM(*mask,n) > 0 )
        {
            if ( (*ptr < valMax) &&  (*ptr > valMin) )
                *ptr = 1;
            else
                *ptr = 0;
        }
    }

    /** ROUND
     *
     * Applies a ROUND (look for the nearest integer) to each array element.
     */
    void selfROUND()
    {
        T* ptr=NULL;
        size_t n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
        *ptr = ROUND(*ptr);
    }

    /** CEILING
     *
     * Applies a CEILING (look for the nearest larger integer) to each
     * array element.
     */
    void selfCEIL()
    {
        T* ptr=NULL;
        size_t n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
        *ptr = CEIL(*ptr);
    }

    /** FLOOR
     *
     * Applies a FLOOR (look for the nearest larger integer) to each
     * array element.
     */
    void selfFLOOR()
    {
        T* ptr=NULL;
        size_t n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
        *ptr = FLOOR(*ptr);
    }

    /** ABS
     *
     * Applies an ABS (absolute value) to each array element.
     */
    void selfABS()
    {
        T* ptr=NULL;
        size_t n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
        *ptr = ABS(*ptr);
    }

    /** MAX
     *
     * Each component of the result is the maximum of the correspoing
     * components of the two input arrays. They must have the same shape, if
     * not an exception is thrown
     */
    friend void MAX(const MultidimArray<T>& v1, const MultidimArray<T>& v2,
                    MultidimArray<T>& result)
    {
        if (!v1.sameShape(v2))
            REPORT_ERROR(ERR_MULTIDIM_SIZE, "MAX: arrays of different shape");

        result.resizeNoCopy(v1);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(result)
        DIRECT_MULTIDIM_ELEM(result,n) = XMIPP_MAX(
                                             DIRECT_MULTIDIM_ELEM(v1,n),
                                             DIRECT_MULTIDIM_ELEM(v2,n));
    }

    /** MIN
     *
     * Each component of the result is the minimum of the correspoing
     * components of the two input arrays. They must have the same shape, if
     * not an exception is thrown
     */
    friend void MIN(const MultidimArray<T>& v1, const MultidimArray<T>& v2,
                    MultidimArray<T>& result)
    {
        if (!v1.sameShape(v2))
            REPORT_ERROR(ERR_MULTIDIM_SIZE, "MIN: arrays of different shape");

        result.resizeNoCopy(v1);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(result)
        DIRECT_MULTIDIM_ELEM(result,n) = XMIPP_MIN(
                                             DIRECT_MULTIDIM_ELEM(v1,n),
                                             DIRECT_MULTIDIM_ELEM(v2,n));
    }

    /** Sqrt.
     *
     * Each component of the result is the square root of the original
     * component.
     */
    void selfSQRT()
    {
        T* ptr=NULL;
        size_t n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
        *ptr = static_cast< T >(sqrt(static_cast< double >(*ptr)));
    }

    /** Sum of matrix values.
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
        size_t n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
        sum += *ptr;
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

        // Unroll the loop
        const size_t unroll=4;
        size_t nmax=(NZYXSIZE(*this)/unroll)*unroll;
        T* ptr = MULTIDIM_ARRAY(*this);
        for (size_t n=0; n<nmax; n+=unroll, ptr+=unroll)
        {
            sum+= (*ptr)*(*ptr);
            T* ptr1=ptr+1;
            sum+= (*ptr1)*(*ptr1);
            T* ptr2=ptr+2;
            sum+= (*ptr2)*(*ptr2);
            T* ptr3=ptr+3;
            sum+= (*ptr3)*(*ptr3);
        }
        // Do the remaining elements
        for (size_t n=nmax; n<NZYXSIZE(*this); ++n, ++ptr)
            sum+=(*ptr)*(*ptr);
        return sum;
    }

    /** Log10.
     *
     * Each component of the result is the log10 of the original components.
     */
    void selfLog10()
    {
        T* ptr=NULL;
        size_t n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
        *ptr = static_cast< T >(log10(static_cast< double >(*ptr)));
    }

    /** Log.
     *
     * Each component of the result is the log of the original components.
     */
    void selfLog()
    {
        T* ptr=NULL;
        size_t n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this,n,ptr)
        *ptr = static_cast< T >(log(static_cast< double >(*ptr)));
    }

    /** Reverse matrix values over X axis, keep in this object.
     *
     * Maybe better with an example:
     *
     *Odd case
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
     *Even case
     * @code
     * slice 0
     * [01 02 03 04         [01 04 03 02
     *  05 06 06 07          05 07 06 06
     *  08 09 10 11          08 11 10 09
     *  12 13 14 15]         12 15 14 13]
     *
     * @endcode
     */
    void selfReverseX()
    {
        size_t xsize=XSIZE(*this);
        size_t halfSizeX = (xsize-2)/2;
        size_t xsize_1=xsize-1;
        for (size_t k = 0; k < ZSIZE(*this); k++)
            for (size_t i = 0; i < YSIZE(*this); i++)
                for (size_t j = 0; j <=halfSizeX; j++)
                {
                    T aux;
                    T& d1=DIRECT_ZYX_ELEM(*this, k, i, j);
                    T& d2=DIRECT_ZYX_ELEM(*this, k, i, xsize_1 - j);
                    SWAP(d1,d2,aux);
                }
        STARTINGX(*this) = -FINISHINGX(*this);
    }

    /** Reverse matrix values over Y axis, keep in this object.
     *
     * Maybe better with an example:
     *
     * @code
     * slice 0
     * [01 02 03          [07 08 09
     *  04 05 06           04 05 06
     *  07 08 09]          01 02 03]
     *
     * @endcode
     *
     */
    void selfReverseY()
    {
        size_t ysize=YSIZE(*this);
        size_t halfSizeY = (ysize-2)/2;
        size_t ysize_1=ysize-1;
        for (size_t k = 0; k < ZSIZE(*this); k++)
            for (size_t i = 0; i <= halfSizeY; i++)
                for (size_t j = 0; j < XSIZE(*this); j++)
                {
                    T aux;
                    T& d1=DIRECT_ZYX_ELEM(*this, k, i, j);
                    T& d2=DIRECT_ZYX_ELEM(*this, k, ysize_1 - i, j);
                    SWAP(d1,d2,aux);
                }
        STARTINGY(*this) = -FINISHINGY(*this);
    }

    /** Reverse matrix values over Z axis, keep result in this object.
     *

     *
     */
    void selfReverseZ()
    {
        size_t zsize=ZSIZE(*this);
        size_t halfSizeZ = (zsize-2)/2;
        size_t zsize_1=zsize-1;
        for (size_t k = 0; k <= halfSizeZ; k++)
            for (size_t i = 0; i <YSIZE(*this); i++)
                for (size_t j = 0; j < XSIZE(*this); j++)
                {
                    T aux;
                    T& d1=DIRECT_ZYX_ELEM(*this, k, i, j);
                    T& d2=DIRECT_ZYX_ELEM(*this, zsize_1- k, i, j);
                    SWAP(d1,d2,aux);
                }
        STARTINGZ(*this) = -FINISHINGZ(*this);
    }

    /** Extracts the 1D profile between two points in a 2D array
     *
     * Given two logical indexes, this function returns samples of the line that
     * joins them. This is done by bilinear interpolation. The number of samples
     * in the line is N.
     */
    void profile(int x0, int y0, int xF, int yF, int N,
                 MultidimArray< double >& profile) const
    {
        checkDimension(2);
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

    /** Show using gnuplot
     *
     * This function uses gnuplot to plot this vector. You must supply the
     * xlabel, ylabel, and title.
     */
    void showWithGnuPlot(const String& xlabel, const String& title)
    {
        checkDimension(1);

        FileName fn_tmp;
        fn_tmp.initRandom(10);
        const char * fnStr = fn_tmp.c_str();
        MultidimArray<T>::write(formatString("PPP%s.txt", fnStr));

        std::ofstream fh_gplot;
        fh_gplot.open(formatString("PPP%s.gpl", fnStr).c_str());
        if (!fh_gplot)
            REPORT_ERROR(ERR_IO_NOTOPEN,
                         formatString("vector::showWithGnuPlot: Cannot open PPP%s.gpl for output", fnStr));
        fh_gplot << "set xlabel \"" + xlabel + "\"\n";
        fh_gplot << "plot \"PPP" + fn_tmp + ".txt\" title \"" + title +
        "\" w l\n";
        fh_gplot << "pause 300 \"\"\n";
        fh_gplot.close();
        system(formatString("(gnuplot PPP%s.gpl; rm PPP%s.txt PPP%s.gpl) &", fnStr, fnStr, fnStr).c_str());
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

        nam = formatString("PPP%s.txt", nam.c_str());
        write(nam);

        system(formatString("xmipp_edit -i %s -remove &", nam.c_str()).c_str());
    }

    /** Write to an ASCII file.
     */
    void write(const FileName& fn) const
    {
        std::ofstream out;
        out.open(fn.c_str(), std::ios::out);
        if (!out)
            REPORT_ERROR(ERR_IO_NOTOPEN,
                         formatString("MultidimArray::write: File %s cannot be opened for output", fn.c_str()));

        out << *this;
        out.close();
    }
    //@}

    /// @name Operators
    /// @{

    /** Assignment.
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
            if (data == NULL || !sameShape(op1))
                resizeNoCopy(op1);
            memcpy(data,op1.data,MULTIDIM_SIZE(op1)*sizeof(T));
        }
        return *this;
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
     *
     * This function is ported to Python as assign.
     */
    MultidimArray<T>& operator=(const Matrix2D<T>& op1)
    {
        resizeNoCopy(MAT_YSIZE(op1), MAT_XSIZE(op1));
        memcpy(data,MATRIX2D_ARRAY(op1), MAT_SIZE(op1)*sizeof(T));

        return *this;
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
    MultidimArray<T> operator-() const
    {
        MultidimArray<T> tmp(*this);
        T* ptr;
        size_t n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(tmp,n,ptr)
        *ptr = -(*ptr);
        return tmp;
    }

    /** Input from input stream.
     *
     * Actual size of the array is used to know how many values must be read.
     *
     * @code
     * v.<3);
     * std::cin >> v;
     * @endcode
     *
     * This function is not ported to Python.
     */
    friend std::istream& operator>>(std::istream& in, MultidimArray<T>& v)
    {
        T* ptr;
        size_t n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(v,n,ptr)
        in >> *ptr;
        return in;
    }

    void copy(Matrix2D<T>& op1) const
    {
        op1.resizeNoCopy(YSIZE(*this), XSIZE(*this));
        memcpy(MATRIX2D_ARRAY(op1), MULTIDIM_ARRAY(*this), MULTIDIM_SIZE(*this)*sizeof(T));
    }

    /** Equality.
     *
     * Returns true if this object has got the same shape (origin and size)
     * than the argument and the same values (within accuracy).
     */
    bool equal(const MultidimArray<T>& op,
               double accuracy = XMIPP_EQUAL_ACCURACY) const
    {
        if (!sameShape(op) || data==NULL || op.data == NULL)
            return false;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(*this)
        {
            if (fabs(DIRECT_MULTIDIM_ELEM(*this,n) -
                     DIRECT_MULTIDIM_ELEM(op,n)) > accuracy)
                return false;
        }
        return true;
    }
    //@}
};



/// @name Functions for all multidimensional arrays
/// @{

/** Conversion from one type to another.
 *
 * If we have an integer array and we need a double one, we can use this
 * function. The conversion is done through a type casting of each element
 * If n >= 0, only the nth volumes will be converted, otherwise all NSIZE volumes
 */
template<typename T1, typename T2>
void typeCast(const MultidimArray<T1>& v1,  MultidimArray<T2>& v2)
{
    if (NZYXSIZE(v1) == 0)
    {
        v2.clear();
        return;
    }

    v2.resizeNoCopy(v1);
    T1* ptr1 = MULTIDIM_ARRAY(v1);

    // Unroll the loop
    const size_t unroll=4;
    size_t nmax=(NZYXSIZE(v1)/unroll)*unroll;
    for (size_t n=0; n<nmax; n+=unroll, ptr1+=unroll)
    {
        DIRECT_MULTIDIM_ELEM(v2,n)   = static_cast< T2 >(*ptr1);
        DIRECT_MULTIDIM_ELEM(v2,n+1) = static_cast< T2 >(*(ptr1+1));
        DIRECT_MULTIDIM_ELEM(v2,n+2) = static_cast< T2 >(*(ptr1+2));
        DIRECT_MULTIDIM_ELEM(v2,n+3) = static_cast< T2 >(*(ptr1+3));
    }
    // Do the remaining elements
    for (size_t n=nmax; n<NZYXSIZE(v1); ++n, ++ptr1)
        DIRECT_MULTIDIM_ELEM(v2,n)   = static_cast< T2 >(*ptr1);
}

template<typename T1>
void typeCastComplex(const MultidimArray<T1>& v1,  MultidimArray<std::complex<double> >& v2)
{
    if (NZYXSIZE(v1) == 0)
    {
        v2.clear();
        return;
    }

    v2.resizeNoCopy(v1);
    T1* ptr1 = MULTIDIM_ARRAY(v1);
    double * ptr2 = (double*) MULTIDIM_ARRAY(v2);

    // Unroll the loop
    const size_t unroll=4;
    size_t nmax=(NZYXSIZE(v1)/unroll)*unroll;
    for (size_t n=0; n<nmax; n+=unroll, ptr1+=unroll)
    {
        *(ptr2++) = static_cast<double>(*ptr1);
        *(ptr2++) = 0;
        *(ptr2++) = static_cast< double >(*(ptr1+1));
        *(ptr2++) = 0;
        *(ptr2++) = static_cast< double >(*(ptr1+2));
        *(ptr2++) = 0;
        *(ptr2++) = static_cast< double >(*(ptr1+3));
        *(ptr2++) = 0;
    }
    // Do the remaining elements
    for (size_t n=nmax; n<NZYXSIZE(v1); ++n, ++ptr1)
    {
        *(ptr2++) = static_cast< double >(*ptr1);
        *(ptr2++) = 0;
    }
}

/** Conversion from one type to another.
 * In some cases, the two types are the same. So a faster way is simply by assignment.
 */
template<typename T1>
void typeCast(const MultidimArray<T1>& v1,  MultidimArray<T1>& v2)
{
    v2=v1;
}

/** Assignment.
 */
template<typename T>
void typeCast(const MultidimArray<T>& v1, Matrix1D<T> &v2)
{
    v2.resizeNoCopy(XSIZE(v1));
    memcpy(&VEC_ELEM(v2,0),&DIRECT_A1D_ELEM(v1,0),MULTIDIM_SIZE(v1)*sizeof(T));
    v2.row=false;
}


/** MultidimArray equality.*/
template<typename T>
bool operator==(const MultidimArray<T>& op1, const MultidimArray<T>& op2)
{
    return op1.equal(op2);
}

/** MultidimArray inequality.*/
template<typename T>
bool operator!=(const MultidimArray<T>& op1, const MultidimArray<T>& op2)
{
    return !(op1==op2);
}

/** Reduce both volumes to a common size.
 *
 * Search the range of logical indexes for which both volumes have got valid
 * values, and cut both to that size, the corresponding origin is automatically
 * computed.
 *
 * @code
 * MultidimArray< double > V1(4, 5, 3);
 * V1.startingX() = -2;
 * V1.startingY() = -2;
 * V1.startingZ() = -2;
 *
 * MultidimArray< double > V2(4, 2, 3);
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


/** Get Sin and Cos of vector x.
 */
void sincos(const MultidimArray<double> &x, MultidimArray<double> &s, MultidimArray<double> &c);

/** Obtains the plane parameters z=p0+p1*x+p2*y.
 */
void planeFit(const MultidimArray<double> &z, const MultidimArray<double> &x, const MultidimArray<double> &y,
              double &p0, double &p1, double &p2);

/*
   mod    Modulus after division.
    mod(x,y) is x - n.*y where n = floor(x./y) if y ~= 0.  If y is not an
    integer and the quotient x./y is within roundoff error of an integer,
    then n is that integer.  The inputs x and y must be real arrays of the
    same size, or real scalars.

    By convention:
       MOD(x, m, 0) is m = x.
       MOD(x,m,x) is m = 0.
       MOD(x,m,y), for x~=y and y~=0, m has the same sign as y.
 */
template <typename T>
void mod(const MultidimArray<T> &x, MultidimArray<T> &m, double y)
{
    m.resizeNoCopy(x);
    double *ptr=NULL;
    double *ptrm=MULTIDIM_ARRAY(m);
    size_t n;
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(x,n,ptr)
    *(ptrm++) = (*ptr) - std::floor((*ptr)/(y))*(y);
}

/** Output to output stream.
 * This function is not ported to Python.
 */
template <typename T>
std::ostream& operator<< (std::ostream& ostrm, const MultidimArray<T>& v)
{
    if (v.xdim == 0)
        ostrm << "NULL Array\n";
    else
        ostrm << std::endl;

    double max_val = ABS(DIRECT_A3D_ELEM(v , 0, 0, 0));

    T* ptr;
    size_t n;
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(v,n,ptr)
    max_val = XMIPP_MAX(max_val, ABS(*ptr));

    int prec = bestPrecision(max_val, 10);

    if (YSIZE(v)==1 && ZSIZE(v)==1)
    {
        for (int j = STARTINGX(v); j <= FINISHINGX(v); j++)
            ostrm << floatToString((double) A3D_ELEM(v, 0, 0, j), 10, prec)
            << std::endl;
    }
    else
    {
        for (size_t l = 0; l < NSIZE(v); l++)
        {
            if (NSIZE(v)>1)
                ostrm << "Image No. " << l << std::endl;
            for (int k = STARTINGZ(v); k <= FINISHINGZ(v); k++)
            {
                if (ZSIZE(v)>1)
                    ostrm << "Slice No. " << k << std::endl;
                for (int i = STARTINGY(v); i <= FINISHINGY(v); i++)
                {
                    for (int j = STARTINGX(v); j <= FINISHINGX(v); j++)
                    {
                        ostrm << floatToString((double) A3D_ELEM(v, k, i, j), 10, prec) << ' ';
                    }
                    ostrm << std::endl;
                }
            }
        }
    }

    return ostrm;
}
//@}

// Specializations cases for complex numbers
template<>
std::ostream& operator<<(std::ostream& ostrm, const MultidimArray< std::complex<double> >& v);
template<>
void MultidimArray< std::complex< double > >::computeDoubleMinMax(double& minval, double& maxval) const;
template<>
void MultidimArray< std::complex< double > >::computeDoubleMinMaxRange(double& minval, double& maxval, size_t pos, size_t size) const;
template<>
void MultidimArray< std::complex< double > >::rangeAdjust(std::complex< double > minF, std::complex< double > maxF);
template<>
double MultidimArray< std::complex< double > >::computeAvg() const;
template<>
void MultidimArray< std::complex< double > >::maxIndex(int &lmax, int& kmax, int& imax, int& jmax) const;
template<>
void MultidimArray<double>::computeAvgStdev(double& avg, double& stddev) const;
template<>
bool operator==(const MultidimArray< std::complex< double > >& op1,
                const MultidimArray< std::complex< double > >& op2);
template<>
double MultidimArray<double>::interpolatedElement2D(double x, double y, double outside_value) const;
template<>
void MultidimArray< std::complex< double > >::getReal(MultidimArray<double> & realImg) const;
template<>
void MultidimArray< std::complex< double > >::getImag(MultidimArray<double> & imagImg) const;

//@}
#endif
