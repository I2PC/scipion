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

#ifndef MULTIDIMARRAY_3D_H
#define MULTIDIMARRAY_3D_H

/*
 * THIS HEADER FILE IS INCLUDED IN multidimensional_array.h
 */


/// @defgroup Volumes Volumes
/// @ingroup MultidimensionalArrays


/** Template class for Xmipp volumes
 */
template<typename T>
class Matrix3D: public MultidimArray<T>
{
public:

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





    /** @defgroup VolumesGeometrical Geometrical Transformations
     * @ingroup MultidimArray
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


};

/** @defgroup VolumesRelated Related functions
 * @ingroup Volumes
 *
 * These functions are not methods of Matrix1D
 */


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
