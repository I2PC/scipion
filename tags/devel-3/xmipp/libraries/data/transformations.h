/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *              Sjors H.W. Scheres
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

#ifndef TRANSFORMATIONS_H
#define TRANSFORMATIONS_H

#include "multidim_array.h"
#include "geometry.h"

#define IS_INV true
#define IS_NOT_INV false
#define DONT_WRAP false
#define WRAP true
#define DBL_EPSILON 1e-50


/// @defgroup GeometricalTransformations Geometrical transformations
/// @ingroup DataLibrary

/** Applies a geometrical transformation.
 * @ingroup GeometricalTransformations
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
 * The result volume (with ndim=1) is resized to the same
 * dimensions as V1 if V2 is empty (0x0) at the beginning, if it
 * is not, ie, if V2 has got some size then only those values in
 * the volume are filled, this is very useful for resizing the
 * volume, then you manually resize the output volume to the
 * desired size and then call this routine.
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

template<typename T>
void applyGeometry(int Splinedegree, 
                   MultidimArray<T>& V2,
                   const MultidimArray<T>& V1,
                   const Matrix2D< double > A, bool inv, 
                   bool wrap, T outside = 0, unsigned long n = 0);


/** Applies a geometrical transformation and overwrites the input matrix.
 * @ingroup GeometricalTransformations
 *
 * The same as the previous function, but input array is overwritten
 */
template<typename T>
void selfApplyGeometry(int Splinedegree, 
                       MultidimArray<T>& V1,
                       const Matrix2D< double > A, bool inv, 
                       bool wrap, T outside = 0, unsigned long n = 0);

//Special cases for complex arrays
void applyGeometry(int Splinedegree, 
                   MultidimArray< std::complex<double> > &V2,
                   const MultidimArray< std::complex<double> > &V1,  
                   const Matrix2D<double> &A, bool inv, 
                   bool wrap, std::complex<double> outside, unsigned long n = 0);
//Special cases for complex arrays
void selfApplyGeometry(int Splinedegree, 
                       MultidimArray< std::complex<double> > &V1,
                       const Matrix2D<double> &A, bool inv, 
                       bool wrap, std::complex<double> outside, unsigned long n = 0);


/** Produce spline coefficients.
 * @ingroup  GeometricalTransformations
 *
 * Create a single image with spline coefficients for the nth image
 *
 */
template<typename T>
void produceSplineCoefficients(int Splinedegree, 
                               MultidimArray< double > &coeffs,
                               const MultidimArray< T > &V1,  
                               unsigned long n = 0);

// Special case for complex arrays
void produceSplineCoefficients(int Splinedegree, 
                               MultidimArray< double > &coeffs,
                               const MultidimArray< std::complex<double> > &V1,  
                               unsigned long n = 0);

/** Produce image from B-spline coefficients.
 * @ingroup  GeometricalTransformations
 * Note that the coeffs and img are only single images!
 */
template<typename T>
void produceImageFromSplineCoefficients(int Splinedegree, 
                                        MultidimArray< T >& img, 
                                        const MultidimArray< double > &coeffs);

/** Rotate an array around a given system axis.
 * @ingroup GeometricalTransformations
 *
 * The rotation angle is in degrees, and the rotational axis is either
 * 'X', 'Y' or 'Z' for 3D arrays and only 'Z' for 2D arrays. An
 * exception is thrown if the axis given is not one of these.
 *
 * @code
 * V2 = V1.rotate(60);
 * @endcode
 */
template<typename T>
void rotate(int Splinedegree, 
            MultidimArray<T>& V2, 
            const MultidimArray<T>& V1, 
            double ang, char axis = 'Z', 
            bool wrap = DONT_WRAP, T outside = 0, unsigned long n = 0);

/** Rotate a 3D array around any axis.
 * @ingroup GeometricalTransformations
 *
 * The rotation angle is in degrees, and the rotational axis is given as a
 * R3 vector. An exception is thrown if the axis is not a R3 vector. The
 * axis needs not to be unitary.
 *
 * @code
 * rotate(3, V2, V1, 60, vectorR3(1, 1, 1));
 * @endcode
 */
template<typename T>
void rotate(int Splinedegree, 
            MultidimArray<T> &V2,
            const MultidimArray<T> &V1, 
            double ang, const Matrix1D< double >& axis, 
            bool wrap = DONT_WRAP, T outside = 0, unsigned long n = 0);

/** Translate a volume.
 * @ingroup GeometricalTransformations
 *
 * The shift is given as a R2 or R3 vector (shift_X, shift_Y, shift_Z) for 2D and 3D arrays, respectively.
 * An exception is thrown if the displacement is not a R3 vector.
 *
 * @code
 * // Displacement of 2 pixels down
 * V2 = V1.translate(vectorR3(0, 0, 2));
 * @endcode
 */
template<typename T>
void translate(int Splinedegree, 
               MultidimArray<T> &V2,
               const MultidimArray<T> &V1, 
               const Matrix1D< double >& v, 
               bool wrap = WRAP, T outside = 0, unsigned long n = 0);

/** Translate center of mass to center
 * @ingroup GeometricalTransformations
 *
 * If the input has very high values, it is better to rescale it to be
 * between 0 and 1.
 */
template<typename T>
void translateCenterOfMassToCenter(int Splinedegree, 
                                   MultidimArray<T> &V2,
                                   const MultidimArray<T> &V1, 
                                   bool wrap = WRAP, unsigned long n = 0);

/** Scales to a new size.
 * @ingroup GeometricalTransformations
 *
 * The volume is scaled (resampled) to fill a new size. It is not the
 * same as "window" in this same class. The size can be larger or smaller
 * than the actual one.
 *
 * @code
 * V2 = V1.scaleToSize(128, 128, 128);
 * @endcode
 */
template<typename T>
void scaleToSize(int Splinedegree, 
                 MultidimArray<T> &V2,
                 const MultidimArray<T> &V1,
                 int Xdim, int Ydim, int Zdim = 1,
                 unsigned long n = 0);

// Special case for complex arrays
void scaleToSize(int Splinedegree, 
                 MultidimArray< std::complex<double> > &V2,
                 const MultidimArray< std::complex<double> > &V1,
                 int Xdim, int Ydim, int Zdim = 1,
                 unsigned long n = 0);

/** Reduce the nth volume by 2 using a BSpline pyramid.
 * @ingroup GeometricalTransformations
 */
template<typename T>
void pyramidReduce(int Splinedegree, 
                   MultidimArray<T> &V2,
                   const MultidimArray<T> &V1,
                   int levels = 1, 
                   unsigned long n = 0);

/** Expand the nth volume by 2 using a BSpline pyramid.
 * @ingroup GeometricalTransformations
 */
template<typename T>
void pyramidExpand(int Splinedegree, 
                   MultidimArray<T> &V2,
                   const MultidimArray<T> &V1,
                   int levels = 1, 
                   unsigned long n = 0);

/** Reduce a set of B-spline coefficients.
 * @ingroup GeometricalTransformations
 *
 * Knowing that V1 is a (single image of) B-spline coefficients, produce the
 * reduced set of B-spline coefficients using the two-scale relationship.
 */
template<typename T>
void reduceBSpline(int Splinedegree, 
                   MultidimArray< double >& V2, 
                   const MultidimArray<T> &V1);
  
/** Expand a set of B-spline coefficients.
 * @ingroup GeometricalTransformations
 *
 * Knowing that V1 is a (single image of) B-spline coefficients, produce the
 * expanded set of B-spline coefficients using the two-scale relationship.
 */
template<typename T>
void expandBSpline3D(int Splinedegree, 
                     MultidimArray< double >& V2, 
                     const MultidimArray<T> &V1);

/** Does a radial average of a 2D/3D image, around the voxel where is the origin.
 * @ingroup GeometricalTransformations
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
void radialAverage(const MultidimArray< T >& m,
                   const Matrix1D< int >& center_of_rot,
                   MultidimArray< T >& radial_mean,
                   MultidimArray< int >& radial_count,
                   const bool& rounding = false,
                   unsigned long n = 0);

/** Creates a rotational matrix (3x3) for images
 * @ingroup GeometricalTransformations
 *
 * The rotation angle is in degrees.
 *
 * @code
 * m = rotation2DMatrix(60);
 * @endcode
 */
Matrix2D< double > rotation2DMatrix(double ang);

/** Creates a rotational matrix (3x3) for images
 * @ingroup GeometricalTransformations
 *
 * The rotation angle is in degrees.
 * m must have been already resized to 3x3
 *
 * @code
 *  rotation2DMatrix(60,m);
 * @endcode
 */
void rotation2DMatrix(double ang, Matrix2D< double > &m);

/** Creates a translational matrix (3x3) for images
 * @ingroup GeometricalTransformations
 *
 * The shift is given as a R2 vector (shift_X, shift_Y). An exception is thrown
 * if the displacement is not a R2 vector.
 *
 * @code
 * // Displacement of 1 pixel to the right
 * m = translation2DMatrix(vectorR2(1, 0));
 * @endcode
 */
Matrix2D< double > translation2DMatrix(const Matrix1D< double > &v);

/** Creates a rotational matrix (4x4) for volumes around system axis
 * @ingroup GeometricalTransformations
 *
 * The rotation angle is in degrees, and the rotational axis is either 'X', 'Y'
 * or 'Z'. An exception is thrown if the axis given is not one of these.
 *
 * The returned matrices are respectively alpha degrees around Z
 *
 * @code
 * [ cos(A) -sin(A)     0   ]
 * [ sin(A)  cos(A)     0   ]
 * [   0       0        1   ]
 * @endcode
 *
 * alpha degrees around Y
 * @code
 * [ cos(A)    0    -sin(A) ]
 * [   0       1       0    ]
 * [ sin(A)    0     cos(A) ]
 * @endcode
 *
 * alpha degrees around X
 * @code
 * [   1       0       0    ]
 * [   0     cos(A) -sin(A) ]
 * [   0     sin(A)  cos(A) ]
 * @endcode
 *
 * @code
 * m = rotation3DMatrix(60, 'X');
 * @endcode
 */
Matrix2D< double > rotation3DMatrix(double ang, char axis);

/** Creates a rotational matrix (4x4) for volumes around any axis
 * @ingroup GeometricalTransformations
 *
 * The rotation angle is in degrees, and the rotational axis is given as a R3
 * vector. An exception is thrown if the axis is not a R3 vector. The axis needs
 * not to be unitary.
 *
 * @code
 * m = rotation3DMatrix(60, vectorR3(1, 1, 1));
 * @endcode
 */
Matrix2D< double > rotation3DMatrix(double ang, const Matrix1D< double >& axis);

/** Matrix which transforms the given axis into Z
 * @ingroup GeometricalTransformations
 *
 * A geometrical transformation matrix (4x4) is returned such that the given
 * axis is rotated until it is aligned with the Z axis. This is very useful in
 * order to produce rotational matrices, for instance, around any axis.
 *
 * @code
 * Matrix2D< double > A = alignWithZ(axis);
 * return A.transpose() * rotation3DMatrix(ang, 'Z') * A;
 * @endcode
 *
 * The returned matrix is such that A*axis=Z, where Z and axis are column
 * vectors.
 */
Matrix2D< double > alignWithZ(const Matrix1D< double >& axis);

/** Creates a translational matrix (4x4) for volumes
 * @ingroup GeometricalTransformations
 *
 * The shift is given as a R3 vector (shift_X, shift_Y, shift_Z). An exception
 * is thrown if the displacement is not a R3 vector.
 *
 * @code
 * // Displacement of 2 pixels down
 * m = translation3DMatrix(vectorR3(0, 0, 2));
 * @endcode
 */
Matrix2D< double > translation3DMatrix(const Matrix1D< double >& v);

/** Creates a scaling matrix (4x4) for volumes
 * @ingroup GeometricalTransformations
 *
 * The scaling factors for the different axis must be given as a vector. So
 * that, XX(sc)=scale for X axis, YY(sc)=...
 */
Matrix2D< double > scale3DMatrix(const Matrix1D< double >& sc);


#endif
