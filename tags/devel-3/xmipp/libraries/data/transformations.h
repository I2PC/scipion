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
void applyGeometry(int SplineDegree, 
                   MultidimArray<T>& V2,
                   const MultidimArray<T>& V1,
                   const Matrix2D< double > A, bool inv, 
                   bool wrap, T outside = 0, unsigned long n = 0)
{

    if (&V1 == &V2)
        REPORT_ERROR(1101,"ApplyGeometry: Input array cannot be the same as output array");

    if ( V1.getDim()==2 && ((A.Xdim() != 3) || (A.Ydim() != 3)) )
        REPORT_ERROR(1102,"ApplyGeometry: 2D transformation matrix is not 3x3");

    if ( V1.getDim()==3 && ((A.Xdim() != 4) || (A.Ydim() != 4)) )
        REPORT_ERROR(1103,"ApplyGeometry: 3D transformation matrix is not 4x4");

    if (A.isIdentity())
    {
        V1.getImage(n, V2);
        return;
    }

    if (XSIZE(V1) == 0)
    {
        V2.clear();
        return;
    }

    MultidimArray<double> Bcoeffs;
    Matrix2D<double> Ainv;
    const Matrix2D<double> * Aptr=&A;
    if (!inv)
    {
        Ainv = A.inv();
        Aptr=&Ainv;
    }
    const Matrix2D<double> &Aref=*Aptr;

    // For scalings the output matrix is resized outside to the final
    // size instead of being resized inside the routine with the
    // same size as the input matrix
    if (XSIZE(V2) == 0)
        V2.resize(1, ZSIZE(V1), YSIZE(V1), XSIZE(V1));

    if (V1.getDim() == 2)
    {
        // 2D transformation

        int m1, n1, m2, n2;
        double x, y, xp, yp;
        double minxp, minyp, maxxp, maxyp;
        int cen_x, cen_y, cen_xp, cen_yp;
        double wx, wy; 
        int Xdim, Ydim;

        if (outside != 0.)
        {
            // Initialise output matrix with value=outside
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(V2)
            {
                DIRECT_A2D_ELEM(V2, i, j) = outside;
            }
        }

        // Find center and limits of image
        cen_y  = (int)(YSIZE(V2) / 2);
        cen_x  = (int)(XSIZE(V2) / 2);
        cen_yp = (int)(YSIZE(V1) / 2);
        cen_xp = (int)(XSIZE(V1) / 2);
        minxp  = -cen_xp;
        minyp  = -cen_yp;
        maxxp  = XSIZE(V1) - cen_xp - 1;
        maxyp  = YSIZE(V1) - cen_yp - 1;
        Xdim   = XSIZE(V1);
        Ydim   = YSIZE(V1);

        if (SplineDegree > 1)
        {
            // Build the B-spline coefficients
            produceSplineCoefficients(SplineDegree, Bcoeffs, V1, n); //Bcoeffs is a single image
            STARTINGX(Bcoeffs) = (int) minxp;
            STARTINGY(Bcoeffs) = (int) minyp;
        }

        // Now we go from the output image to the input image, ie, for any pixel
        // in the output image we calculate which are the corresponding ones in
        // the original image, make an interpolation with them and put this value
        // at the output pixel

#ifdef DEBUG_APPLYGEO
        std::cout << "A\n" << Aref << std::endl
                  << "(cen_x ,cen_y )=(" << cen_x  << "," << cen_y  << ")\n"
                  << "(cen_xp,cen_yp)=(" << cen_xp << "," << cen_yp << ")\n"
                  << "(min_xp,min_yp)=(" << minxp  << "," << minyp  << ")\n"
                  << "(max_xp,max_yp)=(" << maxxp  << "," << maxyp  << ")\n";
#endif

        for (int i = 0; i < YSIZE(V2); i++)
        {
            // Calculate position of the beginning of the row in the output image
            x = -cen_x;
            y = i - cen_y;

            // Calculate this position in the input image according to the
            // geometrical transformation
            // they are related by
            // coords_output(=x,y) = A * coords_input (=xp,yp)
            xp = x * Aref(0, 0) + y * Aref(0, 1) + Aref(0, 2);
            yp = x * Aref(1, 0) + y * Aref(1, 1) + Aref(1, 2);

            for (int j = 0; j < XSIZE(V2); j++)
            {
                bool interp;
                T tmp;

#ifdef DEBUG_APPLYGEO
                std::cout << "Computing (" << i << "," << j << ")\n";
                std::cout << "   (y, x) =(" << y << "," << x << ")\n"
                          << "   before wrapping (y',x')=(" << yp << "," << xp << ") "
                          << std::endl;
#endif
                // If the point is outside the image, apply a periodic extension
                // of the image, what exits by one side enters by the other
                interp = true;
                if (wrap)
                {
                    if (xp < minxp - XMIPP_EQUAL_ACCURACY ||
                        xp > maxxp + XMIPP_EQUAL_ACCURACY)
                        xp = realWRAP(xp, minxp - 0.5, maxxp + 0.5);
                    
                    if (yp < minyp - XMIPP_EQUAL_ACCURACY ||
                        yp > maxyp + XMIPP_EQUAL_ACCURACY)
                        
                        yp = realWRAP(yp, minyp - 0.5, maxyp + 0.5);
                }
                else
                {
                    if (xp < minxp - XMIPP_EQUAL_ACCURACY ||
                        xp > maxxp + XMIPP_EQUAL_ACCURACY)
                        interp = false;
                    
                    if (yp < minyp - XMIPP_EQUAL_ACCURACY ||
                        yp > maxyp + XMIPP_EQUAL_ACCURACY)
                        interp = false;
                }

#ifdef DEBUG_APPLYGEO
                std::cout << "   after wrapping (y',x')=(" << yp << "," << xp << ") "
                          << std::endl;
                std::cout << "   Interp = " << interp << std::endl;
                // The following line sounds dangerous...
                //x++;
#endif

                if (interp)
                {
                    if (SplineDegree==1)
                    {
                        // Linear interpolation

                        // Calculate the integer position in input image, be careful
                        // that it is not the nearest but the one at the top left corner
                        // of the interpolation square. Ie, (0.7,0.7) would give (0,0)
                        // Calculate also weights for point m1+1,n1+1
                        wx = xp + cen_xp;
                        m1 = (int) wx;
                        wx = wx - m1;
                        m2 = m1 + 1;
                        wy = yp + cen_yp;
                        n1 = (int) wy;
                        wy = wy - n1;
                        n2 = n1 + 1;
                    
                        // m2 and n2 can be out by 1 so wrap must be check here
                        if (wrap)
                        {
                            if (m2 >= Xdim)
                                m2 = 0;
                            if (n2 >= Ydim)
                                n2 = 0;
                        }
                    
#ifdef DEBUG_APPLYGEO
                        std::cout << "   From (" << n1 << "," << m1 << ") and ("
                                  << n2 << "," << m2 << ")\n";
                        std::cout << "   wx= " << wx << " wy= " << wy << std::endl;
#endif

                        // Perform interpolation
                        // if wx == 0 means that the rightest point is useless for this
                        // interpolation, and even it might not be defined if m1=xdim-1
                        // The same can be said for wy.
                        tmp  = (T)((1 - wy) * (1 - wx) * DIRECT_NZYX_ELEM(V1, n, 0, n1, m1));
                        
                        if (wx != 0 && m2 < V1.xdim)
                            tmp += (T)((1 - wy) * wx * DIRECT_NZYX_ELEM(V1, n, 0, n1, m2));
                    
                        if (wy != 0 && n2 < V1.ydim)
                        {
                            tmp += (T)(wy * (1 - wx) * DIRECT_NZYX_ELEM(V1, n, 0, n2, m1));
                            
                            if (wx != 0 && m2 < V1.xdim)
                                tmp += (T)(wy * wx * DIRECT_NZYX_ELEM(V1, n, 0, n2, m2));
                        }

                        dAij(V2, i, j) = tmp;
                    }
                    else
                    {
                        // B-spline interpolation

                        dAij(V2, i, j) = (T) Bcoeffs.interpolatedElementBSpline2D(
                            xp, yp, SplineDegree);
                    }
#ifdef DEBUG_APPYGEO
                    std::cout << "   val= " << dAij(V2, i, j) << std::endl;
#endif
                }

                // Compute new point inside input image
                xp += Aref(0, 0);
                yp += Aref(1, 0);
            }
        }
    }
    else
    {
        // 3D transformation

        int m1, n1, o1, m2, n2, o2;
        double x, y, z, xp, yp, zp;
        double minxp, minyp, maxxp, maxyp, minzp, maxzp;
        int cen_x, cen_y, cen_z, cen_xp, cen_yp, cen_zp;
        double wx, wy, wz;

        // Find center of MultidimArray
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

        if (SplineDegree > 1)
        {
            // Build the B-spline coefficients
            produceSplineCoefficients(SplineDegree, Bcoeffs, V1, n); //Bcoeffs is a single image
            STARTINGX(Bcoeffs) = (int) minxp;
            STARTINGY(Bcoeffs) = (int) minyp;
        }

        // Now we go from the output MultidimArray to the input MultidimArray, ie, for any
        // voxel in the output MultidimArray we calculate which are the corresponding
        // ones in the original MultidimArray, make an interpolation with them and put
        // this value at the output voxel

        // V2 is not initialised to 0 because all its pixels are rewritten
        for (int k = 0; k < V2.zdim; k++)
            for (int i = 0; i < V2.ydim; i++)
            {
                // Calculate position of the beginning of the row in the output
                // MultidimArray
                x = -cen_x;
                y = i - cen_y;
                z = k - cen_z;
                
                // Calculate this position in the input image according to the
                // geometrical transformation they are related by
                // coords_output(=x,y) = A * coords_input (=xp,yp)
                xp = x * Aref(0, 0) + y * Aref(0, 1) + z * Aref(0, 2) + Aref(0, 3);
                yp = x * Aref(1, 0) + y * Aref(1, 1) + z * Aref(1, 2) + Aref(1, 3);
                zp = x * Aref(2, 0) + y * Aref(2, 1) + z * Aref(2, 2) + Aref(2, 3);
                
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
                        if (SplineDegree == 1)
                        {

                            // Linear interpolation

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
                            tmp  = (T)((1 - wz) * (1 - wy) * (1 - wx) * DIRECT_NZYX_ELEM(V1, n, o1, n1, m1));
                        
                            if (wx != 0 && m2 < V1.xdim)
                                tmp += (T)((1 - wz) * (1 - wy) * wx * DIRECT_NZYX_ELEM(V1, n, o1, n1, m2));
                        
                            if (wy != 0 && n2 < V1.ydim)
                            {
                                tmp += (T)((1 - wz) * wy * (1 - wx) * DIRECT_NZYX_ELEM(V1, n, o1, n2, m1));
                                if (wx != 0 && m2 < V1.xdim)
                                    tmp += (T)((1 - wz) * wy * wx * DIRECT_NZYX_ELEM(V1, n, o1, n2, m2));
                            }
                        
                            if (wz != 0 && o2 < V1.zdim)
                            {
                                tmp += (T)(wz * (1 - wy) * (1 - wx) * DIRECT_NZYX_ELEM(V1, n, o2, n1, m1));
                                if (wx != 0 && m2 < V1.xdim)
                                    tmp += (T)(wz * (1 - wy) * wx * DIRECT_NZYX_ELEM(V1, n, o2, n1, m2));
                                if (wy != 0 && n2 < V1.ydim)
                                {
                                    tmp += (T)(wz * wy * (1 - wx) * DIRECT_NZYX_ELEM(V1, n, o2, n2, m1));
                                    if (wx != 0 && m2 < V1.xdim)
                                        tmp += (T)(wz * wy * wx * DIRECT_NZYX_ELEM(V1, n, o2, n2, m2));
                                }
                            }

#ifdef DEBUG
                            if (show_debug)
                                std::cout <<
                                    "tmp1=" << DIRECT_NZYX_ELEM(V1, n, o1, n1, m1) << " " 
                                            << (T)((1 - wz) *(1 - wy) *(1 - wx) * DIRECT_NZYX_ELEM(V1, n, o1, n1, m1)) 
                                            << std::endl <<
                                    "tmp2=" << DIRECT_NZYX_ELEM(V1, n, o1, n1, m2) << " " 
                                            << (T)((1 - wz) *(1 - wy) * wx * DIRECT_NZYX_ELEM(V1, n, o1, n1, m2)) 
                                            << std::endl <<
                                    "tmp3=" << DIRECT_NZYX_ELEM(V1, n, o1, n2, m1) << " " 
                                            << (T)((1 - wz) * wy *(1 - wx) * DIRECT_NZYX_ELEM(V1, n, o1, n2, m1)) 
                                            << std::endl <<
                                    "tmp4=" << DIRECT_NZYX_ELEM(V1, n, o1, n2, m2) << " " 
                                            << (T)((1 - wz) * wy * wx * DIRECT_NZYX_ELEM(V1, n, o1, n2, m2)) 
                                            << std::endl <<
                                    "tmp5=" << DIRECT_NZYX_ELEM(V1, n, o2, n1, m1) << " " 
                                            << (T)(wz * (1 - wy) *(1 - wx) * DIRECT_NZYX_ELEM(V1, n, o2, n1, m1))
                                            << std::endl <<
                                    "tmp6=" << DIRECT_NZYX_ELEM(V1, n, o2, n1, m2) << " " 
                                            << (T)(wz * (1 - wy) * wx * DIRECT_NZYX_ELEM(V1, n, o2, n1, m2)) 
                                            << std::endl <<
                                    "tmp7=" << DIRECT_NZYX_ELEM(V1, n, o2, n2, m1) << " " 
                                            << (T)(wz * wy *(1 - wx) * DIRECT_NZYX_ELEM(V1, n, o2, n2, m1)) 
                                            << std::endl <<
                                    "tmp8=" << DIRECT_NZYX_ELEM(V1, n, o2, n2, m2) << " " 
                                            << (T)(wz * wy * wx * DIRECT_NZYX_ELEM(V1, n, o2, n2, m2)) 
                                            << std::endl <<
                                    "tmp= " << tmp << std::endl;
#endif

                            dAkij(V2 , k, i, j) = tmp;
                        }
                        else
                        {
                            // B-spline interpolation

                            dAkij(V2, k, i, j) =
                                (T) Bcoeffs.interpolatedElementBSpline3D(xp, yp, zp,SplineDegree);

                        }
                    }
                    else
                        dAkij(V2, k, i, j) = outside;


                    // Compute new point inside input image
                    xp += Aref(0, 0);
                    yp += Aref(1, 0);
                    zp += Aref(2, 0);
                }
            }
    }

}



/** Applies a geometrical transformation and overwrites the input matrix.
 * @ingroup GeometricalTransformations
 *
 * The same as the previous function, but input array is overwritten
 */
template<typename T>
void selfApplyGeometry(int SplineDegree, 
                       MultidimArray<T>& V1,
                       const Matrix2D< double > A, bool inv, 
                       bool wrap, T outside = 0, unsigned long n = 0)
{
    MultidimArray<T> aux = V1;
    applyGeometry(SplineDegree, V1, aux, A, inv, wrap, outside, n);
}

//Special cases for complex arrays
void applyGeometry(int SplineDegree, 
                   MultidimArray< std::complex<double> > &V2,
                   const MultidimArray< std::complex<double> > &V1,  
                   const Matrix2D<double> &A, bool inv, 
                   bool wrap, std::complex<double> outside, unsigned long n = 0);

//Special cases for complex arrays
void selfApplyGeometry(int SplineDegree, 
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
void produceSplineCoefficients(int SplineDegree, 
                               MultidimArray< double > &coeffs,
                               const MultidimArray< T > &V1,  
                               unsigned long n = 0)
{

    coeffs.initZeros(ZSIZE(V1), YSIZE(V1), XSIZE(V1));
    STARTINGX(coeffs) = STARTINGX(V1);
    STARTINGY(coeffs) = STARTINGY(V1);
    STARTINGZ(coeffs) = STARTINGZ(V1);

    int Status;
    MultidimArray< double > aux;
    typeCast(V1, aux, n); // This will create a single volume!

    ChangeBasisVolume(MULTIDIM_ARRAY(aux), MULTIDIM_ARRAY(coeffs),
                      XSIZE(V1), YSIZE(V1), ZSIZE(V1),
                      CardinalSpline, BasicSpline, SplineDegree,
                      MirrorOffBounds, DBL_EPSILON, &Status);
    if (Status)
        REPORT_ERROR(1200, "Error in produceSplineCoefficients...");

}


// Special case for complex arrays
void produceSplineCoefficients(int SplineDegree, 
                               MultidimArray< double > &coeffs,
                               const MultidimArray< std::complex<double> > &V1,  
                               unsigned long n = 0);

/** Produce image from B-spline coefficients.
 * @ingroup  GeometricalTransformations
 * Note that the coeffs and img are only single images!
 */
void produceImageFromSplineCoefficients(int SplineDegree, 
                                        MultidimArray< double >& img, 
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
void rotate(int SplineDegree, 
            MultidimArray<T>& V2, 
            const MultidimArray<T>& V1, 
            double ang, char axis = 'Z', 
            bool wrap = DONT_WRAP, T outside = 0, unsigned long n = 0)
{
    Matrix2D< double > tmp;
    if (V1.getDim()==2)
    {
        tmp = rotation2DMatrix(ang);
    }
    else if (V1.getDim()==3)
    {
        tmp = rotation3DMatrix(ang, axis);
    }
    else
        REPORT_ERROR(1,"rotate ERROR: rotate only valid for 2D or 3D arrays");

    applyGeometry(SplineDegree, V2, V1, tmp, IS_NOT_INV, wrap, outside, n);

}


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
void translate(int SplineDegree, 
               MultidimArray<T> &V2,
               const MultidimArray<T> &V1, 
               const Matrix1D< double >& v, 
               bool wrap = WRAP, T outside = 0, unsigned long n = 0)
{
    Matrix2D< double > tmp;
    if (V2.getDim()==2)
        tmp = translation2DMatrix(v);
    else if (V2.getDim()==3)
        tmp = translation3DMatrix(v);
    else
        REPORT_ERROR(1,"translate ERROR: translate only valid for 2D or 3D arrays");
    applyGeometry(SplineDegree, V2, V1, tmp, IS_NOT_INV, wrap, outside, n);
}


/** Translate center of mass to center
 * @ingroup GeometricalTransformations
 *
 * If the input has very high values, it is better to rescale it to be
 * between 0 and 1.
 */
template<typename T>
void translateCenterOfMassToCenter(int SplineDegree, 
                                   MultidimArray<T> &V2,
                                   const MultidimArray<T> &V1, 
                                   bool wrap = WRAP, unsigned long n = 0)
{
    V2 = V1;
    V2.setXmippOrigin();
    Matrix1D< double > center;
    V2.centerOfMass(center,NULL,n);
    center *= -1;
    translate(SplineDegree, V2, V1, center, wrap, 0, n);
}

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
void scaleToSize(int SplineDegree, 
                 MultidimArray<T> &V2,
                 const MultidimArray<T> &V1,
                 int Xdim, int Ydim, int Zdim = 1,
                 unsigned long n = 0)
{

    Matrix2D< double > tmp;
    if (V1.getDim()==2)
    {
        tmp.initIdentity(3);
        tmp(0, 0) = (double) Xdim / (double) XSIZE(V1);
        tmp(1, 1) = (double) Ydim / (double) YSIZE(V1);
        V2.resize(1, 1, Ydim, Xdim);
    }
    else if (V1.getDim()==3)
    {
        tmp.initIdentity(4);
        tmp(0, 0) = (double) Xdim / (double) XSIZE(V1);
        tmp(1, 1) = (double) Ydim / (double) YSIZE(V1);
        tmp(2, 2) = (double) Zdim / (double) ZSIZE(V1);
        V2.resize(1, Zdim, Ydim, Xdim);
    }
    else
        REPORT_ERROR(1,"scaleToSize ERROR: scaleToSize only valid for 2D or 3D arrays");

    //FIXME I Dont know why the compiler does not let me use ",0, n" ...
    //applyGeometry(SplineDegree, V2, V1, tmp, IS_NOT_INV, WRAP, 0, n);
    applyGeometry(SplineDegree, V2, V1, tmp, IS_NOT_INV, WRAP, (T)0, n);

}


// Special case for complex arrays
void scaleToSize(int SplineDegree, 
                 MultidimArray< std::complex<double> > &V2,
                 const MultidimArray< std::complex<double> > &V1,
                 int Xdim, int Ydim, int Zdim = 1,
                 unsigned long n = 0);

/** Reduce the nth volume by 2 using a BSpline pyramid.
 * @ingroup GeometricalTransformations
 */
template<typename T>
void pyramidReduce(int SplineDegree, 
                   MultidimArray<T> &V2,
                   const MultidimArray<T> &V1,
                   int levels = 1, 
                   unsigned long n = 0)
{
    MultidimArray< double > coeffs;
    produceSplineCoefficients(SplineDegree, coeffs, V1, n);

    for (int i = 0; i < levels; i++)
    {
        reduceBSpline(3, V2, coeffs);
        coeffs = V2;
    }

    produceImageFromSplineCoefficients(3, V2, coeffs);

}

/** Expand the nth volume by 2 using a BSpline pyramid.
 * @ingroup GeometricalTransformations
 */
template<typename T>
void pyramidExpand(int SplineDegree, 
                   MultidimArray<T> &V2,
                   const MultidimArray<T> &V1,
                   int levels = 1, 
                   unsigned long n = 0)
{
    MultidimArray< double > coeffs;
    produceSplineCoefficients(SplineDegree, coeffs, V1, n);

    for (int i = 0; i < levels; i++)
    {
        expandBSpline(3, V2, coeffs);
        coeffs = V2;
    }

    produceImageFromSplineCoefficients(3, V2, coeffs);

}

/** Reduce a set of B-spline coefficients.
 * @ingroup GeometricalTransformations
 *
 * Knowing that V1 is a (single image of) B-spline coefficients, produce the
 * reduced set of B-spline coefficients using the two-scale relationship.
 */
template<typename T>
void reduceBSpline(int SplineDegree, 
                   MultidimArray< double >& V2, 
                   const MultidimArray<T> &V1)
{
    double g[200]; // Coefficients of the reduce filter
    long ng; // Number of coefficients of the reduce filter
    double h[200]; // Coefficients of the expansion filter
    long nh; // Number of coefficients of the expansion filter
    short IsCentered; // Equal TRUE if the filter is a centered spline

    // Get the filter
    const char *splineType="Centered Spline";
    if (GetPyramidFilter(splineType, SplineDegree,
                         g, &ng, h, &nh, &IsCentered))
        REPORT_ERROR(1, "Unable to load the filter coefficients");

    MultidimArray< double>  aux;
    typeCast(V1, aux);
    if (V1.getDim() == 2)
    {
       if (XSIZE(aux) % 2 != 0 && YSIZE(aux) % 2 != 0)
            aux.resize(YSIZE(aux) - 1, XSIZE(aux) - 1);
        else if (YSIZE(aux) % 2 != 0)
            aux.resize(YSIZE(aux) - 1, XSIZE(aux));
        else if (XSIZE(aux) % 2 != 0)
            aux.resize(YSIZE(aux), XSIZE(aux) - 1);

       V2.resize(YSIZE(aux) / 2, XSIZE(aux) / 2);
       Reduce_2D(MULTIDIM_ARRAY(aux), XSIZE(aux), YSIZE(aux),
                 MULTIDIM_ARRAY(V2), g, ng, IsCentered);
    }
    else if (V1.getDim() == 3)
    {
        if (XSIZE(aux) % 2 != 0 && YSIZE(aux) % 2 != 0 && ZSIZE(aux) % 2 != 0)
            aux.resize(ZSIZE(aux - 1), YSIZE(aux) - 1, XSIZE(aux) - 1);
        else if (XSIZE(aux) % 2 != 0 && YSIZE(aux) % 2 != 0 && ZSIZE(aux) % 2 == 0)
            aux.resize(ZSIZE(aux), YSIZE(aux) - 1, XSIZE(aux) - 1);
        else if (XSIZE(aux) % 2 != 0 && YSIZE(aux) % 2 == 0 && ZSIZE(aux) % 2 != 0)
            aux.resize(ZSIZE(aux) - 1, YSIZE(aux), XSIZE(aux) - 1);
        else if (XSIZE(aux) % 2 != 0 && YSIZE(aux) % 2 == 0 && ZSIZE(aux) % 2 == 0)
            aux.resize(ZSIZE(aux), YSIZE(aux), XSIZE(aux) - 1);
        else if (XSIZE(aux) % 2 == 0 && YSIZE(aux) % 2 != 0 && ZSIZE(aux) % 2 != 0)
            aux.resize(ZSIZE(aux) - 1, YSIZE(aux) - 1, XSIZE(aux));
        else if (XSIZE(aux) % 2 == 0 && YSIZE(aux) % 2 != 0 && ZSIZE(aux) % 2 == 0)
            aux.resize(ZSIZE(aux), YSIZE(aux) - 1, XSIZE(aux));
        else if (XSIZE(aux) % 2 == 0 && YSIZE(aux) % 2 == 0 && ZSIZE(aux) % 2 != 0)
            aux.resize(ZSIZE(aux) - 1, YSIZE(aux), XSIZE(aux));

        V2.resize(ZSIZE(aux) / 2, YSIZE(aux) / 2, XSIZE(aux) / 2);
        Reduce_3D(MULTIDIM_ARRAY(aux), XSIZE(aux), YSIZE(aux), ZSIZE(aux),
                  MULTIDIM_ARRAY(V2), g, ng, IsCentered);
    }
    else
        REPORT_ERROR(1,"reduceBSpline ERROR: only valid for 2D or 3D arrays");


}
  
/** Expand a set of B-spline coefficients.
 * @ingroup GeometricalTransformations
 *
 * Knowing that V1 is a (single image of) B-spline coefficients, produce the
 * expanded set of B-spline coefficients using the two-scale relationship.
 */
template<typename T>
void expandBSpline3D(int SplineDegree, 
                     MultidimArray< double >& V2, 
                     const MultidimArray<T> &V1)
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

    MultidimArray< double > aux;
    typeCast(V1, aux);

    if (V1.getDim() == 2)
    {
        V2.resize(2 * YSIZE(aux), 2 * XSIZE(aux));
        Expand_2D(MULTIDIM_ARRAY(aux), XSIZE(aux), YSIZE(aux),
                  MULTIDIM_ARRAY(V2), h, nh, IsCentered);
    }
    else if (V1.getDim() == 3)
    {
        V2.resize(2 * ZSIZE(aux), 2 * YSIZE(aux), 2 * XSIZE(aux));
        Expand_3D(MULTIDIM_ARRAY(aux), XSIZE(aux), YSIZE(aux), ZSIZE(aux),
                  MULTIDIM_ARRAY(V2), h, nh, IsCentered);
    }
    else
        REPORT_ERROR(1,"expandBSpline ERROR: only valid for 2D or 3D arrays");


}

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
                   unsigned long n = 0)
{
    Matrix1D< double > idx(3);

    // If center_of_rot was written for 2D image
    if (center_of_rot.size() < 3)
        center_of_rot.resize(3);

    // First determine the maximum distance that one should expect, to set the
    // dimension of the radial average vector
    MultidimArray< int > distances(8);

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
    FOR_ALL_ELEMENTS_IN_ARRAY3D(m)
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
        radial_mean(distance) += NZYX_ELEM(m, n, k, i, j);

        // Count the pixel
        radial_count(distance)++;
    }

    // Perform the mean
    FOR_ALL_ELEMENTS_IN_ARRAY1D(radial_mean)
        radial_mean(i) /= (T) radial_count(i);
}


#endif
