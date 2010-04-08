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

#include "transformations.h"

template<typename T>
void applyGeometry(int Splinedegree, 
                   MultidimArray<T>& V2, 
                   const MultidimArray<T>& V1, 
                   const Matrix2D< double > &A, bool inv, 
                   bool wrap, T outside, unsigned long n)
{

    if (&V1 == &V2)
        REPORT_ERROR(1101,"ApplyGeometry: Input array cannot be the same as output array");

    if (ZSIZE(V1) == 1 && ((XSIZE(A) != 3) || (YSIZE(A) != 3)) )
        REPORT_ERROR(1102,"ApplyGeometry: 2D transformation matrix is not 3x3");

    if (ZSIZE(V1) > 1 && ((XSIZE(A) != 4) || (YSIZE(A) != 4)) )
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

    MultidimArray Bcoeffs;
    if (Splinedegree > 1)
    {
        // Build the B-spline coefficients
        //FIXME get produceSplineCoefficients OUTSIDE multidimArray!
        produceSplineCoefficients(Splinedegree, Bcoeffs, V1, n); //Bcoeffs is a single image
        STARTINGX(Bcoeffs) = (int) minxp;
        STARTINGY(Bcoeffs) = (int) minyp;
        STARTINGZ(Bcoeffs) = (int) minzp;
    }

    if (ZSIZE(V1) == 1)
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
            FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(M2)
            {
                DIRECT_MAT_ELEM(M2, i, j) = outside;
            }
        }

        // Find center and limits of image
        cen_y  = (int)(YSIZE(M2) / 2);
        cen_x  = (int)(XSIZE(M2) / 2);
        cen_yp = (int)(YSIZE(M1) / 2);
        cen_xp = (int)(XSIZE(M1) / 2);
        minxp  = -cen_xp;
        minyp  = -cen_yp;
        maxxp  = XSIZE(M1) - cen_xp - 1;
        maxyp  = YSIZE(M1) - cen_yp - 1;
        Xdim   = XSIZE(M1);
        Ydim   = YSIZE(M1);

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

        for (int i = 0; i < YSIZE(M2); i++)
        {
            // Calculate position of the beginning of the row in the output image
            x = -cen_x;
            y = i - cen_y;

            // Calculate this position in the input image according to the
            // geometrical transformation
            // they are related by
            // coords_output(=x,y) = A * coords_input (=xp,yp)
            xp = x * dMij(Aref, 0, 0) + y * dMij(Aref, 0, 1) + dMij(Aref, 0, 2);
            yp = x * dMij(Aref, 1, 0) + y * dMij(Aref, 1, 1) + dMij(Aref, 1, 2);

            for (int j = 0; j < XSIZE(M2); j++)
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
                    if (Splinedegree==1)
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
                        tmp  = (T)((1 - wy) * (1 - wx) * DIRECT_NZYX_ELEM(M1, n, 0, n1, m1));
                        
                        if (wx != 0 && m2 < M1.xdim)
                            tmp += (T)((1 - wy) * wx * DIRECT_NZYX_ELEM(M1, n, 0, n1, m2));
                    
                        if (wy != 0 && n2 < M1.ydim)
                        {
                            tmp += (T)(wy * (1 - wx) * DIRECT_NZYX_ELEM(M1, n, 0, n2, m1));
                            
                            if (wx != 0 && m2 < M1.xdim)
                                tmp += (T)(wy * wx * DIRECT_NZYX_ELEM(M1, n, 0, n2, m2));
                        }

                        dMij(M2, i, j) = tmp;
                    }
                    else
                    {
                        // B-spline interpolation

                        dMij(M2, i, j) = (T) Bcoeffs.interpolatedElementBSpline(
                            xp, yp, Splinedegree);
                    }
#ifdef DEBUG_APPYGEO
                    std::cout << "   val= " << dMij(M2, i, j) << std::endl;
#endif
                }

                // Compute new point inside input image
                xp += dMij(Aref, 0, 0);
                yp += dMij(Aref, 1, 0);
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
                        if (Splinedegree == 1)
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

                            dVkij(V2 , k, i, j) = tmp;
                        }
                        else
                        {
                            // B-spline interpolation

                            dVkij(V2, k, i, j) =
                                (T) Bcoeffs.interpolatedElementBSpline(xp, yp, zp,Splinedegree);

                        }
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

}


// Special case for complex numbers
template <>
void applyGeometry(int Splinedegree, 
                   MultidimArray< std::complex<double> > &V2,
                   const MultidimArray< std::complex<double> > &V1,
                   const Matrix2D<double> &A, bool inv, 
                   bool wrap, std::complex<double> outside, unsigned long n)
{

    if (Splinedegree > 1)
    {
        MultidimArray<double> re, im, rotre, rotim;
        MultidimArray<std::complex<double> > oneImg;
        double outre, outim;
        re.resize(ZSIZE(V1), YSIZE(V1), XSIZE(V1));
        im.resize(ZSIZE(V1), YSIZE(V1), XSIZE(V1));
        outre = outside.real();
        outim = outside.imag();
        V1.getImage(n, oneImg);
        Complex2RealImag(MULTIDIM_ARRAY(oneImg),
                         MULTIDIM_ARRAY(re), MULTIDIM_ARRAY(im),
                         MULTIDIM_SIZE(oneImg));
        applyGeometryBSpline(rotre, A, re, Splinedegree, inv, wrap, outre);
        applyGeometryBSpline(rotim, A, im, Splinedegree, inv, wrap, outim);
        V2.resize(oneImg);
        RealImag2Complex(MULTIDIM_ARRAY(rotre), MULTIDIM_ARRAY(rotim),
                         MULTIDIM_ARRAY(V2), MULTIDIM_SIZE(re));
    }
    else
        applyGeometry(Splinedegree, V2, V1, A, inv, wrap, outside, n);
        

}

template<typename T>
void produceSplineCoefficients(int Splinedegree, 
                               MultidimArray< double > &coeffs,
                               const MultidimArray< T > &V1,  
                               unsigned long n)
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
template<>
void produceSplineCoefficients(int Splinedegree, 
                               MultidimArray< double > &coeffs,
                               const MultidimArray< std::complex<double> > &V1,  
                               unsigned long n)
{
    // TODO Implement
    REPORT_ERROR(222,"Spline coefficients of a complex matrix is not implemented.");
}

template<typename T>
void produceImageFromSplineCoefficients(int Splinedegree, 
                                        MultidimArray< double >& img, 
                                        const MultidimArray< double > &coeffs)
{

    img.initZeros(ZSIZE(coeffs), YSIZE(coeffs), XSIZE(coeffs));
    STARTINGX(img) = STARTINGX(coeffs);
    STARTINGY(img) = STARTINGY(coeffs);
    STARTINGZ(img) = STARTINGZ(coeffs);

    int Status;
    MultidimArray< double > aux;
    typeCast(coeffs, aux); // This will create a single volume!
        
    ChangeBasisVolume(MULTIDIM_ARRAY(aux), MULTIDIM_ARRAY(img),
                      XSIZE(coeffs), YSIZE(coeffs), ZSIZE(coeffs),
                      BasicSpline, CardinalSpline, SplineDegree,
                      MirrorOnBounds, DBL_EPSILON, &Status);
    if (Status)
        REPORT_ERROR(1201, "Error in ImageFromSplineCoefficients...");
    
}

template<typename T>
void rotate(int Splinedegree, 
            MultidimArray<T> &V2,
            const MultidimArray<T> &V1, 
            double ang, char axis, 
            bool wrap, T outside, unsigned long n)
{
    Matrix2D< double > tmp;
    if (ZSIZE(V1) == 1)
    {
        tmp = rotation2DMatrix(ang);
    }
    else
    {
        tmp = rotation3DMatrix(ang, axis);
    }
    
    applyGeometry(Splinedegree, V2, V1, tmp, IS_NOT_INV, wrap, outside, n);

}

template<typename T>
void rotate(int Splinedegree, 
            MultidimArray<T> &V2,
            const MultidimArray<T> &V1, 
            double ang, const Matrix1D< double >& axis, 
            bool wrap, T outside, unsigned long n)
{
    if (ZSIZE(V1)==1)
        REPORT_ERROR(1,"rotate with axis definition is only for 3D arrays");
    else
    {
        Matrix2D< double > tmp = rotation3DMatrix(ang, axis);
        applyGeometry(Splinedegree, V2, V1, tmp, IS_NOT_INV, wrap, outside, n);

    }

}

template<typename T>
void translate(int Splinedegree, 
               MultidimArray<T> &V2,
               const MultidimArray<T> &V1, 
               const Matrix1D< double >& v, 
               bool wrap, T outside, unsigned long n)
{
    Matrix2D< double > tmp;
    if (ZSIZE==1)
        tmp = translation2DMatrix(v);
    else
        tmp = translation3DMatrix(v);
    applyGeometry3D(Splinedegree, V2, V1, tmp, IS_NOT_INV, wrap, outside, n);
}

template<typename T>
void translateCenterOfMassToCenter(int Splinedegree, 
                                   MultidimArray<T> &V2,
                                   const MultidimArray<T> &V1, 
                                   bool wrap, unsigned long n)
{
    V2 = V1;
    V2.setXmippOrigin();
    Matrix1D< double > center;
    V2.centerOfMass(center,NULL,n);
    center *= -1;
    translate(Splinedegree, V2, V1, center, wrap, 0, n);
}

template<typename T>
void scaleToSize(int Splinedegree, 
                 MultidimArray<T> &V2,
                 const MultidimArray<T> &V1,
                 int Xdim, int Ydim, int Zdim,
                 unsigned long n)
{

    Matrix2D< double > tmp;
    if (ZSIZE(V1) == 1)
    {
        tmp.initIdentity(3);
        DIRECT_MAT_ELEM(temp, 0, 0) = (double) Xdim / (double) XSIZE(V1);
        DIRECT_MAT_ELEM(temp, 1, 1) = (double) Ydim / (double) YSIZE(V1);
        V2.resize(1, 1, Ydim, Xdim);
    }
    else
    {
        tmp.initIdentity(4);
        DIRECT_MAT_ELEM(tmp, 0, 0) = (double) Xdim / (double) XSIZE(V1);
        DIRECT_MAT_ELEM(tmp, 1, 1) = (double) Ydim / (double) YSIZE(V1);
        DIRECT_MAT_ELEM(tmp, 2, 2) = (double) Zdim / (double) ZSIZE(V1);
        V2.resize(1, Zdim, Ydim, Xdim);
    }
    applyGeometry(Splinedegree, V2, V1, tmp, IS_NOT_INV, WRAP, n);

}

// Special case for complex arrays
void scaleToSize(int Splinedegree, 
                 MultidimArray< std::complex<double> > &V2,
                 const MultidimArray< std::complex<double> > &V1,
                 int Xdim, int Ydim, int Zdim,
                 unsigned long n)
{
    if (Splinedegree > 1)
    {
        MultidimArray< double > re, im, scre, scim;
        MultidimArray<std::complex<double> > oneImg;

        re.resize(ZSIZE(V1), YSIZE(V1), XSIZE(V1));
        im.resize(ZSIZE(V1), YSIZE(V1), XSIZE(V1));

        V1.getImage(n, oneImg);
        Complex2RealImag(MULTIDIM_ARRAY(oneImg),
                         MULTIDIM_ARRAY(re), MULTIDIM_ARRAY(im),
                         MULTIDIM_SIZE(oneImg));

        re.scaleToSize(Splinedegree, Ydim, Xdim, Zdim, scre);
        im.scaleToSize(Splinedegree, Ydim, Xdim, Zdim, scim);

        result.resize(Zdim, Ydim, Xdim);

        RealImag2Complex(MULTIDIM_ARRAY(re), MULTIDIM_ARRAY(im),
                         MULTIDIM_ARRAY(V2), MULTIDIM_SIZE(re));
    }
    else
        scaleToSize(Splinedegree, V2, V1, Xdim, Ydim, Zdim, n);
    
}



template<typename T>
void pyramidReduce(int Splinedegree, 
                   MultidimArray<T> &V2,
                   const MultidimArray<T> &V1,
                   int levels, 
                   unsigned long n)
{
    MultidimArray< double > coeffs;
    produceSplineCoefficients(Splinedegree, coeffs, V1, n);

    for (int i = 0; i < levels; i++)
    {
        reduceBSpline(3, V2, coeffs);
        coeffs = V2;
    }

    produceImageFromSplineCoefficients(3, V2, coeffs);

}

void pyramidExpand(int Splinedegree, 
                   MultidimArray<T> &V2,
                   const MultidimArray<T> &V1,
                   int levels, 
                   unsigned long n)
{
    MultidimArray< double > coeffs;
    produceSplineCoefficients(Splinedegree, coeffs, V1, n);

    for (int i = 0; i < levels; i++)
    {
        expandBSpline(3, V2, coeffs);
        coeffs = V2;
    }

    produceImageFromSplineCoefficients(3, V2, coeffs);

}

template<typename T>
void reduceBSpline(int Splinedegree, 
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
    if (ZSIZE(V1) == 1)
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
    else
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

template<typename T>
void expandBSpline3D(int Splinedegree, 
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

    if (ZSIZE(V1) == 1)
    {
        V2.resize(2 * YSIZE(aux), 2 * XSIZE(aux));
        Expand_2D(MULTIDIM_ARRAY(aux), XSIZE(aux), YSIZE(aux),
                  MULTIDIM_ARRAY(V2), h, nh, IsCentered);
    }
    else
    {
        V2.resize(2 * ZSIZE(aux), 2 * YSIZE(aux), 2 * XSIZE(aux));
        Expand_3D(MULTIDIM_ARRAY(aux), XSIZE(aux), YSIZE(aux), ZSIZE(aux),
                  MULTIDIM_ARRAY(V2), h, nh, IsCentered);
    }

}

template<typename T>
void radialAverage(const MultidimArray< T >& m,
                   const Matrix1D< int >& center_of_rot,
                   MultidimArray< T >& radial_mean,
                   MultidimArray< int >& radial_count,
                   const bool& rounding,
                   unsigned long n)
{
    Matrix1D< double > idx(3);

    // If center_of_rot was written for 2D image
    if (XSIZE(center_of_rot) < 3)
        center_of_rot.resize(3);

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
        radial_mean(distance) += NZYX_ELEM(m, n, k, i, j);

        // Count the pixel
        radial_count(distance)++;
    }

    // Perform the mean
    FOR_ALL_ELEMENTS_IN_MATRIX1D(radial_mean)
        radial_mean(i) /= (T) radial_count(i);
}

/* Rotation 2D ------------------------------------------------------------- */
Matrix2D<double> rotation2DMatrix(double ang)
{
    Matrix2D<double> result(3, 3);
    rotation2DMatrix(ang, result);
    return result;
}

void rotation2DMatrix(double ang, Matrix2D< double > &result)
{
    double cosine, sine;

    ang = DEG2RAD(ang);
    cosine = cos(ang);
    sine = sin(ang);

    DIRECT_MAT_ELEM(result, 0, 0) = cosine;
    DIRECT_MAT_ELEM(result, 0, 1) = -sine;
    DIRECT_MAT_ELEM(result, 0, 2) = 0;

    DIRECT_MAT_ELEM(result, 1, 0) = sine;
    DIRECT_MAT_ELEM(result, 1, 1) = cosine;
    DIRECT_MAT_ELEM(result, 1, 2) = 0;

    DIRECT_MAT_ELEM(result, 2, 0) = 0;
    DIRECT_MAT_ELEM(result, 2, 1) = 0;
    DIRECT_MAT_ELEM(result, 2, 2) = 1;
}

/* Translation 2D ---------------------------------------------------------- */
Matrix2D<double> translation2DMatrix(const Matrix1D<double> &v)
{
    if (XSIZE(v) != 2)
        REPORT_ERROR(1002, "Translation2D_matrix: vector is not in R2");

    Matrix2D<double> result(3, 3);

    result.initIdentity();
    DIRECT_MAT_ELEM(result, 0, 2) = XX(v);
    DIRECT_MAT_ELEM(result, 1, 2) = YY(v);

    return result;
}

/* Rotation 3D around the system axes -------------------------------------- */
Matrix2D<double> rotation3DMatrix(double ang, char axis)
{
    Matrix2D<double> result(4, 4);
    double cosine, sine;

    ang = DEG2RAD(ang);
    cosine = cos(ang);
    sine = sin(ang);

    result.initZeros();
    DIRECT_MAT_ELEM(result, 3, 3) = 1;
    switch (axis)
    {
    case 'Z':
        DIRECT_MAT_ELEM(result, 0, 0) = cosine;
        DIRECT_MAT_ELEM(result, 0, 1) = -sine;
        DIRECT_MAT_ELEM(result, 1, 0) = sine;
        DIRECT_MAT_ELEM(result, 1, 1) = cosine;
        DIRECT_MAT_ELEM(result, 2, 2) = 1;
        break;
    case 'Y':
        DIRECT_MAT_ELEM(result, 0, 0) = cosine;
        DIRECT_MAT_ELEM(result, 0, 2) = -sine;
        DIRECT_MAT_ELEM(result, 2, 0) = sine;
        DIRECT_MAT_ELEM(result, 2, 2) = cosine;
        DIRECT_MAT_ELEM(result, 1, 1) = 1;
        break;
    case 'X':
        DIRECT_MAT_ELEM(result, 1, 1) = cosine;
        DIRECT_MAT_ELEM(result, 1, 2) = -sine;
        DIRECT_MAT_ELEM(result, 2, 1) = sine;
        DIRECT_MAT_ELEM(result, 2, 2) = cosine;
        DIRECT_MAT_ELEM(result, 0, 0) = 1;
        break;
    default:
        REPORT_ERROR(1105, "rotation3DMatrix: Unknown axis");
    }
    return result;
}

/* Align a vector with Z axis */
Matrix2D<double> alignWithZ(const Matrix1D<double> &axis)
{
    Matrix1D<double>  Axis;
    Matrix2D<double>  A(4, 4);

    if (XSIZE(axis) != 3)
        REPORT_ERROR(1002, "alignWithZ: Axis is not in R3");

    // Copy axis and compute length of the projection on YZ plane
    Axis = axis;
    Axis.selfNormalize();
    double proj_mod = sqrt(YY(Axis) * YY(Axis) + ZZ(Axis) * ZZ(Axis));

    A(3, 3) = 1;
    if (proj_mod > XMIPP_EQUAL_ACCURACY)
    { // proj_mod!=0
        // Build Matrix A, which makes the turning axis coincident with Z
        A(0, 0) = proj_mod;
        A(0, 1) = -XX(Axis) * YY(Axis) / proj_mod;
        A(0, 2) = -XX(Axis) * ZZ(Axis) / proj_mod;
        A(1, 0) = 0;
        A(1, 1) = ZZ(Axis) / proj_mod;
        A(1, 2) = -YY(Axis) / proj_mod;
        A(2, 0) = XX(Axis);
        A(2, 1) = YY(Axis);
        A(2, 2) = ZZ(Axis);
    }
    else
    {
        // I know that the Axis is the X axis
        A(0, 0) = 0;
        A(0, 1) = 0;
        A(0, 2) = -1;
        A(1, 0) = 0;
        A(1, 1) = 1;
        A(1, 2) = 0;
        A(2, 0) = 1;
        A(2, 1) = 0;
        A(2, 2) = 0;
    }
    return A;
}

/* Rotation 3D around any axis -------------------------------------------- */
Matrix2D<double> rotation3DMatrix(double ang, const Matrix1D<double> &axis)
{
    // Compute a matrix which makes the turning axis coincident with Z
    // And turn around this axis
    Matrix2D<double> A = alignWithZ(axis);
    return A.transpose() * rotation3DMatrix(ang, 'Z') * A;
}

/* Translation 3D ---------------------------------------------------------- */
Matrix2D<double> translation3DMatrix(const Matrix1D<double> &v)
{
    if (XSIZE(v) != 3)
        REPORT_ERROR(1002, "Translation3D_matrix: vector is not in R3");

    Matrix2D<double> result(4, 4);

    result.initIdentity();
    DIRECT_MAT_ELEM(result, 0, 3) = XX(v);
    DIRECT_MAT_ELEM(result, 1, 3) = YY(v);
    DIRECT_MAT_ELEM(result, 2, 3) = ZZ(v);

    return result;
}

/* Scale 3D ---------------------------------------------------------------- */
Matrix2D<double> scale3DMatrix(const Matrix1D<double> &sc)
{
    if (XSIZE(sc) != 3)
        REPORT_ERROR(1002, "Scale3D_matrix: vector is not in R3");

    Matrix2D<double> result(4, 4);

    result.initIdentity();
    DIRECT_MAT_ELEM(result, 0, 0) = XX(sc);
    DIRECT_MAT_ELEM(result, 1, 1) = YY(sc);
    DIRECT_MAT_ELEM(result, 2, 2) = ZZ(sc);

    return result;
}

