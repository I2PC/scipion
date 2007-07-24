/***************************************************************************
 *
 * Authors: Sjors H.W. Scheres (scheres@cnb.uam.es)
 *
 *  This code is strongly based on ideas by Pawel Penczek & Zhengfan
 *  Yang as implemented in SPARX at the University of Texas - Houston 
 *  Medical School
 *
 *  see P. A. Penczek, R. Renka, and H. Schomberg,
 *      J. Opt. Soc. Am. _21_, 449 (2004)
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
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/
#ifndef GRIDDING_H
#define GRIDDING_H
#include "funcs.h"
#include "matrix2d.h"
#include "matrix3d.h"
#include "fft.h"

#define GRIDDING_ALPHA 1.75
#define GRIDDING_NPAD 2
#define GRIDDING_K 6

/**  Produces a 2D Fourier-space matrix2d for reverse-gridding
 *  interpolation from a Fourier-space matrix2d
 */
void produceGriddingFourierMatrix2D(const Matrix2D< complex < double> > &in, 
				    Matrix2D< complex < double > > &out,
				    KaiserBessel &kb);

/**  Produces a 2D Fourier-space matrix2d for reverse-gridding
 *  interpolation from a real-space matrix2d 
 *  This real-space matrix2d should have the Xmipp origin set!
 */
void produceGriddingFourierMatrix2D(const Matrix2D< double > &in, 
				    Matrix2D< complex< double > > &out,
				    KaiserBessel &kb);

/**  Produces a 2D real-space matrix2d for reverse-gridding interpolation
 *  from a real-space matrix2d
 */
void produceGriddingMatrix2D(const Matrix2D< double > &in, 
			  Matrix2D< double > &out,
			  KaiserBessel &kb);

/**  Produces a 3D Fourier-space matrix3d for reverse-gridding
 *  interpolation from a Fourier-space matrix3d
 */
void produceGriddingFourierMatrix3D(const Matrix3D< complex < double > > &in, 
				  Matrix3D< complex < double > > &out,
				  KaiserBessel &kb);

/**  Produces a 3D Fourier-space matrix3d for reverse-gridding
 *  interpolation from a real-space matrix3d 
 *  This real-space matrix3d should have the Xmipp origin set!
 */
void produceGriddingFourierMatrix3D(const Matrix3D< double > &in, 
				    Matrix3D< complex< double > > &out,
				    KaiserBessel &kb);

/**  Produces a 3D real-space matrix3d for reverse-gridding
 *  interpolation from a real-space matrix3d
 */
void produceGriddingMatrix3D(const Matrix3D< double > &in, 
			   Matrix3D< double > &out,
			   KaiserBessel &kb);

/** Applies a geometric transformation to a 2D matrix with
 * gridding-based interpolation, knowing that this image has been
 * processed for gridding.
 *
 * A is a 3x3 transformation matrix.
 *
 * To interpolate using gridding you must prepare the image first!
 * An example to apply a gridding-based transformation would be:
 *
 * @code
 * KaiserBessel kb;
 * matrix2D<double> Maux,out;
 * produceGriddingMatrix2D(img(),Maux,kb);
 * Matrix2D<double> A = rotation2DMatrix(63.1);
 * applyGeometryGridding(out,A,Maux,IS_NOT_INV,DONT_WRAP,kb);
 * @endcode
 */
template<typename T>
void applyGeometryGridding(Matrix2D<T> &M2, Matrix2D< double > A, 
			   const Matrix2D<T> &M1, KaiserBessel &kb, 
			   bool inv, bool wrap, T outside = (T) 0)
{
    int m1, n1, m2, n2;
    double x, y, xp, yp;
    double minxp, minyp, maxxp, maxyp;
    int cen_x, cen_y, cen_xp, cen_yp;

    if ((XSIZE(A) != 3) || (YSIZE(A) != 3))
        REPORT_ERROR(1102, "Apply_geom: geometrical transformation is not 3x3");

    if (XSIZE(M1) == 0)
    {
        M2.clear();
        return;
    }

    // Smaller size output image for gridding!
    M2.resize(XSIZE(M1) / GRIDDING_NPAD, YSIZE(M1) / GRIDDING_NPAD);
    M2.setXmippOrigin();

    if (!inv)
        A = A.inv();

    // Find center and limits of image
    cen_y  = (int)(YSIZE(M2) / 2);
    cen_x  = (int)(XSIZE(M2) / 2);
    cen_yp = (int)(YSIZE(M1) / 2);
    cen_xp = (int)(XSIZE(M1) / 2);
    minxp  = -cen_xp;
    minyp  = -cen_yp;
    maxxp  = XSIZE(M1) - cen_xp - 1;
    maxyp  = YSIZE(M1) - cen_yp - 1;

    // Now we go from the output image to the input image, ie, for any pixel
    // in the output image we calculate which are the corresponding ones in
    // the original image, make an interpolation with them and put this value
    // at the output pixel
    for (int i = 0; i < YSIZE(M2); i++)
    {
        // Calculate position of the beginning of the row in the output image
        x = -cen_x;
        y = i - cen_y;

        // Calculate this position in the input image according to the
        // geometrical transformation
        // they are related by
        // coords_output(=x,y) = A * coords_input (=xp,yp)
        xp = x * dMij(A, 0, 0) + y * dMij(A, 0, 1) + dMij(A, 0, 2);
        yp = x * dMij(A, 1, 0) + y * dMij(A, 1, 1) + dMij(A, 1, 2);

        for (int j = 0; j < XSIZE(M2); j++)
        {
            bool interp;

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

            if (interp)
                dMij(M2, i, j) = interpolatedElementGridding(M1,xp,yp,kb);
	    else
		dMij(M2, i, j) = outside;

          // Compute new point inside input image
            xp += dMij(A, 0, 0);
            yp += dMij(A, 1, 0);
        }
    }
}

/** Applies a geometric transformation to a 3D matrix with
 * gridding-based interpolation, knowing that this image has been
 * processed for gridding.
 *
 * A is a 3x3 transformation matrix.
 *
 * To interpolate using gridding you must prepare the matrix3d first!
 * An example to apply a gridding-based transformation would be:
 *
 * @code
 * KaiserBessel kb;
 * matrix3D<double> Maux,out;
 * produceGriddingMatrix3D(vol(),Maux,kb);
 * Matrix2D<double> A = rotation2DMatrix(63.1);
 * applyGeometryGridding(out,A,Maux,IS_NOT_INV,DONT_WRAP,kb);
 * @endcode
 */
template<typename T>
void applyGeometryGridding(Matrix3D<T> &V2, Matrix2D< double > A, 
			   const Matrix3D<T> &V1, KaiserBessel &kb, 
			   bool inv, bool wrap, T outside = (T) 0)
{
    int m1, n1, o1, m2, n2, o2;
    double x, y, z, xp, yp, zp;
    double minxp, minyp, maxxp, maxyp, minzp, maxzp;
    int   cen_x, cen_y, cen_z, cen_xp, cen_yp, cen_zp;

    if ((XSIZE(A) != 4) || (YSIZE(A) != 4))
        REPORT_ERROR(1102, "Apply_geom3D: geometrical transformation is not 4x4");

    if (XSIZE(V1) == 0)
    {
        V2.clear();
        return;
    }

    // Smaller size output image for gridding!
    V2.resize(XSIZE(V1) / GRIDDING_NPAD, YSIZE(V1) / GRIDDING_NPAD, ZSIZE(V1) / GRIDDING_NPAD);
    V2.setXmippOrigin();

    if (!inv)
        A = A.inv();

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
            xp = x * dMij(A, 0, 0) + y * dMij(A, 0, 1) + z * dMij(A, 0, 2)
                 + dMij(A, 0, 3);
            yp = x * dMij(A, 1, 0) + y * dMij(A, 1, 1) + z * dMij(A, 1, 2)
                 + dMij(A, 1, 3);
            zp = x * dMij(A, 2, 0) + y * dMij(A, 2, 1) + z * dMij(A, 2, 2)
                 + dMij(A, 2, 3);

            for (int j = 0; j < V2.xdim; j++)
            {
                bool interp;

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
                    dVkij(V2, k, i, j) = (T) interpolatedElementGridding(V1,xp,yp,zp,kb);
                else
                    dVkij(V2, k, i, j) = outside;

                // Compute new point inside input image
                xp += dMij(A, 0, 0);
                yp += dMij(A, 1, 0);
                zp += dMij(A, 2, 0);
            }
        }
}

/** Interpolates the value of the 2D matrix M at the point (x,y) knowing
 * that this image has been processed for gridding
 *
 * (x,y) are in logical coordinates
 *
 * To interpolate using gridding you must prepare the image first!
 * An example to interpolate an image at (0.5,0.5) using
 * gridding would be:
 *
 * @code
 * KaiserBessel kb;
 * matrix2D<double> Maux;
 * produceGriddingMatrix2D(img(),Maux,kb);
 * interpolated_value = interpolatedElementGridding(Maux,0.5,0.5,kb);
 * @endcode
 */
template <typename T>
T interpolatedElementGridding(const Matrix2D<T> &in, float x, float y, KaiserBessel &kb)
{
    // size of this image:
    int nx = XSIZE(in);
    int ny = YSIZE(in);
    
    // Go from logical to physical coordinates in the small image
    x -= FIRST_XMIPP_INDEX(nx/GRIDDING_NPAD);
    y -= FIRST_XMIPP_INDEX(ny/GRIDDING_NPAD);
    
    // Get the convolution window parameters
    int K = kb.get_window_size();
    int kbmin = -K/2;
    int kbmax = -kbmin;
    int kbc = kbmax+1;
    float pixel =0.0f;
    float w=0.0f;
    
    x = fmod(2*x, float(nx));
    y = fmod(2*y, float(ny));
    int inxold = int(round(x));
    int inyold = int(round(y));
    
    float tablex1 = kb.i0win_tab(x-inxold+3);
    float tablex2 = kb.i0win_tab(x-inxold+2);
    float tablex3 = kb.i0win_tab(x-inxold+1);
    float tablex4 = kb.i0win_tab(x-inxold);
    float tablex5 = kb.i0win_tab(x-inxold-1);
    float tablex6 = kb.i0win_tab(x-inxold-2);
    float tablex7 = kb.i0win_tab(x-inxold-3);

    float tabley1 = kb.i0win_tab(y-inyold+3);
    float tabley2 = kb.i0win_tab(y-inyold+2);
    float tabley3 = kb.i0win_tab(y-inyold+1);
    float tabley4 = kb.i0win_tab(y-inyold);
    float tabley5 = kb.i0win_tab(y-inyold-1);
    float tabley6 = kb.i0win_tab(y-inyold-2);
    float tabley7 = kb.i0win_tab(y-inyold-3); 
	
    int x1, x2, x3, x4, x5, x6, x7, y1, y2, y3, y4, y5, y6, y7;
	
    if ( inxold <= kbc || inxold >=nx-kbc-2 || 
	 inyold <= kbc || inyold >=ny-kbc-2 )  {
	x1 = (inxold-3+nx)%nx;
	x2 = (inxold-2+nx)%nx;
	x3 = (inxold-1+nx)%nx;
	x4 = (inxold  +nx)%nx;
	x5 = (inxold+1+nx)%nx;
	x6 = (inxold+2+nx)%nx;
	x7 = (inxold+3+nx)%nx;

	y1 = ((inyold-3+ny)%ny)*nx;
	y2 = ((inyold-2+ny)%ny)*nx;
	y3 = ((inyold-1+ny)%ny)*nx;
	y4 = ((inyold  +ny)%ny)*nx;
	y5 = ((inyold+1+ny)%ny)*nx;
	y6 = ((inyold+2+ny)%ny)*nx;
	y7 = ((inyold+3+ny)%ny)*nx;
    } else {
	x1 = inxold-3;
	x2 = inxold-2;
	x3 = inxold-1;
	x4 = inxold;
	x5 = inxold+1;
	x6 = inxold+2;
	x7 = inxold+3;

	y1 = (inyold-3)*nx;
	y2 = (inyold-2)*nx;
	y3 = (inyold-1)*nx;
	y4 = inyold*nx;
	y5 = (inyold+1)*nx;
	y6 = (inyold+2)*nx;
	y7 = (inyold+3)*nx;
    }
    
    // The actual convolution
    pixel = ( in.data[x1+y1]*tablex1 + in.data[x2+y1]*tablex2 + in.data[x3+y1]*tablex3 +
	      in.data[x4+y1]*tablex4 + in.data[x5+y1]*tablex5 + in.data[x6+y1]*tablex6 +
	      in.data[x7+y1]*tablex7 ) * tabley1 +

	    ( in.data[x1+y2]*tablex1 + in.data[x2+y2]*tablex2 + in.data[x3+y2]*tablex3 +
	      in.data[x4+y2]*tablex4 + in.data[x5+y2]*tablex5 + in.data[x6+y2]*tablex6 +
	      in.data[x7+y2]*tablex7 ) * tabley2 +

	    ( in.data[x1+y3]*tablex1 + in.data[x2+y3]*tablex2 + in.data[x3+y3]*tablex3 +
	      in.data[x4+y3]*tablex4 + in.data[x5+y3]*tablex5 + in.data[x6+y3]*tablex6 +
	      in.data[x7+y3]*tablex7 ) * tabley3 +

	    ( in.data[x1+y4]*tablex1 + in.data[x2+y4]*tablex2 + in.data[x3+y4]*tablex3 +
	      in.data[x4+y4]*tablex4 + in.data[x5+y4]*tablex5 + in.data[x6+y4]*tablex6 +
	      in.data[x7+y4]*tablex7 ) * tabley4 +

	    ( in.data[x1+y5]*tablex1 + in.data[x2+y5]*tablex2 + in.data[x3+y5]*tablex3 +
	      in.data[x4+y5]*tablex4 + in.data[x5+y5]*tablex5 + in.data[x6+y5]*tablex6 +
	      in.data[x7+y5]*tablex7 ) * tabley5 +

	    ( in.data[x1+y6]*tablex1 + in.data[x2+y6]*tablex2 + in.data[x3+y6]*tablex3 +
	      in.data[x4+y6]*tablex4 + in.data[x5+y6]*tablex5 + in.data[x6+y6]*tablex6 +
	      in.data[x7+y6]*tablex7 ) * tabley6 +
	
	    ( in.data[x1+y7]*tablex1 + in.data[x2+y7]*tablex2 + in.data[x3+y7]*tablex3 +
	      in.data[x4+y7]*tablex4 + in.data[x5+y7]*tablex5 + in.data[x6+y7]*tablex6 +
	      in.data[x7+y7]*tablex7 ) * tabley7;
    
    w = (tablex1+tablex2+tablex3+tablex4+tablex5+tablex6+tablex7) *
	(tabley1+tabley2+tabley3+tabley4+tabley5+tabley6+tabley7);	
    
    return pixel/w;

}

/** Interpolates the value of the 3D matrix M at the point (x,y,z) knowing
 * that this matrix has been processed for gridding
 *
 * (x,y,z) are in logical coordinates
 *
 * To interpolate using gridding you must prepare the matrix first!
 * An example to interpolate an image at (0.5,0.5,0.5) using
 * gridding would be:
 *
 * @code
 * KaiserBessel kb;
 * matrix3D<double> Maux;
 * produceGriddingMatrix3D(vol(),Maux,kb);
 * interpolated_value = interpolatedElementGridding(Maux,0.5,0.5,0.5,kb);
 * @endcode
 */
template <typename T>
T interpolatedElementGridding(const Matrix3D<T> &in, float x, float y, float z, KaiserBessel &kb)
{
    // size of this image:
    int nx = XSIZE(in);
    int ny = YSIZE(in);
    int nz = ZSIZE(in);
    
    // Go from logical to physical coordinates in the small image
    x -= FIRST_XMIPP_INDEX(nx/GRIDDING_NPAD);
    y -= FIRST_XMIPP_INDEX(ny/GRIDDING_NPAD);
    z -= FIRST_XMIPP_INDEX(nz/GRIDDING_NPAD);
    
    // Get the convolution window parameters
    int K = kb.get_window_size();
    int kbmin = -K/2;
    int kbmax = -kbmin;
    int kbc = kbmax+1;
    float pixel =0.0f;
    float w=0.0f;
    
    x = fmod(2*x, float(nx));
    y = fmod(2*y, float(ny));
    z = fmod(2*z, float(nz));
    int inxold = int(round(x));
    int inyold = int(round(y));
    int inzold = int(round(z));
    
    float tablex1 = kb.i0win_tab(x-inxold+3);
    float tablex2 = kb.i0win_tab(x-inxold+2);
    float tablex3 = kb.i0win_tab(x-inxold+1);
    float tablex4 = kb.i0win_tab(x-inxold);
    float tablex5 = kb.i0win_tab(x-inxold-1);
    float tablex6 = kb.i0win_tab(x-inxold-2);
    float tablex7 = kb.i0win_tab(x-inxold-3);

    float tabley1 = kb.i0win_tab(y-inyold+3);
    float tabley2 = kb.i0win_tab(y-inyold+2);
    float tabley3 = kb.i0win_tab(y-inyold+1);
    float tabley4 = kb.i0win_tab(y-inyold);
    float tabley5 = kb.i0win_tab(y-inyold-1);
    float tabley6 = kb.i0win_tab(y-inyold-2);
    float tabley7 = kb.i0win_tab(y-inyold-3); 
	
    float tablez1 = kb.i0win_tab(z-inzold+3);
    float tablez2 = kb.i0win_tab(z-inzold+2);
    float tablez3 = kb.i0win_tab(z-inzold+1);
    float tablez4 = kb.i0win_tab(z-inzold);
    float tablez5 = kb.i0win_tab(z-inzold-1);
    float tablez6 = kb.i0win_tab(z-inzold-2);
    float tablez7 = kb.i0win_tab(z-inzold-3); 

    int x1, x2, x3, x4, x5, x6, x7, y1, y2, y3, y4, y5, y6, y7, z1, z2, z3, z4, z5, z6, z7;
	
    if ( inxold <= kbc || inxold >=nx-kbc-2 || 
	 inyold <= kbc || inyold >=ny-kbc-2 || 
	 inzold <= kbc || inzold >= nz-kbc-2 )  {
	x1 = (inxold-3+nx)%nx;
	x2 = (inxold-2+nx)%nx;
	x3 = (inxold-1+nx)%nx;
	x4 = (inxold  +nx)%nx;
	x5 = (inxold+1+nx)%nx;
	x6 = (inxold+2+nx)%nx;
	x7 = (inxold+3+nx)%nx;

	y1 = ((inyold-3+ny)%ny)*nx;
	y2 = ((inyold-2+ny)%ny)*nx;
	y3 = ((inyold-1+ny)%ny)*nx;
	y4 = ((inyold  +ny)%ny)*nx;
	y5 = ((inyold+1+ny)%ny)*nx;
	y6 = ((inyold+2+ny)%ny)*nx;
	y7 = ((inyold+3+ny)%ny)*nx;

	z1 = ((inzold-3+nz)%nz)*nx*ny;
	z2 = ((inzold-2+nz)%nz)*nx*ny;
	z3 = ((inzold-1+nz)%nz)*nx*ny;
	z4 = ((inzold  +nz)%nz)*nx*ny;
	z5 = ((inzold+1+nz)%nz)*nx*ny;
	z6 = ((inzold+2+nz)%nz)*nx*ny;
	z7 = ((inzold+3+nz)%nz)*nx*ny;
    } else {
	x1 = inxold-3;
	x2 = inxold-2;
	x3 = inxold-1;
	x4 = inxold;
	x5 = inxold+1;
	x6 = inxold+2;
	x7 = inxold+3;

	y1 = (inyold-3)*nx;
	y2 = (inyold-2)*nx;
	y3 = (inyold-1)*nx;
	y4 = inyold*nx;
	y5 = (inyold+1)*nx;
	y6 = (inyold+2)*nx;
	y7 = (inyold+3)*nx;

	z1 = (inzold-3)*nx*ny;
	z2 = (inzold-2)*nx*ny;
	z3 = (inzold-1)*nx*ny;
	z4 = inzold*nx*ny;
	z5 = (inzold+1)*nx*ny;
	z6 = (inzold+2)*nx*ny;		
	z7 = (inzold+3)*nx*ny;
    }
    
    // The actual convolution
    pixel  = ( ( in.data[x1+y1+z1]*tablex1 + in.data[x2+y1+z1]*tablex2 + in.data[x3+y1+z1]*tablex3 +
		 in.data[x4+y1+z1]*tablex4 + in.data[x5+y1+z1]*tablex5 + in.data[x6+y1+z1]*tablex6 +
		 in.data[x7+y1+z1]*tablex7 ) * tabley1 +
	       ( in.data[x1+y2+z1]*tablex1 + in.data[x2+y2+z1]*tablex2 + in.data[x3+y2+z1]*tablex3 +
		 in.data[x4+y2+z1]*tablex4 + in.data[x5+y2+z1]*tablex5 + in.data[x6+y2+z1]*tablex6 +
		 in.data[x7+y2+z1]*tablex7 ) * tabley2 +
	       ( in.data[x1+y3+z1]*tablex1 + in.data[x2+y3+z1]*tablex2 + in.data[x3+y3+z1]*tablex3 +
		 in.data[x4+y3+z1]*tablex4 + in.data[x5+y3+z1]*tablex5 + in.data[x6+y3+z1]*tablex6 +
		 in.data[x7+y3+z1]*tablex7 ) * tabley3 +
	       ( in.data[x1+y4+z1]*tablex1 + in.data[x2+y4+z1]*tablex2 + in.data[x3+y4+z1]*tablex3 +
		 in.data[x4+y4+z1]*tablex4 + in.data[x5+y4+z1]*tablex5 + in.data[x6+y4+z1]*tablex6 +
		 in.data[x7+y4+z1]*tablex7 ) * tabley4 +
	       ( in.data[x1+y5+z1]*tablex1 + in.data[x2+y5+z1]*tablex2 + in.data[x3+y5+z1]*tablex3 +
		 in.data[x4+y5+z1]*tablex4 + in.data[x5+y5+z1]*tablex5 + in.data[x6+y5+z1]*tablex6 +
		 in.data[x7+y5+z1]*tablex7 ) * tabley5 +
	       ( in.data[x1+y6+z1]*tablex1 + in.data[x2+y6+z1]*tablex2 + in.data[x3+y6+z1]*tablex3 +
		 in.data[x4+y6+z1]*tablex4 + in.data[x5+y6+z1]*tablex5 + in.data[x6+y6+z1]*tablex6 +
		 in.data[x7+y6+z1]*tablex7 ) * tabley6 +
	       ( in.data[x1+y7+z1]*tablex1 + in.data[x2+y7+z1]*tablex2 + in.data[x3+y7+z1]*tablex3 +
		 in.data[x4+y7+z1]*tablex4 + in.data[x5+y7+z1]*tablex5 + in.data[x6+y7+z1]*tablex6 +
		 in.data[x7+y7+z1]*tablex7 ) * tabley7 ) *tablez1 +

	     ( ( in.data[x1+y1+z2]*tablex1 + in.data[x2+y1+z2]*tablex2 + in.data[x3+y1+z2]*tablex3 +
		 in.data[x4+y1+z2]*tablex4 + in.data[x5+y1+z2]*tablex5 + in.data[x6+y1+z2]*tablex6 +
		 in.data[x7+y1+z2]*tablex7 ) * tabley1 +
	       ( in.data[x1+y2+z2]*tablex1 + in.data[x2+y2+z2]*tablex2 + in.data[x3+y2+z2]*tablex3 +
		 in.data[x4+y2+z2]*tablex4 + in.data[x5+y2+z2]*tablex5 + in.data[x6+y2+z2]*tablex6 +
		 in.data[x7+y2+z2]*tablex7 ) * tabley2 +
	       ( in.data[x1+y3+z2]*tablex1 + in.data[x2+y3+z2]*tablex2 + in.data[x3+y3+z2]*tablex3 +
		 in.data[x4+y3+z2]*tablex4 + in.data[x5+y3+z2]*tablex5 + in.data[x6+y3+z2]*tablex6 +
		 in.data[x7+y3+z2]*tablex7 ) * tabley3 +
	       ( in.data[x1+y4+z2]*tablex1 + in.data[x2+y4+z2]*tablex2 + in.data[x3+y4+z2]*tablex3 +
		 in.data[x4+y4+z2]*tablex4 + in.data[x5+y4+z2]*tablex5 + in.data[x6+y4+z2]*tablex6 +
		 in.data[x7+y4+z2]*tablex7 ) * tabley4 +
	       ( in.data[x1+y5+z2]*tablex1 + in.data[x2+y5+z2]*tablex2 + in.data[x3+y5+z2]*tablex3 +
		 in.data[x4+y5+z2]*tablex4 + in.data[x5+y5+z2]*tablex5 + in.data[x6+y5+z2]*tablex6 +
		 in.data[x7+y5+z2]*tablex7 ) * tabley5 +
	       ( in.data[x1+y6+z2]*tablex1 + in.data[x2+y6+z2]*tablex2 + in.data[x3+y6+z2]*tablex3 +
		 in.data[x4+y6+z2]*tablex4 + in.data[x5+y6+z2]*tablex5 + in.data[x6+y6+z2]*tablex6 +
		 in.data[x7+y6+z2]*tablex7 ) * tabley6 +
	       ( in.data[x1+y7+z2]*tablex1 + in.data[x2+y7+z2]*tablex2 + in.data[x3+y7+z2]*tablex3 +
		 in.data[x4+y7+z2]*tablex4 + in.data[x5+y7+z2]*tablex5 + in.data[x6+y7+z2]*tablex6 +
		 in.data[x7+y7+z2]*tablex7 ) * tabley7 ) *tablez2 +

       	     ( ( in.data[x1+y1+z3]*tablex1 + in.data[x2+y1+z3]*tablex2 + in.data[x3+y1+z3]*tablex3 +
		 in.data[x4+y1+z3]*tablex4 + in.data[x5+y1+z3]*tablex5 + in.data[x6+y1+z3]*tablex6 +
		 in.data[x7+y1+z3]*tablex7 ) * tabley1 +
	       ( in.data[x1+y2+z3]*tablex1 + in.data[x2+y2+z3]*tablex2 + in.data[x3+y2+z3]*tablex3 +
		 in.data[x4+y2+z3]*tablex4 + in.data[x5+y2+z3]*tablex5 + in.data[x6+y2+z3]*tablex6 +
		 in.data[x7+y2+z3]*tablex7 ) * tabley2 +
	       ( in.data[x1+y3+z3]*tablex1 + in.data[x2+y3+z3]*tablex2 + in.data[x3+y3+z3]*tablex3 +
		 in.data[x4+y3+z3]*tablex4 + in.data[x5+y3+z3]*tablex5 + in.data[x6+y3+z3]*tablex6 +
		 in.data[x7+y3+z3]*tablex7 ) * tabley3 +
	       ( in.data[x1+y4+z3]*tablex1 + in.data[x2+y4+z3]*tablex2 + in.data[x3+y4+z3]*tablex3 +
		 in.data[x4+y4+z3]*tablex4 + in.data[x5+y4+z3]*tablex5 + in.data[x6+y4+z3]*tablex6 +
		 in.data[x7+y4+z3]*tablex7 ) * tabley4 +
	       ( in.data[x1+y5+z3]*tablex1 + in.data[x2+y5+z3]*tablex2 + in.data[x3+y5+z3]*tablex3 +
		 in.data[x4+y5+z3]*tablex4 + in.data[x5+y5+z3]*tablex5 + in.data[x6+y5+z3]*tablex6 +
		 in.data[x7+y5+z3]*tablex7 ) * tabley5 +
	       ( in.data[x1+y6+z3]*tablex1 + in.data[x2+y6+z3]*tablex2 + in.data[x3+y6+z3]*tablex3 +
		 in.data[x4+y6+z3]*tablex4 + in.data[x5+y6+z3]*tablex5 + in.data[x6+y6+z3]*tablex6 +
		 in.data[x7+y6+z3]*tablex7 ) * tabley6 +
	       ( in.data[x1+y7+z3]*tablex1 + in.data[x2+y7+z3]*tablex2 + in.data[x3+y7+z3]*tablex3 +
		 in.data[x4+y7+z3]*tablex4 + in.data[x5+y7+z3]*tablex5 + in.data[x6+y7+z3]*tablex6 +
		 in.data[x7+y7+z3]*tablex7 ) * tabley7 ) *tablez3 +

	     ( ( in.data[x1+y1+z4]*tablex1 + in.data[x2+y1+z4]*tablex2 + in.data[x3+y1+z4]*tablex3 +
		 in.data[x4+y1+z4]*tablex4 + in.data[x5+y1+z4]*tablex5 + in.data[x6+y1+z4]*tablex6 +
		 in.data[x7+y1+z4]*tablex7 ) * tabley1 +
	       ( in.data[x1+y2+z4]*tablex1 + in.data[x2+y2+z4]*tablex2 + in.data[x3+y2+z4]*tablex3 +
		 in.data[x4+y2+z4]*tablex4 + in.data[x5+y2+z4]*tablex5 + in.data[x6+y2+z4]*tablex6 +
		 in.data[x7+y2+z4]*tablex7 ) * tabley2 +
	       ( in.data[x1+y3+z4]*tablex1 + in.data[x2+y3+z4]*tablex2 + in.data[x3+y3+z4]*tablex3 +
		 in.data[x4+y3+z4]*tablex4 + in.data[x5+y3+z4]*tablex5 + in.data[x6+y3+z4]*tablex6 +
		 in.data[x7+y3+z4]*tablex7 ) * tabley3 +
	       ( in.data[x1+y4+z4]*tablex1 + in.data[x2+y4+z4]*tablex2 + in.data[x3+y4+z4]*tablex3 +
		 in.data[x4+y4+z4]*tablex4 + in.data[x5+y4+z4]*tablex5 + in.data[x6+y4+z4]*tablex6 +
		 in.data[x7+y4+z4]*tablex7 ) * tabley4 +
	       ( in.data[x1+y5+z4]*tablex1 + in.data[x2+y5+z4]*tablex2 + in.data[x3+y5+z4]*tablex3 +
		 in.data[x4+y5+z4]*tablex4 + in.data[x5+y5+z4]*tablex5 + in.data[x6+y5+z4]*tablex6 +
		 in.data[x7+y5+z4]*tablex7 ) * tabley5 +
	       ( in.data[x1+y6+z4]*tablex1 + in.data[x2+y6+z4]*tablex2 + in.data[x3+y6+z4]*tablex3 +
		 in.data[x4+y6+z4]*tablex4 + in.data[x5+y6+z4]*tablex5 + in.data[x6+y6+z4]*tablex6 +
		 in.data[x7+y6+z4]*tablex7 ) * tabley6 +
	       ( in.data[x1+y7+z4]*tablex1 + in.data[x2+y7+z4]*tablex2 + in.data[x3+y7+z4]*tablex3 +
		 in.data[x4+y7+z4]*tablex4 + in.data[x5+y7+z4]*tablex5 + in.data[x6+y7+z4]*tablex6 +
		 in.data[x7+y7+z4]*tablex7 ) * tabley7 ) *tablez4 +
	
	     ( ( in.data[x1+y1+z5]*tablex1 + in.data[x2+y1+z5]*tablex2 + in.data[x3+y1+z5]*tablex3 +
		 in.data[x4+y1+z5]*tablex4 + in.data[x5+y1+z5]*tablex5 + in.data[x6+y1+z5]*tablex6 +
		 in.data[x7+y1+z5]*tablex7 ) * tabley1 +
	       ( in.data[x1+y2+z5]*tablex1 + in.data[x2+y2+z5]*tablex2 + in.data[x3+y2+z5]*tablex3 +
		 in.data[x4+y2+z5]*tablex4 + in.data[x5+y2+z5]*tablex5 + in.data[x6+y2+z5]*tablex6 +
		 in.data[x7+y2+z5]*tablex7 ) * tabley2 +
	       ( in.data[x1+y3+z5]*tablex1 + in.data[x2+y3+z5]*tablex2 + in.data[x3+y3+z5]*tablex3 +
		 in.data[x4+y3+z5]*tablex4 + in.data[x5+y3+z5]*tablex5 + in.data[x6+y3+z5]*tablex6 +
		 in.data[x7+y3+z5]*tablex7 ) * tabley3 +
	       ( in.data[x1+y4+z5]*tablex1 + in.data[x2+y4+z5]*tablex2 + in.data[x3+y4+z5]*tablex3 +
		 in.data[x4+y4+z5]*tablex4 + in.data[x5+y4+z5]*tablex5 + in.data[x6+y4+z5]*tablex6 +
		 in.data[x7+y4+z5]*tablex7 ) * tabley4 +
	       ( in.data[x1+y5+z5]*tablex1 + in.data[x2+y5+z5]*tablex2 + in.data[x3+y5+z5]*tablex3 +
		 in.data[x4+y5+z5]*tablex4 + in.data[x5+y5+z5]*tablex5 + in.data[x6+y5+z5]*tablex6 +
		 in.data[x7+y5+z5]*tablex7 ) * tabley5 +
	       ( in.data[x1+y6+z5]*tablex1 + in.data[x2+y6+z5]*tablex2 + in.data[x3+y6+z5]*tablex3 +
		 in.data[x4+y6+z5]*tablex4 + in.data[x5+y6+z5]*tablex5 + in.data[x6+y6+z5]*tablex6 +
		 in.data[x7+y6+z5]*tablex7 ) * tabley6 +
	       ( in.data[x1+y7+z5]*tablex1 + in.data[x2+y7+z5]*tablex2 + in.data[x3+y7+z5]*tablex3 +
		 in.data[x4+y7+z5]*tablex4 + in.data[x5+y7+z5]*tablex5 + in.data[x6+y7+z5]*tablex6 +
		 in.data[x7+y7+z5]*tablex7 ) * tabley7 ) *tablez5 +

	     ( ( in.data[x1+y1+z6]*tablex1 + in.data[x2+y1+z6]*tablex2 + in.data[x3+y1+z6]*tablex3 +
		 in.data[x4+y1+z6]*tablex4 + in.data[x5+y1+z6]*tablex5 + in.data[x6+y1+z6]*tablex6 +
		 in.data[x7+y1+z6]*tablex7 ) * tabley1 +
	       ( in.data[x1+y2+z6]*tablex1 + in.data[x2+y2+z6]*tablex2 + in.data[x3+y2+z6]*tablex3 +
		 in.data[x4+y2+z6]*tablex4 + in.data[x5+y2+z6]*tablex5 + in.data[x6+y2+z6]*tablex6 +
		 in.data[x7+y2+z6]*tablex7 ) * tabley2 +
	       ( in.data[x1+y3+z6]*tablex1 + in.data[x2+y3+z6]*tablex2 + in.data[x3+y3+z6]*tablex3 +
		 in.data[x4+y3+z6]*tablex4 + in.data[x5+y3+z6]*tablex5 + in.data[x6+y3+z6]*tablex6 +
		 in.data[x7+y3+z6]*tablex7 ) * tabley3 +
	       ( in.data[x1+y4+z6]*tablex1 + in.data[x2+y4+z6]*tablex2 + in.data[x3+y4+z6]*tablex3 +
		 in.data[x4+y4+z6]*tablex4 + in.data[x5+y4+z6]*tablex5 + in.data[x6+y4+z6]*tablex6 +
		 in.data[x7+y4+z6]*tablex7 ) * tabley4 +
	       ( in.data[x1+y5+z6]*tablex1 + in.data[x2+y5+z6]*tablex2 + in.data[x3+y5+z6]*tablex3 +
		 in.data[x4+y5+z6]*tablex4 + in.data[x5+y5+z6]*tablex5 + in.data[x6+y5+z6]*tablex6 +
		 in.data[x7+y5+z6]*tablex7 ) * tabley5 +
	       ( in.data[x1+y6+z6]*tablex1 + in.data[x2+y6+z6]*tablex2 + in.data[x3+y6+z6]*tablex3 +
		 in.data[x4+y6+z6]*tablex4 + in.data[x5+y6+z6]*tablex5 + in.data[x6+y6+z6]*tablex6 +
		 in.data[x7+y6+z6]*tablex7 ) * tabley6 +
	       ( in.data[x1+y7+z6]*tablex1 + in.data[x2+y7+z6]*tablex2 + in.data[x3+y7+z6]*tablex3 +
		 in.data[x4+y7+z6]*tablex4 + in.data[x5+y7+z6]*tablex5 + in.data[x6+y7+z6]*tablex6 +
		 in.data[x7+y7+z6]*tablex7 ) * tabley7 ) *tablez6 +

	     ( ( in.data[x1+y1+z7]*tablex1 + in.data[x2+y1+z7]*tablex2 + in.data[x3+y1+z7]*tablex3 +
		 in.data[x4+y1+z7]*tablex4 + in.data[x5+y1+z7]*tablex5 + in.data[x6+y1+z7]*tablex6 +
		 in.data[x7+y1+z7]*tablex7 ) * tabley1 +
	       ( in.data[x1+y2+z7]*tablex1 + in.data[x2+y2+z7]*tablex2 + in.data[x3+y2+z7]*tablex3 +
		 in.data[x4+y2+z7]*tablex4 + in.data[x5+y2+z7]*tablex5 + in.data[x6+y2+z7]*tablex6 +
		 in.data[x7+y2+z7]*tablex7 ) * tabley2 +
	       ( in.data[x1+y3+z7]*tablex1 + in.data[x2+y3+z7]*tablex2 + in.data[x3+y3+z7]*tablex3 +
		 in.data[x4+y3+z7]*tablex4 + in.data[x5+y3+z7]*tablex5 + in.data[x6+y3+z7]*tablex6 +
		 in.data[x7+y3+z7]*tablex7 ) * tabley3 +
	       ( in.data[x1+y4+z7]*tablex1 + in.data[x2+y4+z7]*tablex2 + in.data[x3+y4+z7]*tablex3 +
		 in.data[x4+y4+z7]*tablex4 + in.data[x5+y4+z7]*tablex5 + in.data[x6+y4+z7]*tablex6 +
		 in.data[x7+y4+z7]*tablex7 ) * tabley4 +
	       ( in.data[x1+y5+z7]*tablex1 + in.data[x2+y5+z7]*tablex2 + in.data[x3+y5+z7]*tablex3 +
		 in.data[x4+y5+z7]*tablex4 + in.data[x5+y5+z7]*tablex5 + in.data[x6+y5+z7]*tablex6 +
		 in.data[x7+y5+z7]*tablex7 ) * tabley5 +
	       ( in.data[x1+y6+z7]*tablex1 + in.data[x2+y6+z7]*tablex2 + in.data[x3+y6+z7]*tablex3 +
		 in.data[x4+y6+z7]*tablex4 + in.data[x5+y6+z7]*tablex5 + in.data[x6+y6+z7]*tablex6 +
		 in.data[x7+y6+z7]*tablex7 ) * tabley6 +
	       ( in.data[x1+y7+z7]*tablex1 + in.data[x2+y7+z7]*tablex2 + in.data[x3+y7+z7]*tablex3 +
		 in.data[x4+y7+z7]*tablex4 + in.data[x5+y7+z7]*tablex5 + in.data[x6+y7+z7]*tablex6 +
		 in.data[x7+y7+z7]*tablex7 ) * tabley7 ) *tablez7;
    
    w = (tablex1+tablex2+tablex3+tablex4+tablex5+tablex6+tablex7) *
	(tabley1+tabley2+tabley3+tabley4+tabley5+tabley6+tabley7) *
	(tablez1+tablez2+tablez3+tablez4+tablez5+tablez6+tablez7);	
    
    return pixel/w;

}

#endif
