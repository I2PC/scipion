/***************************************************************************
 *
 * Authors: Sjors H.W. Scheres (scheres@cnb.uam.es)
 *
 *  Part of this code is strongly based on ideas by Pawel Penczek & Zhengfan
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

/// @defgroup Gridding Gridding
/// @ingroup DataLibrary

/// @defgroup ReverseGridding Reverse Gridding
/// @ingroup Gridding

// ***************************************************************************
// ************************ Reverse Gridding *********************************
// ***************************************************************************

/** Produce a complex Matrix2D for reverse gridding from a complex Matrix2D
 * @ingroup ReverseGridding
 *
 *  Produces a 2D Fourier-space matrix2d for reverse-gridding
 *  interpolation from a Fourier-space matrix2d
 */
void produceReverseGriddingFourierMatrix2D(const Matrix2D< std::complex < double> > &in, 
					   Matrix2D< std::complex < double > > &out,
					   KaiserBessel &kb);

/** Produce a complex Matrix2D for reverse gridding from a real Matrix2D
 * @ingroup ReverseGridding
 *
 *  Produces a 2D Fourier-space matrix2d for reverse-gridding
 *  interpolation from a real-space matrix2d 
 *  This real-space matrix2d should have the Xmipp origin set!
 */
void produceReverseGriddingFourierMatrix2D(const Matrix2D< double > &in, 
					   Matrix2D< std::complex< double > > &out,
					   KaiserBessel &kb);

/** Produce a real Matrix2D for reverse gridding from a real Matrix2D
 * @ingroup ReverseGridding
 *
 *  Produces a 2D real-space matrix2d for reverse-gridding interpolation
 *  from a real-space matrix2d
 */
void produceReverseGriddingMatrix2D(const Matrix2D< double > &in, 
				    Matrix2D< double > &out,
				    KaiserBessel &kb);

/** Produce a complex Matrix3D for reverse gridding from a complex Matrix2D
 * @ingroup ReverseGridding
 *
 *  Produces a 3D Fourier-space matrix3d for reverse-gridding
 *  interpolation from a Fourier-space matrix3d
 */
void produceReverseGriddingFourierMatrix3D(const Matrix3D< std::complex < double > > &in, 
					   Matrix3D< std::complex < double > > &out,
					   KaiserBessel &kb);

/** Produce a complex Matrix3D for reverse gridding from a real Matrix2D
 * @ingroup ReverseGridding
 *
 *  Produces a 3D Fourier-space matrix3d for reverse-gridding
 *  interpolation from a real-space matrix3d 
 *  This real-space matrix3d should have the Xmipp origin set!
 */
void produceReverseGriddingFourierMatrix3D(const Matrix3D< double > &in, 
					   Matrix3D< std::complex< double > > &out,
					   KaiserBessel &kb);

/** Produce a real Matrix3D for reverse gridding from a real Matrix2D
 * @ingroup ReverseGridding
 *
 *  Produces a 3D real-space matrix3d for reverse-gridding
 *  interpolation from a real-space matrix3d
 */
void produceReverseGriddingMatrix3D(const Matrix3D< double > &in, 
				    Matrix3D< double > &out,
				    KaiserBessel &kb);

/** Reverse-gridding based interpolation in a 2D matrix
 * @ingroup ReverseGridding
 *
 * Interpolates the value of the 2D matrix M at the point (x,y) knowing
 * that this image has been processed for reverse-gridding
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
 * produceReverseGriddingMatrix2D(img(),Maux,kb);
 * interpolated_value = interpolatedElementReverseGridding(Maux,0.5,0.5,kb);
 * @endcode
 */
template <typename T>
T interpolatedElementReverseGridding(const Matrix2D<T> &in, double x, double y, const KaiserBessel &kb)
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
    T pixel = 0.;
    double w = 0.;
    
    x = fmod(2*x, double(nx));
    y = fmod(2*y, double(ny));
    int inxold = int(ROUND(x));
    int inyold = int(ROUND(y));
    
    double tablex1 = kb.i0win_tab(x-inxold+3);
    double tablex2 = kb.i0win_tab(x-inxold+2);
    double tablex3 = kb.i0win_tab(x-inxold+1);
    double tablex4 = kb.i0win_tab(x-inxold);
    double tablex5 = kb.i0win_tab(x-inxold-1);
    double tablex6 = kb.i0win_tab(x-inxold-2);
    double tablex7 = kb.i0win_tab(x-inxold-3);

    double tabley1 = kb.i0win_tab(y-inyold+3);
    double tabley2 = kb.i0win_tab(y-inyold+2);
    double tabley3 = kb.i0win_tab(y-inyold+1);
    double tabley4 = kb.i0win_tab(y-inyold);
    double tabley5 = kb.i0win_tab(y-inyold-1);
    double tabley6 = kb.i0win_tab(y-inyold-2);
    double tabley7 = kb.i0win_tab(y-inyold-3); 
	
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

/** Reverse-gridding based interpolation in a 3D matrix
 * @ingroup ReverseGridding
 *
 * Interpolates the value of the 3D matrix M at the point (x,y,z) knowing
 * that this matrix has been processed for reverse-gridding
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
 * produceReverseGriddingMatrix3D(vol(),Maux,kb);
 * interpolated_value = interpolatedElementReverseGridding(Maux,0.5,0.5,0.5,kb);
 * @endcode
 */
template <typename T>
T interpolatedElementReverseGridding(const Matrix3D<T> &in, double x, double y, double z, const KaiserBessel &kb)
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
    double pixel =0.;
    double w=0.;
    
    x = fmod(2*x, double(nx));
    y = fmod(2*y, double(ny));
    z = fmod(2*z, double(nz));
    int inxold = int(ROUND(x));
    int inyold = int(ROUND(y));
    int inzold = int(ROUND(z));
    
    double tablex1 = kb.i0win_tab(x-inxold+3);
    double tablex2 = kb.i0win_tab(x-inxold+2);
    double tablex3 = kb.i0win_tab(x-inxold+1);
    double tablex4 = kb.i0win_tab(x-inxold);
    double tablex5 = kb.i0win_tab(x-inxold-1);
    double tablex6 = kb.i0win_tab(x-inxold-2);
    double tablex7 = kb.i0win_tab(x-inxold-3);

    double tabley1 = kb.i0win_tab(y-inyold+3);
    double tabley2 = kb.i0win_tab(y-inyold+2);
    double tabley3 = kb.i0win_tab(y-inyold+1);
    double tabley4 = kb.i0win_tab(y-inyold);
    double tabley5 = kb.i0win_tab(y-inyold-1);
    double tabley6 = kb.i0win_tab(y-inyold-2);
    double tabley7 = kb.i0win_tab(y-inyold-3); 
	
    double tablez1 = kb.i0win_tab(z-inzold+3);
    double tablez2 = kb.i0win_tab(z-inzold+2);
    double tablez3 = kb.i0win_tab(z-inzold+1);
    double tablez4 = kb.i0win_tab(z-inzold);
    double tablez5 = kb.i0win_tab(z-inzold-1);
    double tablez6 = kb.i0win_tab(z-inzold-2);
    double tablez7 = kb.i0win_tab(z-inzold-3); 

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

/// @defgroup ReverseGriddingRelated Reverse Gridding Related Functions
/// @ingroup ReverseGridding

/** Reverse-gridding based 2D projection operation
 * @ingroup ReverseGriddingRelated
 *
 * Extracts a plane from a volume using reverse-gridding interpolation,
 * knowing that this volume has been processed for gridding.
 *
 * rot, tilt and psi and the respective Euler angles
 *
 * To interpolate using gridding you must prepare the volume first!
 * An example to extract a Fourier-plane, i.e. to calculate a
 * real-space projection, would be:
 *
 * @code
 * KaiserBessel kb;
 * matrix3D<std::complex<double> > Faux;
 * TODO!!!!
 * @endcode
 */


/** Reverse-gridding based geometric transformation on a 2D matrix
 * @ingroup ReverseGriddingRelated
 *
 * Applies a geometric transformation to a 2D matrix with
 * reverse-gridding-based interpolation, knowing that this image has been
 * processed for gridding.
 *
 * A is a 3x3 transformation matrix.
 *
 * Note that the output dimensions should be given if a different
 * scale is to be applied.
 *
 * To interpolate using gridding you must prepare the image first!
 * An example to apply a gridding-based transformation would be:
 *
 * @code
 * KaiserBessel kb;
 * matrix2D<double> Maux,out;
 * produceReverseGriddingMatrix2D(img(),Maux,kb);
 * Matrix2D<double> A = rotation2DMatrix(63.1);
 * applyGeometryReverseGridding(out,A,Maux,IS_NOT_INV,DONT_WRAP,kb);
 * @endcode
 */
template<typename T>
void applyGeometryReverseGridding(Matrix2D<T> &M2, Matrix2D< double > A, 
				  const Matrix2D<T> &M1, 
				  const KaiserBessel &kb, bool inv, bool wrap, 
				  int nx = 0, int ny = 0, T outside = (T) 0)
{
    int m1, n1, m2, n2;
    double x, y, xp, yp;
    double minxp, minyp, maxxp, maxyp;
    int cen_x, cen_y;

    if ((XSIZE(A) != 3) || (YSIZE(A) != 3))
        REPORT_ERROR(1102, "Apply_geom: geometrical transformation is not 3x3");

    if (XSIZE(M1) == 0)
    {
        M2.clear();
        return;
    }

    // For scalings the output matrix is explicitly resized to the final size 
    // Otherwise, just take half the size of the gridding image
    if (nx == 0) 
	nx = XSIZE(M1) / GRIDDING_NPAD;
    if (ny == 0) 
	ny = XSIZE(M1) / GRIDDING_NPAD;
    M2.resize(ny, nx);
    M2.setXmippOrigin();

    if (!inv)
        A = A.inv();

    // Find center and limits of image
    cen_y  = (int)(YSIZE(M2) / 2);
    cen_x  = (int)(XSIZE(M2) / 2);
    // Take 2x oversize M1 dims into account for calculating the limits 
    minxp  = FIRST_XMIPP_INDEX(XSIZE(M1)/2);
    minyp  = FIRST_XMIPP_INDEX(YSIZE(M1)/2);
    maxxp  = LAST_XMIPP_INDEX(XSIZE(M1)/2);
    maxyp  = LAST_XMIPP_INDEX(YSIZE(M1)/2);

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
                dMij(M2, i, j) = interpolatedElementReverseGridding(M1,xp,yp,kb);
	    else
		dMij(M2, i, j) = outside;

          // Compute new point inside input image
            xp += dMij(A, 0, 0);
            yp += dMij(A, 1, 0);
        }
    }
}

/** Reverse-gridding based geometric transformation on a 3D matrix
 * @ingroup ReverseGriddingRelated
 *
 * Applies a geometric transformation to a 3D matrix with
 * reverse-gridding-based interpolation, knowing that this image has been
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
 * produceReverseGriddingMatrix3D(vol(),Maux,kb);
 * Matrix2D<double> A = rotation2DMatrix(63.1);
 * applyGeometryReverseGridding(out,A,Maux,IS_NOT_INV,DONT_WRAP,kb);
 * @endcode
 */
template<typename T>
void applyGeometryReverseGridding(Matrix3D<T> &V2, Matrix2D< double > A, 
				  const Matrix3D<T> &V1, const KaiserBessel &kb, 
				  bool inv, bool wrap, 
				  int nx = 0, int ny = 0, int nz = 0, T outside = (T) 0)
{
    int m1, n1, o1, m2, n2, o2;
    double x, y, z, xp, yp, zp;
    double minxp, minyp, maxxp, maxyp, minzp, maxzp;
    int   cen_x, cen_y, cen_z;

    if ((XSIZE(A) != 4) || (YSIZE(A) != 4))
        REPORT_ERROR(1102, "Apply_geom3D: geometrical transformation is not 4x4");

    if (XSIZE(V1) == 0)
    {
        V2.clear();
        return;
    }

    // For scalings the output matrix is explicitly resized to the final size 
    // Otherwise, just take half the size of the gridding image
    if (nx == 0) 
	nx = XSIZE(V1) / GRIDDING_NPAD;
    if (ny == 0) 
	ny = XSIZE(V1) / GRIDDING_NPAD;
    if (nz == 0) 
	nz = ZSIZE(V1) / GRIDDING_NPAD;
    V2.resize(nz, ny, nx);
    V2.setXmippOrigin();

    if (!inv)
        A = A.inv();

    // Find center of Matrix3D
    cen_z = (int)(V2.zdim / 2);
    cen_y = (int)(V2.ydim / 2);
    cen_x = (int)(V2.xdim / 2);
    // Take 2x oversize M1 dims into account for calculating the limits 
    minxp  = FIRST_XMIPP_INDEX(XSIZE(V1)/2);
    minyp  = FIRST_XMIPP_INDEX(YSIZE(V1)/2);
    minzp  = FIRST_XMIPP_INDEX(ZSIZE(V1)/2);
    maxxp  = LAST_XMIPP_INDEX(XSIZE(V1)/2);
    maxyp  = LAST_XMIPP_INDEX(YSIZE(V1)/2);
    maxzp  = LAST_XMIPP_INDEX(ZSIZE(V1)/2);

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
                    dVkij(V2, k, i, j) = (T) interpolatedElementReverseGridding(V1,xp,yp,zp,kb);
                else
                    dVkij(V2, k, i, j) = outside;

                // Compute new point inside input image
                xp += dMij(A, 0, 0);
                yp += dMij(A, 1, 0);
                zp += dMij(A, 2, 0);
            }
        }
}

/// @defgroup ForwardGridding Forward Gridding
/// @ingroup Gridding

// ***************************************************************************
// ************************ Reverse Gridding *********************************
// ***************************************************************************

/** Numerically approximate the voronoi area for a set of 2D points
 * @ingroup ForwardGridding
 *
 *  Calculate the Voronoi area for a set of 2D points. The borders are
 *  set to a rectangle from min_x,min_y to max_x, max_y.
 *  Therefore, border effects should be taken care of in the
 *  generation of the x and y coordinates
 *
 */
void approximateVoronoiArea(std::vector<double> &voronoi_area,
			    const std::vector<double> &xin,
                            const std::vector<double> &yin, 
			    const double oversample = 10.);

/** Interpolate Cartesian coordinates from any irregularly sampled grid
 * @ingroup ForwardGridding
 *
 *  Interpolate Cartesian coordinates from an arbitrarily sampled grid
 *  using (forward) gridding. Note that the voronoi areas of the
 *  irregularly sampled coordinates have to be provided, and that the
 *  output of this routine needs to be passed through
 *  produceForwardGriddingMatrix2D or produceForwardGriddingFourierMatrix2D!
 *  The dimensions of the resulting matrix should be given.
 *
 * @code
 * KaiserBessel kb;
 * std::vector<double> x,y,data,voronoi_area;
 * matrix2D<double> Maux;
 *
 * P.getCartesianCoordinates(x,y,data); // (P is a Polar<double>)
 * approximateVoronoiArea(voronoi_area,x,y);
 * Maux = interpolateCartesianFromArbitrarySampling(64,64,x,y,data,voronoi_area,kb);
 * produceForwardGriddingMatrix2D(Maux,Maux2,kb);
 * @endcode
 *
 */
template <typename T>
Matrix2D<T> interpolateCartesianFromArbitrarySampling(const int xdim, const int ydim,
						      const std::vector<double> &xin, 
						      const std::vector<double> &yin,
						      const std::vector<T> &data, 
						      const std::vector<double> &voronoi_area,
						      const KaiserBessel &kb)
{

    double wx,wy, xx, yy, dx, dy, r, w, sumw;
    Matrix2D<T> result;

    // Oversample result GRIDDING_NPAD times
    result.initZeros(GRIDDING_NPAD*xdim, GRIDDING_NPAD*ydim);
    result.setXmippOrigin();
    
    // 1. Convolution interpolation
    // Loop over all cartesian coordinates: add all sampled points
    // within the window in x or y to the result, multiplied with the
    // corresponding gridding weights
    double window_size = (double)(kb.get_window_size()/2);
    FOR_ALL_ELEMENTS_IN_MATRIX2D(result)
    {
	xx = (double) j;
	yy = (double) i;
	sumw = 0.;
	for (int ii = 0; ii < xin.size(); ii++)
	{
	    dx = GRIDDING_NPAD*xin[ii] - xx;
	    if (ABS(dx) < window_size)
	    { 		
		dy = GRIDDING_NPAD*yin[ii] - yy;
		if (ABS(dy) < window_size)
		{
		    // Gridding weight
		    w = voronoi_area[ii] * kb.i0win_tab(dx) *  kb.i0win_tab(dy);
		    MAT_ELEM(result,i,j) += data[ii] * w;
		    sumw += w;
		}
	    }
	}
	if (sumw > 0.) 
	    MAT_ELEM(result,i,j) /= sumw;
    }

    return result;

}

/** Finish forward gridding after interpolateCartesianFromArbitrarySampling
 * @ingroup ForwardGridding
 *
 *  Produces a 2D real-space matrix2d after having performed 
 *  interpolateCartesianFromArbitrarySampling for real-space coordinates.
 */
void produceForwardGriddingMatrix2D(const Matrix2D< double > &in, 
				    Matrix2D< double > &out,
				    KaiserBessel &kb);

/** Finish forward gridding after interpolateCartesianFromArbitrarySampling
 * @ingroup ForwardGridding
 *
 *  Produces a 2D real-space matrix2d after having performed
 *  interpolateCartesianFromArbitrarySampling for fourier-space coordinates.
 */
void produceForwardGriddingMatrix2D(const Matrix2D< std::complex<double > > &in, 
				    Matrix2D< double > &out,
				    KaiserBessel &kb,
				    bool is_centered = true);

/** Finish forward gridding after interpolateCartesianFromArbitrarySampling
 * @ingroup ForwardGridding
 *
 *  Produces a 2D fourier-space matrix2d after having performed
 *  interpolateCartesianFromArbitrarySampling for fourier-space coordinates.
 */
void produceForwardGriddingFourierMatrix2D(const Matrix2D< std::complex<double > > &in, 
					   Matrix2D< std::complex<double > > &out,
					   KaiserBessel &kb,
					   bool is_centered = true);


#endif
