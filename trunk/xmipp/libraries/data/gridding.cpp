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
#include "gridding.h"

// Prepare a 2D Fourier-space image for gridding
void produceGriddingImage(const Matrix2D< complex< double> > &in, 
			  Matrix2D< complex< double > > &out,
			  KaiserBessel &kb)
{
    // 1. Set up constants and Kaiser-Bessel object
    int N = XSIZE(in) * GRIDDING_NPAD;
    float r = XSIZE(in) / 2.;
    float v = GRIDDING_K/(2.*N);
    kb = KaiserBessel(GRIDDING_ALPHA, GRIDDING_K, r, v , N);

    // 2. iFFT, pad with zeros and divide in real space by a sinhwin
    Matrix2D<double> aux,aux2;
    double wx, wy;
    InverseFourierTransform(in,aux);
    aux.setXmippOrigin();
    aux2.resize(N,N);
    aux2.setXmippOrigin();
    FOR_ALL_ELEMENTS_IN_MATRIX2D(aux) {
	wy=kb.sinhwin(i);
	wx=kb.sinhwin(j);
	MAT_ELEM(aux2,i,j) = MAT_ELEM(aux,i,j) / (wx * wy);
    }
    
    // 3. FFT & (forward) centering
    FourierTransform(aux2,out);
    CenterOriginFFT(out,true);
    
}

// Prepare a 2D real-space image for gridding
void produceGriddingImage(const Matrix2D< double > &in, 
			  Matrix2D< double > &out,
			  KaiserBessel &kb)
{
    // 1. Set up constants and Kaiser-Bessel object
    int N = XSIZE(in) * GRIDDING_NPAD;
    float r = XSIZE(in) / 2.;
    float v = GRIDDING_K/(2.*N);
    kb = KaiserBessel(GRIDDING_ALPHA, GRIDDING_K, r, v , N);
    
    // 2. FFT, pad with zeros and divide in Fourier space by a sinhwin
    Matrix2D<complex<double> > aux,aux2;
    double wx, wy;
    FourierTransform(in,aux);
    CenterFFT(aux,true);
    aux.setXmippOrigin();
    aux2.resize(N,N);
    aux2.setXmippOrigin();
    FOR_ALL_ELEMENTS_IN_MATRIX2D(aux) {
	wy=kb.sinhwin(i);
	wx=kb.sinhwin(j);
	MAT_ELEM(aux2,i,j) = MAT_ELEM(aux,i,j) / (wx * wy);
    }
 
    // 3. (backward) centering & IFFT
    CenterFFT(aux2,false);
    InverseFourierTransform(aux2,out);
    out.setXmippOrigin();

}

// Prepare a 3D Fourier-space volume for gridding
void produceGriddingVolume(const Matrix3D<complex< double> > &in, 
			   Matrix3D< complex< double > > &out,
			   KaiserBessel &kb)
{
    // 1. Set up constants and Kaiser-Bessel object
    int N = XSIZE(in) * GRIDDING_NPAD;
    float r = XSIZE(in) / 2.;
    float v = GRIDDING_K/(2.*N);
    kb = KaiserBessel(GRIDDING_ALPHA, GRIDDING_K, r, v , N);

    // 2. iFFT, pad with zeros and divide in real space by a sinhwin
    Matrix3D<double> aux,aux2;
    double wx, wy, wz;
    InverseFourierTransform(in,aux);
    aux.setXmippOrigin();
    aux2.resize(N,N,N);
    aux2.setXmippOrigin();
    FOR_ALL_ELEMENTS_IN_MATRIX3D(aux) {
	wz=kb.sinhwin(k);
	wy=kb.sinhwin(i);
	wx=kb.sinhwin(j);
	VOL_ELEM(aux2,k,i,j) = VOL_ELEM(aux,k,i,j) / (wx * wy * wz);
    }
    
    // 3. FFT & (forward) centering
    FourierTransform(aux2,out);
    CenterOriginFFT(out,true);
    
}

// Prepare a 3D real-space volume for gridding
void produceGriddingVolume(const Matrix3D< double > &in, 
			   Matrix3D< double > &out,
			   KaiserBessel &kb)
{
    // 1. Set up constants and Kaiser-Bessel object
    int N = XSIZE(in) * GRIDDING_NPAD;
    float r = XSIZE(in) / 2.;
    float v = GRIDDING_K/(2.*N);
    kb = KaiserBessel(GRIDDING_ALPHA, GRIDDING_K, r, v , N);
    
    // 2. Center FFT and divide in Fourier space by a sinhwin
    // I could actually pre-calculate the sinhwin for all images!!
    Matrix3D<complex<double> > aux,aux2;
    double wx, wy, wz;
    FourierTransform(in,aux);
    CenterFFT(aux,true);
    aux.setXmippOrigin();
    aux2.resize(N,N,N);
    aux2.setXmippOrigin();
    FOR_ALL_ELEMENTS_IN_MATRIX3D(aux) {
	wz=kb.sinhwin(k);
	wy=kb.sinhwin(i);
	wx=kb.sinhwin(j);
	VOL_ELEM(aux2,k,i,j) = VOL_ELEM(aux,k,i,j) / (wx * wy * wz);
    }

    // 3. (backward) centering & iFFT
    CenterFFT(aux2,false);
    InverseFourierTransform(aux2,out);
    out.setXmippOrigin();

}

/*
void griddingRotateShiftScale(const Matrix2D<double> &in,
			      Matrix2D<double> &out,
			      float ang, float delx, float dely, 
			      KaiserBessel &kb, float scale_input)
{
	int nx, ny, nxn, nyn;

	if (scale_input == 0.0f) scale_input = 1.0f;
	float  scale = 0.5*scale_input;

	nx = XSIZE(in);
	ny = YSIZE(in);
	nxn = nx/2; nyn = ny/2;
	out.resize(nxn,nyn);

	if(delx >= 0.0f) { delx = fmod(delx, float(nx));} else {delx = -fmod(-delx, float(nx));}
	if(dely >= 0.0f) { dely = fmod(dely, float(ny));} else {dely = -fmod(-dely, float(ny));}
	// center of big image,
	int xc = nxn;
	int ixs = nxn%2;  // extra shift on account of odd-sized images
	int yc = nyn;
	int iys = nyn%2;
	// center of small image
	int xcn = nxn/2;
	int ycn = nyn/2;
	// shifted center for rotation
	float shiftxc = xcn + delx;
	float shiftyc = ycn + dely;
	// bounds if origin at center
	float ymin = -ny/2.0f;
	float xmin = -nx/2.0f;
	float ymax = -ymin;
	float xmax = -xmin;
	if (0 == nx%2) xmax--;
	if (0 == ny%2) ymax--;
	
	// trig
	float cang = cos(ang);
	float sang = sin(ang);
	for (int iy = 0; iy < nyn; iy++) {
		float y = float(iy) - shiftyc;
		float ycang = y*cang/scale + yc;
		float ysang = -y*sang/scale + xc;
		for (int ix = 0; ix < nxn; ix++) {
			float x = float(ix) - shiftxc;
			float xold = x*cang/scale + ysang-ixs;// have to add the fraction on account on odd-sized images for which Fourier zero-padding changes the center location 
			float yold = x*sang/scale + ycang-iys;
			
			xold = xold/2.0;
			yold = yold/2.0;
			cerr<<"xp,yp= "<<xold<<" "<<yold<<" x,y= "<<ix<<" "<<iy<<endl;
			dMij(out,iy,ix) = get_new_pixel_conv2D(nx,ny,xold,yold,in.data,kb);
			
		}
	}

}
*/

