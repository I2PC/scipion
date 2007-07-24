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

void produceGriddingFourierMatrix2D(const Matrix2D< complex< double > > &in, 
				    Matrix2D< complex< double > > &out,
				    KaiserBessel &kb)
{
    Matrix2D<double> aux;
    InverseFourierTransform(in,aux);
    aux.setXmippOrigin();
    produceGriddingFourierMatrix2D(aux,out,kb);
}

void produceGriddingFourierMatrix2D(const Matrix2D< double > &in, 
				    Matrix2D< complex< double > > &out,
				    KaiserBessel &kb)
{
    // 1. Set up constants and Kaiser-Bessel object
    int N = XSIZE(in) * GRIDDING_NPAD;
    float r = XSIZE(in) / 2.;
    float v = GRIDDING_K/(2.*N);
    kb = KaiserBessel(GRIDDING_ALPHA, GRIDDING_K, r, v , N);

    // 2. iFFT, pad with zeros and divide in real space by a sinhwin
    Matrix2D<double> aux2;
    double wx, wy;
    aux2.resize(N,N);
    aux2.setXmippOrigin();
    FOR_ALL_ELEMENTS_IN_MATRIX2D(in) {
	wy=kb.sinhwin(i);
	wx=kb.sinhwin(j);
	MAT_ELEM(aux2,i,j) = MAT_ELEM(in,i,j) / (wx * wy);
    }
    
    // 3. FFT & (forward) centering
    FourierTransform(aux2,out);
    CenterOriginFFT(out,true);
    
}

// Prepare a 2D real-space image for gridding
void produceGriddingMatrix2D(const Matrix2D< double > &in, 
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
void produceGriddingFourierMatrix3D(const Matrix3D< complex< double > > &in, 
				    Matrix3D< complex< double > > &out,
				    KaiserBessel &kb)
{
    Matrix3D<double> aux;
    InverseFourierTransform(in,aux);
    aux.setXmippOrigin();
    produceGriddingFourierMatrix3D(aux,out,kb);
}

void produceGriddingFourierMatrix3D(const Matrix3D< double > &in, 
				    Matrix3D< complex< double > > &out,
				    KaiserBessel &kb)
{
    // 1. Set up constants and Kaiser-Bessel object
    int N = XSIZE(in) * GRIDDING_NPAD;
    float r = XSIZE(in) / 2.;
    float v = GRIDDING_K/(2.*N);
    kb = KaiserBessel(GRIDDING_ALPHA, GRIDDING_K, r, v , N);

    // 2. iFFT, pad with zeros and divide in real space by a sinhwin
    Matrix3D<double> aux2;
    double wx, wy, wz;
    aux2.resize(N,N,N);
    aux2.setXmippOrigin();
    FOR_ALL_ELEMENTS_IN_MATRIX3D(in) {
	wz=kb.sinhwin(k);
	wy=kb.sinhwin(i);
	wx=kb.sinhwin(j);
	VOL_ELEM(aux2,k,i,j) = VOL_ELEM(in,k,i,j) / (wx * wy * wz);
    }
    
    // 3. FFT & (forward) centering
    FourierTransform(aux2,out);
    CenterOriginFFT(out,true);
    
}

// Prepare a 3D real-space volume for gridding
void produceGriddingMatrix3D(const Matrix3D< double > &in, 
			     Matrix3D< double > &out,
			     KaiserBessel &kb)
{
    // 1. Set up constants and Kaiser-Bessel object
    int N = XSIZE(in) * GRIDDING_NPAD;
    float r = XSIZE(in) / 2.;
    float v = GRIDDING_K/(2.*N);
    kb = KaiserBessel(GRIDDING_ALPHA, GRIDDING_K, r, v , N);
    
    // 2. Center FFT and divide in Fourier space by a sinhwin
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

