/***************************************************************************
 * Authors:     Javier Vargas (jvargas@cnb.csic.es)
 *
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

#include "fringe_processing.h"
#include "xmipp_polynomials.h"
#include "xmipp_image.h"
#include "multidim_array.h"
#include "xmipp_funcs.h"
#include "xmipp_fftw.h"

//This member function simulates different types of fringe patterns.
// im is the resultant multidimarray
// type is the type of fringe pattern that we want to generate. Possible values are:
// SIMPLY_OPEN_FRINGES, SIMPLY_CLOSED_FRINGES, COMPLEX_OPEN_FRINGES, COMPLEX_CLOSED_FRINGES
// xdim is the image dimension in x axis
// ydim is the image dimension in y axis
void FringeProcessing::simulPattern(MultidimArray<double> & im, enum FP_TYPE type, int xdim, int ydim, double noiseLevel, const double freq, Matrix1D<int> coefs)
{

    //FIRST WE TRAVEL ALL ELEMENTS fn im
    int ndim = 1;
    int zdim = 1;
    im.resize(ndim,zdim,ydim,xdim,false);
    im.setXmippOrigin();

    double iMaxDim2 = 2./std::max(xdim,ydim);
    PolyZernikes polynom;

    if ( (type == COMPLEX_OPEN_FRINGES) || (type == COMPLEX_CLOSED_FRINGES) )
    {
        polynom.zernikePols(coefs,im);
        im.setXmippOrigin();
    }

    FOR_ALL_ELEMENTS_IN_ARRAY2D(im)
    {
        switch (type)
        {

        case SIMPLY_OPEN_FRINGES :
            {
                A2D_ELEM(im,i,j) = std::cos(j*iMaxDim2*freq)+rnd_gaus(0,noiseLevel);
                break;
            }
        case SIMPLY_CLOSED_FRINGES :
            {
                A2D_ELEM(im,i,j) = std::cos(50*std::exp(-0.5*(std::pow(i*iMaxDim2,2)+std::pow(j*iMaxDim2,2))))+rnd_gaus(0,noiseLevel);
                break;
            }

        case COMPLEX_OPEN_FRINGES :
            {
                A2D_ELEM(im,i,j) = std::cos(j*iMaxDim2*freq+A2D_ELEM(im,i,j))+rnd_gaus(0,noiseLevel);
                break;
            }
        case COMPLEX_CLOSED_FRINGES :
            {
                A2D_ELEM(im,i,j) = std::cos(50*std::exp(-0.5*(std::pow(i*iMaxDim2,2)+std::pow(j*iMaxDim2,2)))+A2D_ELEM(im,i,j))+rnd_gaus(0,noiseLevel);
                break;
            }
        }
    }

    STARTINGX(im)=STARTINGY(im)=0;
}

//Function to simulate some test fringe patterns.
//sd=SPHT(c) computes the quadrture term of c still affected by the
//direction phase factor. Therefore for a real c=b*cos(phi)
//sd=SPHT(c)=i*exp(i*dir)*b*sin(phi)
//Ref: Kieran G. Larkin, Donald J. Bone, and Michael A. Oldfield, "Natural
//demodulation of two-dimensional fringe patterns. I. General background of the spiral phase quadrature transform," J. Opt. Soc. Am. A 18, 1862-1870 (2001)
void FringeProcessing::SPTH(MultidimArray<double> & im, MultidimArray< std::complex <double> > & imProc)
{
	im.setXmippOrigin();
	MultidimArray<std::complex<double> > H, fftIm, imComplex;
	typeCast(im, imComplex);

    // Fourier Transformer
    FourierTransformer ftrans(FFTW_BACKWARD);
    ftrans.FourierTransform(imComplex, fftIm, false);

	H.resizeNoCopy(fftIm);
	H.setXmippOrigin();

	//std::complex<double> compTemp(0, 0);
	//i -> for exterior filas o Y
	//j -> for interior columns o X
    FOR_ALL_ELEMENTS_IN_ARRAY2D(im)
    {
    	A2D_ELEM(H,i,j) = (std::complex<double>(j,i))/std::sqrt(std::pow(i,2)+(std::pow(j,2)));
    	//We subtract the singular value in i=0, j=0
    	if ( (i==0) && (j==0) )
    		A2D_ELEM(H,i,j) = 0;
    }

    CenterFFT(H,true);
    //fftIm *= H;
    ftrans.inverseFourierTransform();
    //Here in the Matlab code there is a complex conjugate s
    //imProc = imComplex;
    imProc = fftIm;
}

