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

#include "xmipp_polynomials.h"
#include "xmipp_image.h"
#include "numerical_tools.h"
#include "mask.h"
#include "xmipp_fft.h"

#define PR(x) std::cout << #x " = " << x << std::endl;
#define ZERNIKE_ORDER(n) ceil((-3+sqrt(9+8*(double)(n)))/2.);

//This function generates the Zernike Polynomials. Please see: Efficient Cartesian representation of Zernike polynomials in computer memory
//SPIE Vol. 3190 pp. 382 to get more details about the implementation
void PolyZernikes::create(const Matrix1D<int> & coef)
{

    Matrix2D<int> * fMatT;

    for (int nZ = 0; nZ < VEC_XSIZE(coef); ++nZ)
    {
        if (VEC_ELEM(coef,nZ) == 0)
            fMatT = NULL;

        else
        {
        	// Note that the paper starts in n=1 and we start in n=0
            int n = ZERNIKE_ORDER(nZ);
            int l = 2*nZ-n*(n+2);
            int m = (n-l)/2;

            Matrix2D<int> testM(n+1,n+1);
            fMatT = new Matrix2D<int>(n+1,n+1);

            int p = (l>0);
            int q = ((n % 2 != 0) ? (abs(l)-1)/2 : (( l > 0 ) ? abs(l)/2-1 : abs(l)/2 ) );
            l = abs(l); //We want the positive value of l
            m = (n-l)/2;

            for (int i = 0; i <= q; ++i)
            {
                double K1=binom(l,2*i+p);
                for (int j = 0; j <= m; ++j)
                {
                    double factor = ( (i+j)%2 ==0 ) ? 1 : -1 ;
                    double K2=factor * K1 * fact(n-j)/(fact(j)*fact(m-j)*fact(n-m-j));
                    for (int k = 0; k <= (m-j); ++k)
                    {
                        int ypow = 2 * (i+k) + p;
                        int xpow = n - 2 * (i+j+k) - p;
                        dMij(*fMatT,xpow,ypow) +=K2 * binom(m-j,k);
                    }
                }
            }
        }

        fMatV.push_back(*fMatT);
    }
};

void PolyZernikes::fit(const Matrix1D<int> & coef, MultidimArray<double> & im)
{
    this->create(coef);

    int xdim = XSIZE(im);
    int ydim = YSIZE(im);
    int numZer = coef.sum();

    int index = 0;
    //Actually polOrder corresponds to the polynomial order +1
    int polOrder=ZERNIKE_ORDER(coef.size());
    CenterFFT(im, true);
    im.setXmippOrigin();

    Matrix2D<double> polValue(polOrder,polOrder);

    // Get the central part of the image
    Mask mask;
    mask.type = BINARY_CIRCULAR_MASK;
    mask.mode = INNER_MASK;
    mask.R1 = XSIZE(im)/2;
    mask.resize(im);
    mask.generate_mask();
    const MultidimArray<int>& mMask=mask.get_binary_mask();

    //First argument means number of images
    //Second argument means number of pixels
    PseudoInverseHelper pseudoInverter;
    Matrix2D<double>& zerMat=pseudoInverter.A;
    zerMat.resizeNoCopy(mMask.sum(), numZer);
    double iMaxDim2 = 2./std::max(xdim,ydim);

    size_t pixel_idx=0;

    pseudoInverter.b.resizeNoCopy(mMask.sum());

    FOR_ALL_ELEMENTS_IN_ARRAY2D(im)
    {
    	if (A2D_ELEM( mMask,i,j)>0)
    	{
			//For one i we swap the different j
			double y=i*iMaxDim2;
			double x=j*iMaxDim2;

			//polValue = [ 0    y   y2    y3   ...
			//             x   xy  xy2    xy3  ...
			//             x2  x2y x2y2   x2y3 ]

			//dMij(polValue,py,px) py es fila, px es columna
			for (int py = 0; py < polOrder; ++py)
			{
				double ypy=std::pow(y,py);
				for (int px = 0; px < polOrder; ++px)
					dMij(polValue,px,py) = ypy*std::pow(x,px);
			}

			int idxZerMat=0;
			Matrix2D<int> *fMat;

			//We generate the representation of the Zernike polynomials
			for (int k=0; k < numZer; ++k)
			{
				fMat = &fMatV[k];

				if (fMat == NULL)
					continue;

				double temp = 0;
				for (int px = 0; px < (*fMat).Xdim(); ++px)
					for (int py = 0; py < (*fMat).Ydim(); ++py)
						temp += dMij(*fMat,py,px)*dMij(polValue,py,px);

				dMij(zerMat,pixel_idx,idxZerMat) = temp;
				++idxZerMat;
			}

			//std::cout << A2D_ELEM(im,i,j) << std::endl;
	    	VEC_ELEM(pseudoInverter.b,pixel_idx)=A2D_ELEM(im,j,i);
	    	++pixel_idx;
		}
    }

    Matrix1D<double> zernikeCoefficients;
    solveLinearSystem(pseudoInverter, zernikeCoefficients);
    fittedCoeffs = zernikeCoefficients;

    for (int i=0; i<zernikeCoefficients.size();i++)
        std::cout << "Zernike Coeff " << i << " : " << VEC_ELEM(fittedCoeffs,i) << std::endl;

    //TEST:
    //VEC_ELEM(fittedCoeffs,0) =     1.2319;
    //VEC_ELEM(fittedCoeffs,1) =    -0.0075;
    //VEC_ELEM(fittedCoeffs,2) =    -0.0081;
    //VEC_ELEM(fittedCoeffs,3) =     0.1824;
    //VEC_ELEM(fittedCoeffs,4) =    -2.0592;
    //VEC_ELEM(fittedCoeffs,5) =    -0.0399;
    //VEC_ELEM(fittedCoeffs,6) =    -0.0005;
    //VEC_ELEM(fittedCoeffs,7) =    -0.0148;
    //VEC_ELEM(fittedCoeffs,8) =    -0.0160;
    //VEC_ELEM(fittedCoeffs,9) =     0.0017;

    reconstructed.resizeNoCopy(im);

    pixel_idx=0;
    double Temp=0;

    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(im)
    if (DIRECT_A2D_ELEM(mMask,i,j)>0)
    {
        for (int k=0; k < numZer; ++k)
        	Temp+=dMij(zerMat,pixel_idx,k)*VEC_ELEM(fittedCoeffs,k);

        DIRECT_A2D_ELEM(reconstructed,i,j)=Temp;
        Temp=0;
    	++pixel_idx;
    }

    pixel_idx=0;


    Image<double> save;
    save().resizeNoCopy(im);
    FileName name;
    save() = reconstructed;
    name.compose("Reconstructed",1,"xmp");
    save.write(name);

    Image<int> save2(mMask);
    FileName name2;
    name2.compose("Mask",1,"xmp");
    save2.write(name2);

//    Image<double> save3;
//    save3() = zerMat;
//    FileName name3;
//    name3.compose("zerMat",1,"xmp");
//    save3.write(name3);
//    //b.mapToFile(name,512*512,1);
    zerMat.write("zerMat.txt");
    pseudoInverter.b.write("b.txt");


}
