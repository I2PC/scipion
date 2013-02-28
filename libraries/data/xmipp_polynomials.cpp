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

    int nMax=(int)VEC_XSIZE(coef);
    for (int nZ = 0; nZ < nMax; ++nZ)
    {
        if (VEC_ELEM(coef,nZ) == 0)
        {
            fMatT = new Matrix2D<int>(1,1);
            dMij(*fMatT,0,0)=0;
            //*fMatT = 0;
        }
        else
        {
            // Note that the paper starts in n=1 and we start in n=0
            int n = (size_t)ZERNIKE_ORDER(nZ);
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
                int K1=binom(l,2*i+p);
                for (int j = 0; j <= m; ++j)
                {
                    int factor = ( (i+j)%2 ==0 ) ? 1 : -1 ;
                    int K2=factor * K1 * fact(n-j)/(fact(j)*fact(m-j)*fact(n-m-j));
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
}

void PolyZernikes::fit(const Matrix1D<int> & coef, MultidimArray<double> & im, MultidimArray<double> &weight,
                       MultidimArray<bool> & ROI, int verbose)
{
    this->create(coef);

    size_t xdim = XSIZE(im);
    size_t ydim = YSIZE(im);
    //int numZer = (size_t)coef.sum();
    int numZer = (size_t)coef.sum();

    //Actually polOrder corresponds to the polynomial order +1
    int polOrder=(int)ZERNIKE_ORDER(coef.size());

    im.setXmippOrigin();

    Matrix2D<double> polValue(polOrder,polOrder);

    //First argument means number of images
    //Second argument means number of pixels
    WeightedLeastSquaresHelper weightedLeastSquaresHelper;
    Matrix2D<double>& zerMat=weightedLeastSquaresHelper.A;

    zerMat.resizeNoCopy((size_t)ROI.sum(), numZer);
    double iMaxDim2 = 2./std::max(xdim,ydim);

    size_t pixel_idx=0;

    weightedLeastSquaresHelper.b.resizeNoCopy((size_t)ROI.sum());
    weightedLeastSquaresHelper.w.resizeNoCopy(weightedLeastSquaresHelper.b);

    FOR_ALL_ELEMENTS_IN_ARRAY2D(im)
    {
        if ( (A2D_ELEM(ROI,i,j)))
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

            Matrix2D<int> *fMat;

            //We generate the representation of the Zernike polynomials
            for (int k=0; k < numZer; ++k)
            {
                fMat = &fMatV[k];

                if (fMat == NULL)
                    continue;

                double temp = 0;
                for (size_t px = 0; px < (*fMat).Xdim(); ++px)
                    for (size_t py = 0; py < (*fMat).Ydim(); ++py)
                        temp += dMij(*fMat,py,px)*dMij(polValue,py,px);

                dMij(zerMat,pixel_idx,k) = temp;
            }

            VEC_ELEM(weightedLeastSquaresHelper.b,pixel_idx)=A2D_ELEM(im,i,j);
            VEC_ELEM(weightedLeastSquaresHelper.w,pixel_idx)=std::abs(A2D_ELEM(weight,i,j));
            ++pixel_idx;
        }
    }

    Matrix1D<double> zernikeCoefficients;
    weightedLeastSquares(weightedLeastSquaresHelper, zernikeCoefficients);
    fittedCoeffs = zernikeCoefficients;

    // Pointer to the image to be fitted
    MultidimArray<double> reconstructed;

    reconstructed.resizeNoCopy(im);
    pixel_idx=0;

    FOR_ALL_ELEMENTS_IN_ARRAY2D(im)
    if (A2D_ELEM(ROI,i,j))
    {
        double temp=0;
        for (int k=0; k < numZer; ++k)
            temp+=dMij(zerMat,pixel_idx,k)*VEC_ELEM(fittedCoeffs,k);

        A2D_ELEM(reconstructed,i,j)=temp;

        if ( fabs(A2D_ELEM(reconstructed,i,j)-A2D_ELEM(im,i,j)) > PI)
            A2D_ELEM(ROI,i,j) = false;

        ++pixel_idx;
    }

    pixel_idx=0;

    if (verbose > 0)
    {
        Image<double> save;
        save()=reconstructed;
        save.write("reconstructedZernikes.xmp");
        ROI.write("ROI.txt");
    }
}

void PolyZernikes::zernikePols(const Matrix1D<int> coef, MultidimArray<double> & im, MultidimArray<bool> & ROI, int verbose)
{

    this->create(coef);

    int polOrder=(int)ZERNIKE_ORDER(coef.size());
    int numZer = coef.size();

    int xdim = XSIZE(im);
    int ydim = YSIZE(im);

    im.setXmippOrigin();

    Matrix2D<double> polValue(polOrder,polOrder);
    double iMaxDim2 = 2./std::max(xdim,ydim);

    double temp = 0;
    FOR_ALL_ELEMENTS_IN_ARRAY2D(im)
    {
        if (A2D_ELEM(ROI,i,j))
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

            Matrix2D<int> *fMat;
            //We generate the representation of the Zernike polynomials

            for (int k=0; k < numZer; ++k)
            {
                fMat = &fMatV[k];

                if ( (dMij(*fMat,0,0) == 0) && MAT_SIZE(*fMat) == 1 )
                    continue;

                for (size_t px = 0; px < (*fMat).Xdim(); ++px)
                    for (size_t py = 0; py < (*fMat).Ydim(); ++py)
                        temp += dMij(*fMat,py,px)*dMij(polValue,py,px)*VEC_ELEM(coef,k);
            }

            A2D_ELEM(im,i,j) = temp;
            temp = 0;
        }
    }

    STARTINGX(im)=STARTINGY(im)=0;

    if (verbose == 1)
    {
        Image<double> save;
        save()=im;
        save.write("PPP1.xmp");
    }
}
