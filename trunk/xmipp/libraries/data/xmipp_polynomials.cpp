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

    //First argument means number of images
    //Second argument means number of pixels
    Matrix2D<double> zerMat(xdim*ydim, numZer);
    double iMaxDim2 = 2./std::max(xdim,ydim);

    int index = 0;
    //Actually polOrder corresponds to the polynomial order +1
    int polOrder=ZERNIKE_ORDER(coef.size());
    im.setXmippOrigin();

    Matrix2D<double> polValue(polOrder,polOrder);
    FOR_ALL_ELEMENTS_IN_ARRAY2D(im)
    {
        //For one i we swap the different j
        double y=i*iMaxDim2;
        double x=j*iMaxDim2;

        int ip=i-FIRST_XMIPP_INDEX(YSIZE(im));
        int jp=j-FIRST_XMIPP_INDEX(XSIZE(im));

        size_t pixel_idx=ip*XSIZE(im)+jp;

        //polValue = [ 0   x   x2    x3  ....
        //          y   yx  yx2   yx3  ...
        //         y2  y2x y3x2  y3x3 ]
        //dMij(polValue,py,px) py es fila, px es columna
        for (int px = 0; px < polOrder; ++px)
        {
            double xpx=std::pow(x,px);
            for (int py = 0; py < polOrder; ++py)
                dMij(polValue,py,px) = xpx*std::pow(y,py);
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
                    temp += dMij(*fMat,px,py)*dMij(polValue,px,py);

            dMij(zerMat,pixel_idx,idxZerMat) = temp;
            ++idxZerMat;
        }
    }

    Image<double> save;
    save().resizeNoCopy(im);
    FileName name;
    for (int k=0; k < numZer; ++k)
    {
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(im)
        DIRECT_A2D_ELEM(save(),i,j)=dMij(zerMat,i*XSIZE(im)+j,k);
        name.compose("PPP",k+1,"xmp");
        save.write(name);
    }

    Matrix1D<double> zernikeCoefficients, image;
    typeCast(im,image);

    //solveNonNegative(zerMat, image.transpose(), zernikeCoefficients);
    PR(zerMat.Xdim());
    PR(image.transpose().size());
    solveBySVD(zerMat, image.transpose(), zernikeCoefficients,1e-6);
    std::cout << zernikeCoefficients << std::endl;

    std::cout << "Ha pasado por fit ";
}
