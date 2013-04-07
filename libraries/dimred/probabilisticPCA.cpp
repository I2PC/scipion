/***************************************************************************
 *
 * Author:    Itziar Benito Ortega     (2013)
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

#include "probabilisticPCA.h"
#include <data/multidim_array.h>

void ProbabilisticPCA::setSpecificParameters(size_t Niters)
{
    this->Niters=Niters;
}

void ProbabilisticPCA::reduceDimensionality()
{
    size_t N=MAT_YSIZE(*X);  // N= number of rows of X
    size_t D=MAT_XSIZE(*X);  // D= number of columns of X

    bool converged=false;
    size_t iter=0;
    double sigma2=rnd_unif()*2;

    // S=cov(X)
    Matrix2D<double> S, W, inW, invM, sigma2invM, Ez, WtX, Wp1, Wp2;
    subtractColumnMeans(*X);
    matrixOperation_AtA(*X,S);
    S/=(double)N;

    W.initRandom(D,outputDim,0,2,RND_UNIFORM);
    matrixOperation_AtA(W,inW);

    MultidimArray <double> Ezz(MAT_YSIZE(*X),outputDim,outputDim);
    while (!converged && iter<=Niters)
    {
        ++iter;

        // Perform E-step
        // inW+=sigma2*I
        for (size_t i=0; i<D; ++i)
        	MAT_ELEM(inW,i,i)+=sigma2;
        inW.inv(invM);

        //Ez=invM*W^t*X
        matrixOperation_AtB(W,*X,WtX);
        matrixOperation_AB(invM,WtX,Ez);

        // Ezz(k,:,:)=sigma2 * invM + Ez(:,i) * Ez(:,i)'
        sigma2invM=invM;
        sigma2invM*=sigma2;
        for (size_t k=0; k<N; ++k)
        {
        	FOR_ALL_ELEMENTS_IN_MATRIX2D(sigma2invM)
        		DIRECT_A3D_ELEM(Ezz,k,i,j)=MAT_ELEM(sigma2invM,i,j)+MAT_ELEM(Ez,i,k)*MAT_ELEM(Ez,j,k);
        }

        // Perform M-step (maximize mapping W)
		Wp1.initZeros(D,outputDim);
		Wp2.initZeros(outputDim,outputDim);
		for (size_t k=0; k<N; ++k)
		{
			FOR_ALL_ELEMENTS_IN_MATRIX2D(Wp1)
				MAT_ELEM(Wp1,i,j)+=MAT_ELEM(*X,k,i)*MAT_ELEM(Ez,k,j);
			FOR_ALL_ELEMENTS_IN_MATRIX2D(Wp2)
				MAT_ELEM(Wp2,i,j)+=DIRECT_A3D_ELEM(Ezz,k,i,j);
		}
    }
}
