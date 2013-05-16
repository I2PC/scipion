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

Matrix2D<double> eye(size_t dim, double factor){
	Matrix2D<double> res(dim,dim);
	res.initZeros(dim,dim);
	FOR_ALL_ELEMENTS_IN_MATRIX2D(res)
		if(i==j) MAT_ELEM(res,i,j)=factor;
	return res;
}

void ProbabilisticPCA::reduceDimensionality()
{
    size_t N=MAT_YSIZE(*X);  // N= number of rows of X
    size_t D=MAT_XSIZE(*X);  // D= number of columns of X

    bool converged=false;
    size_t iter=0;
    double sigma2=rnd_unif()*2;
    double Q=MAXDOUBLE, oldQ;



    Matrix2D<double> S, W, inW, invM, sigma2invM, Ez, WtX, Wp1, Wp2, invWp2, Xt,invWt;
    subtractColumnMeans(*X);
    matrixOperation_AtA(*X,S);
    S/=(double)N;


    W.initRandom(D,outputDim,0,2,RND_UNIFORM);


    Xt=X->transpose();


    matrixOperation_AtA(W,inW);


    MultidimArray <double> Ezz(N,outputDim,outputDim);
    while (!converged && iter<=Niters)
    {
        ++iter;


        // Perform E-step

        for (size_t i=0; i<outputDim; ++i)
        	MAT_ELEM(inW,i,i)+=sigma2;
        inW.inv(invM);

        matrixOperation_AtB(W,Xt,WtX);

        matrixOperation_AB(invM,WtX,Ez);


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
				MAT_ELEM(Wp1,i,j)+=MAT_ELEM(Xt,i,k)*MAT_ELEM(Ez,j,k);
			FOR_ALL_ELEMENTS_IN_MATRIX2D(Wp2)
				MAT_ELEM(Wp2,i,j)+=DIRECT_A3D_ELEM(Ezz,k,i,j);
		}


		Wp2.inv(invWp2);
		matrixOperation_AB(Wp1,invWp2,W);
		matrixOperation_AtA(W,inW);


		Matrix1D<double> normX(N,false);
		normX.initZeros(N);
		FOR_ALL_ELEMENTS_IN_MATRIX2D(Xt)
			VEC_ELEM(normX,j)+=pow(MAT_ELEM(Xt,i,j),2);

		double sigma2_nueva=0;


		Matrix1D<double> colEz(MAT_YSIZE(Ez)),colX(N);
		for (size_t k=0; k<N; ++k){

			Ez.getCol(k,colEz);

			(Xt).getCol(k,colX);

			Matrix1D<double> EzWt(MAT_YSIZE(W));
			EzWt.initZeros(MAT_YSIZE(W));
			FOR_ALL_ELEMENTS_IN_MATRIX2D(W)
				VEC_ELEM(EzWt,i) += VEC_ELEM(colEz,j)*MAT_ELEM(W,i,j);

			double EzWtX=0;
			FOR_ALL_ELEMENTS_IN_MATRIX1D(EzWt)
				EzWtX += VEC_ELEM(EzWt,i)*VEC_ELEM(colX,i);

			Matrix2D<double> Ezz_k;
			Ezz_k.initZeros(outputDim,outputDim);
			FOR_ALL_ELEMENTS_IN_MATRIX2D(Ezz_k)
				MAT_ELEM(Ezz_k,i,j) = DIRECT_A3D_ELEM(Ezz,k,i,j);
			Matrix2D<double> Ezz_kinW;
			matrixOperation_AB(Ezz_k,inW,Ezz_kinW);
			double t = Ezz_kinW.trace();

			sigma2_nueva += VEC_ELEM(normX,k) - 2 * EzWtX + t;
		}
		double ND;
		ND = (double) N * (double) D;
		sigma2_nueva = (1 / (ND)) * sigma2_nueva;


		//Compute likelihood of new model
		oldQ = Q;


		if (iter > 1){

			Matrix2D<double> invC;
			Matrix2D<double> WinvM,WinvMWt,sigmaI,dimI,sigma_1I,WtSDI,WtSDIW,dimI_WtSDIW;

			sigma_1I=eye(D,1/sigma2_nueva);
			matrixOperation_AB(W,invM,WinvM);
			matrixOperation_AB(WinvM,W.transpose(),WinvMWt);
			WinvMWt = (1/sigma2_nueva)*WinvMWt;
			invC = sigma_1I - WinvMWt;

			double detC;
			sigmaI=eye(D,sigma2_nueva);
			dimI=eye(outputDim,1);

			matrixOperation_AtB(W,sigma_1I,WtSDI);
			matrixOperation_AB(WtSDI,W,WtSDIW);
			dimI_WtSDIW = dimI + WtSDIW;

			detC = sigmaI.det() * dimI_WtSDIW.det();

			Matrix2D<double> invCS;
			matrixOperation_AB(invC,S,invCS);
			if(iter==2) {
				std::cout << "N= " << N << "\n";
				std::cout << "-N*0.5= " << N*(-0.5) << "\n";
			}
			Q = (N*(-0.5)) * (D * log (2*PI) + log(detC) + invCS.trace());
		}


		// Stop condition to detect convergence
		// Must not apply to the first iteration, because then it will end inmediately
		if(iter>2)
		if (abs(oldQ-Q) < 0.001){
			converged=true;
		}

		sigma2=sigma2_nueva;

    }

    //mapping.M = (inW \ W')';
    Matrix2D<double> mappingM,Wt2(MAT_YSIZE(W),MAT_XSIZE(W)),auxW,inversa_inW;
    Wt2=W.transpose();
    inW.inv(inversa_inW);


    matrixOperation_AB(inversa_inW,Wt2,mappingM);
    mappingM=mappingM.transpose();

    //mappedX, lo guardamos en Y
    matrixOperation_AB(*X,mappingM,Y);
    Y.write("dimred/Y.txt");


}
