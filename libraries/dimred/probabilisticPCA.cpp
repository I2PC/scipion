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

double trace(Matrix2D<double> m)
{
	double res=0;
	FOR_ALL_ELEMENTS_IN_MATRIX2D(m)
		if(i==j) res+= MAT_ELEM(m,i,j);
	return res;
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
    std::cout << sigma2 << "\n";
    double Q=MAXDOUBLE, oldQ;


    // S=cov(X);
    Matrix2D<double> S, W, inW, invM, sigma2invM, Ez, WtX, Wp1, Wp2, invWp2, Xt,invWt;
    subtractColumnMeans(*X);
    matrixOperation_AtA(*X,S);
    S/=(double)N;
    S.write("dimred/S.txt");

    W.initRandom(D,outputDim,0,2,RND_UNIFORM);
    W.write("dimred/W.txt");

    Xt=X->transpose();
    // ----- Hasta aqui ------

    //Multiplica At por A -> Wt * W y lo guarda en inW
    matrixOperation_AtA(W,inW);
    inW.write("dimred/inW.txt");

    MultidimArray <double> Ezz(N,outputDim,outputDim);
    while (!converged && iter<=Niters)
    {
        ++iter;
        std::cout << iter << "\n";

        // Perform E-step
        // inW+=sigma2*I
        for (size_t i=0; i<outputDim; ++i)
        	MAT_ELEM(inW,i,i)+=sigma2;
        inW.inv(invM);

        //Ez=invM*W^t*X

        matrixOperation_AtB(W,Xt,WtX);
        std::cout << "WtX: " << MAT_YSIZE(WtX) << "x" << MAT_XSIZE(WtX) << "\n";

        matrixOperation_AB(invM,WtX,Ez);
        std::cout << "Ez: " << MAT_YSIZE(Ez) << "x" << MAT_XSIZE(Ez) << "\n";

        // Ezz(k,:,:)=sigma2 * invM + Ez(:,k) * Ez(:,k)'
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

		// W = Wp1 / Wp2;
		Wp2.inv(invWp2);
		matrixOperation_AB(Wp1,invWp2,W);
		W.write("dimred/W_nueva.txt");
		matrixOperation_AtA(W,inW);


		// normX = sum(X .^ 2, 1); sum(A,1) suma los elementos de cada columna de A
		Matrix1D<double> normX(N,false);
		normX.initZeros(N);
		FOR_ALL_ELEMENTS_IN_MATRIX2D(Xt)
			VEC_ELEM(normX,j)+=pow(MAT_ELEM(Xt,i,j),2);

		double sigma2_nueva=0;

		//sigma2 = sigma2 + (normX(i) - 2 * Ez(:,i)' * W' * X(:,i) + trace(Ezz(:,:,i) * inW));
		Matrix1D<double> colEz(MAT_YSIZE(Ez)),colX(N);
		for (size_t k=0; k<N; ++k){
			//Ez(:,i)
			Ez.getCol(k,colEz);
			//X(:,i)
			(Xt).getCol(k,colX);
			//Ez(:,i)*W'
			Matrix1D<double> EzWt(MAT_YSIZE(W));
			EzWt.initZeros(MAT_YSIZE(W));
			FOR_ALL_ELEMENTS_IN_MATRIX2D(W)
				VEC_ELEM(EzWt,i) += VEC_ELEM(colEz,j)*MAT_ELEM(W,i,j);
			//Ez(:,i)' * W' * X(:,i)
			double EzWtX=0;
			FOR_ALL_ELEMENTS_IN_MATRIX1D(EzWt)
				EzWtX += VEC_ELEM(EzWt,i)*VEC_ELEM(colX,i);
			//trace(Ezz(:,:,i) * inW));
			Matrix2D<double> Ezz_k;
			Ezz_k.initZeros(outputDim,outputDim);
			FOR_ALL_ELEMENTS_IN_MATRIX2D(Ezz_k)
				MAT_ELEM(Ezz_k,i,j) = DIRECT_A3D_ELEM(Ezz,k,i,j);
			Matrix2D<double> Ezz_kinW;
			matrixOperation_AB(Ezz_k,inW,Ezz_kinW);
			double t = trace(Ezz_kinW);
			//sigma2 = sigma2 + (normX(i) - 2 * Ez(:,i)' * W' * X(:,i) + trace(Ezz(:,:,i) * inW));
			sigma2_nueva += VEC_ELEM(normX,k) - 2 * EzWtX + t;
		}
		sigma2_nueva = (1 / (N * D)) * sigma2_nueva;

		//Compute likelihood of new model
		oldQ = Q;
		std::cout << "oldQ= " << oldQ << "\n";
		if (iter > 1){
			//invC = ((1 / sigma2) * eye(D)) - ((1 / sigma2) * W * invM * W');
			Matrix2D<double> invC;//, aux1, aux2;
			Matrix2D<double> WinvM,WinvMWt,sigmaI,dimI,sigma_1I,WtSDI,WtSDIW,dimI_WtSDIW;

			sigma_1I=eye(D,1/sigma2_nueva);
			matrixOperation_AB(W,invM,WinvM);
			matrixOperation_AB(WinvM,W.transpose(),WinvMWt);
			WinvMWt = (1/sigma2_nueva)*WinvMWt;
			invC = sigma_1I - WinvMWt;
//			matrixOperation_AB(W,invM,aux1);
//			matrixOperation_AB(aux1,W.transpose(),aux2);
//			aux2 = (1/sigma2_nueva)*aux2;
//			invC = eye(D,1/sigma2_nueva) - aux2;

			//detC = det(sigma2 * eye(D)) * det(eye(no_dims) + W' * ((sigma2 .^ -1) * eye(D)) * W);
			double detC;
			sigmaI=eye(D,sigma2_nueva);
			dimI=eye(outputDim,1);

			matrixOperation_AtB(W,sigma_1I,WtSDI);
			matrixOperation_AB(WtSDI,W,WtSDIW);
			dimI_WtSDIW = dimI + WtSDIW;
			std::cout << "sigmaI: " << MAT_YSIZE(sigmaI) << "x" << MAT_XSIZE(sigmaI) << "\n";
			std::cout << "dimI_WtSDIW: " << MAT_YSIZE(dimI_WtSDIW) << "x" << MAT_XSIZE(dimI_WtSDIW) << "\n";
			std::cout << "sigmaI:\n" << sigmaI << std::endl;
			std::cout << "dimI_WtSDIW:\n" << dimI_WtSDIW << std::endl;

			detC = sigmaI.det() * dimI_WtSDIW.det();
//			matrixOperation_AB(eye(D,1/sigma2_nueva),W,aux1);
//			matrixOperation_AB(W.transpose(),aux1,aux2);
//			aux1=(eye(outputDim,1)+aux2);
//			detC = (eye(D,sigma2_nueva).det()) * aux1.det();

			//Q = (-n / 2) * (D * log(2 * pi) + log(detC) + trace(invC * S));
			Matrix2D<double> invCS;
			matrixOperation_AB(invC,S,invCS);
			Q = (-N/2) * (D * log (2*PI) + log(detC) + trace(invCS));
		}

		// Stop condition to detect convergence
		// Must not apply to the first iteration, because then it will end inmediately
		if(iter!=1)
		if (abs(oldQ-Q) < 0.001){
			converged=true;
		}
		sigma2=sigma2_nueva;

    }
    // ------ Desde aqui ------
    //mapping.M = (inW \ W')';
    Matrix2D<double> mappingM,Wt2(MAT_YSIZE(W),MAT_XSIZE(W)),auxW,inversa_inW;
    Wt2=W.transpose();
    W.write("dimred/W_ultima.txt");
    inW.inv(inversa_inW);

    inversa_inW.write("dimred/inversa_inW.txt");
    //Wt2.inv(invWt);
    matrixOperation_AB(inversa_inW,Wt2,mappingM);
    mappingM=mappingM.transpose();

    //mappedX, lo guardamos en Y
    matrixOperation_AB(*X,mappingM,Y);
    Y.write("dimred/Y.txt");


}
