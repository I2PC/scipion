/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.csic.es (2013)
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

#include "nca.h"
#include <data/numerical_tools.h>

void NeighbourhoodCA::setLabels(const Matrix1D<unsigned char> &labels)
{
	this->labels=labels;
}

void NeighbourhoodCA::setSpecificParameters(double lambda, int K)
{
	this->lambda=lambda;
	this->K=K;
}

double NeighbourhoodCA::objectiveFunction()
{
	// Project X onto the subspace given by A, Y=X*A
	matrixOperation_AB(*X,A,Y);

	// Compute the distance between the nearest neighbours of each observation
	// in this new subspace
	D2Y.resizeNoCopy(MAT_YSIZE(Y),K);
	FOR_ALL_ELEMENTS_IN_MATRIX2D(D2Y)
	{
		double aux=0.;
		int actualj=MAT_ELEM(idx,i,j);
		for (int l=0; l<MAT_XSIZE(Y); ++l)
		{
			double diff=MAT_ELEM(Y,i,l)-MAT_ELEM(Y,actualj,l);
			aux+=diff*diff;
		}
		MAT_ELEM(D2Y,i,j)=exp(-aux);
	}
	D2Y.rowSum(D2YRowSum);

	// Compute F
	double F=0.;
	for (int i=0; i<MAT_YSIZE(Y); ++i)
	{
		double sumClassi=0.;
		int labeli=VEC_ELEM(labels,i);
		for (int j=0; j<K; ++j)
		{
			int actualj=MAT_ELEM(idx,i,j);
			if (VEC_ELEM(labels,actualj)==labeli)
				sumClassi+=MAT_ELEM(D2Y,i,j);
		}
		if (fabs(VEC_ELEM(D2YRowSum,i))>1e-16)
			F+=sumClassi/VEC_ELEM(D2YRowSum,i);
	}

	// Compute regularization
	double sumA2=0.0;
	FOR_ALL_ELEMENTS_IN_MATRIX2D(A)
		sumA2+=MAT_ELEM(A,i,j)*MAT_ELEM(A,i,j);
	sumA2/=MAT_XSIZE(A)*MAT_YSIZE(A);

	return -(F-lambda*sumA2);
}

double ncaObjectiveFuntion(double *p, void *prm)
{
	NeighbourhoodCA *nca=(NeighbourhoodCA *)prm;
	Matrix2D<double> &A=nca->A;
	memcpy(&MAT_ELEM(A,0,0),&(p[1]),MAT_XSIZE(A)*MAT_YSIZE(A)*sizeof(double));
	double c=nca->objectiveFunction();
	return c;
}

void NeighbourhoodCA::reduceDimensionality()
{
	Matrix2D<double> D2;
	subtractColumnMeans(*X);
	kNearestNeighbours(*X, K, idx, D2, distance, false);

	size_t d=MAT_XSIZE(*X);
    A.initGaussian(d,outputDim,0,0.01);

    // Optimize A
	Matrix1D<double> pA(MAT_XSIZE(A)*MAT_YSIZE(A)), steps(MAT_XSIZE(A)*MAT_YSIZE(A));
	steps.initConstant(1);
	memcpy(&VEC_ELEM(pA,0),&MAT_ELEM(A,0,0),VEC_XSIZE(pA)*sizeof(double));
	double goal;
	int iter;
	powellOptimizer(pA,1,VEC_XSIZE(pA),&ncaObjectiveFuntion,this,0.001,goal,iter,steps,true);
	memcpy(&MAT_ELEM(A,0,0),&VEC_ELEM(pA,0),VEC_XSIZE(pA)*sizeof(double));

	// Reduction
	matrixOperation_AB(*X,A,Y);
}
