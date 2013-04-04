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


void NeighbourhoodCA::setSpecificParameters(double lambda)
{
	this->lambda=lambda;
}

//double NeighbourhoodCA::objectiveFunction(Matrix1D<double> x){
double NeighbourhoodCA::objectiveFunction(Matrix2D<double> &A, Matrix2D<double> &X, Matrix1D<double> &labels)
{
	matrixOperation_AB(X,A,Y); // Y= X*A
	Matrix2D<double> D2;
	computeDistance(Y,D2,distance,false);
	computeSimilarityMatrix(D2,1);

	// Compute the sum of rows, normalize and compute F
	for (size_t i=0; i<MAT_XSIZE(D2); ++i)
		MAT_ELEM(D2,i,i)=0.;
	Matrix1D<double> rowSum;
	D2.rowSum(rowSum);
	FOR_ALL_ELEMENTS_IN_MATRIX1D(rowSum)
		VEC_ELEM(rowSum,i)=1.0/VEC_ELEM(rowSum,i);
	const double eps=1e-10;
	double F=0.;
	FOR_ALL_ELEMENTS_IN_MATRIX2D(D2)
	{
		MAT_ELEM(D2,i,j)*=VEC_ELEM(rowSum,i);
		if (MAT_ELEM(D2,i,j)<eps)
			MAT_ELEM(D2,i,j)=eps;
		F+=MAT_ELEM(D2,i,j);
	}

	double sumA2=0;
	FOR_ALL_ELEMENTS_IN_MATRIX2D(A)
		sumA2+=MAT_ELEM(A,i,j)*MAT_ELEM(A,i,j);
	double avgA2=sumA2/(MAT_XSIZE(A)*MAT_YSIZE(A));
	return -(F-lambda*avgA2);
}

void NeighbourhoodCA::reduceDimensionality()
{
	// See drtoolbox/techniques/nca.m
	size_t d=MAT_XSIZE(*X)-1;

	Matrix2D<double> cur_X;
	Matrix1D<double> cur_labels;
	X->getRow(0,cur_labels);
	cur_X=*X;
	eraseFirstColumn(cur_X);
	subtractColumnMeans(cur_X);
	//cur_X.write("dimred/MatrixcurX.txt");

    Matrix2D<double> A;
    A.initGaussian(d,outputDim,0,0.01);

    //Debug
	A(0,0)=0.001;
	A(1,0)=-0.0098;
	A(0,1)=-0.0014;
	A(1,1)=-0.0026;
	double F=objectiveFunction(A,cur_X,cur_labels);
	std::cout<<F<<std::endl;

	size_t max_iter = 200;
}
