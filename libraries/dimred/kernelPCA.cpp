/***************************************************************************
 *
 * Authors:    Albrecht Wolf
 *
 ***************************************************************************/

#include "kernelPCA.h"

void KernelPCA::setSpecificParameters(double sigma) {
	this->sigma = sigma;
}

void KernelPCA::reduceDimensionality() {
	// Compute Gram matrix
	Matrix2D<double> D2;
	computeDistance(*X, D2, distance, false);
	computeSimilarityMatrix(D2,sigma);

	// Normalize it
	Matrix1D<double> mean_i;
	D2.rowSum(mean_i);
	mean_i/=MAT_XSIZE(D2);
	double mean=mean_i.computeMean();
	FOR_ALL_ELEMENTS_IN_MATRIX2D(D2)
		MAT_ELEM(D2,i,j)+=mean-VEC_ELEM(mean_i,i)-VEC_ELEM(mean_i,j);

	// Compute the largest eigenvalues
	Matrix1D<double> lambda;
	firstEigs(D2,outputDim,lambda,Y);

	// Readjust variances
	FOR_ALL_ELEMENTS_IN_MATRIX1D(lambda)
		VEC_ELEM(lambda,i)=sqrt(VEC_ELEM(lambda,i));
	FOR_ALL_ELEMENTS_IN_MATRIX2D(Y)
		MAT_ELEM(Y,i,j)*=VEC_ELEM(lambda,j);
}
