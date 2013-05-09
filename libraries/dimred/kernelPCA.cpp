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
	computeDistance(*X, D2, distance, false);
	computeSimilarityMatrix(D2,sigma);
	normMatrix();							//Z=1/det(Z)*Z

	getGreatestEigen();			//compute Eigen and get greates Eigen

	Y.initZeros(MAT_YSIZE(U_dR), MAT_XSIZE(U_dR));
	calcMapping();						//compute Mapping
}

void KernelPCA::normMatrix() {
	Matrix2D<double> A;	//A=sums of rows
	A.initZeros(1, MAT_YSIZE(*X));
	double a = 0;

	int i, j;
	for (i = 0; i < MAT_YSIZE(*X); i++) {
		for (j = 0; j < MAT_YSIZE(*X); j++) {
			MAT_ELEM(A,0,i)+=MAT_ELEM(D2,j,i);	//sum of vector z_i (row)
		}
		MAT_ELEM(A,0,i)/=MAT_YSIZE(*X); //sum_row/(dimension N)

		a+=MAT_ELEM(A,0,i);

	}
	a /= MAT_YSIZE(*X);

	for (i = 0; i < MAT_YSIZE(*X); i++) {
		for (j = 0; j < MAT_YSIZE(*X); j++) {
			MAT_ELEM(D2,j,i)=(MAT_ELEM(D2,j,i)-MAT_ELEM(A,0,j)-MAT_ELEM(A,0,i));
			MAT_ELEM(D2,j,i)=MAT_ELEM(D2,j,i)+a;
		}
	}
}

void KernelPCA::getGreatestEigen() {
	//std::cout << "Entering getGreatestEigen" << std::endl;
	//calc Eigen-Value/Vector
	D2.eigs(U, W, V, index);

	U_dR.initZeros(MAT_YSIZE(*X), outputDim);
	W_dR.initZeros(outputDim);

	//search for and get i-th greates eigenvalue and associated eigenvalue
	for (size_t i = 0; i < outputDim; i++) {
		for (int j = 0; j < MAT_YSIZE(*X); j++) {
			if (i == VEC_ELEM(index,j)){
			VEC_ELEM(W_dR,i)=VEC_ELEM(W,j);

			for(int n=0;n<MAT_YSIZE(*X);n++)
			{	MAT_ELEM(U_dR,n,i)=MAT_ELEM(U,n,j);}

			j=MAT_YSIZE(*X);
		}}}
}

void KernelPCA::calcMapping() {
	Matrix1D<double> sqrtL(outputDim);
	for (size_t i = 0; i < outputDim; i++)
		VEC_ELEM(sqrtL,i)=sqrt(VEC_ELEM(W_dR,i));

	for (size_t i = 0; i < MAT_YSIZE(U_dR); i++)
	{
		MAT_ELEM(Y,i,0)=(MAT_ELEM(U_dR,i,0)*VEC_ELEM(sqrtL,0));
		MAT_ELEM(Y,i,1)=(MAT_ELEM(U_dR,i,1)*VEC_ELEM(sqrtL,1));
	}
}
