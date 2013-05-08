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
	euclideanNorm();						//d=(v_1-v_2)^2
	gramMatrix();						    //k(v_1,v_2)=exp(-1/2*(d/sigma)^2)
	normMatrix();							//Z=1/det(Z)*Z

	getGreatestEigen();			//compute Eigen and get greates Eigen

	Y.initZeros(MAT_YSIZE(U_dR), MAT_XSIZE(U_dR));
	calcMapping();						//compute Mapping
}

void KernelPCA::euclideanNorm() {
	//std::cout << "Entering euclideanNorm" << std::endl;
	Z.initZeros(MAT_YSIZE(*X),MAT_YSIZE(*X));
	int i, j, m;
	for (j = 0; j < MAT_YSIZE(*X); j++) {
		for (i = 0; i < MAT_YSIZE(*X); i++) {
			for (m = 0; m < MAT_XSIZE(*X); m++) {
				double diff = MAT_ELEM(*X,i,m)-MAT_ELEM(*X,j,m);
				MAT_ELEM(Z,i,j)+=diff*diff;
			}}}
}

void KernelPCA::gramMatrix() {
	//std::cout << "Entering gramMatrix" << std::endl;
	FOR_ALL_ELEMENTS_IN_MATRIX2D(Z)
		MAT_ELEM(Z,i,j)=exp(-MAT_ELEM(Z,i,j)/(2*sigma*sigma));

	}

void KernelPCA::normMatrix() {
	Matrix2D<double> A;	//A=sums of rows
	A.initZeros(1, MAT_YSIZE(*X));
	double a = 0;

	int i, j;
	for (i = 0; i < MAT_YSIZE(*X); i++) {
		for (j = 0; j < MAT_YSIZE(*X); j++) {
			MAT_ELEM(A,0,i)+=MAT_ELEM(Z,j,i);	//sum of vector z_i (row)
		}
		MAT_ELEM(A,0,i)/=MAT_YSIZE(*X); //sum_row/(dimension N)

		a+=MAT_ELEM(A,0,i);

	}
	a /= MAT_YSIZE(*X);

	for (i = 0; i < MAT_YSIZE(*X); i++) {
		for (j = 0; j < MAT_YSIZE(*X); j++) {
			MAT_ELEM(Z,j,i)=(MAT_ELEM(Z,j,i)-MAT_ELEM(A,0,j)-MAT_ELEM(A,0,i));
			MAT_ELEM(Z,j,i)=MAT_ELEM(Z,j,i)+a;
		}
	}
}

void KernelPCA::getGreatestEigen() {
	//std::cout << "Entering getGreatestEigen" << std::endl;
	//calc Eigen-Value/Vector
	Z.eigs(U, W, V, index);

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
