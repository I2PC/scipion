/***************************************************************************
 *
 * Authors:    Albrecht Wolf
 *
 ***************************************************************************/

#include "kernelPCA.h"

void KernelPCA::setSpecificParameters(double sigma)
{
    this->sigma=sigma;
}

void KernelPCA::reduceDimensionality()
{
    Matrix2D<double> Z;
    computeDistance(*X,Z,distance,false);
    computeSimilarityMatrix(Z,sigma);

    // Normalize the graph
    Matrix1D<double> zi;
    Z.computeRowMeans(zi);
    double z=zi.computeMean();
    FOR_ALL_ELEMENTS_IN_MATRIX2D(Z)
    	MAT_ELEM(Z,i,j)+=z-VEC_ELEM(zi,i)-VEC_ELEM(zi,j);

    // Get the largest eigenvalues
    Matrix2D<double> U;
	Matrix1D<double> W;
    firstEigs(Z, outputDim, W, U);

    // Compute mapping
    Y.initZeros(MAT_YSIZE(*X),outputDim);
    Matrix1D<double> sqrtL(outputDim);
    for(size_t i=0;i<outputDim; i++)
        VEC_ELEM(sqrtL,i)=sqrt(VEC_ELEM(W,i));

    FOR_ALL_ELEMENTS_IN_MATRIX2D(Y)
        MAT_ELEM(Y,i,j)=(MAT_ELEM(U,i,j)*VEC_ELEM(sqrtL,j));
}
