/***************************************************************************
 *
 * Authors:    Albrecht Wolf
 *
 ***************************************************************************/
#ifndef _KERNEL_PCA
#define _KERNEL_PCA

#include <data/matrix2d.h>
#include <data/matrix1d.h>
#include "dimred_tools.h"

/**@defgroup KernelPCA Kernel PCA
   @ingroup DimRedLibrary */
//@{
/** Class for making a Kernel PCA dimensionality reduction */
class KernelPCA: public DimRedAlgorithm
{
public:
	double sigma;
	int N,M;

	Matrix2D<double> Z, U, V, U_dR;
	Matrix1D<double> W, W_dR, e_W_dR, invsqrtL;
	Matrix1D<int> index;
public:
	/// Set specific parameters
	void setSpecificParameters(double sigma=1.0);

	/// Reduce dimensionality
	void reduceDimensionality();

	/// Euclidean Norm
	void euclideanNorm();

	/// compute gram Matrix
	void gramMatrix();

	/// Normalize Matrix
	void normMatrix();

	/// Eigendecomposition and largest Eigen-Value/-Vector
	void getGreatestEigen();

	// calc Mapping Values
	void calcMapping();
};
//@}
#endif
