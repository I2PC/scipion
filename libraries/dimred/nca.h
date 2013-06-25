/*
 * NCA.h
 *
 *  Created on: Feb 25, 2013
 *      Author: xmipp
 */

#ifndef NCA_H_
#define NCA_H_

#include <data/matrix2d.h>
#include <data/matrix1d.h>
#include "dimred_tools.h"

/**@defgroup Neighbourhood Component Analysis
   @ingroup DimRedLibrary */
//@{
/** Class for making a NCA dimensionality reduction */

class NeighbourhoodCA: public DimRedAlgorithm
{
public:
	/// Weight factor for regularization
	double lambda;

	/// Labels
	Matrix1D<unsigned char> labels;

	/// Number of neighbours
	int K;
public:
	// Auxiliary variables for evaluating the objective function
	Matrix2D<double> D2Y, A;
	Matrix1D<double> D2YRowSum;

	// Index of the nearest neighbors to each one of the observations
	Matrix2D<int> idx;
public:
	/// Set labels
	void setLabels(const Matrix1D<unsigned char> &labels);

	/// Set specific parameters
	void setSpecificParameters(double lambda=0., int K=12);

	double objectiveFunction();

	/// Reduce dimensionality
	void reduceDimensionality();
};
//@}

#endif /* NCA_H_ */
