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
public:
	// Auxiliary variables for evaluating the objective function
	Matrix2D<double> A, Y, P;
public:
	/// Set specific parameters
	void setSpecificParameters(double lambda=0.);

	double objectiveFunction(Matrix2D<double> &A, Matrix2D<double> &cur_X, Matrix1D<double> &cur_labels);

	/// Reduce dimensionality
	void reduceDimensionality();
};
//@}

#endif /* NCA_H_ */
