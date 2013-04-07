/***************************************************************************
 *
 * Authors:    Francisco Sanz Encinas      franciscosanz89@gmail.com (2013)
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

#include "gplvm.h"

void GPLVM::setSpecificParameters(double sigma)
{
	this->sigma=sigma;
}

double GPLVM::objectiveFunction(const Matrix2D<double> &Y)
{
	computeDistance(*X,D2,distance);
	computeSimilarityMatrix(D2,sigma);
	return 0;
}

void GPLVM::reduceDimensionality()
{
	subtractColumnMeans(*X);

	// Pretreat the data with PCA
	Matrix2D<double> C, M;
	matrixOperation_AAt(*X,C);
	Matrix1D<double> lambda;
	firstEigs(C, outputDim, lambda, M);

	Matrix2D<double> Y=*X*M;

	double c = objectiveFunction(Y);
}
