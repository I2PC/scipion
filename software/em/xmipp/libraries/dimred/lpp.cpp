/***************************************************************************
 *
 * Authors:    Oscar Hueros Fernandez      o.hueros@gmail.com (2013)
 *
 * Universidad San Pablo CEU
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

#include "lpp.h"

void LPP::setSpecificParameters(int k, double sigma)
{
	this->k=k;
	this->sigma=sigma;
}
/** Reduce dimensionality method based on the Locality Preserving Projections (LPP) algorithm.
 *  These are linear projective maps that arise by solving a variational problem
 *  that optimally preserves the neighborhood structure of the data set.
 */
void LPP::reduceDimensionality()
{
	// Compute the distance to the k nearest neighbors
	Matrix2D<double> D2;
	computeDistanceToNeighbours(*X, k, D2, distance, false);

	// Compute similarity matrix
	computeSimilarityMatrix(D2,sigma,true,true);

	// Compute graph laplacian
	Matrix2D<double> L;
	computeGraphLaplacian(D2,L);

	Matrix2D<double> DP, LP;
	matrixOperation_XtAX_symmetric(*X,D2,DP);
	matrixOperation_XtAX_symmetric(*X,L,LP);

	// Compute eigenvalues and eigenvectors resolving the generalized eigenvector problem
	Matrix2D<double> Peigvec, eigvector;
	Matrix1D<double> Deigval;
	generalizedEigs(LP, DP, Deigval, Peigvec);

	// Sort the eigenvalues in ascending order
	Matrix1D<int> idx;
	Deigval.indexSort(idx);

	// Sort the eigenvalues in descending order and get the smallest eigenvectors
	eigvector.resizeNoCopy(MAT_YSIZE(Peigvec),outputDim);
	for (size_t j=0; j<outputDim; ++j)
	{
		int idxj=VEC_ELEM(idx,j)-1;
		for (int i=0; i<MAT_YSIZE(Peigvec); ++i)
			MAT_ELEM(eigvector,i,j)=MAT_ELEM(Peigvec,i,idxj);
	}

	// Compute the result of the reduce dimensionality method
	Y=*X * eigvector;
	if (fnMapping!="")
		eigvector.write(fnMapping);
}
