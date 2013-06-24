/***************************************************************************
 *
 * Authors:    Javier Gamas (2013)
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

#include "npe.h"

void NPE::setSpecificParameters(int k)
{
	this->k=k;
}

void NPE::reduceDimensionality()
{
	if (MAT_YSIZE(*X) > MAT_XSIZE(*X))
			REPORT_ERROR(ERR_MATRIX_DIM, "Number of samples should be higher than number of dimensions.");

	int n = MAT_XSIZE(*X);
	int d = MAT_YSIZE(*X);

	Matrix2D<double> D;
	Matrix2D<int> idx, Eye;
	kNearestNeighbours(*X,k,idx,D,distance,false);

	Matrix2D<double> W(n,k), Xi, C;
	for (int i=0; i<n; ++i)
	{
		extractNearestNeighbours(*X, idx, i, Xi);
		subtractColumnMeans(Xi); //Xi = z
		matrixOperation_AtA(Xi,C);
	}
	if(k>d)
		double tol= 1e-5;
	else
		double tol= 0;

	Eye.initIdentity(k);

	Matrix2D<double> DP, WP;
	matrixOperation_XtAX_symmetric(*X,D,WP);
	matrixOperation_AtA(*X, DP);

	//Solve eigenvector problem
	Matrix2D<double> Peigvec, eigvector;
	Matrix1D<double> Deigval;
	generalizedEigs(WP,DP,Deigval,Peigvec);

	//Sort eigenvalues in descending order
	Matrix1D<int> idx2;
	Deigval.indexSort(idx2);

	eigvector.resizeNoCopy(MAT_YSIZE(Peigvec),outputDim);
	for(int j =0;j<outputDim;++j){
		int idxj=VEC_ELEM(idx2,j)-1;
		for(int i=0;i<MAT_YSIZE(Peigvec);++i)
			MAT_ELEM(eigvector, i, j)=MAT_ELEM(Peigvec,i,idxj);
	}

	//Compute results
	Y=*X*eigvector;
}
