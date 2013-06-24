/***************************************************************************
 *
 * Authors:    Carlos Oscar Sanchez Sorzano      coss@cnb.csic.es (2013)
 *             Alejandro Gomez Rodriguez (2013)
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

#include "spe.h"
#include <data/numerical_tools.h>

void SPE::setSpecificParameters(bool global, int k)
{
	this->global = true;
	this->k = 1;
}

void SPE::reduceDimensionality()
{
	double lambda = 1;  // lambda
	int s = 100; // number of updates per iteration
	int max_iter = 20000 + round(0.04 * MAT_YSIZE(*X) * MAT_YSIZE(*X)); // number of iterations
	double tol = 1e-5;
	int n = MAT_YSIZE(*X);
	Matrix2D<double> R;
	Matrix1D<double> D, Rt, W;
	MultidimArray<int> J;
	Matrix1D<int> ind1,ind2;

	ind1.resize(s);
	ind2.resize(s);
	Rt.initZeros(VEC_XSIZE(ind1));
	D.initZeros(VEC_XSIZE(ind1));
	W.initZeros(VEC_XSIZE(ind1));

	Y.initRandom(MAT_YSIZE(*X),outputDim,0,1);

	// Compute distance all vs. all
	computeDistance(*X, R, distance, true);
	R /= R.computeMax() / sqrt(2);

	if (global)
		max_iter *= 3;
	int ind1_0=0, ind1_F=s;
	int ind2_F=std::min(2*s,n-1);
	for(int nn=0;nn<max_iter;nn++)
	{
		// Compute random points
		randomPermutation(n,J);

		// Get the ind1 and ind2 indexes vectors
		if(n>s){
			memcpy(&VEC_ELEM(ind1, 0),&A1D_ELEM(J,ind1_0),ind1_F*sizeof(int));
			memcpy(&VEC_ELEM(ind2, 0),&A1D_ELEM(J,ind1_F),ind1_F*sizeof(int));
		}
		else{
			memcpy(&VEC_ELEM(ind1, 0),&A1D_ELEM(J,0),ind2_F*sizeof(int));
			memcpy(&VEC_ELEM(ind2, 0),&A1D_ELEM(J,ind1_F),ind2_F*sizeof(int));
		}

		// Compute distances between points in embedded space
		computeRandomPointsDistance(Y,D,ind1,ind2,distance,true);

		// Get corresponding distances in real space
		for (int l=0;l<VEC_XSIZE(ind1);l++)
			VEC_ELEM(Rt,l)= R(VEC_ELEM(ind1,l),VEC_ELEM(ind2,l));

		// Compute (Rt-D) / (D + tol)
		FOR_ALL_ELEMENTS_IN_MATRIX1D(Rt)
			VEC_ELEM(W,i) = (VEC_ELEM(Rt,i)-VEC_ELEM(D,i))/(VEC_ELEM(D,i)+tol);

		// Get the index to update locations
		for(int ii=0;ii<s;ii++)
		{
			int i1 = VEC_ELEM(ind1,ii);
			int i2 = VEC_ELEM(ind2,ii);

		    // Update locations
			double factor=lambda * 1/2;
			for (size_t j = 0; j< outputDim; ++j){
				double aux= MAT_ELEM(Y, i1, j) + factor*VEC_ELEM(W,ii)*(MAT_ELEM(Y, i1,j)-MAT_ELEM(Y, i2,j));
				MAT_ELEM(Y, i2, j) += factor*VEC_ELEM(W,ii)*(MAT_ELEM(Y,i2,j)-MAT_ELEM(Y,i1,j));
				MAT_ELEM(Y,i1,j) = aux;
			}
		}
		// Update lambda
		lambda = lambda - (lambda / max_iter);
	}
}
