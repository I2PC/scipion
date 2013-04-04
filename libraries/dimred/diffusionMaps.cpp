/***************************************************************************
 *
 * Authors:    Carlos Oscar Sanchez Sorzano      coss@cnb.csic.es (2013)
 * 			   Tomas Bolgiani La Placa
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

#include "diffusionMaps.h"

void DiffusionMaps::setSpecificParameters(double t, double sigma)
{
	this->t=t;
	this->sigma=sigma;
}

void DiffusionMaps::reduceDimensionality()
{
	//Normalize data (between 0 and 1)
	normalizeColumnsBetween0and1(*X);

	// Compute Gaussian Kernel Matrix.
	// First, compute the distance of all vs all.
	Matrix2D<double> L2distance;
	computeDistance(*X,L2distance,distance,false);

	// Now use this distance to build a similarity matrix
	double auxOperator=-1/(2.0*(sigma*sigma));
	FOR_ALL_ELEMENTS_IN_MATRIX2D(L2distance)
		MAT_ELEM(L2distance,i,j)=exp(MAT_ELEM(L2distance,i,j)*auxOperator);

	// Normalize L2distance to be a stochastic matrix
	Matrix1D<double> p;
	L2distance.rowSum(p);
	if (t!=1.)
		FOR_ALL_ELEMENTS_IN_MATRIX2D(L2distance)
				MAT_ELEM(L2distance,i,j)/=pow(VEC_ELEM(p,i)*VEC_ELEM(p,j),t);

	// Normalize L2distance again
	L2distance.rowSum(p);
	FOR_ALL_ELEMENTS_IN_MATRIX1D(p)
		VEC_ELEM(p,i)=sqrt(VEC_ELEM(p,i));
	FOR_ALL_ELEMENTS_IN_MATRIX2D(L2distance)
			MAT_ELEM(L2distance,i,j)/=VEC_ELEM(p,i)*VEC_ELEM(p,j);

	// L2distance=U*S*V^t
	Matrix2D<double> U,V;
	Matrix1D<double> S;
	svdcmp(L2distance,U,S,V);

	// Get columns 1 to outputDim of U as output
	// normalzied by the first element in its row
	Y.resizeNoCopy(MAT_YSIZE(U),outputDim);
	for (size_t i=0;i<MAT_YSIZE(U);++i)
	{
		double iK=1/MAT_ELEM(U,i,0);
		for (size_t j=1;j<=outputDim;++j)
		    MAT_ELEM(Y,i,j-1)=MAT_ELEM(U,i,j)*iK;
	}
}
