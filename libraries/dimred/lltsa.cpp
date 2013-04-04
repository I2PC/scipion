/***************************************************************************
 *
 * Authors:    Pablo Perez Tirador            p.perez44@usp.ceu.es (2013)
 *
 * Escuela Politecnica Superior, Universidad San Pablo - CEU
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

#include "lltsa.h"

void LLTSA::setSpecificParameters(int neighbourNumber)
{
	this->neighbourNumber=neighbourNumber;
}

void LLTSA::reduceDimensionality()
{
	Matrix2D<double> D;
	Matrix2D<int> ni;

	Y.resizeNoCopy(MAT_YSIZE(*X),outputDim);

	// Step 1 - Make columns in original data zero-mean
	subtractColumnMeans(*X);

	// Step 2 - Find nearest neighbours
	kNearestNeighbours(*X, neighbourNumber, ni, D);

	// Step 3 - Compute local information matrix
	Matrix2D<double> Xi,W,Vi,Si,Gi, Bi, B;
	Matrix1D<int> Ji;
	Matrix1D<double> s;
	B.initZeros(MAT_YSIZE(*X),MAT_YSIZE(*X));
	double iNeighbourNumber=1.0/sqrt(neighbourNumber);
	for(size_t n=0; n<MAT_YSIZE(*X); n++)
	{
        extractNearestNeighbours(*X, ni, n, Xi);
		subtractColumnMeans(Xi);
        matrixOperation_AAt(Xi, W); // W=X*X^t
        schur(W, Vi, Si);           // W=Vi*Si*Vi^t

		Si.getDiagonal(s);
		s.indexSort(Ji);

		if (VEC_XSIZE(Ji)<outputDim)
			outputDim=VEC_XSIZE(Ji);

        // Approximate tangent space
		Gi.resizeNoCopy(MAT_YSIZE(Vi),outputDim+1);
		Gi.setConstantCol(0,iNeighbourNumber);
		FOR_ALL_ELEMENTS_IN_MATRIX1D(Ji) // Copy the columns of Vi pointed by Ji
		{
			size_t idxi=VEC_ELEM(Ji,i)-1;
			size_t i_1=i+1;
			for (int j=0; j<MAT_YSIZE(Vi); ++j)
				MAT_ELEM(Gi,i_1,j)=MAT_ELEM(Vi,idxi,j);
		}
		matrixOperation_AAt(Gi,Bi);
		matrixOperation_IminusA(Bi);

		for (int j1=0; j1<neighbourNumber; ++j1)
		{
			int ni1=MAT_ELEM(ni,n,j1);
			for (int j2=j1; j2<neighbourNumber; ++j2)
			{
				int ni2=MAT_ELEM(ni,n,j2);
				MAT_ELEM(B,ni1,ni2)+=MAT_ELEM(Bi,j1,j2);
				MAT_ELEM(B,ni2,ni1)=MAT_ELEM(B,ni1,ni2);
			}
		}
	}
}
