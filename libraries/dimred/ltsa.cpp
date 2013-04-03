/***************************************************************************
 *
 * Authors:    Sergio Calvo Gonzalez       sergiocg90@gmail.com (2013)
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

#include "ltsa.h"


void LTSA::setSpecificParameters(int k)
{
    this->k = k;
}

void computeWeightsVector (const Matrix2D<double> &A, Matrix1D<int> &weightVector)
{
    weightVector.initZeros(MAT_XSIZE(A));
    for (size_t i = 0; i < MAT_XSIZE(A); ++i)
    {
        double diagelem = MAT_ELEM(A, i, i);
        for (size_t j = 0; j < MAT_XSIZE(A); ++j)
            if (((i != j) && (MAT_ELEM(A, j, j) > diagelem)) ||  ((MAT_ELEM(A, j, j) == diagelem) && (j < i)))
                VEC_ELEM(weightVector,i)++;
    }
}

void getLessWeightNColumns (const Matrix2D<double> &A, const Matrix1D<int> &weightVector, Matrix2D<double> &B)
{
    size_t outputDim = MAT_XSIZE(B) - 1;
    for (size_t index = 0; index < outputDim; ++index)
        FOR_ALL_ELEMENTS_IN_MATRIX1D(weightVector)
            if (VEC_ELEM(weightVector,i) == index)
            {
                for (size_t j = 0; j < MAT_YSIZE(A); ++j)
                    MAT_ELEM(B, j, index + 1) = MAT_ELEM(A, j, i);
                break;
            }
}

void LTSA::reduceDimensionality()
{
    subtractColumnMeans(*X);

    size_t n = MAT_YSIZE(*X);
    Matrix2D<int> ni;
    Matrix2D<double> D, Xi(MAT_XSIZE(ni), MAT_XSIZE(*X)), B(n, n), W, Vi, Vi2, Si, Gi;
    kNearestNeighbours(*X, k, ni, D);

    B.initIdentity();
    Matrix1D<int> weightVector;
    for (size_t iLoop = 0; iLoop < n; ++iLoop)
    {
        extractNearestNeighbours(*X, ni, iLoop, Xi);
        subtractColumnMeans(Xi);

        matrixOperation_AtA(Xi, W); // W=X^t*X
        schur(W, Vi, Si);           // W=Vi*Si*Vi^t

        if (MAT_XSIZE(Si) < outputDim)
        {
            outputDim = MAT_XSIZE(Si);
            std::cout<<"Target dimensionality reduced to "<<outputDim;
        }

        computeWeightsVector(Si, weightVector);

        Vi2.resizeNoCopy(MAT_YSIZE(Vi), outputDim + 1);
        Vi2.setConstantCol(0, 1/sqrt(k)); //Vi2(0,:)=1/sqrt(k)
        getLessWeightNColumns(Vi, weightVector, Vi2);

        matrixOperation_AtA(Vi2, Gi); // Gi=Vi2^t*Vi2
        matrixOperation_IminusA(Gi);  // Gi=I-Gi

        // Compute partial B with correlation matrix Gi
        FOR_ALL_ELEMENTS_IN_MATRIX2D(Gi)
            MAT_ELEM(B, MAT_ELEM(ni,iLoop,i), MAT_ELEM(ni,iLoop,j)) += MAT_ELEM(Gi, i, j);
        MAT_ELEM(B, iLoop, iLoop)-=1;
    }

    Matrix1D<double> DEigs;
    lastEigs(B,outputDim+1,DEigs,Y);

    if (MAT_XSIZE(Y) < outputDim + 1)
        outputDim = MAT_XSIZE(Y) - 1;

    eraseFirstColumn(Y);
}
