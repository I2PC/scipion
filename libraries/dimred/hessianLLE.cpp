/***************************************************************************
 *
 * Authors:    Emma Sesmero      emmasesmero@gmail.com (2013)
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

#include "hessianLLE.h"

void HessianLLE::setSpecificParameters(int kNeighbours)
{
    this->kNeighbours = kNeighbours;
}

void HessianLLE::reduceDimensionality()
{
    Matrix2D<int> neighboursMatrix;
    Matrix2D<double> distanceNeighboursMatrix;

    kNearestNeighbours(*X, kNeighbours, neighboursMatrix, distanceNeighboursMatrix);

    size_t sizeY = MAT_YSIZE(*X);
    size_t dp = outputDim * (outputDim+1)/2;
    Matrix2D<double> weightMatrix, thisX, U, V, Vpr, Yi, Yi_complete, Yt, R, Pii;
    Matrix1D<double> D, vector;

    weightMatrix.initZeros(dp*sizeY,sizeY);

    for(size_t index=0; index<MAT_YSIZE(*X);++index)
    {
        extractNearestNeighbours(*X, neighboursMatrix, index, thisX);
        subtractColumnMeans(thisX);
        thisX = thisX.transpose();
        svdcmp(thisX, U, D, Vpr); // thisX = U * D * Vpr^t

        // Copy the first columns of Vpr onto V
        V.resizeNoCopy(MAT_YSIZE(Vpr),outputDim);
        for (size_t y=0; y<MAT_YSIZE(V); ++y)
        	memcpy(&MAT_ELEM(V,y,0),&MAT_ELEM(Vpr,y,0),outputDim*sizeof(double));

        //Basically, the above is applying PCA to the neighborhood of Xi.
        //The PCA mapping that is found (and that is contained in V) is an
        //approximation for the tangent space at Xi.

        //Build Hessian estimator
        buildYiHessianEstimator(V,Yi,outputDim,dp);
        completeYt(V, Yi, Yt);
        orthogonalizeColumnsGramSchmidt(Yt);

        //Get the transpose of the last columns
        size_t indexExtra = outputDim+1;
        size_t Ydim = MAT_XSIZE(Yt)-indexExtra;
        Pii.resizeNoCopy(Ydim,MAT_YSIZE(Yt));
        FOR_ALL_ELEMENTS_IN_MATRIX2D(Pii)
        	MAT_ELEM(Pii,i,j) = MAT_ELEM(Yt,j,indexExtra+i);

        //Double check weights sum to 1
        for(size_t j=0; j<dp; j++)
        {
        	Pii.getRow(j,vector);
        	double sum = vector.sum();
        	if(sum > 0.0001)
        		vector*=1.0/sum;

        	//Fill weight matrix
          	for(int k = 0; k<kNeighbours; k++){
        		size_t neighbourElem = MAT_ELEM(neighboursMatrix,index,k);
        		MAT_ELEM(weightMatrix,index*dp+j,neighbourElem) = VEC_ELEM(vector,k);
          	}
        }
    }
  	Matrix2D<double> G;
  	matrixOperation_AtA(weightMatrix,G);

  	Matrix1D<double> v;
  	eigsBetween(G,1,outputDim,v,Y);
  	Y*=sqrt(sizeY);
}

void HessianLLE::completeYt(const Matrix2D<double> &V,
		const Matrix2D<double> &Yi, Matrix2D<double> &Yt_complete)
{
    size_t Xdim = 1+MAT_XSIZE(V)+MAT_XSIZE(Yi);
    size_t Ydim = MAT_YSIZE(Yi);
    Yt_complete.resizeNoCopy(Ydim, Xdim);

    for (size_t i=0; i<Ydim; ++i)
    {
    	MAT_ELEM(Yt_complete,i,0)=1.;
    	memcpy(&MAT_ELEM(Yt_complete,i,1),             &MAT_ELEM(V,i,0), MAT_XSIZE(V)*sizeof(double));
    	memcpy(&MAT_ELEM(Yt_complete,i,MAT_XSIZE(V)+1),&MAT_ELEM(Yi,i,0),MAT_XSIZE(Yi)*sizeof(double));
    }
}

void HessianLLE::buildYiHessianEstimator(const Matrix2D<double> &V,
		Matrix2D<double> &Yi, size_t no_dim, size_t dp)
{
    Matrix1D<double> startp;
    Matrix1D<double> vector;

    size_t ct = 0;
    Yi.resizeNoCopy(MAT_YSIZE(V),dp);

    for(size_t mm=0; mm<no_dim; mm++)
    {
        V.getCol(mm,startp);

        size_t length = no_dim-mm;
        size_t indle=mm;
        for(size_t nn=0; nn<length; nn++)
        {
        	V.getCol(indle,vector);
            size_t column = ct+nn;
            for(size_t element = 0; element<MAT_YSIZE(V); element++)
                MAT_ELEM(Yi, element, column) = VEC_ELEM(startp, element)*VEC_ELEM(vector, element);
            ++indle;
        }
        ct += length;
    }
}
