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

// QR =====================================================================
int QR(Matrix2D<double> *F)
{
  size_t jQ=0;
  Matrix1D<double> qj1, qj2;
  int iBlockMax=MAT_XSIZE(*F)/4;

  for (size_t j1=0; j1<MAT_YSIZE(*F); j1++)
  {
    F->getRow(j1,qj1);
    // Project twice in the already established subspace
    // One projection should be enough but Gram-Schmidt suffers
    // from numerical problems
    for (int it=0; it<2; it++)
    {
      for (size_t j2=0; j2<jQ; j2++)
      {
        F->getRow(j2,qj2);

        // Compute dot product
        double s12=qj1.dotProduct(qj2);

        // Subtract the part of qj2 from qj1
        double *ptr1=&VEC_ELEM(qj1,0);
        const double *ptr2=&VEC_ELEM(qj2,0);
        for (int i=0; i<iBlockMax; i++)
        {
          (*ptr1++)-=s12*(*ptr2++);
          (*ptr1++)-=s12*(*ptr2++);
          (*ptr1++)-=s12*(*ptr2++);
          (*ptr1++)-=s12*(*ptr2++);
        }
        for (size_t i=iBlockMax*4; i<MAT_XSIZE(*F); ++i)
        (*ptr1++)-=s12*(*ptr2++);
      }
    }

    // Keep qj1 in Q if it has enough norm
    double Rii=qj1.module();
    if (Rii>1e-14)
    {
      // Make qj1 to be unitary and store in Q
      qj1/=Rii;
      F->setRow(jQ++,qj1);
    }
  }
  return jQ;
}

void HessianLLE::setSpecificParameters(int kNeighbours)
{
    this->kNeighbours = kNeighbours;
}

void HessianLLE::reduceDimensionality()
{
    Matrix2D<int> neighboursMatrix;
    Matrix2D<double> distanceNeighboursMatrix;
    kNearestNeighbours(*X, kNeighbours, neighboursMatrix, distanceNeighboursMatrix);

    size_t sizeX = MAT_XSIZE(*X);
    size_t dp = outputDim * (outputDim+1)/2;
    Matrix2D<double> weightMatrix(dp*sizeX,sizeX), thisX, U, V, Vpr, Yi, Yi_complete, Yt, R, Pii;
    Matrix1D<double> D;

    std::cout<<"Building Hessian estimator for neighboring points...\n";
    for(size_t i=0; i<MAT_YSIZE(*X);++i)
    {
        extractNearestNeighbours(*X, neighboursMatrix, i, thisX);
        subtractColumnMeans(thisX);
        svdcmp(thisX, U, D, Vpr); // thisX = U * D * Vpr^t

        if (MAT_XSIZE(Vpr)<outputDim)
        {
            outputDim = MAT_XSIZE(Vpr);
            dp = outputDim * (outputDim+1)/2;
            std::cout<<"Target dimensionality reduced to "<<outputDim<<"\n";
        }

        // Copy the first columns of Vpr onto V
        V.resizeNoCopy(MAT_YSIZE(Vpr),outputDim);
        for (size_t i=0; i<MAT_YSIZE(V); ++i)
        	memcpy(&MAT_ELEM(V,i,0),&MAT_ELEM(Vpr,i,0),outputDim*sizeof(double));

        //Basically, the above is applying PCA to the neighborhood of Xi.
        //The PCA mapping that is found (and that is contained in V) is an
        //approximation for the tangent space at Xi.

        //Build Hessian estimator
        buildHessianEstimator(V,Yi,outputDim,dp);

        completeYt(V, Yi, Yt);
        QR(&Yt);

        size_t indexExtra = outputDim+1;
        size_t Ydim =MAT_XSIZE(Yt)-indexExtra;
        Pii.resizeNoCopy(Ydim,MAT_YSIZE(Yt));
        FOR_ALL_ELEMENTS_IN_MATRIX2D(Pii)
        	MAT_ELEM(Pii,i,j)=MAT_ELEM(Yt,indexExtra+j,i);
    }
}

void HessianLLE::completeYt(const Matrix2D<double> &V,
		const Matrix2D<double> &Yt, Matrix2D<double> &Yt_complete)
{
    size_t Xdim = 1+MAT_XSIZE(V)+MAT_XSIZE(Yt);
    size_t Ydim = MAT_YSIZE(Yt);
    Yt_complete.resizeNoCopy(Ydim, Xdim);

    for (size_t i=0; i<Ydim; ++i)
    {
    	MAT_ELEM(Yt,i,0)=1.;
    	memcpy(&MAT_ELEM(Yt,i,1),             &MAT_ELEM(V,i,0), MAT_XSIZE(V)*sizeof(double));
    	memcpy(&MAT_ELEM(Yt,i,MAT_XSIZE(V)+1),&MAT_ELEM(Yt,i,0),MAT_XSIZE(Yt)*sizeof(double));
    }
}

void HessianLLE::buildHessianEstimator(const Matrix2D<double> &V,
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
