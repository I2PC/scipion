/***************************************************************************
 *
 * Authors:    Clara M. Mejias Perez      clrmejias7@gmail.com (2013)
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

#include "chartingmanifold.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

//Inter MoFA using information from EM algorithm in MPPCA
void intermfa(const Matrix2D<double> &X, Matrix2D<double> &LX, Matrix2D<double> &MX, Matrix1D<double> &PX, Matrix2D<double> &R,
		      MultidimArray<double> &Z, int no_analyzers)
{
    //Initialize some parameters
    int N=MAT_YSIZE(X);   /*column size*/
    R.initZeros(no_analyzers,N);

    //Estimate the cluster centers based on input (= E-step)
    Matrix1D<double> pii=1/PX;
}

// mmpca runs EM algorithm and computes local factor analyzers
void mppca (const Matrix2D<double> &X, int no_dims, int no_analyzers, double tol, int maxiter, double minstd,
            Matrix2D<double> &LX, Matrix2D<double> &MX, Matrix1D<double> &PX)
{
    // Initialize some variables:
    int D=MAT_YSIZE(X);       /* row size */
    int N=MAT_XSIZE(X);       /* column size */

    // Randomly initialize factor analyzers
    Matrix1D<double> MM;
    X.computeRowMeans(MM);
    Matrix2D<double> SS,CC,RND1,MUL,X2;
    matrixOperation_AtA(X,SS);
    cholesky(SS,CC);

    RND1.initGaussian(D, no_analyzers);
    MUL=RND1*CC.transpose();
    MultidimArray<double> RND2;

    //Initialize variables
    double lik= -1e-38;

    MultidimArray<double> CZZ, ZZ;
    CZZ.initZeros(no_dims,no_dims,no_analyzers);
    ZZ.initZeros(no_dims,N,no_analyzers);

    Matrix2D<double> iltpil,COVZ,M1,MEANZ, R;
    R.initZeros(no_analyzers, N);
    for (int i=0; i<maxiter; i++)
    {
        //E step
        for (int k=0; k< no_analyzers; ++k)
        	;

        // Update likelihood of estimation
        float oldlik=lik;
        if (fabs(oldlik - lik)<tol)
            break;

        // M step
    }
}

void ChartingManifold::setSpecificParameters(int max_iterations, int no_analyzers)
{
    this->max_iterations=max_iterations;
    this->no_analyzers=no_analyzers;
}

// m=Z(:,i,j)*Z(:,i,j)^t
void multiplyZijZijt(const MultidimArray<double>& Z, size_t i, size_t j, Matrix2D<double> &m)
{
    m.resizeNoCopy(ZSIZE(Z),ZSIZE(Z));
    for (size_t k1=0; k1<ZSIZE(Z); ++k1)
        for (size_t k2=0; k2<ZSIZE(Z); ++k2)
            MAT_ELEM(m,k1,k2)=DIRECT_A3D_ELEM(Z,k1,i,j)*DIRECT_A3D_ELEM(Z,k2,i,j);
}

void ChartingManifold::reduceDimensionality()
{
    Matrix1D<double> XrowMeans, PX;
    X->computeRowMeans(XrowMeans);

    //Initialize some parameters
    double tol=1e-10;
    double min_std=1e-3;

    //Construct MFA model on data
    Matrix2D<double> MX,LX,R;
    mppca(*X,outputDim,no_analyzers,tol,max_iterations,min_std,LX,MX,PX);
    MultidimArray<double> Z;

    // Construct blockdiagonal matrix D
    Matrix2D<double> D, Ds, TRA;
    D.initZeros((outputDim+1)*no_analyzers,(outputDim+1)*no_analyzers);
    for (int i=0;i<no_analyzers;++i)
    {
        Ds.initZeros(outputDim+1,outputDim+1);
        for (size_t j=1; j<MAT_XSIZE(*X); ++j)
        {
            multiplyZijZijt(Z,i,j,TRA);
            double Rij=MAT_ELEM(R,i,j);
            FOR_ALL_ELEMENTS_IN_MATRIX2D(Ds)
            MAT_ELEM(Ds,i,j)+=Rij*MAT_ELEM(TRA,i,j);
        }
    }

    // Construct responsibility weighted local representation matrix U
    R.resize(1,no_analyzers, MAT_YSIZE(*X));
}
