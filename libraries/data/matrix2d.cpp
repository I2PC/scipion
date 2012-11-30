/***************************************************************************
 *
 * Authors:     Sjors H.W. Scheres (scheres@cnb.csic.es)
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

#include <algorithm>
#include "matrix2d.h"

/* Cholesky decomposition -------------------------------------------------- */
void cholesky(const Matrix2D<double> &M, Matrix2D<double> &L)
{
	L=M;
	Matrix1D<double> p;
	p.initZeros(MAT_XSIZE(M));
	choldc(L.adaptForNumericalRecipes2(), MAT_XSIZE(M), p.adaptForNumericalRecipes());
	FOR_ALL_ELEMENTS_IN_MATRIX2D(L)
	if (i==j)
		MAT_ELEM(L,i,j)=VEC_ELEM(p,i);
	else if (i<j)
		MAT_ELEM(L,i,j)=0.0;
}

/* Interface to numerical recipes: svbksb ---------------------------------- */
void svbksb(Matrix2D<double> &u, Matrix1D<double> &w, Matrix2D<double> &v,
            Matrix1D<double> &b, Matrix1D<double> &x)
{
    // Call to the numerical recipes routine. Results will be stored in X
    svbksb(u.adaptForNumericalRecipes2(),
           w.adaptForNumericalRecipes(),
           v.adaptForNumericalRecipes2(),
           u.mdimy, u.mdimx,
           b.adaptForNumericalRecipes(),
           x.adaptForNumericalRecipes());
}

// Solve linear systems ---------------------------------------------------
void solveLinearSystem(PseudoInverseHelper &h, Matrix1D<double> &result)
{
	Matrix2D<double> &A=h.A;
	Matrix1D<double> &b=h.b;
	Matrix2D<double> &AtA=h.AtA;
	Matrix2D<double> &AtAinv=h.AtAinv;
	Matrix1D<double> &Atb=h.Atb;

	// Compute AtA and Atb
	int I=MAT_YSIZE(A);
	int J=MAT_XSIZE(A);
	AtA.initZeros(J,J);
	Atb.initZeros(J);
	for (int i=0; i<J; ++i)
	{
		for (int j=0; j<J; ++j)
		{
			double AtA_ij=0;
			for (int k=0; k<I; ++k)
				AtA_ij+=MAT_ELEM(A,k,i)*MAT_ELEM(A,k,j);
			MAT_ELEM(AtA,i,j)=AtA_ij;
		}
		double Atb_i=0;
		for (int k=0; k<I; ++k)
			Atb_i+=MAT_ELEM(A,k,i)*VEC_ELEM(b,k);
		VEC_ELEM(Atb,i)=Atb_i;
	}

	// Compute the inverse of AtA
	AtA.inv(AtAinv);

	// Now multiply by Atb
	result.initZeros(J);
	FOR_ALL_ELEMENTS_IN_MATRIX2D(AtAinv)
		VEC_ELEM(result,i)+=MAT_ELEM(AtAinv,i,j)*VEC_ELEM(Atb,j);
}

// Solve linear systems ---------------------------------------------------
void weightedLeastSquares(WeightedLeastSquaresHelper &h, Matrix1D<double> &result)
{
	Matrix2D<double> &A=h.A;
	Matrix1D<double> &b=h.b;
	Matrix1D<double> &w=h.w;

	// See http://en.wikipedia.org/wiki/Least_squares#Weighted_least_squares
	FOR_ALL_ELEMENTS_IN_MATRIX1D(w)
	{
		double wii=sqrt(VEC_ELEM(w,i));
		VEC_ELEM(b,i)*=wii;
		for (int j=0; j<MAT_XSIZE(A); ++j)
			MAT_ELEM(A,i,j)*=wii;
	}
	solveLinearSystem(h,result);
}

// Solve linear system with RANSAC ----------------------------------------
void ransacWeightedLeastSquares(WeightedLeastSquaresHelper &h, Matrix1D<double> &result,
		double tol, int Niter, double outlierFraction)
{
	int N=MAT_YSIZE(h.A); // Number of equations
	int M=MAT_XSIZE(h.A); // Number of unknowns

	// Initialize a vector with all equation indexes
	std::vector<int> eqIdx;
	eqIdx.reserve(N);
	for (int n=0; n<N; ++N)
		eqIdx.push_back(n);
	int *eqIdxPtr=&eqIdx[0];

	// Resize a WLS helper for solving the MxM equation systems
	PseudoInverseHelper haux;
	haux.A.resizeNoCopy(M,M);
	haux.b.resizeNoCopy(M);
	Matrix2D<double> &A=haux.A;
	Matrix1D<double> &b=haux.b;

	// Solve Niter randomly chosen equation systems
	double bestError=1e38;
	const int Mfloat=M*sizeof(float);
	int minNewM=(int)((1.0-outlierFraction)*N-M);
	if (minNewM<0)
		minNewM=0;

	Matrix1D<double> resultAux;
	Matrix1D<int> idxIn(N);
	WeightedLeastSquaresHelper haux2;
	for (int it=0; it<Niter; ++it)
	{
        idxIn.initZeros();

        // Randomly select M equations
        std::random_shuffle(eqIdx.begin(), eqIdx.end());

        // Select the equation system
        for (int i=0; i<M; ++i)
        {
        	int idx=eqIdxPtr[i];
        	memcpy(&MAT_ELEM(A,i,0),&MAT_ELEM(h.A,idx,0),Mfloat);
        	VEC_ELEM(b,i)=VEC_ELEM(h.b,idx);
        	VEC_ELEM(idxIn,idx)=1;
        }

        // Solve the equation system
        // We use LS because the weight of some of the equations might be too low
        // and then the system is ill conditioned
        solveLinearSystem(haux, resultAux);

        // Study the residuals of the rest
        int newM=0;
        for (int i=M+1; i<N; ++i)
        {
        	int idx=eqIdxPtr[i];
        	double bp=0;
        	for (int j=0; j<M; ++j)
        		bp+=MAT_ELEM(h.A,idx,j)*VEC_ELEM(resultAux,j);
        	if (fabs(bp-VEC_ELEM(h.b,idx))<tol)
        	{
        		VEC_ELEM(idxIn,idx)=1;
        		++newM;
        	}
        }

        // If the model represent more points
        if (newM>minNewM)
        {
        	Matrix2D<double> &A2=haux2.A;
        	Matrix1D<double> &b2=haux2.b;
        	Matrix1D<double> &w2=haux2.w;
        	A2.resizeNoCopy(M+newM,M);
        	b2.resizeNoCopy(M+newM);
        	w2.resizeNoCopy(M+newM);

            // Select the equation system
        	int targeti=0;
            for (int i=0; i<N; ++i)
            	if (VEC_ELEM(idxIn,i))
            	{
					memcpy(&MAT_ELEM(A2,targeti,0),&MAT_ELEM(h.A,i,0),Mfloat);
					VEC_ELEM(b2,targeti)=VEC_ELEM(h.b,i);
					VEC_ELEM(w2,targeti)=VEC_ELEM(h.w,i);
					++targeti;
            	}

            // Solve it with WLS
            weightedLeastSquares(haux2, resultAux);

            // Compute the mean error
            double err=0;
            for (int i=0; i<M+newM; ++i)
			{
				double bp=0;
				for (int j=0; j<M; ++j)
					bp+=MAT_ELEM(A2,i,j)*VEC_ELEM(resultAux,j);
				err+=fabs(VEC_ELEM(b2,i)-bp)*VEC_ELEM(w2,i);
			}
            err/=(M+newM);
            if (err<bestError)
            {
            	bestError=err;
            	result=resultAux;
            }
        }
	}
}
