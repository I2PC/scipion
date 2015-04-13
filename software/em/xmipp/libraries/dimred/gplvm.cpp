/***************************************************************************
 *
 * Authors:    Francisco Sanz Encinas      franciscosanz89@gmail.com (2013)
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

#include "gplvm.h"
#include <data/numerical_tools.h>

void GPLVM::setSpecificParameters(double sigma)
{
	this->sigma=sigma;
}

double GPLVM::objectiveFunction()
{
	sumY2.initZeros(MAT_YSIZE(Y));
	FOR_ALL_ELEMENTS_IN_MATRIX2D(Y)
		VEC_ELEM(sumY2,i)+=MAT_ELEM(Y,i,j)*MAT_ELEM(Y,i,j);

	double aux=1.0/(2*sigma*sigma);
	K.resizeNoCopy(MAT_YSIZE(Y),MAT_YSIZE(Y));
	FOR_ALL_ELEMENTS_IN_MATRIX2D(K)
	{
		double aux2=0;
		for (int k=0; k<MAT_XSIZE(Y); k++)
			aux2+=MAT_ELEM(Y,i,k)*MAT_ELEM(Y,j,k);
        MAT_ELEM(K,i,j)=exp((2*aux2-VEC_ELEM(sumY2,i)-VEC_ELEM(sumY2,j))*aux);
	}

	matrixOperation_AAt(*X,tmp);
	tmp= K.inv() * tmp;

	double d=MAT_YSIZE(*X);
	double n=MAT_XSIZE(*X);
	return(-(d*n/2) * log( 2 * M_PI) - (n/2) * log(K.det()+2.2e-308) - 0.5 * tmp.trace());
}

double gplvmObjectiveFuntion(double *p, void *prm)
{
	GPLVM *gpvlm=(GPLVM *)prm;
	Matrix2D<double> &Y=gpvlm->Y;
	memcpy(&MAT_ELEM(Y,0,0),&(p[1]),MAT_XSIZE(Y)*MAT_YSIZE(Y)*sizeof(double));
	double c=gpvlm->objectiveFunction();
	return c;
}

void GPLVM::reduceDimensionality()
{
	// Compute the first guess
	PCA::reduceDimensionality();

	// Refine
	Matrix1D<double> pY(MAT_XSIZE(Y)*MAT_YSIZE(Y)), steps(MAT_XSIZE(Y)*MAT_YSIZE(Y));
	steps.initConstant(1);
	memcpy(&VEC_ELEM(pY,0),&MAT_ELEM(Y,0,0),VEC_XSIZE(pY)*sizeof(double));
	double goal;
	int iter;
	powellOptimizer(pY,1,VEC_XSIZE(pY),&gplvmObjectiveFuntion,this,0.01,goal,iter,steps,true);
	memcpy(&MAT_ELEM(Y,0,0),&VEC_ELEM(pY,0),VEC_XSIZE(pY)*sizeof(double));
}
