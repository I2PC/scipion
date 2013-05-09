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

#include "lltsaSCG.h"

void eraseLastNColumns(Matrix2D<double> &A, int N)
{
	Matrix2D<double> Ap;
	Ap.resize(MAT_YSIZE(A),MAT_XSIZE(A)-N);
    for (size_t i = 0; i < MAT_YSIZE(A); ++i)
    	memcpy(&MAT_ELEM(Ap,i,0),&MAT_ELEM(A,i,0),MAT_XSIZE(Ap)*sizeof(double));
    A=Ap;
}

void LLTSASCG::reduceDimensionality()
{
	Matrix2D<double> B, XtBX, XtX;
	computeAlignmentMatrix(B);

    Matrix1D<double> DEigs;
    matrixOperation_AtB(*X, B, XtBX);
    matrixOperation_AtA(*X, XtX);
    XtBX = XtBX * (*X);
    generalizedEigs(XtBX, XtX, DEigs, A);
    eraseLastNColumns(A, MAT_XSIZE(A) - outputDim);
    Y = *X * A;
}
