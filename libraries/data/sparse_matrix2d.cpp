/***************************************************************************
 * Authors:     Javier Vargas (jvargas@cnb.csic.es)
 *
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


#include "sparse_matrix2d.h"


// Sparse matrices --------------------------------------------------------
SparseMatrix2D::SparseMatrix2D(){
	N 		= 0;
}

SparseMatrix2D::SparseMatrix2D(std::vector<SparseElement> &_elements, int _Nelements)
{
	//First of all, we sort the elements by rows and then by columns
	std::sort(_elements.begin(),_elements.end());

	size_t ln = _elements.size();
	N      = _Nelements;

	values.resizeNoCopy(ln);
	jIdx.resizeNoCopy(ln);
	iIdx.resizeNoCopy(N);

	int actualRow = -1;
	int i         =  0; // Iterator for the vectors "values" and "jIdx"

	for (size_t k=0 ; k< ln ; ++k )// Iterator for vector of SparseElemets
	{
		if(_elements.at(k).value != 0.0) // Searching that there isn't any zero value
		{
			DIRECT_MULTIDIM_ELEM(values,i) = _elements.at(k).value;
			DIRECT_MULTIDIM_ELEM(jIdx,i)   = _elements.at(k).j +1;

			int rse = _elements.at(k).i;
			while( rse > actualRow )
			{
				actualRow++;
				if( rse == actualRow )
					// The first element in row number "x" is in position iIdx(x-1)
					DIRECT_MULTIDIM_ELEM(iIdx,actualRow) = i+1;

				else
					// If there isn't any nonzero value on this row, the value in this vector is 0
					DIRECT_MULTIDIM_ELEM(iIdx,actualRow) = 0;
			}
			++i;
		}
	}
}

/*
 * Fills the matrix form a vector of SparseElements
 */
void SparseMatrix2D::sparseMatrix2DFromVector(std::vector<SparseElement> &_elements)
{
	std::sort(_elements.begin(),_elements.end());

	size_t ln = _elements.size();
	N = _elements.at(ln-1).j;
	values.resizeNoCopy(ln);
	jIdx.resizeNoCopy(ln);
	iIdx.resizeNoCopy(N);

	int actualRow = -1;
	int i         =  0;

	for (size_t k=0 ; k< ln ; ++k )
	{
		if(_elements.at(k).value != 0.0) // Seaching that there isn't any zero value
		{
			DIRECT_MULTIDIM_ELEM(values,i) = _elements.at(k).value;
			DIRECT_MULTIDIM_ELEM(jIdx,i)   = _elements.at(k).j +1;

			int rse = _elements.at(k).i;
			while( rse > actualRow )
			{
				actualRow++;
				if( rse == actualRow )
					DIRECT_MULTIDIM_ELEM(iIdx,actualRow) = i+1;

				else
					DIRECT_MULTIDIM_ELEM(iIdx,actualRow) = 0;
			}
			++i;
		}
	}
	N = actualRow+1;
}

SparseMatrix2D &SparseMatrix2D::operator =(const SparseMatrix2D &X)
{
	if (this!=&X)
	{
		N=X.N;
		values=X.values;
		iIdx=X.iIdx;
		jIdx=X.jIdx;
	}
	return *this;
}

/**
 * It computes y <- this*x
 */
void SparseMatrix2D::multMv(double* x, double* y)
{
	int col, rowEnd, rowBeg;
	double val;
	memset(y,0,N*sizeof(double));

	size_t nnz=XSIZE(values);
    for(int i = 0; i< N; i++)
    {
    	// We get the gap where are the elements of the row i in value's vector
		rowBeg = DIRECT_MULTIDIM_ELEM(iIdx,i) -1;
		if( i != N-1 )
			rowEnd = DIRECT_MULTIDIM_ELEM(iIdx,i+1) -1;
		else
			rowEnd = nnz;

		val = 0.0;
		for(int j = rowBeg; j < rowEnd ; j++)
		{
			col = DIRECT_MULTIDIM_ELEM(jIdx,j) -1;// Column with a nonzero element in this row of the matrix
			val += DIRECT_MULTIDIM_ELEM(values,j) * x[col];
		}
		y[i] = val;
    }
}

/*
 * It shows the sparse matrix as a full matrix. If the sparse matrix is real big, you shoudn't use it
 * */
std::ostream & operator << (std::ostream &out, const SparseMatrix2D &X)
{
	int N=X.nrows();
	for(int i =0 ; i< N ; i++)
	{
		for(int j =0; j<N ; j++)
			out << X.getElemIJ(i, j) << "\t";
		out << std::endl;
	}
	return out;
}

/**
 * Get element
 */
double SparseMatrix2D::getElemIJ(int row, int col) const
{
	int rowBeg = DIRECT_MULTIDIM_ELEM(iIdx,row) -1;
	int rowEnd;

	if( row != N-1 )
		rowEnd = DIRECT_MULTIDIM_ELEM(iIdx,row+1) -1;
	else
		rowEnd = XSIZE(values);

	// If there is a non-zero element, the column is in jIdx
	for(int i = rowBeg; i < rowEnd ; i++)
		if( DIRECT_MULTIDIM_ELEM(jIdx,i)-1 == col )
			return DIRECT_MULTIDIM_ELEM(values,i);
	return 0.0;
}

/// Computes y=SparseMatrixThis*SparseMatrixX
void SparseMatrix2D::multMM(const SparseMatrix2D &X, SparseMatrix2D &Y)
{
	size_t nnz=XSIZE(X.values);

	Y.N = X.N;
 	Y.values.initZeros(nnz);
	Y.jIdx.initZeros(nnz);
	Y.iIdx.initZeros(N);

	size_t nnzY = 0;
	int Nelems  = X.N;
	for(int row = 0; row < Nelems ; ++row){
		int firstElemRow = 0;
		for(int col = 0; col < Nelems ; ++col){
			double aij = 0.0;
			for(int k = 0; k < Nelems ; ++k ){
				double val=getElemIJ(row, k);
				if(val != 0.0){ // If GetElemIJ(row, k) is 0, the product is 0.
					double val2=X.getElemIJ(k,col);
					aij += val * val2;
				}
			}

			if(aij != 0.0){ // If there is a nonzero element in (row,col) position, we include that in the matrix
				if(firstElemRow == 0)
					firstElemRow = nnzY +1;
				DIRECT_MULTIDIM_ELEM(Y.values,nnzY) = aij;
				DIRECT_MULTIDIM_ELEM(Y.jIdx,nnzY)   = col+1;
				++nnzY;
			}
		}
		// Indicates where is the first element of the row
		DIRECT_MULTIDIM_ELEM(Y.iIdx,row) = firstElemRow;
	}
	Y.values.resize(nnzY);
	Y.jIdx.resize(nnzY);
}

/// Computes y=SparseMatrixThis*SparseMatrix
/*
 * The matrix D is diagonal.
 *
 * | a11	a12		a13|   | d11	 0		 0	|	| a11*d1	a12*d1		a13*d1|
 * | a21	a22		a23| * |  0		d22		 0	| = | a21*d2	a22*d2		a23*d2|
 * | a31	a32		a33|   |  0		 0		a33 |	| a31*d3	a32*d3		a33*d3|
 *
 *	Vectors iIdx and jIdx are the same as the non diagonal matrix
 *	Vector Values is the same as the non diagonal matrix multiply each row by each d_row
 *
 *	The values are obtain from the last to the beginning.
 * *
 */
void SparseMatrix2D::multMMDiagonal(const MultidimArray<double> &D, SparseMatrix2D &Y)
{
	int actualRow = N-1;

	Y=*this;

	// We see where is the first element in the next row
	int until = DIRECT_MULTIDIM_ELEM(iIdx,actualRow) -1;

	// Obtain the value that we have to multiply by the first row.
	double dx = DIRECT_MULTIDIM_ELEM(D,actualRow);
	size_t nnz=XSIZE(values);
	for(int i = nnz-1; i >= 0 ; --i)
	{
		// Yij = Aij * Dii / Y = A => Yij = Yij * Dii
		DIRECT_MULTIDIM_ELEM(Y.values,i) *= dx;

		if(until == i && i > 0){
			--actualRow;
			// Search for a row with nonzero elements
			while( DIRECT_MULTIDIM_ELEM(iIdx,actualRow) == 0 )
				--actualRow;

			until = DIRECT_MULTIDIM_ELEM(iIdx,actualRow) -1;
			dx    = DIRECT_MULTIDIM_ELEM(D,actualRow);
		}
	}
}

/*
 * Load matrix */
void SparseMatrix2D::loadMatrix(const FileName &fn)
{
	std::ifstream fhIn;
    fhIn.open(fn.c_str());
    if (!fhIn)
       REPORT_ERROR(ERR_IO_NOTEXIST,fn);

	double dobVecSize, auxDob;
	int vectorSize = 0;
	fhIn >> dobVecSize;
	vectorSize = (int)dobVecSize ;

	std::vector<SparseElement> elems(vectorSize);
	for(int i =0; i< vectorSize; ++i){
		fhIn >> auxDob;
		elems.at(i).i = (int) auxDob -1;

		fhIn >> auxDob;
		elems.at(i).j = (int) auxDob -1;

		fhIn >> elems.at(i).value;
	}
	sparseMatrix2DFromVector(elems);
}

