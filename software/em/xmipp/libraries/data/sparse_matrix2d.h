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

#ifndef SPARSE_MATRIX2D_H_
#define SPARSE_MATRIX2D_H_

#include "multidim_array.h"

/** @ingroup Matrices
 */
//@{

/** Square, sparse matrices.
 *
 * To create a Sparse Matrix we need a vector with SparseElements that contains
 * a value, its row position and its column position.
 *
 * We store that information with a Compressed Row Storage (CRS).
 * This method stores a vector with the non-zero values (values), their column position (jIdx)
 * and another vector with the position on values's vector of the element
 * which is the first non-zero elements of each row (iIdx). The value 0 on iIdx vector means that
 * there aren't nonzero values in that row.
 * iIdx and jIdx have the first element in position 1. Value 0 is for no elements on a row.
 *
 * i.e.:

 Vector values
 _4_, 3, _1_, 2, _1_, 4, _3_

 Vector jIdX
  1,  2,  1,  2,  2,  4,  4

 Vector iIdx
 1, 3, 5, 7

 Matriz A:
 _4_   3  0  0
 _1_   2  0  0
  0  _1_  0  4
  0   0  0 _3_
 */
class SparseMatrix2D
{
public:
    /// The matrix is of size NxN
    int N;

    /// List of i positions
    MultidimArray<int>    iIdx;
    /// List of j positions
    MultidimArray<int>    jIdx;
    /// List of values
    MultidimArray<double> values;
public:
    /// Y size of the matrix
    int nrows() const
    {
        return N;
    }

    /// X size of the matrix
    int ncols() const
    {
        return N;
    }

    /// Empty constructor
    SparseMatrix2D();

    /** Constructor from a set of i,j indexes and their corresponding values.
     * N is the total dimension of the square, sparse matrix.
     */
    SparseMatrix2D(std::vector<SparseElement> &_elements, int _Nelements);

    /** Assig operator *this=X */
    SparseMatrix2D &operator =(const SparseMatrix2D &X);

    /// Fill the sparse matrix A with the elements of the vector.
    void sparseMatrix2DFromVector(std::vector<SparseElement> &_elements);

    /** Computes y=this*x
     * y and x are vectors of size Nx1
     */
    void multMv(double* x, double* y);

    /// Computes Y=this*X
    void multMM(const SparseMatrix2D &X, SparseMatrix2D &Y);

    /** Computes Y=this*D where D is a diagonal matrix.
     * It is assumed that the size of this sparse matrix is NxN and that the length of y is N.
    */
    void multMMDiagonal(const MultidimArray<double> &D, SparseMatrix2D &Y);

    /// Shows the dense Matrix associated
    friend std::ostream & operator << (std::ostream &out, const SparseMatrix2D &X);

    /** Returns the value on position (row,col).
     * The matrix goes from 0 to N-1
     */
    double getElemIJ(int row, int col) const;

    /** Loads a sparse matrix from a file.
     * Loads a SparseMatrix2D from a file which has the following format:
     *
     * sizeOfMatrix
     * 1 1  value
     * 2 1 value
     * . . .
     * i j value
     * . . .
     * last_i last_j last_value
     */
    void loadMatrix(const FileName &fn);
};
//@}

#endif /* SPARSE_MATRIX2D_H_ */
