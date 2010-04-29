/*

CONDOR 1.06 - COnstrained, Non-linear, Direct, parallel Optimization
              using trust Region method for high-computing load,
              noisy functions
Copyright (C) 2004 Frank Vanden Berghen

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation version 2
of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

If you want to include this tools in any commercial product,
you can contact the author at fvandenb@iridia.ulb.ac.be

*/

#ifndef _INCLUDE_MATRIX_H
#define _INCLUDE_MATRIX_H

#include "Vector.h"
#include "VectorInt.h"

#ifndef NOMATRIXTRIANGLE
#include "MatrixTriangle.h"
#endif

class Matrix
{
protected:
    typedef struct MatrixDataTag
    {
        char **columnNames;
        int nLine ,nColumn, extColumn, extLine;
        int ref_count;
        double **p;
    } MatrixData;
    MatrixData *d;

    void init(int _nLine, int _nColumn, int _extLine, int _extColumn, MatrixData* d=NULL);
    void setExtSize(int _extLine, int _extColumn);
    void destroyCurrentBuffer();

public:

// creation & management of Matrix:
    Matrix(int _ligne=0,int _nColumn=0);
    Matrix(int _ligne,int _nColumn, int _extLine,int _extColumn);
    Matrix(const char *filename, char ascii=0);
    Matrix(Vector a, Vector b);  // a * b^T
    Matrix(Vector a);
    void save(char *filename,char ascii);
    void save(FILE *f,char ascii);
    void updateSave(char *saveFileName); // only binary
    void extendLine();
    void setNLine(int _nLine);
    void extendColumn();
    void setNColumn(int _nColumn);
    void setSize(int _nLine,int _nColumn);
    void exactshape();
    void print();
    void setColNames(char **c, int nc=0);

// allow shallow copy:
    ~Matrix();
    Matrix(const Matrix &A);
    Matrix& operator=( const Matrix& A );
    Matrix clone();
    void copyFrom(Matrix a);

// accessor method
    inline bool operator==( const Matrix& A )
    {
        return (A.d==d);
    }
    inline int nLine()
    {
        return d->nLine;
    };
    inline int nColumn()
    {
        return d->nColumn;
    };
    inline char **getColumnNames()
    {
        return d->columnNames;
    }
    inline double *operator [](int i)
    {
        return d->p[i];
    };
    inline operator double**() const
    {
        return d->p;
    };
    Vector getLine(int i, int n=0, int startCol=0);
    void getLine(int i, Vector r, int n=0, int startCol=0);
    Vector getColumn(int i, int n=0);
    void getColumn(int i, Vector r, int n=0);
    void getSubMatrix(Matrix R, int startL, int StartC, int nl=0, int nc=0);
    void setLine(int i, Vector v, int n=0);
    void setLines(int indexDest, Matrix Source, int indexSource=0, int number=0);
    void swapLines(int i, int j);
    int lineIndex(Vector r, int nn=0);
    void merge(Matrix m, int eliminateDoubles=1);
    void append(Vector tmp);

// simple math tools:
    void zero();
    void diagonal(double d);
    Matrix multiply(Matrix B);
    void multiplyInPlace(double d);
    void multiply(Vector R, Vector v);  // result in R
    void transposeAndMultiply(Vector R, Vector a);// result in R
    void multiply(Matrix R, Matrix a); // result in R
    void transposeAndMultiply(Matrix R, Matrix a);// result in R
    void multiplyByTranspose(Matrix R, Matrix a); // result in R
    void multiplyByDiagonalMatrix(Vector v);
    Vector multiply(Vector v);
    void addInPlace(Matrix B);
    void addMultiplyInPlace(double d, Matrix B);
    void addUnityInPlace(double d);
    void transposeInPlace();
    Matrix transpose();
    void transpose(Matrix trans);
    double scalarProduct(int nl, Vector v);

#ifndef NOMATRIXTRIANGLE
    Matrix(MatrixTriangle A, char bTranspose=0);
    bool cholesky(MatrixTriangle matL, double lambda=0, double *lambdaCorrection=NULL);
    void choleskySolveInPlace(Vector b);
    void QR(Matrix Q=Matrix::emptyMatrix, MatrixTriangle R=MatrixTriangle::emptyMatrixTriangle,
            VectorInt permutation=VectorInt::emptyVectorInt);
#endif

    double frobeniusNorm();
    double LnftyNorm();
    double euclidianNorm(int i);
    Vector getMaxColumn();

// default return matrix in case of problem in a function
    static Matrix emptyMatrix;
};

#endif
