/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE LSMatrxA.h
   Function template for the matrix

                          | T -I          |
                          |-I  T -I       |
                      A = |   -I  T       |
                          |        ...  -I|
                          |           -I T|

   derived from the standard central difference discretization of the
   2-dimensional Laplacian on the unit square with zero Dirichlet 
   boundary conditions.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef LSMATRXA_H
#define LSMATRXA_H

#include <math.h>

template<class ARFLOAT, class ARINT>
void SymmetricMatrixA(ARINT nx, ARINT& n, ARINT& nnz, ARFLOAT* &A, 
                      ARINT* &irow, ARINT* &pcol, char uplo = 'L')

{

  // Defining internal variables.

  ARINT    i, j;
  ARFLOAT  h2, df, dd;

  // Defining constants.

  h2  = 1.0/(ARFLOAT(nx+1)*ARFLOAT(nx+1));
  dd  = 4.0/h2;
  df  = -1.0/h2;

  // Defining the number of columns and nonzero elements of matrix.

  n   = nx*nx;
  nnz = 3*n-2*nx;

  // Creating output vectors.

  A    = new ARFLOAT[nnz];
  irow = new ARINT[nnz];
  pcol = new ARINT[n+1];

  // Defining  matrix A.

  pcol[0] = 0;
  i       = 0;

  if (uplo == 'U') {

    for (j = 0; j < n; j++) {
      if (j >= nx) {
        A[i] = df;   irow[i++] = j-nx;
      }
      if ((j%nx) != 0) {
        A[i] = df;   irow[i++] = j-1;
      }
      A[i] = dd;     irow[i++] = j;
      pcol[j+1] = i;
    }

  }
  else {

    for (j = 0; j < n; j++) {
      A[i] = dd;     irow[i++] = j;
      if (((j+1)%nx) != 0) {
        A[i] = df;   irow[i++] = j+1;
      }
      if (j < n-nx) {
        A[i] = df;   irow[i++] = j+nx;
      }
      pcol[j+1] = i;
    }

  }

} // SymmetricMatrixA.

#endif // LSMATRXA_H

