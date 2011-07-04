/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE DSMatrxB.h
   Function template for the one dimensional discrete Laplacian
   on the interval [0, 1], with zero Dirichlet boundary conditions.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef DSMATRXB_H
#define DSMATRXB_H

#include <math.h>

template<class ARFLOAT, class ARINT>
void DenseMatrixB(ARINT n, ARFLOAT* &A, char uplo = 'L')
{

  // Defining internal variables.

  ARINT    i, j;
  ARFLOAT  h, df, dd;

  // Defining constants.

  h  = ARFLOAT(n+1);
  dd = 2.0*h;
  df = -h;

  // Creating matrix A.

  A  = new ARFLOAT[(n*n+n)/2];

  if (uplo == 'L') {

    for (i=0, j=0; i<n; j+=(n-(i++))) {

      A[j] = dd;
      if (i != (n-1)) {
        A[j+1] = df;
      }
    
    }

  }
  else { // uplo = "U"

    for (i=0, j=0; i<n; j+=(++i)) {

      A[j+i] = dd;
      if (i != 0) {
        A[j+i-1] = df;
      }
    
    }  

  }

} // DenseMatrixB.

#endif // DSMATRXB_H

