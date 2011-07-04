/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE DSMatrxC.h
   Function template for the mass matrix formed by using piecewise
   linear elements on the interval [0, 1].

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef DSMATRXC_H
#define DSMATRXC_H

#include <math.h>

template<class ARFLOAT, class ARINT>
void DenseMatrixC(ARINT n, ARFLOAT* &A, char uplo='L')

{

  // Defining internal variables.

  ARINT    i, j;
  ARFLOAT  h, df, dd;

  // Defining constants.

  h  = 1.0/ARFLOAT(n+1);
  dd = (4.0/6.0)*h;
  df = (1.0/6.0)*h;

  // Creating output vector A.

  A   = new ARFLOAT[(n*n+n)/2];

  if (uplo == 'L') {

    for (i=0, j=0; i<n; j+=(n-(i++))) {

      A[j] = dd;
      if (i != (n-1)) {
        A[j+1] = df;
      }
    
    }

  }
  else {

    for (i=0, j=0; i<n; j+=(++i)) {

      A[j+i] = dd;
      if (i != 0) {
        A[j+i-1] = df;
      }
    
    }  

  }
    
} // DenseMatrixC.

#endif // DSMATRXC_H

