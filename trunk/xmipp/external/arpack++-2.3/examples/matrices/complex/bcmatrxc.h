/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE BCMatrxC.h
   Function template for the mass matrix formed by using piecewise 
   linear elements on [0, 1].

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef BCMATRXC_H
#define BCMATRXC_H

#include "arcomp.h"

template<class ARFLOAT, class ARINT>
void CompMatrixC(ARINT n, ARINT& nL, ARINT& nU, arcomplex<ARFLOAT>* &A)
{

  ARINT              i;
  arcomplex<ARFLOAT> diag, sub;

  // Defining constants.

  sub  = arcomplex<ARFLOAT>(1.0/(ARFLOAT)(n+1), 0.0);
  diag = arcomplex<ARFLOAT>(4.0/(ARFLOAT)(n+1), 0.0);

  // Defining the lower and upper bandwidth.

  nL = 1;
  nU = 1;

  // Creating output vector A.

  A  = new arcomplex<ARFLOAT>[3*n];
  for (i=1; i<(3*n); i+=3) {
    if (i-1)   A[i-1] = sub;
               A[i]   = diag;
    if (n-i-1) A[i+1] = sub;
  }

} // CompMatrixC.

#endif // BCMATRXC_H

