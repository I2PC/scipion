/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE DCMatrxC.h
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

#ifndef DCMATRXC_H
#define DCMATRXC_H

#include "arcomp.h"

template<class ARFLOAT, class ARINT>
void CompMatrixC(ARINT n, arcomplex<ARFLOAT>* &A)
{

  ARINT              i, j;
  arcomplex<ARFLOAT> diag, sub;

  // Defining constants.

  sub  = arcomplex<ARFLOAT>(1.0/(ARFLOAT)(n+1), 0.0);
  diag = arcomplex<ARFLOAT>(4.0/(ARFLOAT)(n+1), 0.0);

  // Creating output vector A.

  A  = new arcomplex<ARFLOAT>[n*n];

  for (i=0; i<n*n; i++) A[i] = arcomplex<ARFLOAT>(0.0, 0.0);

  for (i=0, j=0; i<n; i++, j+=n+1) {
    if (i)      A[j-1]  = sub;
                A[j]    = diag;
    if (n-i-1)  A[j+1]  = sub;
  }

} // CompMatrixC.

#endif // DCMATRXC_H

