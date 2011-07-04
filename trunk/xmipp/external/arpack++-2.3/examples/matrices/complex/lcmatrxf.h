/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE LCMatrxF.h
   Function template for a tridiagonal complex matrix.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef LCMATRXF_H
#define LCMATRXF_H

#include "arcomp.h"

template<class ARFLOAT, class ARINT>
void CompMatrixF(ARINT n, ARINT& nnz, arcomplex<ARFLOAT>* &A, 
                 ARINT* &irow, ARINT* &pcol)
{

  ARINT              i, j;
  arcomplex<ARFLOAT> h, dd, ds;

  // Defining constants.

  const arcomplex<ARFLOAT> one(1.0, 0.0);
  const arcomplex<ARFLOAT> four(4.0, 0.0);

  h  = one/arcomplex<ARFLOAT>((ARFLOAT)(n+1),0.0);
  dd = four*h;
  ds = one*h;

  // Defining the number of nonzero matrix elements.

  nnz = 3*n-2;

  // Creating output vectors.

  A    = new arcomplex<ARFLOAT>[nnz];
  irow = new ARINT[nnz];
  pcol = new ARINT[n+1];

  // Filling A, irow and pcol.

  pcol[0] = 0;
  j       = 0;

  for (i=0; i!=n; i++) {

    if (i) {
      irow[j] = i-1;
      A[j++]  = ds;
    }

    irow[j]   = i;
    A[j++]    = dd;

    if (i!=(n-1)) {
      irow[j] = i+1;
      A[j++]  = ds;
    }

    pcol[i+1] = j;

  }

} // CompMatrixF.


#endif // LCMATRXF_H

