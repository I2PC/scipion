/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE LCMatrxE.h
   Function template for the stiffness matrix formed by using
   piecewise linear elements on [0,1].

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef LCMATRXE_H
#define LCMATRXE_H

#include "arcomp.h"

template<class ARFLOAT, class ARINT>
void CompMatrixE(ARINT n, arcomplex<ARFLOAT> rho, ARINT& nnz, 
                 arcomplex<ARFLOAT>* &A, ARINT* &irow, ARINT* &pcol)
{

  ARINT              i, j;
  arcomplex<ARFLOAT> dd, dl, du, s, h;

  // Defining constants.

  const arcomplex<ARFLOAT> one( 1.0, 0.0);
  const arcomplex<ARFLOAT> two( 2.0, 0.0);

  h  = one/arcomplex<ARFLOAT>((ARFLOAT)(n+1),0.0);
  s  = rho/two;
  dd = two/h;
  dl = -(one/h) - s;
  du = -(one/h) + s;

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
      A[j++]  = du;
    }

    irow[j]   = i;
    A[j++]    = dd;

    if (i!=(n-1)) {
      irow[j] = i+1;
      A[j++]  = dl;
    }

    pcol[i+1] = j;

  }

} //  CompMatrixE.


#endif // LCMATRXE_H


