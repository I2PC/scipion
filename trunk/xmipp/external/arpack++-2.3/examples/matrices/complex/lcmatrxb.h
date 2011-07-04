/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE LCMatrxB.h
   Function template for the tridiagonal matrix derived from 
   the standard central difference of the 1-d convection diffusion 
   operator u" + rho*u' on the interval [0, 1] with zero
   Dirichlet boundary conditions.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef LCMATRXB_H
#define LCMATRXB_H

#include "arcomp.h"

template<class ARFLOAT, class ARINT>
void CompMatrixB(ARINT n, arcomplex<ARFLOAT> rho, ARINT& nnz, 
                 arcomplex<ARFLOAT>* &A, ARINT* &irow, ARINT* &pcol) 
{

  ARINT              i, j;
  arcomplex<ARFLOAT> dd, dl, du, s, h, h2;

  // Defining constants.

  const arcomplex<ARFLOAT> one( 1.0, 0.0);
  const arcomplex<ARFLOAT> two( 2.0, 0.0);

  h  = one/arcomplex<ARFLOAT>((ARFLOAT)(n+1),0.0);
  h2 = h*h;
  s  = rho/two;
  dd = two/h2;
  dl = -(one/h2) - s/h;
  du = -(one/h2) + s/h;

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

} //  CompMatrixB.


#endif // LCMATRXB_H

