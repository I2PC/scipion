/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE DCMatrxB.h
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

#ifndef DCMATRXB_H
#define DCMATRXB_H

#include "arcomp.h"

template<class ARFLOAT, class ARINT>
void CompMatrixB(ARINT n, arcomplex<ARFLOAT> rho, arcomplex<ARFLOAT>* &A)
{

  ARINT              i, j;
  arcomplex<ARFLOAT> dd, dl, du, h;

  // Defining constants.

  const arcomplex<ARFLOAT> one( 1.0, 0.0);
  const arcomplex<ARFLOAT> two( 2.0, 0.0);
  const arcomplex<ARFLOAT> half( 0.5, 0.0);

  h  = one/arcomplex<ARFLOAT>((ARFLOAT)(n+1),0.0);
  dd = two/h;
  dl = -(one/h) - half*rho;
  du = -(one/h) + half*rho;

  // Creating output vector A.

  A  = new arcomplex<ARFLOAT>[n*n];

  for (i=0; i<n*n; i++) A[i] = arcomplex<ARFLOAT>(0.0, 0.0);

  for (i=0, j=0; i<n; i++, j+=n+1) {
    if (i)      A[j-1]  = du;
                A[j]    = dd;
    if (n-i-1)  A[j+1]  = dl;
  }

} //  CompMatrixB.


#endif // DCMATRXB_H

