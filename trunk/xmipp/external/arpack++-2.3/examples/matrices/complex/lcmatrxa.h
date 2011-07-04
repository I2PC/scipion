/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE LCMatrxA.h
   Function template for the nx*nx by nx*nx block tridiagonal matrix

                              | T -I          |
                              |-I  T -I       |
                         OP = |   -I  T       |
                              |        ...  -I|
                              |           -I T|

   derived from the standard central difference discretization
   of the 2 dimensional convection-diffusion operator
                       (Laplacian u) + rho*(du/dx) 
   on a unit square with zero boundary conditions.
   T is a nx by nx tridiagonal matrix with DD on the diagonal,
   DL on the subdiagonal, and DU on the superdiagonal.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef LCMATRXA_H
#define LCMATRXA_H

#include "arcomp.h"
#include "blas1c.h"

template<class ARFLOAT, class ARINT>
void CompMatrixA(ARINT nx, ARINT& nnz, arcomplex<ARFLOAT>* &A, 
                 ARINT* &irow, ARINT* &pcol)
{

  ARINT              i, j, k, id;
  arcomplex<ARFLOAT> h, h2, dd, dl, du, f;

  // Defining constants.

  const arcomplex<ARFLOAT> half(0.5,0.0);
  const arcomplex<ARFLOAT> one(1.0,0.0);
  const arcomplex<ARFLOAT> four(4.0,0.0);
  const arcomplex<ARFLOAT> rho(1.0e2,0.0);

  h   = one/arcomplex<ARFLOAT>((ARFLOAT)(nx+1),0);
  h2  = h*h;
  f   = -(one/h2);
  dd  = four/h2;
  dl  = f - half*rho/h;
  du  = f + half*rho/h;

  // Defining the number of nonzero matrix elements.

  nnz = (5*nx-4)*nx;

  // Creating output vectors.

  A    = new arcomplex<ARFLOAT>[nnz];
  irow = new ARINT[nnz];
  pcol = new ARINT[nx*nx+1];

  // Filling A, irow and pcol.

  pcol[0] = 0;
  j       = 0;
  id      = 0;

  for (k=0; k!=nx; k++) {
    for (i=0; i!=nx; i++) {

      if (k) {
        irow[j] = id-nx;
        A[j++]  = f; 
      }

      if (i) {
        irow[j] = id-1;
        A[j++]  = du;
      }

      irow[j]   = id;
      A[j++]    = dd;

      if (i!=(nx-1)) {
        irow[j] = id+1;
        A[j++]  = dl;
      }

      if (k!=(nx-1)) {
        irow[j] = id+nx;
        A[j++]  = f;
      }

      pcol[++id]= j;
    }
  }     
  
} // CompMatrixA.


#endif // LCMATRXA_H
