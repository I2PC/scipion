/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE BCMatrxA.h
   Function template for the nx*nx by nx*nx block tridiagonal matrix

                              | T -I          |
                              |-I  T -I       |
                         OP = |   -I  T       |
                              |        ...  -I|
                              |           -I T|

   derived from the standard central difference discretization
   of the 2 dimensional convection-diffusion operator
                      -(Laplacian u) + rho*(du/dx) 
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

#ifndef BCMATRXA_H
#define BCMATRXA_H

#include "arcomp.h"
#include "blas1c.h"

template<class ARFLOAT, class ARINT>
void CompMatrixA(ARINT nx, ARINT& n, ARINT& nL, ARINT& nU, 
                 arcomplex<ARFLOAT>* &A)
{

  ARINT              i, j, lda;
  arcomplex<ARFLOAT> h, h2, dd, dl, du, df;

  // Defining constants.

  const arcomplex<ARFLOAT> half(0.5,0.0);
  const arcomplex<ARFLOAT> one(1.0,0.0);
  const arcomplex<ARFLOAT> four(4.0,0.0);
  const arcomplex<ARFLOAT> rho(1.0e2,0.0);

  h   = one/arcomplex<ARFLOAT>((ARFLOAT)(nx+1),0);
  h2  = h*h;
  df  = -(one/h2);
  dd  = four/h2;
  dl  = df - half*rho/h;
  du  = df + half*rho/h;

  // Defining the dimension of the problem, 
  // the upper and the lower bandwidth.

  n  = nx*nx;
  nL = nx;
  nU = nx;

  // Creating output vector A.

  lda = nL+nU+1;
  A    = new arcomplex<ARFLOAT>[n*lda];

  // Filling A with zeros.

  for (j=0; j<(n*lda); j++ ) A[j] = arcomplex<ARFLOAT>(0.0, 0.0);

  // Creating matrix A.

  for (i=0, j=nU; i<n; i++, j+=lda) {
    if (i>=nx)     A[j-nx] = df;
    if (i%nx)      A[j-1]  = du;
                   A[j]    = dd;
    if ((i+1)%nx)  A[j+1]  = dl;
    if (i<(n-nx))  A[j+nx] = df;
  }
    
} // CompMatrixA.


#endif // BCMATRXA_H
