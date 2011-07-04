/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE DNMatrxA.h
   Function template for the matrix

                          | T -I          |
                          |-I  T -I       |
                      A = |   -I  T       |
                          |        ...  -I|
                          |           -I T|

   derived from the standard central difference discretization of the
   2-dimensional convection-diffusion operator (Laplacian u) + rho*(du/dx)
   on a unit square with zero Dirichlet boundary conditions.
   When rho*h/2 <= 1, the discrete convection-diffusion operator has real
   eigenvalues.  When rho*h/2 > 1, it has COMPLEX eigenvalues.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef DNMATRXA_H
#define DNMATRXA_H

#include <math.h>

template<class ARFLOAT, class ARINT>
void DenseMatrixA(ARINT nx, ARFLOAT rho, ARINT& n, ARFLOAT* &A)

{

  // Defining internal variables.

  ARINT    i, j;
  ARFLOAT  h, h2, df;
  ARFLOAT  dd, dl, du;

  // Defining constants.

  h   = 1.0/ARFLOAT(nx+1);
  h2  = h*h;
  dd  = 4.0/h2;
  df  = -1.0/h2;
  dl  = df - 5.0e-1*rho/h;
  du  = df + 5.0e-1*rho/h;

  // Defining the dimension of A.

  n   = nx*nx;

  // Creating matrix A.

  A   = new ARFLOAT[n*n];

  for (i=0; i<n*n; i++) A[i]=0.0;

  for (i=0, j=0; i<n; i++, j+=n+1) {

      if (i>=nx)    A[j-nx] = df;
      if (i%nx)     A[j-1]  = du;
                    A[j]    = dd;
      if ((i+1)%nx) A[j+1]  = dl;
      if (i<(n-nx)) A[j+nx] = df;
    
    }
   
} // DenseMatrixA.

#endif // DNMATRXA_H

