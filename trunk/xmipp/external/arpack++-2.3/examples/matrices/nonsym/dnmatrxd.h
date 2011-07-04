/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE DNMatrxD.h
   Function template for the block tridiagonal matrix A, whose
   diagonal blocks are tridiagonal with 4 on the diagonal, 
   -1-rho*h/2 on the subdiagonal and -1+rho*h/2 on the superdiagonal.  
   Each off-diagonal block of A is an identity matrix.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef DNMATRXD_H
#define DNMATRXD_H

template<class ARFLOAT, class ARINT>
void DenseMatrixD(ARINT nx, ARFLOAT rho, ARINT& n, ARFLOAT* &A)

{

  ARINT   i, j;
  ARFLOAT dd, dl, du, h;

  // Defining constants.

  h  = 1.0/ARFLOAT(nx+1);
  dd = 4.0;
  dl = -1.0 - 0.5*rho*h;
  du = -1.0 + 0.5*rho*h;

  // Defining the dimension of A.

  n   = nx*nx;

  // Creating output vector A.

  A  = new ARFLOAT[n*n];

  for (i=0; i<n*n; i++) A[i]=0.0;

  for (i=0, j=0; i<n; i++, j+=n+1) {
    if (i>=nx)     A[j-nx] = -1.0;
    if (i%nx)      A[j-1]  = du;
                   A[j]    = dd;
    if ((i+1)%nx)  A[j+1]  = dl;
    if (i<(n-nx))  A[j+nx] = -1.0;
  }
    
} //  DenseMatrixD.

#endif // DNMATRXD_H

