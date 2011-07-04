/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE BNMatrxD.h
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

#ifndef BNMATRXD_H
#define BNMATRXD_H

template<class ARFLOAT, class ARINT>
void BandMatrixD(ARINT nx, ARFLOAT rho, ARINT& n, ARINT& nL, ARINT& nU, 
                 ARFLOAT* &A)

{

  ARINT   i, j, lda;
  ARFLOAT dd, dl, du, h;

  // Defining constants.

  h  = 1.0/ARFLOAT(nx+1);
  dd = 4.0;
  dl = -1.0 - 0.5*rho*h;
  du = -1.0 + 0.5*rho*h;

  // Defining the number of columns and the upper and lower bandwidth.

  n   = nx*nx;
  nL  = nx;
  nU  = nx;

  // Creating output vector A.

  lda = nL+nU+1;
  A   = new ARFLOAT[n*lda];

  // Filling A with zeros.

  for (j=0; j<(n*lda); j++ ) A[j] = (ARFLOAT)0.0;

  // Creating matrix A.

  for (i=0, j=nU; i<n; i++, j+=lda) {
    if (i>=nx)     A[j-nx] = -1.0;
    if (i%nx)      A[j-1]  = du;
                   A[j]    = dd;
    if ((i+1)%nx)  A[j+1]  = dl;
    if (i<(n-nx))  A[j+nx] = -1.0;
  }
    
} //  BandMatrixD.

#endif // BNMATRXD_H

