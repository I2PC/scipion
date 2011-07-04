/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE DNMatrxE.h
   Function template for the tridiagonal matrix with 4 on the main
   diagonal and 1 on the subdiagonal and superdiagonal.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef DNMATRXE_H
#define DNMATRXE_H

template<class ARFLOAT, class ARINT>
void DenseMatrixE(ARINT n, ARFLOAT* &A)

{

  ARINT i, j;

  // Creating output vector A.

  A  = new ARFLOAT[n*n];

  for (i=0; i<n*n; i++) A[i]=0.0;

  for (i=0, j=0; i<n; i++, j+=n+1) {
    if (i)     A[j-1] = 1.0;
               A[j]   = 4.0;
    if (n-i-1) A[j+1] = 1.0;
  }
   
} //  DenseMatrixE.

#endif // DNMATRXE_H
