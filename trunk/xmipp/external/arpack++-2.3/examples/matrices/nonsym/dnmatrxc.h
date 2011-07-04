/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE DNMatrxC.h
   Function template for the mass matrix formed by using piecewise 
   linear elements on [0, 1].

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef DNMATRXC_H
#define DNMATRXC_H

template<class ARFLOAT, class ARINT>
void MassMatrix(ARINT n, ARFLOAT* &A)
{

  ARINT   i, j;
  ARFLOAT diag, sub;

  // Defining constants.

  sub  = 1.0/6.0/ARFLOAT(n+1);
  diag = 4.0/6.0/ARFLOAT(n+1);

  // Creating output vector A.

  A  = new ARFLOAT[n*n];

  for (i=0; i<n*n; i++) A[i]=0.0;

  for (i=0, j=0; i<n; i++, j+=n+1) {
    if (i)     A[j-1] = sub;
               A[j]   = diag;
    if (n-i-1) A[j+1] = sub;
  }

} // MassMatrix.

#endif // DNMATRXC_H
