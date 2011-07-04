/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE BNMatrxC.h
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

#ifndef BNMATRXC_H
#define BNMATRXC_H

template<class ARFLOAT, class ARINT>
void MassMatrix(ARINT n, ARINT& nL, ARINT& nU, ARFLOAT* &A)
{

  ARINT   i;
  ARFLOAT diag, sub;

  // Defining constants.

  sub  = 1.0/6.0/ARFLOAT(n+1);
  diag = 4.0/6.0/ARFLOAT(n+1);

  // Defining the lower and upper bandwidth.

  nL = 1;
  nU = 1;

  // Creating output vector A.

  A  = new ARFLOAT[3*n];
  for (i=1; i<(3*n); i+=3) {
    if (i-1)   A[i-1] = sub;
               A[i]   = diag;
    if (n-i-1) A[i+1] = sub;
  }

} // MassMatrix.

#endif // BNMATRXC_H

