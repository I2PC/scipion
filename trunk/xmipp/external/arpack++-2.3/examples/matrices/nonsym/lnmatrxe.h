/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE LNMatrxE.h
   Function template that generates a nonsymmetric 
   tridiagonal matrix with 2 on the main diagonal,
   3 on the superdiagonal and -2 on the subdiagonal.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef LNMATRXE_H
#define LNMATRXE_H


template<class ARFLOAT, class ARINT>
void NonSymMatrixE(ARINT n, ARINT& nnz, ARFLOAT* &A, ARINT* &irow, 
                   ARINT* &pcol)
{

  ARINT   i, j;

  // Defining constants.
  
  const ARFLOAT three = 3.0;
  const ARFLOAT two   = 2.0;

  // Defining the number of nonzero matrix elements.

  nnz = 3*n-2;

  // Creating output vectors.

  A    = new ARFLOAT[nnz];
  irow = new ARINT[nnz];
  pcol = new ARINT[n+1];

  // Filling A, irow and pcol.

  pcol[0] = 0;
  j = 0;
  for (i=0; i!=n; i++) {
    if (i != 0) {
      irow[j] = i-1;
      A[j++]  = three;
    }
    irow[j] = i;
    A[j++]  = two;
    if (i != (n-1)) {
      irow[j] = i+1;
      A[j++]  = -two;
    }
    pcol[i+1] = j;
  }

} // NonSymMatrixE.


#endif // LNMATRXE_H



