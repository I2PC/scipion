/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE LNMatrxV.h
   Function template that generates a  sparse (2n x n) rectangular matrix.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef LNMATRXV_H
#define LNMATRXV_H


template<class ARFLOAT, class ARINT>
void RectangularMatrix(ARINT n, ARINT& m, ARINT& nnz, ARFLOAT* &A, 
                       ARINT* &irow, ARINT* &pcol)
{

  ARINT   i, j;
  ARFLOAT dd, dl, du;

  // Defining constants.

  dl = 1.0;
  dd = 4.0;
  du = 2.0;

  // Defining the number of nonzero elements and lines of the matrix.

  nnz =  n*6-2;
  m   =  n*2;

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
      A[j++]  = du;
    }

    irow[j] = i;
    A[j++]  = dd;

    irow[j] = i+1;
    A[j++]  = dl;

    irow[j] = i+n-1;
    A[j++]  = dl;

    irow[j] = i+n;
    A[j++]  = dd;

    if (i != (n-1)) {
      irow[j] = i+n+1;
      A[j++]  = du;
    }

    pcol[i+1] = j;

  }

} // Rectangular matrix.

#endif // LNMATRXV_H

