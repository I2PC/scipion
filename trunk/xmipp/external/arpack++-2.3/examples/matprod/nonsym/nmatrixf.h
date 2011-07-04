/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE NMatrixF.h
   Class template for the symetric tridiagonal matrix with 4
   on the main diagonal and -1 on the other diagonals.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef NMATRIXF_H
#define NMATRIXF_H

#include "matprod.h"
#include "blas1c.h"
#include "lapackc.h"

template<class T>
class NonSymMatrixF: public MatrixWithProduct<T> {

 public:

  void MultMv(T* v, T* w);
  // Matrix vector multiplication w <- M*v.

  NonSymMatrixF(int nx): MatrixWithProduct<T>(nx) { }
  // Constructor.

}; // NonSymMatrixF.


template<class T>
void NonSymMatrixF<T>::MultMv(T* v, T* w)
{

  int j;

  const T one  = 1.0;
  const T four = 4.0;

  w[0] = four*v[0] + one*v[1];
  for (j=1; j<this->ncols()-1; j++) {
    w[j] = one*v[j-1] + four*v[j] + one*v[j+1];
  }
  w[this->ncols()-1] = one*v[this->ncols()-2] + four*v[this->ncols()-1];

} // MultMv.


#endif // NMATRIXF_H

