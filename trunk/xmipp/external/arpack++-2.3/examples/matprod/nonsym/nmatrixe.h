/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE NMatrixE.h
   Class template for the tridiagonal matrix with 2 on the main
   diagonal, -2 on the subdiagonal and 3 on the superdiagonal.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef NMATRIXE_H
#define NMATRIXE_H

#include "matprod.h"

template<class T>
class NonSymMatrixE: public MatrixWithProduct<T> {

 public:

  void MultMv(T* v, T* w);
  // Matrix vector multiplication w <- A*v.

  NonSymMatrixE(int nx): MatrixWithProduct<T>(nx) { }
  // Constructor.

}; // NonSymMatrixE.

template<class T>
void NonSymMatrixE<T>::MultMv(T* v, T* w)
{

  int j;

  const T three = 3.0;
  const T two   = 2.0;

  w[0] = two*v[0] + three*v[1];
  for (j=1; j<this->ncols()-1; j++) {
    w[j] = -two*v[j-1] + two*v[j] + three*v[j+1];
  }
  w[this->ncols()-1] = -two*v[this->ncols()-2] + two*v[this->ncols()-1];

} //  MultMv.


#endif // NMATRIXE_H



