/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE CMatrixF.h
   Class template for the symmetric tridiagonal matrix with 4.0 on
   the main diagonal, and 1.0 on the off-diagonals.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef CMATRIXF_H
#define CMATRIXF_H

#include "arcomp.h"
#include "matprod.h"
#include "blas1c.h"
#include "lapackc.h"

template<class T>
class ComplexMatrixF: public MatrixWithProduct<arcomplex<T> > {

 public:

  void MultMv(arcomplex<T> *v, arcomplex<T> *w);
  // Matrix vector multiplication w <- M*v.

  ComplexMatrixF(int nx): MatrixWithProduct<arcomplex<T> >(nx) { }
  // Constructor

}; // ComplexMatrixF.


template<class T>
void ComplexMatrixF<T>::MultMv(arcomplex<T> *v, arcomplex<T> *w)
{

  int          j;
  arcomplex<T> h;

  const arcomplex<T> one(1.0, 0.0);
  const arcomplex<T> four(4.0, 0.0);

  w[0] = four*v[0] + one*v[1];
  for (j=1; j<this->ncols()-1; j++) {
    w[j] = one*v[j-1] + four*v[j] + one*v[j+1];
  }
  w[this->ncols()-1] = one*v[this->ncols()-2] + four*v[this->ncols()-1];

  h = one/arcomplex<T>((this->ncols()+1),0.0);
  scal(this->ncols(), h, w, 1);

} // MultMv.


#endif // CMATRIXF_H

