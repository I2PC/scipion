/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE CMatrixC.h
   Class template for the stiffness matrix formed by using
   piecewise linear elements on [0,1].

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef CMATRIXC_H
#define CMATRIXC_H

#include "arcomp.h"
#include "matprod.h"

template<class T>
class ComplexMatrixC: public MatrixWithProduct<arcomplex<T> > {

 public:

  void MultMv(arcomplex<T>* v, arcomplex<T>* w);
  // Matrix vector multiplication w <- A*v.

  ComplexMatrixC(int nx): MatrixWithProduct<arcomplex<T> >(nx) { }
  // Constructor

}; // ComplexMatrixC.


template<class T>
void ComplexMatrixC<T>::MultMv(arcomplex<T>* v, arcomplex<T>* w)
{

  int          j;
  arcomplex<T> dd, dl, du, s, h;

  const arcomplex<T> one( 1.0, 0.0);
  const arcomplex<T> two( 2.0, 0.0);
  const arcomplex<T> rho(10.0, 0.0);

  h  = one/arcomplex<T>((this->ncols()+1),0.0);
  s  = rho/two;
  dd = two/h;
  dl = -one/h - s;
  du = -one/h + s;

  w[0] = dd*v[0] + du*v[1];
  for (j=1; j<this->ncols()-1; j++) {
    w[j] = dl*v[j-1] + dd*v[j] + du*v[j+1];
  }
  w[this->ncols()-1] = dl*v[this->ncols()-2] + dd*v[this->ncols()-1];

} //  MultMv.


#endif // CMATRIXC_H




