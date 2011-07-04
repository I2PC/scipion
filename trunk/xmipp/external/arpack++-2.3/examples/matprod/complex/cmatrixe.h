/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE CMatrixE.h
   Complex matrix template used to illustrate how to solve
   a generalized eigenvalue problem in shift and invert mode.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef CMATRIXE_H
#define CMATRIXE_H

#include "arcomp.h"
#include "matprod.h"

template<class T>
class ComplexMatrixE: public MatrixWithProduct<arcomplex<T> > {

 private:

  arcomplex<T> rho; 

 public:

  void MultMv(arcomplex<T>* v, arcomplex<T>* w);
  // Matrix vector multiplication w <- M*v.

  ComplexMatrixE(int nx, arcomplex<T> rhop);
  // Constructor

}; // ComplexMatrixE.


template<class T>
void ComplexMatrixE<T>::MultMv(arcomplex<T>* v, arcomplex<T>* w)
{

  int          j;
  arcomplex<T> dd, dl, du, s, h;

  const arcomplex<T> one( 1.0, 0.0);
  const arcomplex<T> two( 2.0, 0.0);

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


template<class T>
ComplexMatrixE<T>::
ComplexMatrixE(int nx, arcomplex<T> rhop): MatrixWithProduct<arcomplex<T> >(nx) 
{ 

  rho = rhop;

} // Constructor

#endif // CMATRIXE_H


