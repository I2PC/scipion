/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE NMatrixC.h
   Class template for the stiffness matrix obtained from the 
   finite element discretization of the 1-dimensional convection-
   diffusion operator
                        d^2u/dx^2 + rho*(du/dx)
   on the interval [0,1], with zero Dirichlet boundary conditions,
   using linear elements.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef NMATRIXC_H
#define NMATRIXC_H

#include "matprod.h"

template<class T>
class NonSymMatrixC: public MatrixWithProduct<T> {

 public:

  void MultMv(T* v, T* w);
  // Computes the matrix vector multiplication w <- A*v.

  NonSymMatrixC(int nx): MatrixWithProduct<T>(nx) { }
  // Constructor

}; // NonSymMatrixC.


template<class T>
void NonSymMatrixC<T>::MultMv(T* v, T* w)
{

  int  j;
  T    dd, dl, du, s, h;

  const T one = 1.0;
  const T two = 2.0;
  const T rho = 10.0;

  h  = one/T(this->ncols()+1);
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


#endif // NMATRIXC_H


