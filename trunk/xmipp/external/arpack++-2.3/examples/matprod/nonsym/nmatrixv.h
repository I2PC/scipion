/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE NMatrixV.h
   Class template for the matrix derived from the simplest
   finite difference discretization of the integral operator

          f(s) = integral[K(s,t) x(t) dt].

   A is a discretization of the 2-dimensional kernel K(s,t)dt 
   where

          k(s,t) = s(t-1),  if 0 <= s <= t <= 1.
                   t(s-1),  if 0 <= t <= s <= 1.

   Thus A is a m by n matrix with entries

          A(i,j) = k*(s(i))*(t(j)-1), if i <= j;
                   k*(t(j))*(s(i)-1), if i >  j,

   where s(i) = i/(m+1), t(j) = j/(n+1) and k = 1/(n+1).       

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef NMATRIXV_H
#define NMATRIXV_H

#include "matprod.h"
#include "blas1c.h"
#include "lapackc.h"

template<class T>
class NonSymMatrixV: public MatrixWithProduct<T> {

 public:

  void MultMv(T* v, T* w);

  void MultMtv(T* v, T* w);

  void MultOPv(T* v, T* w);

  NonSymMatrixV(int mv, int nv): MatrixWithProduct<T>(mv, nv) { }

  virtual ~NonSymMatrixV() { }

}; // NonSymMatrixV.


template<class T>
void NonSymMatrixV<T>::MultMv(T* v, T* w)
/*
  Performs the matrix-vector multiplication w <- A*v.
*/

{

  int  i,j;
  T    h,k,s,t,z;

  h  = 1.0/T(this->nrows()+1);
  k  = 1.0/T(this->ncols()+1);

  for (i=0; i!=this->nrows(); i++) w[i] = 0.0;

  t = 0.0;
  for (j=0; j!=this->ncols(); j++) {
    t += k;
    s  = 0.0;
    z  = v[j];
    for (i=0; i<=j; i++) {
      s    += h;
      w[i] += k*s*(t-1.0)*z;
    }
    for (i=j+1; i!=this->nrows(); i++) {
      s    += h;
      w[i] += k*t*(s-1.0)*z;
    }
  }

  return;

} //  MultMv.


template<class T>
void NonSymMatrixV<T>::MultMtv(T* v, T* w)
/*
  Performs the matrix-vector multiplication w <- A'*v.
*/

{

  int  i,j;
  T    h,k,s,t,z;

  h  = 1.0/T(this->nrows()+1);
  k  = 1.0/T(this->ncols()+1);

  t = 0.0;
  for (j=0; j!=this->ncols(); j++) {
    t += k;
    z  = 0.0;
    s  = 0.0;
    for (i=0; i<=j; i++) {
      s += h;
      z += k*s*(t-1.0)*v[i];
    }
    for (i=j+1; i!=this->nrows(); i++) {
      s += h;
      z += k*t*(s-1.0)*v[i];
    }
    w[j] = z;
  }

  return;

} //  MultMtv.


template<class T>
void NonSymMatrixV<T>::MultOPv(T* v, T* w)
/*
  Performs the matrix-vector product w <- (A'*A)*v.
*/

{

  T* t = new T[this->nrows()];

  MultMv(v,t);
  MultMtv(t,w);

  delete[] t;

} // MultOPv


#endif // NMATRIXV_H
