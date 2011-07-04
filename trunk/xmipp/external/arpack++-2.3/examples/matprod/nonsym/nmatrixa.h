/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE NMatrixA.h
   Class template for the matrix

                          | T -I          |
                          |-I  T -I       |
                      A = |   -I  T       |
                          |        ...  -I|
                          |           -I T|

   derived from the standard central difference discretization of the
   2-dimensional convection-diffusion operator (Laplacian u) + rho*(du/dx)
   on a unit square with zero boundary condition.
   When rho*h/2 <= 1, the discrete convection-diffusion operator has real
   eigenvalues.  When rho*h/2 > 1, it has COMPLEX eigenvalues.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef NMATRIXA_H
#define NMATRIXA_H

#include "matprod.h"

template<class T>
class NonSymMatrixA: public MatrixWithProduct<T> {

 private:

  int nx;

  void MultTv(T* x, T* y);

 public:

  void MultMv(T* v, T* w);

  NonSymMatrixA(int nxval): MatrixWithProduct<T>(nxval*nxval) { nx = nxval; }

}; // NonSymMatrixA.


template<class T>
void NonSymMatrixA<T>::MultTv(T* x, T* y)
/*
  Computes the matrix vector multiplication y <- T*x
  where T is a nx by nx tridiagonal matrix with DD on the
  diagonal, DL on the subdiagonal, and DU on the superdiagonal.
*/

{

  int  j;
  T    h, h2, dd, dl, du;

  const T one = 1.0;
  const T rho = 0.0;

  h   = one/T(nx+1);
  h2  = h*h;
  dd  = 4.0/h2;
  dl  = -one/h2 - 5.0e-1*rho/h;
  du  = -one/h2 + 5.0e-1*rho/h;

  y[0] = dd*x[0] + du*x[1];
  for (j = 1; j<nx-1; j++) {
    y[j] = dl*x[j-1] + dd*x[j] + du*x[j+1];
  }
  y[nx-1] =  dl*x[nx-2] + dd*x[nx-1];

  return;

} // MultTv


template<class T>
void NonSymMatrixA<T>::MultMv(T* v, T* w)
/*
  Computes w <- M*v.
*/

{

  int  j, lo;
  T    h2;

  const T one = 1.0;

  h2 = one/T((nx+1)*(nx+1));

  MultTv(v,w);
  axpy(nx, -one/h2, &v[nx], 1, w, 1);

  for (j = 2; j<=nx-1; j++) {
    lo = (j-1)*nx;
    MultTv(&v[lo], &w[lo]);
    axpy(nx, -one/h2, &v[lo-nx], 1, &w[lo], 1);
    axpy(nx, -one/h2, &v[lo+nx], 1, &w[lo], 1);
  }

  lo = (nx-1)*nx;
  MultTv(&v[lo], &w[lo]);
  axpy(nx, -one/h2, &v[lo-nx], 1, &w[lo], 1);

  return;

} //  MultMv


#endif // NMATRIXA_H








