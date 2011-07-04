/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE CMatrixA.h
   Class template for the nx*nx by nx*nx block tridiagonal matrix

                              | T -I          |
                              |-I  T -I       |
                         OP = |   -I  T       |
                              |        ...  -I|
                              |           -I T|

   derived from the standard central difference discretization
   of the 2 dimensional convection-diffusion operator
                       (Laplacian u) + rho*(du/dx) 
   on a unit square with zero boundary condition.
   T is a nx by nx tridiagonal matrix with DD on the diagonal,
   DL on the subdiagonal, and DU on the superdiagonal.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef CMATRIXA_H
#define CMATRIXA_H

#include "arcomp.h"
#include "blas1c.h"
#include "arcomp.h"
#include "matprod.h"

template<class T>
class CompMatrixA: public MatrixWithProduct<arcomplex<T> > {

 private:

  int nx;

  void MultTv(arcomplex<T>* x, arcomplex<T>* y);

 public:

  void MultMv(arcomplex<T>* v, arcomplex<T>* w);

  CompMatrixA(int nxval);

}; // CompMatrixA.


template<class T>
void CompMatrixA<T>::MultTv(arcomplex<T>* x, arcomplex<T>* y)
/*
  Computes the matrix vector multiplication y <- T*x
  where T is a nx by nx tridiagonal matrix with DD on the
  diagonal, DL on the subdiagonal, and DU on the superdiagonal.
*/

{

  int          j;
  arcomplex<T> h, h2, dd, dl, du;

  const arcomplex<T> half(0.5,0.0);
  const arcomplex<T> one(1.0,0.0);
  const arcomplex<T> four(4.0,0.0);
  const arcomplex<T> rho(1.0e2,0.0);

  h   = one/arcomplex<T>(nx+1,0);
  h2  = h*h;
  dd  = four/h2;
  dl  = -one/h2 - half*rho/h;
  du  = -one/h2 + half*rho/h;

  y[0] = dd*x[0] + du*x[1];
  for (j = 1; j<nx-1; j++) {
    y[j] = dl*x[j-1] + dd*x[j] + du*x[j+1];
  }
  y[nx-1] =  dl*x[nx-2] + dd*x[nx-1];

  return;

} // MultTv


template<class T>
void CompMatrixA<T>::MultMv(arcomplex<T>* v, arcomplex<T>* w)
/*
  Matrix-vector subroutine. Computes w <- M*v.
*/

{

  int          j, lo;
  arcomplex<T> h2;

  const arcomplex<T> one(1.0,0.0);

  h2 = one/arcomplex<T>((nx+1)*(nx+1),0.0);

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


template<class T>
CompMatrixA<T>::CompMatrixA(int nxval):
  MatrixWithProduct<arcomplex<T> >(nxval*nxval)
/*
  Constructor.
*/

{

  nx = nxval;

} // constructor.


#endif // CMATRIXA_H
