/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE SMatrixA.h
   Class template for the 2-dimensional discrete Laplacian on unit
   square with zero Dirichlet boundary conditions. A is the nx*nx 
   by nx*nx block tridiagonal matrix

                              | T -I          |
                              |-I  T -I       |
                          A = |   -I  T       |
                              |        ...  -I|
                              |           -I T|

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef SMATRIXA_H
#define SMATRIXA_H

#include "matprod.h"

template<class ART>
class SymMatrixA: public MatrixWithProduct<ART> {

 private:

  int nx;

  void MultTv(ART* x, ART* y);

 public:

  void MultMv(ART* v, ART* w);

  SymMatrixA(int nxval): MatrixWithProduct<ART>(nxval*nxval) { nx = nxval; }

}; // SymMatrixA.


template<class ART>
void SymMatrixA<ART>::MultTv(ART* x, ART* y)
/*
  Computes the matrix vector multiplication y<---T*x
  where T is a nx by nx tridiagonal matrix with DD on the
  diagonal, DL on the subdiagonal, and DU on the superdiagonal.
*/

{

  int  j;
  ART  dd, dl, du;

  const ART one  = 1.0;
  const ART four = 4.0;

  dd  = four;
  dl  = -one;
  du  = -one;

  y[0] =  dd*x[0] + du*x[1];
  for (j = 1; j <= nx-2; j++) {
    y[j] = dl*x[j-1] + dd*x[j] + du*x[j+1];
  }
  y[nx-1] =  dl*x[nx-2] + dd*x[nx-1];

  return;

} // MultTv


template<class ART>
void SymMatrixA<ART>::MultMv(ART* v, ART* w)
/*
  Computes w <- A*v.
*/

{

  int  j, lo, n2;
  ART  h2;

  const ART one = 1.0;

  MultTv(&v[0],&w[0]);
  axpy(nx, -one, &v[nx], 1, &w[0], 1);

  for (j = 2; j<=nx-1; j++) {
    lo = (j-1)*nx;
    MultTv(&v[lo], &w[lo]);
    axpy(nx, -one, &v[lo-nx], 1, &w[lo], 1);
    axpy(nx, -one, &v[lo+nx], 1, &w[lo], 1);
  }

  lo = (nx-1)*nx;
  MultTv(&v[lo], &w[lo]);
  axpy(nx, -one, &v[lo-nx], 1, &w[lo], 1);

  // Scaling the vector w by (1/h^2), where h is the mesh size.

  n2 = nx*nx;
  h2 = ART((nx+1)*(nx+1));
  scal(n2, h2, w, 1L);

  return;

} //  MultMv.


#endif // SMATRIXA_H
