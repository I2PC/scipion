/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE SMatrixC.h
   Class template for the one dimensional discrete Laplacian on
   the interval [0,1] with zero Dirichlet boundary conditions.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef SMATRIXC_H
#define SMATRIXC_H

#include "matprod.h"
#include "blas1c.h"
#include "lapackc.h"

template<class ART>
class SymMatrixC: public MatrixWithProduct<ART> {

 public:

  void MultMv(ART* v, ART* w);

  SymMatrixC(int nv): MatrixWithProduct<ART>(nv) { }

  virtual ~SymMatrixC() { }

}; // SymMatrixC.


template<class ART>
void SymMatrixC<ART>::MultMv(ART* v, ART* w)
/*
  Matrix-vector product.
*/

{

  int  j;
  ART  h2;

  const ART two = 2.0;

  w[0] =  two*v[0] - v[1];
  for (j=1; j<this->ncols()-1; j++) {
    w[j] = - v[j-1] + two*v[j] - v[j+1];
  }
  w[this->ncols()-1] = - v[this->ncols()-2] + two*v[this->ncols()-1];

  // Scaling the vector w by (1 / h).

  h2 = ART(this->ncols()+1);
  scal(this->ncols(), h2, w, 1L);

  return;

} //  MultMv.


#endif // SMATRIXC_H

