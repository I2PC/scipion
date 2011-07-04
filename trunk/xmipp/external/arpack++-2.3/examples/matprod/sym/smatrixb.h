/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE SMatrixB.h
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

#ifndef SMATRIXB_H
#define SMATRIXB_H

#include "matprod.h"
#include "blas1c.h"
#include "lapackc.h"

template<class ART>
class SymMatrixB: public MatrixWithProduct<ART> {

 private:

  ART  shift;
  ART  *Ad, *Adl, *Adu, *Adu2;
  int  *ipiv;
  int  decsize;

  void FactorDataDeallocate();

 public:

  void FactorOP();

  void MultMv(ART* v, ART* w);

  void MultOPv(ART* v, ART* w);

  SymMatrixB(int nv);

  SymMatrixB(int nv, ART shiftv);

  virtual ~SymMatrixB();

}; // SymMatrixB.


template<class ART>
inline void SymMatrixB<ART>::FactorDataDeallocate()
// Eliminates the data structure used on matrix factorization.

{

  delete[] Ad;
  delete[] Adl;
  delete[] Adu;
  delete[] Adu2;
  delete[] ipiv;

} // FactorDataDeallocate.


template<class ART>
void SymMatrixB<ART>::FactorOP()
/*
  Factors (M-shift*I).
*/

{

  int  i, ierr;
  ART  h2;

  const ART one = 1.0;
  const ART two = 2.0;

  if (decsize != this->ncols()) {
    decsize = this->ncols();
    FactorDataDeallocate();
    Ad   = new ART[this->ncols()];
    Adl  = new ART[this->ncols()];
    Adu  = new ART[this->ncols()];
    Adu2 = new ART[this->ncols()];
    ipiv = new int[this->ncols()];
  }

  h2 = ART((this->ncols()+1)*(this->ncols()+1));

  for (i=0; i<this->ncols(); i++) {
    Ad[i]  = two*h2 - shift;
    Adl[i] = -one*h2;
  }

  copy(this->ncols(), Adl, 1, Adu, 1);
  gttrf(this->ncols(), Adl, Ad, Adu, Adu2, ipiv, ierr);

} // FactorOP.


template<class ART>
void SymMatrixB<ART>::MultMv(ART* v, ART* w)
/*
  Matrix-vector multiplication w <- M*v.
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

  // Scaling the vector w by (1 / h^2).

  h2 = ART((this->ncols()+1)*(this->ncols()+1));
  scal(this->ncols(), h2, w, 1L);

  return;

} //  MultMv.


template<class ART>
void SymMatrixB<ART>::MultOPv(ART* v, ART* w)
/*
  Matrix-vector multiplication w <- inv(M-shift*I)*v.
*/

{

  int  ierr;
  char *type = "N";

  copy(this->ncols(), v, 1, w, 1);
  gttrs(type, this->ncols(), 1, Adl, Ad, Adu, Adu2, ipiv, w, this->ncols(), ierr);

} // MultOPv


template<class ART>
inline SymMatrixB<ART>::SymMatrixB(int nval): MatrixWithProduct<ART>(nval)
// Constructor

{

  decsize = 0;
  Ad      = 0;
  Adl     = 0;
  Adu     = 0;
  Adu2    = 0;
  ipiv    = 0;
  shift   = 0.0;

} // Constructor.


template<class ART>
inline SymMatrixB<ART>::
SymMatrixB(int nv, ART shiftv): MatrixWithProduct<ART>(nv)
// Constructor with shift.

{

  decsize = 0;
  Ad      = 0;
  Adl     = 0;
  Adu     = 0;
  Adu2    = 0;
  ipiv    = 0;
  shift   = shiftv;
  FactorOP();

} // Constructor with shift.


template<class ART>
inline SymMatrixB<ART>::~SymMatrixB()
// Destructor

{

  FactorDataDeallocate();

} // Destructor.


#endif // SMATRIXB_H

