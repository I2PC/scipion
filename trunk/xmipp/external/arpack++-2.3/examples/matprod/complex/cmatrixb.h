/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE CMatrixB.h
   Class template for the tridiagonal matrix derived from 
   the standard central difference of the 1-d convection diffusion 
   operator u" + rho*u' on the interval [0, 1] with zero
   Dirichlet boundary conditions.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef CMATRIXB_H
#define CMATRIXB_H

#include "arcomp.h"
#include "matprod.h"
#include "blas1c.h"
#include "lapackc.h"

template<class T>
class CompMatrixB: public MatrixWithProduct<arcomplex<T> > {

 private:

  arcomplex<T> rho;
  arcomplex<T> shift;
  arcomplex<T> *Ad, *Adl, *Adu, *Adu2;
  int          *ipiv;
  int          decsize;

  void FactorDataDeallocate();

 public:

  void FactorOP();

  void MultMv(arcomplex<T>* v, arcomplex<T>* w);

  void MultOPv(arcomplex<T>* v, arcomplex<T>* w);

  CompMatrixB(int nv, arcomplex<T> rhov);

  CompMatrixB(int nv, arcomplex<T> shiftv, arcomplex<T> rhov);

  virtual ~CompMatrixB();

}; // CompMatrixB.


template<class T>
inline void CompMatrixB<T>::FactorDataDeallocate()
/*
  Eliminates the data structure used on matrix factorization.
*/

{

  delete[] Ad;
  delete[] Adl;
  delete[] Adu;
  delete[] Adu2;
  delete[] ipiv;

} // FactorDataDeallocate.


template<class T>
void CompMatrixB<T>::FactorOP()
/*
  Factors (M-shift*I).
*/

{

  int          j, ierr;
  arcomplex<T> h, h2, s, s1, s2, s3;

  const arcomplex<T> one(1.0, 0.0);
  const arcomplex<T> two(2.0, 0.0);

  if (decsize != this->ncols()) {
    decsize = this->ncols();
    FactorDataDeallocate();
    Ad   = new arcomplex<T>[this->ncols()];
    Adl  = new arcomplex<T>[this->ncols()];
    Adu  = new arcomplex<T>[this->ncols()];
    Adu2 = new arcomplex<T>[this->ncols()];
    ipiv = new int[this->ncols()];
  }

  h  = one/arcomplex<T>((this->ncols()+1),0.0);
  h2 = h*h;
  s  = rho/two;
  s1 = -one/h2 - s/h;
  s2 = two/h2 - shift;
  s3 = -one/h2 + s/h;

  for (j=0; j<this->ncols()-1; j++) {
    Adl[j] = s1;
    Ad[j]  = s2;
    Adu[j] = s3;
  }
  Ad[this->ncols()-1]  = s2;

  gttrf(this->ncols(), Adl, Ad, Adu, Adu2, ipiv, ierr);

} // FactorOP.


template<class T>
void CompMatrixB<T>::MultMv(arcomplex<T>* v, arcomplex<T>* w)
/*
  Computes the matrix-vector multiplication w <- A*v.
*/

{

  int          j;
  arcomplex<T> dd, dl, du, s, h, h2;

  const arcomplex<T> one( 1.0, 0.0);
  const arcomplex<T> two( 2.0, 0.0);

  h  = one/arcomplex<T>((this->ncols()+1),0.0);
  h2 = h*h;
  s  = rho/two;
  dd = two/h2;
  dl = -one/h2 - s/h;
  du = -one/h2 + s/h;

  w[0] = dd*v[0] + du*v[1];
  for (j=1; j<this->ncols()-1; j++) {
    w[j] = dl*v[j-1] + dd*v[j] + du*v[j+1];
  }
  w[this->ncols()-1] = dl*v[this->ncols()-2] + dd*v[this->ncols()-1];

} //  MultMv.


template<class T>
void CompMatrixB<T>::MultOPv(arcomplex<T>* v, arcomplex<T>* w)
/*
  Computes the matrix-vector product w <- inv(M-shift*I)*v.
*/

{

  int  ierr;
  char *type = "N";

  copy(this->ncols(), v, 1, w, 1);
  gttrs(type, this->ncols(), 1, Adl, Ad, Adu, Adu2, ipiv, w, this->ncols(), ierr);

} // MultOPv.


template<class T>
inline CompMatrixB<T>::CompMatrixB(int nval, arcomplex<T> rhov):
  MatrixWithProduct<arcomplex<T> >(nval)
/*
  Constructor
*/

{

  decsize = 0;
  Ad      = 0;
  Adl     = 0;
  Adu     = 0;
  Adu2    = 0;
  ipiv    = 0;
  shift   = 0.0;
  rho     = rhov;

} // Constructor.


template<class T>
inline CompMatrixB<T>::
CompMatrixB(int nv, arcomplex<T> shiftv, arcomplex<T> rhov):
  MatrixWithProduct<arcomplex<T> >(nv)
/*
  Constructor with shift
*/

{

  decsize = 0;
  Ad      = 0;
  Adl     = 0;
  Adu     = 0;
  Adu2    = 0;
  ipiv    = 0;
  shift   = shiftv;
  rho     = rhov;
  FactorOP();

} // Constructor with shift.


template<class T>
inline CompMatrixB<T>::~CompMatrixB()
/*
  Destructor
*/

{

  FactorDataDeallocate();

} // Destructor.


#endif // CMATRIXB_H

