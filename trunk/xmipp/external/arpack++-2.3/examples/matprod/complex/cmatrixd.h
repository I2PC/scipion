/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE CMatrixD.h
   Class template for the mass matrix formed by using
   piecewise linear elements on [0,1].

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef CMATRIXD_H
#define CMATRIXD_H

#include "arcomp.h"
#include "matprod.h"
#include "blas1c.h"
#include "lapackc.h"

template<class T>
class ComplexMatrixD: public MatrixWithProduct<arcomplex<T> > {

 private:

  arcomplex<T> *Ad, *Adl, *Adu, *Adu2;
  int          *ipiv;
  int          decsize;

  void FactorDataDeallocate();
  // Eliminates the data structure used on matrix factorization.

 public:

  void FactorM();
  // Factors B, a symmetric tridiagonal matrix with 4 on the
  // main diagonal and 1 on the off-diagonals.

  void MultMv(arcomplex<T> *v, arcomplex<T> *w);
  // Performs the matrix-vector multiplication w <- M*v

  void SolveM(arcomplex<T> *v);
  // Solves M.w = v. v is overwritten with vector w.

  ComplexMatrixD(int nx);
  // Constructor

  virtual ~ComplexMatrixD();
  // Destructor

}; // ComplexMatrixD.


template<class T>
inline void ComplexMatrixD<T>::FactorDataDeallocate()
{

  if (Ad != 0) {
    delete[] Ad;
    delete[] Adl;
    delete[] Adu;
    delete[] Adu2;
    delete[] ipiv;
  }

} // FactorDataDeallocate


template<class T>
void ComplexMatrixD<T>::FactorM()
{

  int          j, ierr;
  arcomplex<T> h;

  const arcomplex<T> one(1.0, 0.0);
  const arcomplex<T> four(4.0, 0.0);

  if (decsize != this->ncols()) {
    decsize = this->ncols();
    FactorDataDeallocate();
    Ad   = new arcomplex<T>[decsize];
    Adl  = new arcomplex<T>[decsize];
    Adu  = new arcomplex<T>[decsize];
    Adu2 = new arcomplex<T>[decsize];
    ipiv = new int[decsize];
  }

  h  = one/arcomplex<T>((this->ncols()+1),0.0);

  for (j=0; j<this->ncols()-1; j++) {
    Adl[j] = one*h;
    Ad[j]  = four*h;
    Adu[j] = one*h;
  }
  Ad[this->ncols()-1]  = four*h;

  gttrf(this->ncols(), Adl, Ad, Adu, Adu2, ipiv, ierr);

} // FactorM.


template<class T>
void ComplexMatrixD<T>::MultMv(arcomplex<T> *v, arcomplex<T> *w)
{

  int          j;
  arcomplex<T> h;

  const arcomplex<T> one(1.0, 0.0);
  const arcomplex<T> four(4.0, 0.0);

  w[0] = four*v[0] + one*v[1];
  for (j=1; j<this->ncols()-1; j++) {
    w[j] = one*v[j-1] + four*v[j] + one*v[j+1];
  }
  w[this->ncols()-1] = one*v[this->ncols()-2] + four*v[this->ncols()-1];

  h = one/arcomplex<T>((this->ncols()+1),0.0);
  scal(this->ncols(), h, w, 1);

} // MultMv.


template<class T>
inline void ComplexMatrixD<T>::SolveM(arcomplex<T> *v)
{

  int  ierr;
  char *type = "N";

  gttrs(type, this->ncols(), 1, Adl, Ad, Adu, Adu2, ipiv, v, this->ncols(), ierr);

} // SolveM.


template<class T>
inline ComplexMatrixD<T>::
ComplexMatrixD(int nx): MatrixWithProduct<arcomplex<T> >(nx)
{

  decsize = 0;
  Ad      = 0;
  Adl     = 0;
  Adu     = 0;
  Adu2    = 0;
  ipiv    = 0;

} // Constructor.


template<class T>
inline ComplexMatrixD<T>::~ComplexMatrixD()
{

  FactorDataDeallocate();

} // Destructor.


#endif // CMATRIXD_H

