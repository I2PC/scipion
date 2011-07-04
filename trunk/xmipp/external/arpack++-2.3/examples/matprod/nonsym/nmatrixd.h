/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE NMatrixD.h
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

#ifndef NMATRIXD_H
#define NMATRIXD_H

#include "matprod.h"
#include "blas1c.h"
#include "lapackc.h"

template<class T>
class NonSymMatrixD: public MatrixWithProduct<T> {

 private:

  T    *Md, *Me;
  int  decsize;

  void FactorDataDeallocate();
  // Eliminates of the data structure used on matrix factorization.

 public:

  void FactorM();
  // Factors M. 

  void MultMv(T* v, T* w);
  // Computes w <- M*v.

  void SolveM(T* v);
  // Solves M.w = v. v is overwritten with vector w.

  NonSymMatrixD(int nx);
  // Constructor.

  virtual ~NonSymMatrixD();
  // Destructor.

}; // NonSymMatrixD.

template<class T>
inline void NonSymMatrixD<T>::FactorDataDeallocate()
{

  delete[] Md;
  delete[] Me;

} // FactorDataDeallocate


template<class T>
void NonSymMatrixD<T>::FactorM()
/*
  Factors M, the mass matrix formed by using
  piecewise linear elements on [0, 1].
*/

{

  int  i, ierr;
  T    h;

  const T one  = 1.0;
  const T four = 4.0;

  if (decsize != this->ncols()) {
    decsize = this->ncols();
    FactorDataDeallocate();
    Md   = new T[this->ncols()];
    Me   = new T[this->ncols()];
  }

  h  = one/T(this->ncols()+1);

  for (i=0; i<this->ncols()-1; i++) {
    Md[i] = four*h;
    Me[i] = h;
  }
  Md[this->ncols()-1] = four*h;

  pttrf(this->ncols(), Md, Me, ierr);

} // FactorM.


template<class T>
void NonSymMatrixD<T>::MultMv(T* v, T* w)
{

  int  j;
  T    h;

  const T one  = 1.0;
  const T four = 4.0;

  w[0] = four*v[0] + one*v[1];
  for (j=1; j<this->ncols()-1; j++) {
    w[j] = one*v[j-1] + four*v[j] + one*v[j+1];
  }
  w[this->ncols()-1] = one*v[this->ncols()-2] + four*v[this->ncols()-1];

  h = one/T(this->ncols()+1);
  scal(this->ncols(), h, w, 1);

} // MultMv.


template<class T>
inline void NonSymMatrixD<T>::SolveM(T* v)
{

  int  ierr;

  pttrs(this->ncols(), 1, Md, Me, v, this->ncols(), ierr);

} // SolveM.


template<class T>
inline NonSymMatrixD<T>::NonSymMatrixD(int nx): MatrixWithProduct<T>(nx)
{

  decsize = 0;
  Md      = 0;
  Me      = 0;

} // Constructor.


template<class T>
inline NonSymMatrixD<T>::~NonSymMatrixD()
{

  FactorDataDeallocate();

} // Destructor.


#endif // NMATRIXD_H

