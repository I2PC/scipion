/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE CGenPrbB.h
   Very simple template class intended to illustrate how to
   use ARPACK++ to find some few eigenvalues and eigenvectors
   of complex generalized problems in shift and invert mode.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef CGENPRBB_H
#define CGENPRBB_H

#include "arcomp.h"
#include "blas1c.h"
#include "lapackc.h"
#include "cmatrixe.h"
#include "cmatrixf.h"


template <class T>
class ComplexGenProblemB {

 private:

  int          n, decsize;
  int          *ipiv;
  arcomplex<T> rho;
  arcomplex<T> shift;
  arcomplex<T> *Ad, *Adl, *Adu, *Adu2;

  void FactorDataDeallocate();
  // Eliminates the data structure used on matrix factorization.

  void FactorOP();
  // Factors (A-shift*B).

 public:

  ComplexMatrixE<T> A;
  ComplexMatrixF<T> B;

  void MultAv(arcomplex<T> *v, arcomplex<T> *w) { A.MultMv(v,w); }
  // Performs the matrix-vector multiplication w <- A*v.

  void MultBv(arcomplex<T> *v, arcomplex<T> *w) { B.MultMv(v,w); }
  // Performs the matrix-vector multiplication w <- B*v.

  void MultOPv(arcomplex<T>* v, arcomplex<T>* w);
  // Performs the matrix-vector multiplication w <- inv(A-shift*B)*v.

  ComplexGenProblemB(int nx, arcomplex<T> rhop, arcomplex<T> shiftp);
  // Constructor.

  ~ComplexGenProblemB();
  // Destructor

}; // struct ComplexGenProblemA.


template<class T>
inline void ComplexGenProblemB<T>::FactorDataDeallocate()
{

  delete[] Ad;
  delete[] Adl;
  delete[] Adu;
  delete[] Adu2;
  delete[] ipiv;

} // FactorDataDeallocate.


template<class T>
void ComplexGenProblemB<T>::FactorOP()
{

  int          j, ierr;
  arcomplex<T> h, s, s1, s2, s3;

  const arcomplex<T> one(1.0, 0.0);
  const arcomplex<T> two(2.0, 0.0);
  const arcomplex<T> four(4.0, 0.0);

  if (decsize != n) {
    decsize = n;
    FactorDataDeallocate();
    Ad   = new arcomplex<T>[n];
    Adl  = new arcomplex<T>[n];
    Adu  = new arcomplex<T>[n];
    Adu2 = new arcomplex<T>[n];
    ipiv = new int[n];
  }

  h  = one/arcomplex<T>((n+1),0.0);
  s  = rho/two;
  s1 = -one/h - s - shift*h;
  s2 = two/h - four*shift*h;
  s3 = -one/h + s - shift*h;

  for (j=0; j<n-1; j++) {
    Adl[j] = s1;
    Ad[j]  = s2;
    Adu[j] = s3;
  }
  Ad[n-1]  = s2;

  gttrf(n, Adl, Ad, Adu, Adu2, ipiv, ierr);

} // FactorOP.


template<class T>
void ComplexGenProblemB<T>::MultOPv(arcomplex<T>* v, arcomplex<T>* w)
{

  int  ierr;
  char *type = "N";

  copy(n, v, 1, w, 1);
  gttrs(type, n, 1, Adl, Ad, Adu, Adu2, ipiv, w, n, ierr);

} // MultOPv.


template<class T>
ComplexGenProblemB<T>::
ComplexGenProblemB(int nx, arcomplex<T> rhop, arcomplex<T> shiftp): A(nx,rhop), B(nx)
{

  rho     = rhop;
  shift   = shiftp;
  decsize = 0;
  Ad      = 0;
  Adl     = 0;
  Adu     = 0;
  Adu2    = 0;
  ipiv    = 0;
  n       = A.ncols();
  FactorOP();

} // Constructor.


template<class T>
ComplexGenProblemB<T>::~ComplexGenProblemB()
{

  FactorDataDeallocate();

} // Destructor.


#endif // CGENPRBB_H

