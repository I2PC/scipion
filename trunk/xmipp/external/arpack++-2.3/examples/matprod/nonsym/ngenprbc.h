/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE NGenPrbC.h
   Very simple template class intended to illustrate how to use
   ARPACK++ to find some few eigenvalues and eigenvectors of
   nonsymmetric generalized problems in complex shift and invert mode.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef NGENPRBC_H
#define NGENPRBC_H

#include "arcomp.h"
#include "blas1c.h"
#include "lapackc.h"
#include "nmatrixe.h"
#include "nmatrixf.h"

template<class T>
class NonSymGenProblemC {

 private:

  int          n, decsize;
  int          *ipiv;
  T            shiftR, shiftI;
  arcomplex<T> *Ad, *Adl, *Adu, *Adu2, *ctemp;

  void FactorDataDeallocate();
  // Eliminates the data structure used on matrix factorization.

  void FactorOP();
  // Factors (A-shift*B).
 
 public:

  NonSymMatrixE<T> A;
  NonSymMatrixF<T> B;

  void MultAv(T* v, T* w) { A.MultMv(v,w); }
  // Matrix-vector multiplication w <- A*v.

  void MultBv(T* v, T* w) { B.MultMv(v,w); }
  // Matrix-vector multiplication w <- B*v.

  void MultOPvRe(T* v, T* w);
  // Matrix-vector multiplication w <- Real{inv(M-shift*B)}*v.

  void MultOPvIm(T* v, T* w);
  // Matrix-vector multiplication w <- Imag{inv(M-shift*B)}*v.

  NonSymGenProblemC(int nx, T shiftRp, T shiftIp = 0.0);
  // Constructor.

  ~NonSymGenProblemC();
  // Destructor.

}; // NonSymGenProblemC


template<class T>
inline void NonSymGenProblemC<T>::FactorDataDeallocate()
{

  if (decsize) {
    delete[] Ad;
    delete[] Adl;
    delete[] Adu;
    delete[] Adu2;
    delete[] ctemp;
    delete[] ipiv;
  }

} // FactorDataDeallocate.


template<class T>
void NonSymGenProblemC<T>::FactorOP()
{

  int          j, ierr;

  arcomplex<T> c1(-2.0-shiftR, -shiftI);
  arcomplex<T> c2( 2.0-4.0*shiftR, -4.0*shiftI);
  arcomplex<T> c3( 3.0-shiftR, -shiftI);

  if (decsize != n) {
    decsize = n;
    FactorDataDeallocate();
    Ad    = new arcomplex<T>[n];
    Adl   = new arcomplex<T>[n];
    Adu   = new arcomplex<T>[n];
    Adu2  = new arcomplex<T>[n];
    ctemp = new arcomplex<T>[n];
    ipiv  = new int[n];
  }

  for (j=0; j<n-1; j++) {
    Adl[j] = c1;
    Ad[j]  = c2;
    Adu[j] = c3;
  }
  Ad[n-1]  = c2;

  gttrf(n, Adl, Ad, Adu, Adu2, ipiv, ierr);

} // FactorOP.


template<class T>
void NonSymGenProblemC<T>::MultOPvRe(T* v, T* w)
{

  int   i, ierr;
  char* type = "N";

  for (i=0; i<n; i++) {
    ctemp[i] = arcomplex<T>(v[i],0.0);
  }
  gttrs(type, n, 1, Adl, Ad, Adu, Adu2, ipiv, ctemp, n, ierr);
  for (i=0; i<n; i++) {
    w[i]=real(ctemp[i]);
  }

} // MultOPvRe


template<class T>
void NonSymGenProblemC<T>::MultOPvIm(T* v, T* w)
{

  int   i, ierr;
  char* type = "N";

  for (i=0; i<n; i++) {
    ctemp[i] = arcomplex<T>(v[i],0.0);
  }
  gttrs(type, n, 1, Adl, Ad, Adu, Adu2, ipiv, ctemp, n, ierr);
  for (i=0; i<n; i++) {
    w[i]=imag(ctemp[i]);
  }

} // MultOPvIm


template<class T>
inline NonSymGenProblemC<T>::
NonSymGenProblemC(int nx, T shiftRp, T shiftIp): A(nx), B(nx)
{

  shiftR  = shiftRp;
  shiftI  = shiftIp;
  decsize = 0;
  Ad      = 0;
  Adl     = 0;
  Adu     = 0;
  Adu2    = 0;
  ipiv    = 0;
  ctemp   = 0;
  n       = A.ncols();
  FactorOP();

} // Constructor.


template<class T>
inline NonSymGenProblemC<T>::~NonSymGenProblemC()
{

  FactorDataDeallocate();

} // Destructor.


#endif // NGENPRBC_H

