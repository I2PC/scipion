/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE NGenPrbB.h
   Very simple template class intended to illustrate how to
   use ARPACK++ to find some few eigenvalues and eigenvectors
   of nonsymmetric generalized problems in shift and invert mode.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef NGENPRBB_H
#define NGENPRBB_H

#include "blas1c.h"
#include "lapackc.h"
#include "nmatrixc.h"
#include "nmatrixd.h"

template<class T>
class NonSymGenProblemB {

 private:

  int  n, decsize;
  int  *ipiv;
  T    rho;
  T    shiftR, shiftI;
  T    *Ad, *Adl, *Adu, *Adu2;

  void FactorDataDeallocate();
  // Eliminates the data structure used on matrix factorization.

  void FactorOP();
  // Factors (A-shift*B).

 public:

  NonSymMatrixC<T> A;
  NonSymMatrixD<T> B;

  void MultAv(T* v, T* w) { A.MultMv(v,w); }
  // Matrix-vector multiplication w <- A*v, where A is a n by n
  // nonsymmetric tridiagonal matrix derived from the central difference
  // discretization of the 1-dimensional convection-diffusion operator
  // on the interval [0,1] with zero Dirichlet boundary condition.

  void MultBv(T* v, T* w) { B.MultMv(v,w); }
  // Matrix-vector multiplication w <- B*v, where B is the mass matrix
  // formed by using piecewise linear elements on [0,1].

  void MultOPv(T* v, T* w);
  // Matrix-vector multiplication w <- inv(A-shift*B)*v.

  NonSymGenProblemB(int nx, T rhop, T shiftRp, T shiftIp = 0.0);
  // Constructor.

  ~NonSymGenProblemB();
  // Destructor

}; // NonSymGenProblemB


template<class T>
inline void NonSymGenProblemB<T>::FactorDataDeallocate()
{

  delete[] Ad;
  delete[] Adl;
  delete[] Adu;
  delete[] Adu2;
  delete[] ipiv;

} // FactorDataDeallocate.


template<class T>
void NonSymGenProblemB<T>::FactorOP()
{

  int  j, ierr;
  T    h, s, s1, s2, s3;

  const T one = 1.0;
  const T two = 2.0;

  if (decsize != n) {
    decsize = n;
    FactorDataDeallocate();
    Ad   = new T[n];
    Adl  = new T[n];
    Adu  = new T[n];
    Adu2 = new T[n];
    ipiv = new int[n];
  }

  h  = one/T(n+1);
  s  = rho/two;
  s1 = -one/h - s - shiftR*h;
  s2 = two/h - 4.0*shiftR*h;
  s3 = -one/h + s - shiftR*h;

  for (j=0; j<n-1; j++) {
    Adl[j] = s1;
    Ad[j]  = s2;
    Adu[j] = s3;
  }
  Ad[n-1]  = s2;

  gttrf(n, Adl, Ad, Adu, Adu2, ipiv, ierr);

} // FactorOP.


template<class T>
void NonSymGenProblemB<T>::MultOPv(T* v, T* w)
{

  int  ierr;
  char *type = "N";

  copy(n, v, 1, w, 1);
  gttrs(type, n, 1, Adl, Ad, Adu, Adu2, ipiv, w, n, ierr);

} // MultOPv.


template<class T>
inline NonSymGenProblemB<T>::
NonSymGenProblemB(int nx, T rhop, T shiftRp, T shiftIp): A(nx), B(nx)
{

  rho     = rhop;
  shiftR  = shiftRp;
  shiftI  = shiftIp;
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
inline NonSymGenProblemB<T>::~NonSymGenProblemB()
{

  FactorDataDeallocate();

} // Destructor.


#endif // NGENPRBB_H

