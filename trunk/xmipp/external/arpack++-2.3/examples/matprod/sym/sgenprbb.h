/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE SGenPrbB.h
   Very simple template class intended to illustrate how to 
   use ARPACK++ to find some few eigenvalues and eigenvectors 
   of symmetric generalized problems in shift and invert,
   buckling and Cayley modes.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef SGENPRBB_H
#define SGENPRBB_H

#include "blas1c.h"
#include "lapackc.h"
#include "smatrixc.h"
#include "smatrixd.h"

template<class ART>
class SymGenProblemB {

 private:

  int  n, decsize;
  int  *ipiv;
  ART  shift;
  ART  *Ad, *Adl, *Adu, *Adu2;

  void FactorDataDeallocate();
  // Eliminates the data structure used on matrix factorization.

  void FactorOP();
  // Factors (A-shift*B).

 public:

  SymMatrixC<ART> A;
  SymMatrixD<ART> B;

  void MultAv(ART* v, ART* w) { A.MultMv(v,w); }
  // Matrix-vector product w <- A*v,
  // where the matrix is the 1 dimensional discrete Laplacian on
  // the interval [0,1] with zero Dirichlet boundary conditions.

  void MultBv(ART* v, ART* w) { B.MultMv(v,w); }
  // Matrix-vector product w <- B*v,
  // where the matrix is the 1 dimensional mass matrix
  // on the interval [0,1].

  void MultOPv(ART* v, ART* w);
  // Matrix-vector product w <- inv(A-shift*B)*v.

  SymGenProblemB(int nx, ART shiftp);
  // Constructor.

  ~SymGenProblemB();
  // Destructor

}; // SymGenProblemB


template<class ART>
inline void SymGenProblemB<ART>::FactorDataDeallocate()
{

  delete[] Ad;
  delete[] Adl;
  delete[] Adu;
  delete[] Adu2;
  delete[] ipiv;

} // FactorDataDeallocate.


template<class ART>
void SymGenProblemB<ART>::FactorOP()
{

  int  i, ierr;
  ART  h, r1, r2;

  const ART one  = 1.0;
  const ART two  = 2.0;
  const ART four = 4.0;
  const ART six  = 6.0;

  if (decsize != n) {
    decsize = n;
    FactorDataDeallocate();
    Ad   = new ART[n];
    Adl  = new ART[n];
    Adu  = new ART[n];
    Adu2 = new ART[n];
    ipiv = new int[n];
  }

  h  = one/ART(n+1);
  r1 =  two/h - shift*h*four/six;
  r2 = -one/h - shift*h*one/six;

  for (i=0; i<n; i++) {
    Ad[i]  = r1;
    Adl[i] = r2;
  }

  copy(n, Adl, 1, Adu, 1);
  gttrf(n, Adl, Ad, Adu, Adu2, ipiv, ierr);

} // FactorOP.


template<class ART>
void SymGenProblemB<ART>::MultOPv(ART* v, ART* w)
{

  int  ierr;
  char *type = "N";

  copy(n, v, 1, w, 1);
  gttrs(type, n, 1, Adl, Ad, Adu, Adu2, ipiv, w, n, ierr);

} // MultOPv.


template<class ART>
inline SymGenProblemB<ART>::SymGenProblemB(int nx, ART shiftp): A(nx), B(nx)
{

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


template<class ART>
inline SymGenProblemB<ART>::~SymGenProblemB()
{

  FactorDataDeallocate();

} // Destructor.


#endif // SGENPRBB_H

