/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE NGenPrbA.h
   Very simple template class intended to illustrate how to
   use ARPACK++ to find some few eigenvalues and eigenvectors
   of nonsymmetric generalized problems in regular mode.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef NGENPRBA_H
#define NGENPRBA_H

#include "blas1c.h"
#include "nmatrixc.h"
#include "nmatrixd.h"


template <class T>
struct NonSymGenProblemA {

  NonSymMatrixC<T> A;
  NonSymMatrixD<T> B;

  void MultOPv(T* v, T* w)
  // Matrix vector subroutine where the matrix is inv(B)*A.
  {

    A.MultMv(v, w);
    B.SolveM(w);

  } // MultOPv.

  NonSymGenProblemA(int nx): A(nx), B(nx) 
  // Constructor.
  {
    B.FactorM();
  }

}; // struct NonSymGenProblemA.

#endif // NGENPRBA_H
