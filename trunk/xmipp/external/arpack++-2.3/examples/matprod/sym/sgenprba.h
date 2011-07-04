/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE SGenPrbA.h
   Very simple template class intended to illustrate how to
   use ARPACK++ to find some few eigenvalues and eigenvectors
   of symmetric generalized problems in regular mode.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef SGENPRBA_H
#define SGENPRBA_H

#include "blas1c.h"
#include "smatrixc.h"
#include "smatrixd.h"


template <class ART>
struct SymGenProblemA {

  SymMatrixC<ART> A;
  SymMatrixD<ART> B;

  void MultOPv(ART* v, ART* w)
  // Matrix-vector multiplication w <- inv(B)*A*v.
  {
    A.MultMv(v, w);
    copy(A.ncols(), w, 1, v, 1);
    B.SolveM(w);
  }

  SymGenProblemA(int nx): A(nx), B(nx) 
  // Constructor.
  {
    B.FactorM();
  }

}; // struct SymGenProblemA.

#endif // SGENPRBA_H

