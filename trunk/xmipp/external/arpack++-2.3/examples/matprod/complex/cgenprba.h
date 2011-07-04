/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE CGenPrbA.h
   Very simple template class intended to illustrate how to
   use ARPACK++ to find some few eigenvalues and eigenvectors
   of complex generalized problems in regular mode.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef CGENPRBA_H
#define CGENPRBA_H

#include "arcomp.h"
#include "blas1c.h"
#include "cmatrixc.h"
#include "cmatrixd.h"


template <class T>
struct ComplexGenProblemA {

  ComplexMatrixC<T> A;
  ComplexMatrixD<T> B;

  void MultOPv(arcomplex<T>* v, arcomplex<T>* w)
  // Matrix vector subroutine where the matrix is inv(B)*A.
  {

    A.MultMv(v, w);
    B.SolveM(w);

  } // MultOPv.

  ComplexGenProblemA(int nx): A(nx), B(nx) 
  // Constructor.
  {
    B.FactorM();
  }

}; // struct ComplexGenProblemA.

#endif // CGENPRBA_H

