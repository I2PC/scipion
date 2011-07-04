/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE SymShft.cc.
   Example program that illustrates how to solve a real symmetric
   standard eigenvalue problem in shift and invert mode using the
   ARSymStdEig class.

   1) Problem description:

      In this example we try to solve A*x = x*lambda in shift and
      invert mode, where A is derived from the central difference
      discretization of the 1-dimensional Laplacian on [0,1] with
      zero Dirichlet boundary conditions.

   2) Data structure used to represent matrix A:

      When using ARSymStdEig, the user is required to provide a
      class that contains a member function which computes the
      matrix-vector product w = OPv, where OP = inv[A - sigma*I].
      In this example, this class is called SymMatrixB, and 
      MultOPv is the function.

   3) Included header files:

      File             Contents
      -----------      -------------------------------------------
      smatrixb.h       The SymMatrixB class definition.
      arssym.h         The ARSymStdEig class definition.
      symsol.h         The Solution function.

   4) ARPACK Authors:

      Richard Lehoucq
      Kristyn Maschhoff
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#include "arssym.h"
#include "smatrixb.h"
#include "symsol.h"


template<class T>
void Test(T type)
{

  // Creating a symmetric matrix.

  SymMatrixB<T> A(100,0.0); // n = 100, shift = 0.0.

  // Defining what we need: the four eigenvectors of B nearest to 0.0.
  // A.MultOPv is the function that performs the product w <- OPv.

  ARSymStdEig<T, SymMatrixB<T> >
  dprob(A.ncols(), 4, &A, &SymMatrixB<T>::MultOPv, 0.0);

  // Finding eigenvalues and eigenvectors.

  dprob.FindEigenvectors();

  // Printing solution.

  Solution(A, dprob);

} // Test.


int main()
{

  // Solving a double precision problem with n = 100.

  Test((double)0.0);

  // Solving a single precision problem with n = 100.

  Test((float)0.0);

} // main

