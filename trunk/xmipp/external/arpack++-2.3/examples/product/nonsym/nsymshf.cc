/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE NSymShf.cc.
   Example program that illustrates how to solve a real
   nonsymmetric standard eigenvalue problem in shift and invert
   mode using the ARNonSymStdEig class.

   1) Problem description:

      In this example we try to solve A*x = x*lambda in regular mode,
      where A is derived from the centered difference discretization
      of the 1-dimensional convection-diffusion operator
                        (d^2u / dx^2) + rho*(du/dx)
      on the interval [0,1] with zero Dirichlet boundary condition.
      The shift sigma is a real number.

   2) Data structure used to represent matrix A:

      When using ARNonSymStdEig, the user is required to provide a
      class that contains a member function which computes the
      matrix-vector product w = OPv, where OP = inv[A - sigma*I].
      In this example, this class is called NonSymMatrixB, and 
      MultOPv is the function.

   3) Included header files:

      File             Contents
      -----------      -------------------------------------------
      nmatrixb.h       The NonSymMatrixB class definition.
      arsnsym.h        The ARNonSymStdEig class definition.
      nsymsol.h        The Solution function.

   4) ARPACK Authors:

      Richard Lehoucq
      Kristyn Maschhoff
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#include "nmatrixb.h"
#include "nsymsol.h"
#include "arsnsym.h"


template<class T>
void Test(T type)
{

  // Creating a nonsymmetric matrix.

  NonSymMatrixB<T> A(100, 1.0, 10.0); // n = 100, shift = 1, rho = 10.

  // Defining what we need: the four eigenvectors of A nearest to 1.0.
  // A.MultOPv is the function that performs the product w <- OPv.

  ARNonSymStdEig<T, NonSymMatrixB<T> >
  dprob(A.ncols(), 4, &A, &NonSymMatrixB<T>::MultOPv, 1.0);

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

