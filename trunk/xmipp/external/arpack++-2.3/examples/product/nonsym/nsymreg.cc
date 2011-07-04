/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE NSymReg.cc.
   Example program that illustrates how to solve a real
   nonsymmetric standard eigenvalue problem in regular mode
   using the ARNonSymStdEig class.

   1) Problem description:

      In this example we try to solve A*x = x*lambda in regular mode,
      where A is obtained from the standard central difference
      discretization of the convection-diffusion operator
                     (Laplacian u) + rho*(du / dx)
      on the unit square [0,1]x[0,1] with zero Dirichlet boundary
      conditions.

   2) Data structure used to represent matrix A:

      When using ARNonSymStdEig, the user is required to provide a
      class that contains a member function which computes the
      the matrix-vector product w = Av. In this example, this class 
      is called NonSymMatrixA, and MultMv is the function.

   3) Included header files:

      File             Contents
      -----------      -------------------------------------------
      nmatrixa.h       The NonSymMatrixA class definition.
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

#include "nmatrixa.h"
#include "nsymsol.h"
#include "arsnsym.h"


template<class T>
void Test(T type)
{

  // Creating a nonsymetric matrix.

  NonSymMatrixA<T> A(10); // n = 10*10.

  // Defining what we need: the four eigenvectors of A with largest magnitude.
  // A.MultMv is the function that performs the product w <- A.v.

  ARNonSymStdEig<T, NonSymMatrixA<T> >
  dprob(A.ncols(), 4, &A, &NonSymMatrixA<T>::MultMv);

  // Finding eigenvalues and eigenvectors.

  dprob.FindEigenvectors();

  // Printing solution.

  Solution(A, dprob);

} // Test.


int main()
{

  // Solving a single precision problem with n = 100.

  Test((float)0.0);

  // Solving a double precision problem with n = 100.

  Test((double)0.0);

} // main


