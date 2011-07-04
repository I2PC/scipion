/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE SymReg.cc.
   Example program that illustrates how to solve a real symmetric
   standard eigenvalue problem in regular mode using the
   ARSymStdEig class.

   1) Problem description:

      In this example we try to solve A*x = x*lambda in regular mode,
      where A is derived from the central difference discretization
      of the 2-dimensional Laplacian on the unit square [0,1]x[0,1]
      with zero Dirichlet boundary conditions.

   2) Data structure used to represent matrix A:

      When using ARSymStdEig, the user is required to provide a
      class that contains a member function which computes the
      the matrix-vector product w = Av. In this example, this class 
      is called SymMatrixA, and MultMv is the function.

   3) Included header files:

      File             Contents
      -----------      -------------------------------------------
      smatrixa.h       The SymMatrixA class definition.
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
#include "smatrixa.h"
#include "symsol.h"

template<class T>
void Test(T type)
{

  // Creating a symmetric matrix.

  SymMatrixA<T> A(10L); // n = 10*10.

  // Defining what we need: the four eigenvectors of A with smallest magnitude.
  // A.MultMv is the function that performs the product w <- A.v.

  ARSymStdEig<T, SymMatrixA<T> >
    dprob(A.ncols(), 4L, &A, &SymMatrixA<T>::MultMv, "SM");

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

} // main.


