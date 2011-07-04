/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE RSymReg.cc.
   Example program that illustrates how to solve a real symmetric
   standard eigenvalue problem in regular mode using the
   ARrcSymStdEig class.

   1) Problem description:

      In this example we try to solve A*x = x*lambda in regular mode,
      where A is derived from the central difference discretization
      of the 2-dimensional Laplacian on the unit square [0,1]x[0,1]
      with zero Dirichlet boundary conditions.

   2) Data structure used to represent matrix A:

      ARrcSymStdEig is a class that requires the user to provide a
      way to perform the matrix-vector product w = Av. In this
      example a class called SymMatrixA was created with this purpose.
      SymMatrixA contains a member function, MultMv(v,w), that takes a
      vector v and returns the product Av in w.

   3) The reverse communication interface:

      This example uses the reverse communication interface, which
      means that the desired eigenvalues cannot be obtained directly
      from an ARPACK++ class.
      Here, the overall process of finding eigenvalues by using the
      Arnoldi method is splitted into two parts. In the first, a
      sequence of calls to a function called TakeStep is combined
      with matrix-vector products in order to find an Arnoldi basis.
      In the second part, an ARPACK++ function like FindEigenvectors
      (or EigenValVectors) is used to extract eigenvalues and
      eigenvectors.

   4) Included header files:

      File             Contents
      -----------      -------------------------------------------
      smatrixa.h       The SymMatrixA class definition.
      arrssym.h        The ARrcSymStdEig class definition.
      rsymsol.h        The Solution function.

   5) ARPACK Authors:

      Richard Lehoucq
      Kristyn Maschhoff
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#include "arrssym.h"
#include "smatrixa.h"
#include "rsymsol.h"


template<class T>
void Test(T type)
{

  // Defining a matrix.

  SymMatrixA<T> A(10); // n = 10*10 is the dimension of the problem.

  // Creating a symmetric eigenvalue problem and defining what we need:
  // the four eigenvectors of A with smallest magnitude.

  ARrcSymStdEig<T> prob(A.ncols(), 4L, "SM");

  // Finding an Arnoldi basis.

  while (!prob.ArnoldiBasisFound()) {

    // Calling ARPACK FORTRAN code. Almost all work needed to
    // find an Arnoldi basis is performed by TakeStep.

    prob.TakeStep();

    if ((prob.GetIdo() == 1)||(prob.GetIdo() == -1)) {

      // Performing matrix-vector multiplication.
      // In regular mode, w = Av must be performed whenever
      // GetIdo is equal to 1 or -1. GetVector supplies a pointer
      // to the input vector, v, and PutVector a pointer to the
      // output vector, w.

      A.MultMv(prob.GetVector(), prob.PutVector());

    }

  }

  // Finding eigenvalues and eigenvectors.

  prob.FindEigenvectors();

  // Printing solution.

  Solution(prob);

} // Test


int main()
{

  // Solving a double precision problem with n = 100.

  Test((double)0.0);

  // Solving a single precision problem with n = 100.

  Test((float)0.0);

} // main

