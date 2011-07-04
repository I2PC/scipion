/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE RSVD.cc.
   Example program that illustrates how to use ARPACK++ 
   to compute the singular value decomposition of a real
   nonsymmetric matrix using the ARrcNonSymStdEig class.

   1) Problem description:

      In this example we try to solve (A'*A)*x = x*lambda in 
      regular mode (where A is a m by n real matrix) in order
      to find a few of the largest singular values of A and 
      the corresponding right singular vectors. 
      In this example, m is greater than n. To find the SVD
      of a matrix with m < n, the user must change the roles
      of A and A' in the above equation.

   2) Data structure used to represent matrix A:

      ARrcNonSymStdEig is a class that requires the user to provide
      a way to perform the matrix-vector product w = (A'*A)v. In this
      example a class called NonSymMatrixA was created with this 
      purpose. NonSymMatrixA contains a member function, MultOPv(v,w), 
      that takes a vector v and returns the product (A'*A)*v in w.

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
      nmatrixv.h       The NonSymMatrixV class definition.
      arrsnsym.h       The ARrcNonSymStdEig class definition.
      rnsymvsl.h       The Solution function.

   5) ARPACK Authors:

      Richard Lehoucq
      Kristyn Maschhoff
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#include "nmatrixv.h"
#include "rnsymvsl.h"
#include "arrsnsym.h"


template<class T>
void Test(T type)
{

  // Defining a nonsymetric matrix.

  NonSymMatrixV<T> A(500, 100); // m = 500, n = 100.

  // Creating a nonsymmetric eigenvalue problem and defining what we need:
  // the four eigenvectors of (A'*A) with largest magnitude.

  ARrcNonSymStdEig<T> prob(A.ncols(), 4L);

  // Finding an Arnoldi basis.

  while (!prob.ArnoldiBasisFound()) {

    // Calling ARPACK FORTRAN code. Almost all work needed to
    // find an Arnoldi basis is performed by TakeStep.

    prob.TakeStep();

    if ((prob.GetIdo() == 1)||(prob.GetIdo() == -1)) {

      // Performing matrix-vector multiplication.
      // In regular mode, w = (A'*A)*v must be performed whenever
      // GetIdo is equal to 1 or -1. GetVector supplies a pointer
      // to the input vector, v, and PutVector a pointer to the
      // output vector, w.

      A.MultOPv(prob.GetVector(), prob.PutVector());

    }

  }

  // Finding eigenvalues and eigenvectors.

  prob.FindEigenvectors();

  // Printing singular values.

  Solution(prob);

} // Test.


int main()
{

  // Solving a double precision problem with n = 100.

  Test((float)0.0);

  // Solving a single precision problem with n = 100.

  Test((double)0.0);

} // main

