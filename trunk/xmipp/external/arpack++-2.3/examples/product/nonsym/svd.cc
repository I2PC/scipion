/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE SVD.cc.
   Example program that illustrates how to use ARPACK++ 
   to compute the singular value decomposition of a real
   nonsymmetric matrix using the ARNonSymStdEig class.

   1) Problem description:

      In this example we try to solve (A'*A)*x = x*lambda in 
      regular mode (where A is a m by n real matrix) in order
      to find a few of the largest singular values of A and 
      the corresponding right singular vectors. 
      In this example, m is greater than n. To find the SVD
      of a matrix with m < n, the user must change the roles
      of A and A' in the above equation.

   2) Data structure used to represent matrix A:

      When using ARNonSymStdEig, the user is required to provide a
      class that contains a member function which computes the
      matrix-vector product w = (A'*A)*v. In this example, this 
      class is called NonSymMatrixV, and MultOPv is the function.

   3) Included header files:

      File             Contents
      -----------      -------------------------------------------
      nmatrixv.h       The NonSymMatrixV class definition.
      arsnsym.h        The ARNonSymStdEig class definition.
      nsymvsol.h       The Solution function.

   4) ARPACK Authors:

      Richard Lehoucq
      Kristyn Maschhoff
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#include "nmatrixv.h"
#include "nsymvsol.h"
#include "arsnsym.h"


template<class T>
void Test(T type)
{

  // Creating a nonsymetric matrix.

  NonSymMatrixV<T> A(500, 100); // m = 500, n = 100.

  // Defining what we need: the four eigenvectors of A'*A with largest 
  // magnitude. A.MultOPv is the function that performs the product 
  // w <- (A'*A).v.

  ARNonSymStdEig<T, NonSymMatrixV<T> >
  dprob(A.ncols(), 4, &A, &NonSymMatrixV<T>::MultOPv);

  // Finding eigenvalues and eigenvectors.

  dprob.FindEigenvectors();

  // Printing singular values, sigma, and the residual norm 
  // ||A*v - sigma*u||, where v are the right and u the left 
  // singular vectors corresponding to sigma.

  Solution(A, dprob);

} // Test.


int main()
{

  // Solving a double precision problem with n = 100.

  Test((float)0.0);

  // Solving a single precision problem with n = 100.

  Test((double)0.0);

} // main


