/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE NSymGSCR.cc.
   Example program that illustrates how to solve a real
   nonsymmetric generalized eigenvalue problem in complex shift
   and invert mode (taking the real part of OP*x) using the 
   ARNonSymGenEig class.

   1) Problem description:

      In this example we try to solve A*x = B*x*lambda in shift and
      invert mode, where A is the tridiagonal matrix with 2 on the
      diagonal, -2 on the subdiagonal and 3 on the superdiagonal, and
      B is the tridiagonal matrix with 4 on the diagonal and 1 on the
      off-diagonals.
      The shift sigma is a complex number.

   2) Data structure used to represent matrix A:

      When using ARNonSymGenEig, the user is required to provide some 
      classes that contain as member functions the matrix-vector
      products w = OP*Bv = real{inv(A-sigma*B)}*B*v, w = A*v and 
      w = B*v. In this example, NonSymGenProblemC is a class that 
      contains three member functions, MultOPvRe, MultAv and MultBv. 
      The first takes a vector v and returns the product OPv. The 
      second performs the product Av. The last one performs Bv.

   3) Included header files:

      File             Contents
      -----------      -------------------------------------------
      ngenprbc.h       The NonSymGenProblemC class definition.
      argnsym.h        The ARNonSymGenEig class definition.
      nsymgsol.h       The Solution function.

   4) ARPACK Authors:

      Richard Lehoucq
      Kristyn Maschhoff
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#include "ngenprbc.h"
#include "nsymgsol.h"
#include "argnsym.h"


template<class T>
void Test(T type)
{

  // Creating Eig A*x = lambda*B*x.

  NonSymGenProblemC<T> P(100, 0.4, 0.6); // n = 100, sigma = (0.4, 0.6).

  // Defining what we need: the four eigenvectors nearest to sigma.
  // P.MultOPvRe is the function that performs the product w <- OPv.
  // P.MultAv is the function that performs the product w <- Av.
  // P.MultBv is the function that performs the product w <- Bv.
  // 'R' is passed as a parameter to indicate that we will use only 
  // the real part of OP*v. 

  ARNonSymGenEig<T, NonSymGenProblemC<T>, NonSymGenProblemC<T> >
    dprob(P.A.ncols(), 4, &P, &NonSymGenProblemC<T>::MultOPvRe,
          &P, &NonSymGenProblemC<T>::MultAv, &P,
          &NonSymGenProblemC<T>::MultBv, 'R', 0.4, 0.6);

  // Finding eigenvalues and eigenvectors.

  dprob.FindEigenvectors();

  // Printing solution.

  Solution(P.A, P.B, dprob);

} // Test.


int main()
{

  // Solving a single precision problem with n = 100.

#ifndef __SUNPRO_CC

  Test((float)0.0);

#endif 

  // Solving a double precision problem with n = 100.

  Test((double)0.0);

} // main







