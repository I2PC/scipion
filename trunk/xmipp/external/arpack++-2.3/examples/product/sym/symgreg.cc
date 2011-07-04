/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE SymGReg.cc
   Example program that illustrates how to solve a real symmetric
   generalized eigenvalue problem in regular mode using the
   ARSymGenEig class.

   1) Problem description:

      In this example we try to solve A*x = B*x*lambda in regular mode,
      where A and B are obtained from the finite element discretization
      of the 1-dimensional discrete Laplacian
                                  d^2u / dx^2
      on the interval [0,1] with zero Dirichlet boundary conditions
      using piecewise linear elements.

   2) Data structure used to represent matrices A and B:

      When using ARSymGenEig, the user is required to provide some
      classes that contain as member functions the matrix-vector
      products w = OPv = inv(B)*A*v and w = B*v. 
      In this example, SymGenProblemA is a class that contains a 
      member function, MultOPv(v,w), that takes a vector v and returns 
      the product OPv in w. It also contains an object, B, that stores
      matrix B data. The product Bv is performed by MultMv, a member
      function of B.

   3) Included header files:

      File             Contents
      -----------      -------------------------------------------
      sgenprba.h       The SymGenProblemA class definition.
      argsym.h         The ARSymGenEig class definition.
      symgsol.h        The Solution function.

   5) ARPACK Authors:

      Richard Lehoucq
      Kristyn Maschhoff
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#include "sgenprba.h"
#include "symgsol.h"
#include "argsym.h"

template<class T>
void Test(T type)
{


  // Creating Eig A*x = lambda*B*x.

  SymGenProblemA<T> P(100); // n = 100.

  // Defining what we need: the four eigenvectors with largest magnitude.
  // P.MultOPv is the function that performs the product w <- OPv.
  // P.B.MultMv is the function that performs the product w <- Bv.

  ARSymGenEig<T, SymGenProblemA<T>, SymMatrixD<T> >
    dprob(P.A.ncols(), 4, &P, &SymGenProblemA<T>::MultOPv,
          &P.B, &SymMatrixD<T>::MultMv);

  // Finding eigenvalues and eigenvectors.

  dprob.FindEigenvectors();

  // Printing solution.

  Solution(P.A, P.B, dprob);

} // Test.


int main()
{

  // Solving a double precision problem with n = 100.

  Test((double)0.0);

  // Solving a single precision problem with n = 100.

  Test((float)0.0);

} // main


