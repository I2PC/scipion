/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE SymGShft.cc.
   Example program that illustrates how to solve a real symmetric
   generalized eigenvalue problem in shift and invert mode using
   the ARSymGenEig class.

   1) Problem description:

      In this example we try to solve A*x = B*x*lambda in shift and
      invert mode, where A and B are obtained from the finite element
      discretrization of the 1-dimensional discrete Laplacian
                                d^2u / dx^2
      on the interval [0,1] with zero Dirichlet boundary conditions
      using piecewise linear elements.

   2) Data structure used to represent matrices A and B:

      When using ARSymGenEig, the user is required to provide some 
      classes that contain as member functions the matrix-vector
      products w = OP*Bv = inv(A-sigma*B)*B*v and w = B*v. 
      In this example, SymGenProblemB is a class that contains two 
      member functions, MultOPv and MultBv. The first takes a vector v 
      and returns the product OPv. The second performs the product Bv.

   3) Included header files:

      File             Contents
      -----------      -------------------------------------------
      sgenprbb.h       The SymGenProblemB class definition.
      argsym.h         The ARSymGenEig class definition.
      symgsol.h        The Solution function.

   4) ARPACK Authors:

      Richard Lehoucq
      Kristyn Maschhoff
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#include "sgenprbb.h"
#include "symgsol.h"
#include "argsym.h"


template<class T>
void Test(T type)
{

  // Creating a symmetric generalized problem with n = 100.

  SymGenProblemB<T> P(100, 0.0);

  // Defining what we need: the four eigenvectors nearest to 0.0.
  // P.MultOPv is the function that performs the product w <- OPv.
  // P.MultBv is the function that performs the product w <- Bv.
  // 'S' is passed as a parameter to indicate that we will use the 
  // shift and invert mode. 

  ARSymGenEig<T, SymGenProblemB<T>, SymGenProblemB<T> >
    dprob('S', P.A.ncols(), 4L, &P, &SymGenProblemB<T>::MultOPv,
          &P, &SymGenProblemB<T>::MultBv, 0.0);

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


