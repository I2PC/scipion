/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE NSymGShf.cc.
   Example program that illustrates how to solve a real
   nonsymmetric generalized eigenvalue problem in real shift
   and invert mode using the ARNonSymGenEig class.

   1) Problem description:

      In this example we try to solve A*x = B*x*lambda in shift and
      invert mode, where A and B are derived from the finite element
      discretization of the 1-dimensional convection-diffusion operator
                         (d^2u / dx^2) + rho*(du/dx)
      on the interval [0,1] with zero Dirichlet boundary conditions
      using linear elements.
      The shift sigma is a real number.

   2) Data structure used to represent matrix A:

      When using ARNonSymGenEig, the user is required to provide some 
      classes that contain as member functions the matrix-vector
      products w = OP*Bv = inv(A-sigma*B)*B*v and w = B*v. 
      In this example, NonSymGenProblemB is a class that contains two 
      member functions, MultOPv and MultBv. The first takes a vector v 
      and returns the product OPv. The second performs the product Bv.

   3) Included header files:

      File             Contents
      -----------      -------------------------------------------
      ngenprbb.h       The NonSymGenProblemB class definition.
      argnsym.h        The ARNonSymGenEig class definition.
      nsymgsol.h        The Solution function.

   4) ARPACK Authors:

      Richard Lehoucq
      Kristyn Maschhoff
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#include "ngenprbb.h"
#include "nsymgsol.h"
#include "argnsym.h"


template<class T>
void Test(T type)
{

  // Creating Eig A*x = lambda*B*x.

  NonSymGenProblemB<T> P(100, 10, 1.0); // n = 100, rho = 10, sigma = 1.

  // Defining what we need: the four eigenvectors nearest to 1.0.
  // P.MultOPv is the function that performs the product w <- OPv.
  // P.MultBv is the function that performs the product w <- Bv.

  ARNonSymGenEig<T, NonSymGenProblemB<T>, NonSymGenProblemB<T> >
    dprob(P.A.ncols(), 4, &P, &NonSymGenProblemB<T>::MultOPv,
          &P, &NonSymGenProblemB<T>::MultBv, 1.0);

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


