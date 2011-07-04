/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE CompGShf.cc.
   Example program that illustrates how to solve a complex
   generalized eigenvalue problem in shift and invert mode 
   using the ARCompGenEig class.

   1) Problem description:

      In this example we try to solve A*x = B*x*lambda in shift and
      invert mode, where A and B are derived from a finite element 
      discretization of a 1-dimensional convection-diffusion operator
                         (d^2u/dx^2) + rho*(du/dx)
      on the interval [0,1] with zero boundary conditions using 
      piecewise linear elements.

   2) Data structure used to represent matrix A:

      When using ARCompGenEig, the user is required to provide some 
      classes that contain as member functions the matrix-vector
      products w = OP*Bv = inv(A-sigma*B)*B*v and w = B*v. 
      In this example, ComplexGenProblemB is a class that contains two 
      member functions, MultOPv and MultBv. The first takes a vector v 
      and returns the product OPv. The second performs the product Bv.

   3) Included header files:

      File             Contents
      -----------      -------------------------------------------
      cgenprbb.h       The ComplexGenProblemB class definition.
      argcomp.h        The ARCompGenEig class definition.
      compgsl.h        The Solution function.
      arcomp.h         The "arcomplex" (complex) type definition.

   4) ARPACK Authors:

      Richard Lehoucq
      Kristyn Maschhoff
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#include "arcomp.h"
#include "cgenprbb.h"
#include "compgsol.h"
#include "argcomp.h"


template<class T>
void Test(T type)
{

  // Creating a complex problem with n = 100, rho = 10, sigma = 1.

  ComplexGenProblemB<T>
    P(100, arcomplex<T>(10.0,0.0), arcomplex<T>(1.0,0.0));

  // Defining what we need: the four eigenvectors nearest to sigma.
  // P.MultOPv is the function that performs the product w <- OPv.
  // P.MultBv is the function that performs the product w <- Bv.

  ARCompGenEig<T, ComplexGenProblemB<T>, ComplexGenProblemB<T> >
    dprob(P.A.ncols(), 4, &P, &ComplexGenProblemB<T>::MultOPv, &P,
          &ComplexGenProblemB<T>::MultBv, arcomplex<T>(1.0,0.0));

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


